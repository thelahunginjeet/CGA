#!/usr/bin/env python

import unittest, os, time, copy
from scipy import mean,log
from scipy.linalg import inv
import numpy as MATH
from numpy.random import rand as urand
from random import choice as rchoice
from CGAPreprocessing import Utilities
from CGAGenerator import CGAGenerator
from CGAParameters import CGAParameters
from CGALogging import Subject, DataLogger, SqliteLogger

# TODO : 
#  - write the simulation as a generator
#  - write a more functional Chromosome container class
#  - decorators to compute co-dependent fitness functions


class CGAChromosome(object):
    """CGA chromosomes are a container class that stores both the unit of selection (a tree in this problem) 
    and their fitness (multiple fitness criteria are stored in the list fitnessVals). Right now this is very 
    basic, but helps with bookkeeping in CGASimulation.  __cmp__ is overloaded to allow fitness comparisons 
    and population sorting based on fitness.  More functionality could be added as desired."""
    def __init__(self, tree=None, fitnessVals=None):
        self.tree = tree
        self.fitnessVals = fitnessVals
    
    def copy(self):
        return CGAChromosome(self.tree.copy(), copy.copy(self.fitnessVals))
    
    """
    def __cmp__(self, other):
        assert self.fitness is not None
        # add in some handling for NaN, -Inf, and Inf to automatically put them at the bottom
        return cmp(self.fitness, other.fitness)
    """
    

class CGASimulation(Subject):
    """Main class that does the simulation, records results, evaluates trees, etc."""
    def __init__(self, databaseFile, pdbFile, cgap=CGAParameters()):
        super(CGASimulation,self).__init__()
        # database with protein information
        if not os.path.exists(pdbFile):
            raise IOError, 'something is wrong with your pdb file; check yourself . . .'
        database = Utilities.readDatabase(databaseFile)
        self.indices = database['indices']
        self.singleFrequencies = database['singleFrequencies']
        self.jointFrequencies = database['jointFrequencies']
        self.prepare_data()
        # parameters - perform checks on the parameters
        self.cgap = cgap
        self.prepare_parameters()
        # distances and protein diameter
        self.distances = Utilities.calculateAtomicDistances(pdbFile)
        self.proteinDiameter = mean(self.distances.values()) - min(self.distances.values())
        self.proteinMinimum = min(self.distances.values())
        self.proteinSum = sum(self.distances.values())
        # will be a list of CGAChromosome objects
        self.population = []
        # running simulation time, not total time (that is a parameter)
        self.time = 0
        # 'private' dictionary of data shared by fitness functions
        self._sharedData = {}
    
    
    def prepare_parameters(self):
        """Prepares checks on the CGAParameters object to ensure all the information necessary to
        perform a simulation are contained therein."""
        # correct forestSize to be even
        if MATH.mod(self.cgap.treep['forestSize'],2):
            print 'Adjusting forest size to be even.'
            self.cgap.treep['forestSize'] = self.cgap.treep['forestSize'] + 1    
        # check mutation probabilities - only one kind of mutation is performed per offspring, 
        #    so pG + pP + pM + pC + pHC <= 1. Throw a warning if the sum is > 1
        if sum(self.cgap.mutation.values()) > 1.0:
            print 'WARNING : Sum of mutation parameters > 1; unexpected behavior may result!'
        # functions to evaluate and whether they are involved in selection
        self.fitStrings = []
        self.selectOn = []
        for k in self.cgap.fitness:
            if hasattr(self,'calculate_'+k):
                self.fitStrings.append(k)
                self.selectOn.append(self.cgap.fitness[k])
                
        
    def prepare_data(self):
        """Array-izes the input frequency data."""
        # change input frequencies from the db to 20X20 arrays (including single frequencies)
        for i in range(0,len(self.indices)):
            indxi = self.indices[i]
            self.singleFrequencies[indxi] = MATH.tile(MATH.asarray(self.singleFrequencies[indxi]),(1,20))
            for j in range(i+1,len(self.indices)):
                indxj = self.indices[j]
                self.jointFrequencies[(indxi,indxj)] = MATH.asarray(self.jointFrequencies[(indxi,indxj)])
    
    
    def initialize_tree(self):
        """Creates a single random tree using the parameters in self.treeGenDict.  This functionality
        has been separated from self.populate() in order to make headless chicken crossover (in which we
        crossover with a randomly generated tree) easier."""
        r = self.cgap.treep['r']
        if self.cgap.treep['treetype'] == 'exponential':
            if self.cgap.treep['p'] < 1.0:
                p = self.cgap.treep['p']
            else:
                # ensures consistent parameter values
                p = 1/self.cgap.treep['p']
            return CGAGenerator.expgenerate(p,r)
        elif self.cgap.treep['treetype'] == 'fixed':
            if self.cgap.treep['p'] > 1:
                treeSize = self.cgap.treep['p']
            else:
                treeSize = MATH.int(1.0/self.cgap.treep['p'])
            return CGAGenerator.generate(treeSize,r)
        else:
            raise TypeError, 'Unknown tree type %s' % self.cgap.treep['treetype']
        
        
    def populate(self):
        """Populate the forest of function trees.  You can choose to either use probabilistic tree
        generation (the number of nodes will be exponentially distributed) or trees with a fixed 
        number of non-terminal nodes.  The treeGenDict determines which method will be used, and
        incompatible parameters result in default behavior.
        """
        for i in range(self.cgap.treep['forestSize']):
            tree = self.initialize_tree()
            fitVals = self.evaluate_fitness(tree)
            self.population.append(CGAChromosome(tree,fitVals))


    def advance(self):
        """Step forward one step in time."""
        # log data BEFORE manipulating the population
        if MATH.remainder(MATH.int(self.time),self.cgap.timing['sampGen']) == 0:
            self.notify(time=self.time)
        # first select a round of parents
        parents = [self.select_parent(method=self.cgap.selection['method']) for x in xrange(self.cgap.treep['forestSize'])]
        offspring = list()
        for i in xrange(0, self.cgap.treep['forestSize']):
            offspring.append(self.mate(rchoice(parents),rchoice(parents)))
        # overwrite current forest
        self.population = offspring
        # advance time
        self.time += 1
        
            
    def mate(self, parentOne, parentTwo):
        """Accepts two parents and returns ONE offspring; if the offspring is unchanged from one of
        the parents, the fitness values are not re-evaluated, just copied.  Only one category of 
        mutation is performed on a single offspring."""
        # one parent will form the basis for the new offspring; the other is just there for 
        #    potential crossover (by default the mother will become the offspring)
        mother, father  = parentOne.copy(), parentTwo.copy()
        # fitEval will be set to True if any mutations occur that make the offspring different from
        #    its copied parent; this way we avoid unnecessary fitness evaluations
        fitEval = False
        # only one category of mutation is allowed per mating event
        r = urand()
        mutP = self.cgap.mutation
        if r < mutP['pC']: # parental crossover
            fitEval = True
            mNode = rchoice(mother.tree.getNodes())
            fNode = rchoice(father.tree.getNodes())
            CGAGenerator.single_crossover(mNode,fNode)
        elif r < mutP['pC'] + mutP['pHC']: # headless chicken crossover
            fitEval = True
            mNode = rchoice(mother.tree.getNodes())
            rNode = rchoice(self.initialize_tree().getNodes())
            CGAGenerator.single_crossover(mNode,rNode)
        elif r < mutP['pC'] + mutP['pHC'] + mutP['pM']: # point mutation (uses pM/node for mutation prob.)
            fitEval = True
            for n in mother.tree.getNodes():
                if urand() < mutP['pM']:
                    CGAGenerator.point_mutate(mother.tree, n)
        elif r < mutP['pC'] + mutP['pHC'] + mutP['pM'] + mutP['pP']: # pruning - guaranteed to do one pruning operation
            fitEval = True
            mNode = rchoice(mother.tree.getNodes())
            CGAGenerator.prune(mother.tree,mNode)
        elif r < mutP['pC'] + mutP['pHC'] + mutP['pM'] + mutP['pP'] + mutP['pG']: # growth - guaranteed to do one growth op.
            fitEval = True
            mTerm = rchoice(mother.tree.getTermini())
            CGAGenerator.grow(mother.tree,mTerm)
        else: # offspring will just be a copy
            pass
        if fitEval:
            mother.fitnessVals = self.evaluate_fitness(mother.tree)
        return mother
         
    
    def select_parent(self,method,**kwargs):
        """A dispatcher that implements a variety of selection methods and returns a single parent
        to place in the pool for the next generation.
            Current allowed methods (all strings):
                -'tournament' : requires a parameter k. two individuals are chosen at random; if 
                    rand < k, the fitter individual is selected.  if rand > k, the less fit one is.
        Regardless of selected method, the parent is returned as a CGAChromosome object."""
        mstring = self.cgap.selection['method'] + '_selection'
        if hasattr(self, mstring):
            parent = getattr(self,mstring)(**kwargs)
        else:
            parent = rchoice(self.population) # pick at random if there's a problem
        return parent
    
        
    def pareto_tournament_selection(self, k=0.75):
        '''Psuedo-Pareto tournament selection for multi-objective offspring.  Offspring that are
        better in a majority (2/3) of the fitness terms win the tournaments. For a single fitness
        function (to select on), this defaults to regular tournament selection.'''
        pOne, pTwo = rchoice(self.population), rchoice(self.population)
        f1 = MATH.array(pOne.fitnessVals)[self.selectOn]
        f2 = MATH.array(pTwo.fitnessVals)[self.selectOn]
        if urand() < k:
            return pOne if MATH.sign((f1-f2)).sum() > 0 else pTwo
        else:
            return pOne if MATH.sign((f1-f2)).sum() < 0 else pTwo
    
    
    def evaluate_fitness(self,tree,**kwargs):
        """A dispatcher to do fitness evaluations; the functions to be evaluated are
        stored as strings in self.fitStrings."""
        # clear the shared dict from the last calculation
        self._sharedData = {}
        fitVals = []
        for fs in self.fitStrings:
            fsToEval = 'calculate_'+fs
            fval = getattr(self,fsToEval)(tree,**kwargs)
            fitVals.append(fval)
        return fitVals 
    
    
    def calculate_wij(self,tree):
        """Almost all conceivable fitness functions compute a matrix of scores, so that 
        functionality is coded here."""
        weights = {}
        pi = [x for x in tree.getTermini() if x.string == 'p_i']
        pj = [x for x in tree.getTermini() if x.string == 'p_j']
        pij = [x for x in tree.getTermini() if x.string == 'p_ij']
        for i in self.indices:
            for j in [x for x in self.indices if x > i]:
                map(lambda x : x.replaceData(self.singleFrequencies[i]), pi)
                map(lambda x : x.replaceData(self.singleFrequencies[j]), pj)
                map(lambda x : x.replaceData(self.jointFrequencies[(i, j)]), pij)
                weights[(i, j)] = tree.getFunction()
        return weights

    
    def calculate_parsimony(self,tree):
        """Simplest possible parsimony term; the negative of the tree size."""
        if self._sharedData.has_key('nNodes'):
            return -1.0*self._sharedData['nNodes']
        nNodes = len(tree.getNodes())
        self._sharedData['nNodes'] = nNodes
        return -1.0*nNodes

    
    def calculate_finitewts(self,tree):
        """Computes the fraction of weights which are finite (not NaN, Inf, or -Inf)."""
        if self._sharedData.has_key('wij'):
            wij = self._sharedData['wij']
        else:
            wij = self.calculate_wij(tree)
            self._sharedData['wij'] = wij
        if self._sharedData.has_key('nFinite'):
            nFinite = self._sharedData['nFinite']
        else:
            nFinite = 0
            for i,j in wij:
                if (i,j) in self.distances:
                    if MATH.isfinite(wij[(i,j)]):
                        nFinite += 1
            self._sharedData['nFinite'] = nFinite
        return MATH.float(nFinite)/len(wij)
    
    # FIX THIS
    def calculate_topNm1_finitewts(self,tree):
        """Computes the fraction of weights which are finite (not NaN, Inf, or -Inf),
        restricted to just the largest N-1 weights."""
        if self._sharedData.has_key('wij'):
            wij = self._sharedData['wij']
        else:
            wij = self.calculate_wij(tree)
            self._sharedData['wij'] = wij
        finitewts = 0.0
        N = len(self.indices)
        wkeys,wvals = zip(*sorted(wij.iteritems(), key = lambda (k,v): (v,k))[-1:-N:-1])
        for i in range(0,len(wkeys)):
            if wkeys[i] in self.distances:
                if MATH.isfinite(wvals[i]):
                    finitewts += 1.0
        return finitewts/len(wvals)
        
    
    def calculate_distance_matrix(self,tree):
        """Accepts an input tree (member of the forest) and evaluates its fitness, defined
        as a transformed direct fit to the dimensionless distance matrix.  Specifically:
            fitness = -0.5*sum((a*w_ij + b - d_ij)^2)
        This maps fitness to [0,1], and allows an arbitrary linear rescaling of the scores.
        The optimal rescaling is a simple linear algebra subproblem.
        """
        if self._sharedData.has_key('wij'):
            wij = self._sharedData['wij']
        else:
            wij = self.calculate_wij(tree)
            self._sharedData['wij'] = wij
        # accumulators allow calculation of best linear transformation
        nFinite = 0
        SDW = 0.0
        W = []
        D = []
        for i,j in wij:
            if (i,j) in self.distances:
                if MATH.isfinite(wij[(i,j)]):
                    W.append(wij[(i,j)])
                    D.append((self.distances[(i,j)] - self.proteinMinimum)/self.proteinDiameter)
                    SDW += ((self.distances[(i,j)] - self.proteinMinimum)/self.proteinDiameter)*wij[(i,j)]
                    nFinite += 1
        # linear algebra to get transformation
        SW2 = sum([x**2 for x in W])
        SW = sum([x for x in W])
        SD = sum(D)
        try:
            rescale = MATH.dot(inv(MATH.array([[SW2,SW],[SW,len(W)]])),MATH.array([[SDW],[SD]])).flatten()
            if not MATH.isfinite(rescale).all():
                rescale = MATH.array([1.0,0.0])
        except:
            rescale = MATH.array([1.0,0.0]) # just give up if there's any numerical weirdness
        # apply transformation to get weights
        resid = MATH.array([rescale[0]*W[i] + rescale[1] - D[i] for i in range(0,len(W))])
        # residuals might be empty - if so, there's not a single finite weight
        if len(resid) == 0:
            fitness = -1.0*MATH.inf
        else:
            fitness = -0.5*((resid*resid).sum())
        if not self._sharedData.has_key('nFinite'):
            self._sharedData['nFinite'] = nFinite
        return fitness
        
    
    
    def calculate_weighted_accuracy(self,tree):
        """Accepts an input tree (member of the forest) and evaluates its fitness, currently 
        defined as:
            fitness = -1.0*sum(ws_ij*d_ij)/sum(ws_ij),
        where d_ij is the dimensionless matrix of (positive) distances, and ws_ij = w_ij + min(min(w_ij),0).
        The -1.0 multiplier insures better outcomes = increasing fitness.
        """
        nFinite = 0
        if self._sharedData.has_key('wij'):
            wij = self._sharedData['wij']
        else:
            wij = self.calculate_wij(tree)
            self._sharedData['wij'] = wij
        try:
            minw = MATH.min([x for x in wij.values() if MATH.isfinite(x)])
        except ValueError:
            minw = 0.0 # all weights are infinite or non-numeric
            fitness = MATH.inf
            return fitness
        minw = minw if minw < 0 else 0.0 # no min shift necessary if all weights positive
        accuracy = []
        normalization = minw*len(self.distances)
        for i,j in wij:
            if (i,j) in self.distances:
                if MATH.isfinite(wij[(i,j)]):
                    value = (wij[(i,j)] + minw)*((self.distances[(i,j)] - self.proteinMinimum)/self.proteinDiameter)
                    accuracy.append(value)
                    normalization += wij[(i,j)]
                    nFinite += 1
        # normalization should be a meaningful number
        fitness = -1.0*sum(accuracy)/normalization
        if not self._sharedData.has_key('nFinite'):
            self._sharedData['nFinite'] = nFinite
        return fitness
    
    
    def calculate_topNm1_weighted_accuracy(self,tree):
        """Accepts an input tree (member of the forest) and evaluates its fitness, currently 
        defined as:
            fitness = -1.0*sum(ws_ij*d_ij)/sum(ws_ij),
        where d_ij is the dimensionless matrix of (positive) distances, and ws_ij = w_ij + min(min(w_ij),0).
        The -1.0 multiplier insures better outcomes = increasing fitness.  This function is only computed
        using the N-1 largest scores, a degree of truncation typically used for these algorithms .
        """
        if self._sharedData.has_key('wij'):
            wij = self._sharedData['wij']
        else:
            wij = self.calculate_wij(tree)
            self._sharedData['wij'] = wij
        N = len(self.indices)
        # ignore everything but the top N weights
        wkeys,wvals = zip(*sorted(wij.iteritems(), key = lambda (k,v): (v,k))[-1:-N:-1])
        try:
            minw = MATH.min([x for x in wvals if MATH.isfinite(x)])
        except ValueError:
            minw = 0.0         # there are NO finite topN weights; no need to keep calculating
            fitness = MATH.inf # bad value
            return fitness
        # not all weights are nonsense
        minw = minw if minw < 0 else 0.0
        accuracy = []
        normalization = minw*len(self.distances)
        for i in range(0,len(wkeys)):
            if wkeys[i] in self.distances:
                if MATH.isfinite(wvals[i]):
                    value = (wvals[i] + minw)*((self.distances[wkeys[i]] - self.proteinMinimum)/self.proteinDiameter)
                    accuracy.append(value)
                    normalization += wvals[i]
        # normalization should be a meaningful number
        fitness = -1.0*sum(accuracy)/normalization
        return fitness
   
    

class CGASimulationTests(unittest.TestCase):
    def setUp(self):
        # created with default parameters
        self.cgap = CGAParameters()
        self.cgap.set('timing',timeSteps=10)
        self.cgap.set('treep',forestSize=6)
        self.cgap.set(fitness={'weighted_accuracy':True,'parsimony':True,'finitewts':False})
        self.mySimulation = CGASimulation('../tests/pdz_test.db', '../tests/1iu0.pdb',cgap=self.cgap)
        self.mySimulation.populate()
        #self.sqliteLogger = SqliteLogger('../tests')
        #self.mySimulation.attach(self.sqliteLogger)
    
    def testKnownTrees(self):
        print "\n\n----- calculating fitness of known trees -----"
        # create trees - MI, OMES, whatever - that have a special form and
        #    check their fitness
        treenames = ['MI','OMES']
        for t in treenames:
            tree = CGAGenerator.generate_special_tree(t)
            print '%s : %s' % (t,tree.getString())
            print self.mySimulation.evaluate_fitness(tree)
       
    def testAdvance(self):
        print "\n\n----- testing one-step advancement -----"
        print 'Before advancement:'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitnessVals) for x in self.mySimulation.population]
        self.mySimulation.advance()
        print 'After advancement (one step):'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitnessVals) for x in self.mySimulation.population]
 

    def testMultiAdvance(self):
        print "\n\n----- testing advancement of the tree over many steps -----"
        print 'Before advancement:'
        print 'Pop. size : ', len(self.mySimulation.population)
        print 'Fitness order : ', self.mySimulation.fitStrings
        n = 10
        print [(x.tree.getString(),x.fitnessVals) for x in self.mySimulation.population]
        for i in range(n):
            self.mySimulation.advance()
            print [(x.tree.getString(),x.fitnessVals) for x in self.mySimulation.population]
        print 'After advancement (%d steps):' % n
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitnessVals) for x in self.mySimulation.population]

       
    def testPopulation(self):
        print "\n\n----- testing population generation -----"
        print "Population size : ", len(self.mySimulation.population)
        for x in self.mySimulation.population:
            x.tree()
            print 'String rep : %s' % x.tree.getString()
            
    def testShape(self):
        print "\n\n----- testing shape conversion of input data -----"
        print self.mySimulation.indices
        print 'P_7(:,1)  : ', self.mySimulation.singleFrequencies[7][:,0]
        print 'sum(P_7(:,1)) : ', sum(self.mySimulation.singleFrequencies[7][:,0])
        print 'Dimensions of P_i :          %d X %d' % (self.mySimulation.singleFrequencies[7].shape)
        print 'Dimensions of P_ij :         %d X %d' % (self.mySimulation.jointFrequencies[(7,20)].shape)
        print 'Shape of entropy(P_7) :      %d X %d' % ((self.mySimulation.singleFrequencies[7]*log(self.mySimulation.singleFrequencies[7])).shape)
        print 'Shape of entropy(P_{7,20}) : %d X %d' % ((self.mySimulation.jointFrequencies[(7,20)]*log(self.mySimulation.jointFrequencies[(7,20)])).shape)
        print 'Distance(7,20) : ' , self.mySimulation.distances[(7, 20)]
        
    def testSelectParent(self):
        print"\n\n----- testing parental selection -----"
        t1 = time.clock()
        parent = self.mySimulation.select_parent(method='pareto_tournament')
        print 'Parent selected: %s, %s' % (parent.tree.getString(),parent.fitnessVals)
        print 'Elapsed time : %f sec' % (time.clock()-t1)
 
    """
    # bring this back in later when we actually are trying to optimize
    def testEvaluateFitness(self):
        print "\n----- testing and timing fitness evaluation -----"
        # double loop is to not count eval() timing against the old version
        nRepeats = 1
        for tree in self.mySimulation.population:
            tree()
            print 'String : ', tree.string
            # non-optimized version
            t1 = time.clock()
            for tree in self.mySimulation.population:
                for i in xrange(0,nRepeats):
                    self.mySimulation.evaluate_fitness(tree)
            t2 = time.clock()
            print 'Elapsed clock time (basic) : %f seconds' %(t2-t1)
            # optimized version
            t1 = time.clock()
            for tree in self.mySimulation.population:
                for i in xrange(0,nRepeats):
                    self.mySimulation.evaluate_fitness_iter(tree)
            t2 = time.clock()
            print 'Elapsed clock time (optimized) : %f seconds' %(t2-t1)
    """
        
if __name__ == '__main__':
    unittest.main()