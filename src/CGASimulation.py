#!/usr/bin/env python

import unittest, os, time
from scipy import mean,log
from scipy.linalg import inv
import numpy as MATH
from numpy.random import rand as urand
from random import choice as rchoice
from CGAPreprocessing import Utilities
from CGAGenerator import CGAGenerator
from CGALogging import Subject, DataLogger, SqliteLogger

# TODO : 
#  - write the simulation as a generator
#  - undo hard-coding of parameters passed into selection functions
#  - dangerous for CGAChromosome to do all comparisons based on fitness?

class CGAChromosome(object):
    """CGA chromosomes are a container class that stores both the unit of selection (a tree in this problem) 
    and their fitness (multiple fitness criteria are stored). Right now this is very basic, but helps with 
    bookkeeping in CGASimulation.  __cmp__ is overloaded to allow fitness comparisons and population sorting 
    based on fitness.  More functionality could be added as desired."""
    def __init__(self, tree=None, fitness=MATH.nan, parsimony=MATH.nan, finitewts=MATH.nan):
        self.tree = tree
        self.fitness = fitness
        self.parsimony = parsimony
        self.finitewts = finitewts
    
    def copy(self):
        return CGAChromosome(self.tree.copy(), self.fitness, self.parsimony, self.finitewts)

    def __cmp__(self, other):
        assert self.fitness is not None
        # add in some handling for NaN, -Inf, and Inf to automatically put them at the bottom
        return cmp(self.fitness, other.fitness)


class CGASimulation(Subject):
    """Main class that does the simulation, records results, evaluates trees, etc."""
    def __init__(self, databaseFile, pdbFile, treeGenDict = {'treetype':'fixed','p':5,'r':0.6}, fitnessMethod='distance_matrix', selectionMethod='pareto_tournament',elitism=True,forestSize=50, timeSteps=100, sampGen=1, pG = 0.025, pP = 0.025, pHC = 0.1, pM = 0.05, pC = 0.7):
        super(CGASimulation,self).__init__()
        if not os.path.exists(pdbFile):
            raise IOError, 'something is wrong with your pdb file; check yourself . . .'
        database = Utilities.readDatabase(databaseFile)
        self.indices = database['indices']
        self.singleFrequencies = database['singleFrequencies']
        self.jointFrequencies = database['jointFrequencies']
        self.prepare_data()
        if MATH.mod(forestSize,2):
            self.forestSize = forestSize + 1
            print 'Forest size has been adjusted to be even.'
        else:
            self.forestSize = forestSize
        self.treeGenDict = treeGenDict
        self.fitnessMethod = fitnessMethod
        self.selectionMethod = selectionMethod
        # elitism can essentially be grafted onto any selection method - ~1/10 of the
        #    population are advanced as copies, without mutation, into the next
        #    round
        self.elitism = elitism
        # ensures eliteN neq 0 for popSize > 10, and is at least 1
        if self.elitism:
            self.eliteN = max(MATH.int(self.forestSize/10.0),1)
        else:
            self.eliteN = 0
        # distances and protein diameter do not change
        self.distances = Utilities.calculateAtomicDistances(pdbFile)
        self.proteinDiameter = mean(self.distances.values()) - min(self.distances.values())
        self.proteinMinimum = min(self.distances.values())
        self.proteinSum = sum(self.distances.values())
        # will be a list of CGAChromosome objects
        self.population = []
        # mutation probabilities - only one kind of mutation is performed per offspring, so
        #    pG + pP + pM + pC + pHC <= 1
        # the leftover probability will just advance the individual without mutation
        self.pG = pG   # growth of terminal nodes (probability is per-terminus)
        self.pP = pP   # pruning any non-terminal node (probability is per-node)
        self.pM = pM   # mutation probability (per node)
        self.pC = pC   # crossover probability (per event - only one allowed per mating)
        self.pHC = pHC # "headless chicken" crossover probability
        self.sampGen = sampGen # write to database once every sampGen generations
        self.time = 0
    
        
    def prepare_data(self):
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
        r = self.treeGenDict['r']
        if self.treeGenDict['treetype'] == 'exponential':
            if self.treeGenDict['p'] < 1.0:
                p = self.treeGenDict['p']
            else:
                # ensures consistent parameter values
                p = 1/self.treeGenDict['p']
            return CGAGenerator.expgenerate(p,r)
        elif self.treeGenDict['treetype'] == 'fixed':
            if self.treeGenDict['p'] > 1:
                treeSize = self.treeGenDict['p']
            else:
                treeSize = MATH.int(1.0/self.treeGenDict['p'])
            return CGAGenerator.generate(treeSize,r)
        else:
            raise TypeError, 'Unknown tree type %s' % self.treeGenDict['treetype']
        
        
    def populate(self):
        """Populate the forest of function trees.  You can choose to either use probabilistic tree
        generation (the number of nodes will be exponentially distributed) or trees with a fixed 
        number of non-terminal nodes.  The treeGenDict determines which method will be used, and
        incompatible parameters result in default behavior.
        """
        for i in range(self.forestSize):
            tree = self.initialize_tree()
            fitness,parsimony,finitewts = self.evaluate_fitness(tree)
            self.population.append(CGAChromosome(tree,fitness,parsimony,finitewts))


    def advance(self):
        """Step forward one step in time."""
        # log data BEFORE manipulating the population
        if MATH.remainder(MATH.int(self.time),self.sampGen) == 0:
            self.notify(time=self.time)
        # first select a round of parents
        parents = [self.select_parent(method=self.selectionMethod) for x in xrange(self.forestSize)]
        offspring = list()
        if not self.elitism:
            for i in xrange(0, self.forestSize):
                offspring.append(self.mate(rchoice(parents),rchoice(parents)))
        else:
            fitvals = MATH.asarray([x.fitness for x in self.population])
            toCopy = fitvals.argsort()[-self.eliteN:]
            offspring += [self.population[x].copy() for x in toCopy]
            # now mate to create the remaining population
            for i in xrange(0, self.forestSize-self.eliteN):
                offspring.append(self.mate(rchoice(parents), rchoice(parents)))
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
        if r < self.pC: # parental crossover
            fitEval = True
            mNode = rchoice(mother.tree.getNodes())
            fNode = rchoice(father.tree.getNodes())
            CGAGenerator.single_crossover(mNode,fNode)
        elif r < self.pC + self.pHC: # headless chicken crossover
            fitEval = True
            mNode = rchoice(mother.tree.getNodes())
            rNode = rchoice(self.initialize_tree().getNodes())
            CGAGenerator.single_crossover(mNode,rNode)
        elif r < self.pC + self.pHC + self.pM: # point mutation (uses pM/node for mutation prob.)
            fitEval = True
            for n in mother.tree.getNodes():
                if urand() < self.pM:
                    CGAGenerator.point_mutate(mother.tree, n)
                    fitEval = True
        elif r < self.pC + self.pHC + self.pM + self.pP: # pruning - guaranteed to do one pruning operation
            fitEval = True
            mNode = rchoice(mother.tree.getNodes())
            CGAGenerator.prune(mother.tree,mNode)
        elif r < self.pC + self.pHC + self.pM + self.pP + self.pG: # growth - guaranteed to do one growth op.
            fitEval = True
            mTerm = rchoice(mother.tree.getTermini())
            CGAGenerator.grow(mother.tree,mTerm)
        else: # offspring will just be a copy
            pass
        if fitEval:
            mother.fitness,mother.parsimony,mother.finitewts = self.evaluate_fitness(mother.tree)
        return mother
         
    
    def select_parent(self,method,**kwargs):
        """A dispatcher that implements a variety of selection methods and returns a single parent
        to place in the pool for the next generation.
            Current allowed methods (all strings):
                -'tournament' : requires a parameter k. two individuals are chosen at random; if 
                    rand < k, the fitter individual is selected.  if rand > k, the less fit one is.
        Regardless of selected method, the parent is returned as a CGAChromosome object."""
        mstring = method + '_selection'
        if hasattr(self, mstring):
            parent = getattr(self,mstring)(**kwargs)
        else:
            parent = rchoice(self.population) # pick at random if there's a problem
        return parent
    
    
    def tournament_selection(self, k=0.75):
        pOne, pTwo = rchoice(self.population), rchoice(self.population)
        # modified to avoid using __cmp__ method of chromosome
        if urand() < k:
            return pOne if pOne > pTwo else pTwo
        else:
            return pOne if pOne < pTwo else pTwo
        
    
    def pareto_tournament_selection(self, k=0.75):
        '''Psuedo-Pareto tournament selection for multi-objective offspring.'''
        pOne, pTwo = rchoice(self.population), rchoice(self.population)
        f1 = MATH.array([pOne.fitness,pOne.parsimony,pOne.finitewts])
        f2 = MATH.array([pTwo.fitness,pTwo.parsimony,pTwo.finitewts])
        if urand() < k:
            return pOne if MATH.sign((f1-f2)).sum() > 0 else pTwo
        else:
            return pOne if MATH.sign((f1-f2)).sum() < 0 else pTwo
    
    
    def evaluate_fitness(self,tree,**kwargs):
        """A dispatcher to do fitness evaluations; the different fitness definitions are stored
        in the fitnessMethod attribute and as different member functions.
        """
        fstring = 'evaluate_fitness_'+self.fitnessMethod
        if hasattr(self, fstring):
            fitness,finitewts = getattr(self,fstring)(tree,**kwargs)
        else:
            fitness,finitewts = -1.0*MATH.inf,0.0
        parsimony = -1.0*len(tree.getNodes())
        return fitness,parsimony,finitewts    
    
    
    def compute_wij(self,tree):
        """Almost all conceivable fitness functions compute a matrix of scores, so that 
        functionality is coded here.
        """
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
    
    
    def evaluate_fitness_distance_matrix(self,tree):
        """Accepts an input tree (member of the forest) and evaluates its fitness, defined
        as a transformed direct fit to the dimensionless distance matrix.  Specifically:
            fitness = -0.5*sum((a*w_ij + b - d_ij)^2)
        This maps fitness to [0,1], and allows an arbitrary linear rescaling of the scores.
        The optimal rescaling is a simple linear algebra subproblem.
        """
        # calculate the weights
        weights = self.compute_wij(tree)
        finitewts = 0.0
        # accumulators allow calculation of best linear transformation
        SDW = 0.0
        W = []
        D = []
        for i,j in weights:
            if (i,j) in self.distances:
                if MATH.isfinite(weights[(i,j)]):
                    W.append(weights[(i,j)])
                    D.append((self.distances[(i,j)] - self.proteinMinimum)/self.proteinDiameter)
                    SDW += ((self.distances[(i,j)] - self.proteinMinimum)/self.proteinDiameter)*weights[(i,j)]
                    finitewts += 1.0
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
        return fitness,finitewts/len(weights)
        
    
    
    def evaluate_fitness_weighted_accuracy(self,tree):
        """Accepts an input tree (member of the forest) and evaluates its fitness, currently 
        defined as:
            fitness = -1.0*sum(ws_ij*d_ij)/sum(ws_ij),
        where d_ij is the dimensionless matrix of (positive) distances, and ws_ij = w_ij + min(min(w_ij),0).
        The -1.0 multiplier insures better outcomes = increasing fitness.
        """
        # compute the weights
        weights = self.compute_wij(tree)
        finitewts = 0.0
        try:
            minw = MATH.min([x for x in weights.values() if MATH.isfinite(x)])
        except ValueError:
            minw = 0.0 # if there are NO finite weights, don't bother
        minw = minw if minw < 0 else 0.0 # no min shift necessary if all weights positive
        accuracy = []
        normalization = minw*len(self.distances)
        for i,j in weights:
            if (i,j) in self.distances:
                if MATH.isfinite(weights[(i,j)]):
                    value = (weights[(i,j)] + minw)*((self.distances[(i,j)] - self.proteinMinimum)/self.proteinDiameter)
                    accuracy.append(value)
                    normalization += weights[(i,j)]
                    finitewts += 1.0
        # normalization might be zero, if there are no finite weights
        if normalization > 0.0:
            fitness = -1.0*sum(accuracy)/normalization
        else:
            fitness = -1.0*MATH.inf
        return fitness,finitewts/len(weights)
                
    
    # potentially deprecated?
    def calculate_accuracy(self,weights):
        """Calculates the accuracy from an input dictionary of weights, keyed on the same indices as the
        distances matrix."""
        accuracy = []
        normalization = 0.0
        for i,j in weights:
            if (i,j) in self.distances:
                value = weights[(i,j)]*((self.distances[(i,j)] - self.proteinMinimum)/self.proteinDiameter)
                accuracy.append(value)
                normalization += weights[(i,j)]
        return 1 - sum(accuracy)/normalization
   
    

class CGASimulationTests(unittest.TestCase):
    def setUp(self):
        self.mySimulation = CGASimulation('../tests/pdz_test.db', '../tests/1iu0.pdb',forestSize=6,fitnessMethod='distance_matrix',selectionMethod='pareto_tournament',treeGenDict={'treetype':'fixed','p':5,'r':0.5})
        self.mySimulation.populate()
        # create and attach a DataLogger
        self.dataLogger = DataLogger()
        self.sqliteLogger = SqliteLogger('../tests')
        self.mySimulation.attach(self.dataLogger)
        self.mySimulation.attach(self.sqliteLogger)
    
    def testKnownTrees(self):
        print "\n\n----- calculating fitness of known trees -----"
        # create trees - MI, OMES, whatever - that have a special form and
        #    check their fitness
        treenames = ['MI','OMES']
        for t in treenames:
            tree = CGAGenerator.generate_special_tree(t)
            print '%s : %s' % (t,tree.getString())
            print 'Weighted accuracy fitness : ',self.mySimulation.evaluate_fitness_weighted_accuracy(tree)
            print 'Distance matrix fitness : ',self.mySimulation.evaluate_fitness_distance_matrix(tree)
           
    #def testCGAChromosome(self):
    #    print "\n\n----- testing comparison operator overloads in CGAChromosome() -----"
    #    c1 = CGAChromosome(tree=None,fitness=1.0)
    #    c2 = CGAChromosome(tree=None,fitness=0.0)
    #    if c1 > c2:
    #        print 'Fitness %f > %f' %(c1.fitness,c2.fitness)
    #    else:
    #        print 'Fitness %f <= %f' %(c1.fitness,c2.fitness)
    #    print 'Maximum fitness : ',MATH.max([c1,c2]).fitness
    #    print 'Minimum fitness : ',MATH.min([c2,c2]).fitness
    """
    def testAdvance(self):
        print "\n\n----- testing one-step advancement -----"
        print 'Before advancement:'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitness) for x in self.mySimulation.population]
        self.mySimulation.advance()
        print 'After advancement (one step):'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitness) for x in self.mySimulation.population]
    """    
    def testMultiAdvance(self):
        print "\n\n----- testing advancement of the tree over many steps -----"
        print 'Before advancement:'
        print 'Pop. size : ', len(self.mySimulation.population)
        n = 10
        print [(x.tree.getString(),x.fitness) for x in self.mySimulation.population]
        for i in range(n):
            print 'Max fitness: %e' % MATH.max([x.fitness for x in self.mySimulation.population])
            self.mySimulation.advance()
        print 'After advancement (%d steps):' % n
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitness) for x in self.mySimulation.population]
    """   
    def testDataLogging(self):
        print "\n\n----- testing data logging -----"
        print 'Advancing twice.'
        self.mySimulation.advance()
        self.mySimulation.advance()
        print 'Logged Data:'
        for k in self.dataLogger.data.keys():
            print k,self.dataLogger.data[k]
            
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
        parent = self.mySimulation.select_parent(method='tournament')
        print 'Parent selected: %s, %f' % (parent.tree.getString(),parent.fitness)
        print 'Elapsed time : %f sec' % (time.clock()-t1)
 

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