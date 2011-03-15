#!/usr/bin/env python

import unittest, os, time, copy
from scipy import mean,log
import numpy as MATH
import random
from CGAPreprocessing import Utilities
from CGAStructures import AlgorithmTree
import CGAGenerator
from CGALogging import Subject,Observer,DataLogger

# TODO : 
#  - write the simulation as a generator
#  - undo hard-coding of parameters passed into selection functions
#  - dangerous for CGAChromosome to do all comparisons based on fitness?

class CGAChromosome(object):
    """CGA chromosomes are a container class that stores both the unit of selection (a tree in this problem) 
    and their fitness. Right now this is very basic, but helps with bookkeeping in CGASimulation.  __cmp__
    is overloaded to allow fitness comparisons and population sorting based on fitness.  More functionality 
    could be added as desired."""
    def __init__(self,tree=None,fitness=MATH.nan):
        self.tree = tree
        self.fitness = fitness
    
    def __cmp__(self,other):
        if hasattr(other,'fitness'):
            return cmp(self.fitness,other.fitness)


class CGASimulation(Subject):
    """Main class that does the simulation, records results, evaluates trees, etc."""
    def __init__(self, databaseFile, pdbFile, forestSize=100, timeSteps=100, pG = 0.01, pP = 0.01, pM = 0.01, pC = 0.05):
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
        # distances and protein diameter do not change
        self.distances = Utilities.calculateAtomicDistances(pdbFile)
        self.proteinDiameter = mean(self.distances.values()) - min(self.distances.values())
        self.proteinMinimum = min(self.distances.values())
        self.cgag = CGAGenerator.CGAGenerator()
        # will be a list of CGAChromosome objects
        self.population = []
        # mutation probabilities - some are per node, some are global (per tree)
        #    pG : probability of growth
        #    pP : probability of pruning
        #    pM : probability of mutation
        #    pC : probability of crossover (between two trees)
        self.pG = pG
        self.pP = pP
        self.pM = pM
        self.pC = pC
        self.time = 0
    
        
    def prepare_data(self):
        # change input frequencies from the db to 20X20 arrays (including single frequencies)
        for i in range(0,len(self.indices)):
            indxi = self.indices[i]
            self.singleFrequencies[indxi] = MATH.tile(MATH.asarray(self.singleFrequencies[indxi]),(1,20))
            for j in range(i+1,len(self.indices)):
                indxj = self.indices[j]
                self.jointFrequencies[(indxi,indxj)] = MATH.asarray(self.jointFrequencies[(indxi,indxj)])
        
    def populate(self,treetype='exponential',p=0.65,treeSize=10):
        """Populate the forest of function trees.  You can choose to either use probabilistic tree
        generation (the number of nodes will be exponentially distributed) or trees with a fixed 
        number of non-terminal nodes.
            Input:
                treetype : str, optional, can be 'exponential' or 'fixed'
                p        : float, optional, controls number of nodes in exponential trees
                treeSize : int, optional, number of nodes for fixed trees
                popSize  : int, optional, forest size
        """
        if treetype == 'exponential':
            for i in range(self.forestSize):
                tree = self.cgag.expgenerate(p)
                fitness = self.evaluate_fitness(tree)
                self.population.append(CGAChromosome(tree,fitness))
        elif treetype == 'fixed':
            for i in range(self.forestSize):
                tree = self.cgag.generate(treeSize)
                fitness = self.evaluate_fitness(tree)
                self.population.append(CGAChromosome(tree,fitness))
        else:
            raise TypeError, 'Unknown tree type : %s' % treetype

    def advance(self):
        """Step forward one step in time recording.  Before advancing, """
        # first select a round of parents
        parents = [self.select_parent(method='tournament') for x in range(0,self.forestSize)]
        # now obtain the offspring, two at a time
        offspring = []
        for i in range(0,self.forestSize,2):
            offspring += self.mate(random.choice(parents),random.choice(parents))
        # overwrite current forest
        self.population = offspring
        # a few things we want to save
        minN = MATH.min([len(k.tree.nodes) for k in self.population])
        maxN = MATH.max([len(k.tree.nodes) for k in self.population])
        maxFit = MATH.nanmax([k.fitness for k in self.population])
        wellFormed = len([k.fitness for k in self.population if ~MATH.isnan(k.fitness) and ~MATH.isinf(k.fitness)])/MATH.float64(self.forestSize)
        meanFit = MATH.mean([k.fitness for k in self.population if ~MATH.isnan(k.fitness) and ~MATH.isinf(k.fitness)])
        self.notify(time=self.time,minSize=minN,maxSize=maxN,maxFit=maxFit,wellFormed=wellFormed,meanFit=meanFit)
        self.time += 1
        
    def mate(self,parentOne,parentTwo):
        """Accepts two parents and returns two offspring; if the offspring is unchanged from one of
        the parents, the fitness is not re-evaluated, just copied."""
        # offspring begin as copies of their parents
        offOne = copy.copy(parentOne)
        offTwo = copy.copy(parentTwo)
        # fitEval[i] will be set to True if any mutations occur that make offspring i different from
        #    parent i or j; this way we avoid unnecessary fitness evaluations
        fitEval = [False,False]
        # CROSSOVER
        # sexual - check for a crossover first.  Allow only one per mating event?
        if MATH.random.rand() < self.pC:
            fitEval = [True,True]
            # pick the nodes (roots won't crossover)
            nodeOne = random.choice(offOne.tree.nodes)
            nodeTwo = random.choice(offTwo.tree.nodes)
            self.cgag.single_crossover(offOne.tree,nodeOne,offTwo.tree,nodeTwo)
        # POINT MUTATION
        # now check for point mutations, in both trees
        for n in offOne.tree.nodes:
            if MATH.random.rand() < self.pM:
                self.cgag.point_mutate(offOne.tree,n)
                fitEval[0] = True
        for n in offTwo.tree.nodes:
            if MATH.random.rand() < self.pM:
                self.cgag.point_mutate(offTwo.tree,n)
                fitEval[1] = True
        # GROWTH/EXTENSION
        for t in offOne.tree.termini:
            if MATH.random.rand() < self.pG:
                self.cgag.grow(offOne.tree,t)
                fitEval[0] = True
        for t in offTwo.tree.termini:
            if MATH.random.rand() < self.pG:
                self.cgag.grow(offTwo.tree,t)
                fitEval[1] = True
        # PRUNING
        for n in offOne.tree.nodes:
            if MATH.random.rand() < self.pP:
                self.cgag.prune(offOne.tree,n)
                fitEval[0] = True
            if MATH.random.rand() < self.pP:
                self.cgag.prune(offOne.tree,n)
                fitEval[1] = True
        # compute fitnesses, if they have changed
        if fitEval[0]:
            offOne.fitness = self.evaluate_fitness(offOne.tree)
        if fitEval[1]:
            offTwo.fitness = self.evaluate_fitness(offTwo.tree)
        return offOne,offTwo
         
    
    def select_parent(self,method,**kwargs):
        """A dispatcher that implements a variety of selection methods and returns a single parent
        to place in the pool for the next generation.
            Current allowed methods (all strings):
                -'tournament' : requires a parameter k. two individuals are chosen at random; if 
                    rand < k, the fitter individual is selected.  if rand > k, the less fit one is.
        Regardless of selected method, the parent is returned as a CGAChromosome object."""
        mstring = method+'_selection'
        if hasattr(self,method):
            parent = getattr(self,method)(**kwargs)
        else:
            # to avoid the simulation grinding to a total halt, just pick parents at random
            #    if there's a problem
            parent = random.choice(self.population)
        return parent
    
    
    def tournament_selection(self,k=0.75):
        pOne,pTwo = (random.choice(self.population),random.choice(self.population))
        if MATH.random.rand() < k:
            return pOne if pOne > pTwo else pTwo
        else:
            return pOne if pOne < pTwo else pTwo
        
    
    def evaluate_fitness(self,tree):
        """Accepts an input tree (member of the forest) and evaluates its fitness, currently 
        defined as:
            fitness = 1 + accuracy
        Accuracy is computed in a different function, in order to allow easy swapping in of
        different definitions."""
        weights = {}
        pi = [x for x in tree.termini if x.string == 'p_i']
        pj = [x for x in tree.termini if x.string == 'p_j']
        pij = [x for x in tree.termini if x.string == 'p_ij']
        for i in self.indices:
            for j in [x for x in self.indices if x > i]:
                map(lambda x : x.replaceData(self.singleFrequencies[i]), pi)
                map(lambda x : x.replaceData(self.singleFrequencies[j]), pj)
                map(lambda x : x.replaceData(self.jointFrequencies[(i, j)]), pij)
                tree()
                weights[(i, j)] = tree.function
        fitness = 1 + self.calculate_accuracy(weights)
        return fitness
    

    def evaluate_fitness_iter(self,tree):
        """Accepts an input tree (member of the forest) and evaluates its fitness, currently 
        defined as:
            fitness = 1 + accuracy
        Accuracy is computed in a different function, in order to allow easy swapping in of
        different definitions."""
        weights = {}
        pi = [x for x in tree.termini if x.string == 'p_i']
        pj = [x for x in tree.termini if x.string == 'p_j']
        pij = [x for x in tree.termini if x.string == 'p_ij']
        for i in self.indices:
            for j in [x for x in self.indices if x > i]:
                map(lambda x : x.replaceData(self.singleFrequencies[i]), pi)
                map(lambda x : x.replaceData(self.singleFrequencies[j]), pj)
                map(lambda x : x.replaceData(self.jointFrequencies[(i, j)]), pij)
                tree()
                weights[(i, j)] = tree.function
        fitness = 1 + self.calculate_accuracy(weights)
        return fitness
                

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
        self.mySimulation = CGASimulation('../tests/pdz_test.db', '../tests/1iu0.pdb',forestSize=6)
        self.mySimulation.populate(treetype='fixed',treeSize=5)
        # create and attach a DataLogger
        self.dataLogger = DataLogger()
        self.mySimulation.attach(self.dataLogger)
        
    def testCGAChromosome(self):
        print "\n\n----- testing comparison operator overloads in CGAChromosome() -----"
        c1 = CGAChromosome(tree=None,fitness=1.0)
        c2 = CGAChromosome(tree=None,fitness=0.0)
        if c1 > c2:
            print 'Fitness %f > %f' %(c1.fitness,c2.fitness)
        else:
            print 'Fitness %f <= %f' %(c1.fitness,c2.fitness)
        print 'Maximum fitness : ',MATH.max([c1,c2]).fitness
        print 'Minimum fitness : ',MATH.min([c2,c2]).fitness
        
    def testAdvance(self):
        print "\n\n----- testing one-step advancement -----"
        print 'Before advancement:'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.string,x.fitness) for x in self.mySimulation.population]
        self.mySimulation.advance()
        print 'After advancement:'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.string,x.fitness) for x in self.mySimulation.population]
        
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
        for x in self.mySimulation.population:
            x.tree()
            print 'String rep : %s' % x.tree.string

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
        print 'Parent selected: %s, %f' % (parent.tree.string,parent.fitness)
        print 'Elapsed time : %f sec' % (time.clock()-t1)
 
        
    # bring this back in later when we actually are trying to optimize
    """def testEvaluateFitness(self):
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
            print 'Elapsed clock time (optimized) : %f seconds' %(t2-t1)"""
    
        
if __name__ == '__main__':
    unittest.main()