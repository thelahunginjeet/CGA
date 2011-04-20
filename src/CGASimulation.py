#!/usr/bin/env python

import unittest, os, time
from scipy import mean,log
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
    and their fitness. Right now this is very basic, but helps with bookkeeping in CGASimulation.  __cmp__
    is overloaded to allow fitness comparisons and population sorting based on fitness.  More functionality 
    could be added as desired."""
    def __init__(self, tree=None, fitness=MATH.nan):
        self.tree = tree
        self.fitness = fitness
    
    def copy(self):
        return CGAChromosome(self.tree.copy(), self.fitness)

    def __cmp__(self, other):
        assert self.fitness is not None
        # add in some handling for NaN, -Inf, and Inf to automatically put them at the bottom
        return cmp(self.fitness, other.fitness)


class CGASimulation(Subject):
    """Main class that does the simulation, records results, evaluates trees, etc."""
    def __init__(self, databaseFile, pdbFile, treeGenDict = {'treetype':'fixed','p':5,'r':0.6}, selectionMethod='tournament',forestSize=100, timeSteps=100, sampGen=10, pG = 0.01, pP = 0.01, pM = 0.01, pC = 0.05):
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
        self.selectionMethod = selectionMethod
        # distances and protein diameter do not change
        self.distances = Utilities.calculateAtomicDistances(pdbFile)
        self.proteinDiameter = mean(self.distances.values()) - min(self.distances.values())
        self.proteinMinimum = min(self.distances.values())
        # will be a list of CGAChromosome objects
        self.population = []
        # mutation probabilities - some are per node, some are global (per tree)
        self.pG = pG # growth of terminal nodes
        self.pP = pP # pruning any non-terminal node
        self.pM = pM # mutation probability
        self.pC = pC # crossover probability
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
    
        
    def populate(self):
        """Populate the forest of function trees.  You can choose to either use probabilistic tree
        generation (the number of nodes will be exponentially distributed) or trees with a fixed 
        number of non-terminal nodes.  The treeGenDict determines which method will be used, and
        incompatible parameters result in default behavior.
        """
        r = self.treeGenDict['r']
        if self.treeGenDict['treetype'] == 'exponential':
            if self.treeGenDict['p'] < 1.0:
                p = self.treeGenDict['p']
            else:
                # use a default value; parameters are inconsistent
                p = 0.65
            for i in xrange(self.forestSize):
                tree = CGAGenerator.expgenerate(p,r)
                fitness = self.evaluate_fitness(tree)
                self.population.append(CGAChromosome(tree, fitness))
        elif self.treeGenDict['treetype'] == 'fixed':
            if self.treeGenDict['p'] > 1:
                treeSize = self.treeGenDict['p']
            else:
                # default value; inconsistent parameters
                treeSize = 5
            for i in xrange(self.forestSize):
                tree = CGAGenerator.generate(treeSize,r)
                fitness = self.evaluate_fitness(tree)
                self.population.append(CGAChromosome(tree, fitness))
        else:
            raise TypeError, 'Unknown tree type : %s' % treetype


    def advance(self):
        """Step forward one step in time recording."""
        # log data BEFORE manipulating the population
        if MATH.remainder(MATH.int(self.time),self.sampGen) == 0:
            self.notify(time=self.time)
        # first select a round of parents
        parents = [self.select_parent(method=self.selectionMethod) for x in xrange(self.forestSize)]
        # now obtain the offspring, two at a time
        offspring = []
        for i in xrange(0, self.forestSize, 2):
            offspring += self.mate(rchoice(parents), rchoice(parents))
        # overwrite current forest
        self.population = offspring
        # advance time
        self.time += 1
        
    def mate(self, parentOne, parentTwo):
        """Accepts two parents and returns two offspring; if the offspring is unchanged from one of
        the parents, the fitness is not re-evaluated, just copied."""
        # offspring begin as copies of their parents
        offOne, offTwo = parentOne.copy(), parentTwo.copy()
        # fitEval[i] will be set to True if any mutations occur that make offspring i different from
        #    parent i or j; this way we avoid unnecessary fitness evaluations
        fitEval = [False, False]
        # NONCONSERVATIVE (GROWTH, PRUNING) FIRST
        # GROWTH/EXTENSION
        for t in offOne.tree.getTermini():
            if urand() < self.pG:
                CGAGenerator.grow(offOne.tree, t)
                fitEval[0] = True
        for t in offTwo.tree.getTermini():
            if urand() < self.pG:
                CGAGenerator.grow(offTwo.tree, t)
                fitEval[1] = True
        # PRUNING - SHOULD HAVE KEPT LOCAL COPIES, OR NOT?
        for n in offOne.tree.getNodes():
            if urand() < self.pP:
                CGAGenerator.prune(offOne.tree, n)
                fitEval[0] = True
        for n in offTwo.tree.getNodes():
            if urand() < self.pP:
                CGAGenerator.prune(offOne.tree, n)
                fitEval[1] = True
        # CROSSOVER - ONLY ONE PER MATING EVENT IS ALLOWED
        if urand() < self.pC:
            fitEval = [True,True]
            # pick the nodes (roots won't crossover)
            nodeOne = rchoice(offOne.tree.getNodes())
            nodeTwo = rchoice(offTwo.tree.getNodes())
            CGAGenerator.single_crossover(nodeOne, nodeTwo)
        # POINT MUTATION 
        for n in offOne.tree.getNodes():
            if urand() < self.pM:
                CGAGenerator.point_mutate(offOne.tree, n)
                fitEval[0] = True
        for n in offTwo.tree.getNodes():
            if urand() < self.pM:
                CGAGenerator.point_mutate(offTwo.tree, n)
                fitEval[1] = True
        # compute fitnesses, if they have changed, and update the trees
        if fitEval[0]:
            offOne.fitness = self.evaluate_fitness(offOne.tree)
        if fitEval[1]:
            offTwo.fitness = self.evaluate_fitness(offTwo.tree)
        return offOne, offTwo
         
    
    def select_parent(self,method,**kwargs):
        """A dispatcher that implements a variety of selection methods and returns a single parent
        to place in the pool for the next generation.
            Current allowed methods (all strings):
                -'tournament' : requires a parameter k. two individuals are chosen at random; if 
                    rand < k, the fitter individual is selected.  if rand > k, the less fit one is.
        Regardless of selected method, the parent is returned as a CGAChromosome object."""
        mstring = method + '_selection'
        if hasattr(self, method):
            parent = getattr(self,method)(**kwargs)
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
        
    
    def evaluate_fitness(self,tree):
        """Accepts an input tree (member of the forest) and evaluates its fitness, currently 
        defined as:
            fitness = 1 + accuracy
        Accuracy is computed in a different function, in order to allow easy swapping in of
        different definitions."""
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
        self.mySimulation = CGASimulation('../tests/pdz_test.db', '../tests/1iu0.pdb',forestSize=6,selectionMethod='tournament',treeGenDict={'treetype':'fixed','p':5,'r':0.5})
        self.mySimulation.populate()
        # create and attach a DataLogger
        self.dataLogger = DataLogger()
        self.sqliteLogger = SqliteLogger('../tests')
        self.mySimulation.attach(self.dataLogger)
        self.mySimulation.attach(self.sqliteLogger)
        
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
   
    def testAdvance(self):
        print "\n\n----- testing one-step advancement -----"
        print 'Before advancement:'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitness) for x in self.mySimulation.population]
        self.mySimulation.advance()
        print 'After advancement (one step):'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitness) for x in self.mySimulation.population]
        
    def testMultiAdvance(self):
        print "\n\n----- testing advancement of the tree over many steps -----"
        print 'Before advancement:'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitness) for x in self.mySimulation.population]
        for i in range(100):
            self.mySimulation.advance()
        print 'After advancement (100 steps):'
        print 'Pop. size : ', len(self.mySimulation.population)
        print [(x.tree.getString(),x.fitness) for x in self.mySimulation.population]
   
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