#!/usr/bin/env python

import unittest, os, time
from scipy import mean,log
import numpy as MATH
import random
from CGAPreprocessing import Utilities
import CGAGenerator

# TODO : 
#  - write the simulation as a generator
#  - undo hard-coding of parameters passed into selection functions
#  - add mutations (right now offspring are just reshuffled, with duplicates)


class CGASimulation(object):
    """Main class that does the simulation, records results, evaluates trees, etc."""
    def __init__(self, databaseFile, pdbFile, forestSize=100, timeSteps=100, pG = 0.01, pP = 0.01, pM = 0.01, pC = 0.01):
        if not os.path.exists(pdbFile):
            raise IOError, 'something is wrong with your pdb file; check yourself . . .'
        database = Utilities.readDatabase(databaseFile)
        self.indices = database['indices']
        self.singleFrequencies = database['singleFrequencies']
        self.jointFrequencies = database['jointFrequencies']
        self.prepare_data()
        # distances and protein diameter do not change
        self.distances = Utilities.calculateAtomicDistances(pdbFile)
        self.proteinDiameter = mean(self.distances.values()) - min(self.distances.values())
        self.proteinMinimum = min(self.distances.values())
        self.cgag = CGAGenerator.CGAGenerator()
        # parameters specific to the simulation; fitness evaluations are expensive, so we 
        #    keep track of them for any children which do not differ from their parents
        self.population = []
        self.fitnesses = []
        # mutation probabilities - some are per node, some are global (per tree)
        #    pG : probability of growth
        #    pP : probability of pruning
        #    pM : probability of mutation
        #    pC : probability of crossover (between two trees)
        self.pG = pG
        self.time = 0
        self.forestSize = forestSize
    
        
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
                self.population.append(self.cgag.expGenerate(p))
        elif treetype == 'fixed':
            for i in range(self.forestSize):
                self.population.append(self.cgag.generate(treeSize))
        else:
            raise TypeError, 'Unknown tree type : %s' % treetype
        # evaluate the fitness of the initial forest
        self.fitnesses = [self.evaluate_fitness(t) for t in self.population]
    
    def advance(self):
        """Step forward one step in time recording"""
        # first select a round of parents
        parents = []
        parentfit = []
        for i in range(0,self.forestSize):
            pt,pfit = self.select_parent(method='tournament')
            parents.append(pt)
            parentfit.append(pfit)
        # now the offspring, two at a time
        offspring = []
        ofit = []
        for i in range(0,self.forestSize,2):
            indices = MATH.random.randint(self.forestSize,size=2)
            offspring += [parents[indices[0]],parents[indices[1]]]
            ofit += [parentfit[indices[0]],parentfit[indices[1]]]
        # overwrite forest/fitness
        self.population = offspring
        self.fitnesses = ofit
        
    
    def select_parent(self,method,**kwargs):
        """A dispatcher that implements a variety of selection methods and returns a single parent
        to place in the pool for the next generation.
            Current allowed methods (all strings):
                -'tournament' : requires a parameter k. two individuals are chosen at random; if 
                    rand < k, the fitter individual is selected.  if rand > k, the less fit one is.
        Regardless of selected method, the parent and its fitness are returned (as a tuple)"""
        mstring = method+'_selection'
        if hasattr(self,method):
            ptree,pfit = getattr(self,method)(**kwargs)
        else:
            # to avoid the simulation grinding to a total halt, just pick parents at random
            #    if there's a problem
            pindx = MATH.random.randint(self.forestSize)
            pt = self.population[pindx]
            pfit = self.fitnesses[pindx]
        return pt,pfit
    
    
    def tournament_selection(self,k=0.75):
        indices = MATH.random.randint(self.forestSize,size=2)
        if MATH.random.rand() < k:
            pindx = MATH.argmax(indices)
        else:
            pindx = MATH.argmin(indices)
        return self.population[pindx],self.fitnesses[pindx]
        
    
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
        self.mySimulation = CGASimulation('../tests/pdz_test.db', '../tests/1iu0.pdb',forestSize=5)
        self.mySimulation.populate(treetype='fixed',treeSize=5)
        
    def testAdvance(self):
        print "\n----- testing one-step advancement -----"
        print 'Before advancement:'
        print [x.string for x in self.mySimulation.population]
        print self.mySimulation.fitnesses
        self.mySimulation.advance()
        print 'After advancement:'
        print [x.string for x in self.mySimulation.population]
        print self.mySimulation.fitnesses

    def testPopulation(self):
        print "\n----- testing population generation -----"
        for x in self.mySimulation.population:
            x()
            print 'String rep : %s' % x.string

    def testShape(self):
        print "\n----- testing shape conversion of input data -----"
        print self.mySimulation.indices
        print 'P_7(:,1)  : ', self.mySimulation.singleFrequencies[7][:,0]
        print 'sum(P_7(:,1)) : ', sum(self.mySimulation.singleFrequencies[7][:,0])
        print 'Dimensions of P_i :          %d X %d' % (self.mySimulation.singleFrequencies[7].shape)
        print 'Dimensions of P_ij :         %d X %d' % (self.mySimulation.jointFrequencies[(7,20)].shape)
        print 'Shape of entropy(P_7) :      %d X %d' % ((self.mySimulation.singleFrequencies[7]*log(self.mySimulation.singleFrequencies[7])).shape)
        print 'Shape of entropy(P_{7,20}) : %d X %d' % ((self.mySimulation.jointFrequencies[(7,20)]*log(self.mySimulation.jointFrequencies[(7,20)])).shape)
        print 'Distance(7,20) : ' , self.mySimulation.distances[(7, 20)]
        
    def testSelectParent(self):
        print"\n----- testing parental selection -----"
        t1 = time.clock()
        ptree,pfit = self.mySimulation.select_parent(method='tournament')
        print 'Parent selected: %s, %f' % (ptree.string,pfit)
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