#!/usr/bin/env python

import unittest, os
from scipy import mean,log
import numpy as MATH
from CGAPreprocessing import Utilities
import CGAGenerator

# TODO : 
#  - write the simulation as a generator


class CGASimulation(object):
    """Main class that does the simulation, records results, evaluates trees, etc."""
    def __init__(self, databaseFile, pdbFile, forestSize=100, timeSteps=100):
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
        # parameters specific to the simulation
        self.population = []
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
    
    def advance(self):
        """Step forward one step in time recording"""
        pass
    
    def evaluate_fitness(self, tree):
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
        print "\tfitness : ", fitness
        
                
    def calculate_accuracy(self, weights):
        """Calculates the accuracy from an input dictionary of weights, keyed on the same indices as the
        distances matrix.  The weights are trees evaluated at a specific (i, j) corresponding to two indices
        from the protein."""
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
        self.mySimulation = CGASimulation('../tests/pdz_test.db', '../tests/1iu0.pdb',forestSize=10)
        self.mySimulation.populate(treetype='fixed',treeSize=5)

    def testPopulation(self):
        print "\n----- testing population generation -----"
        for x in self.mySimulation.population:
            x()
            print 'String rep : %s, Value : %f' %(x.string,x.function)

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
        
    def testEvaluateFitness(self):
        print "\n----- testing fitness evaluation -----"
        for tree in self.mySimulation.population:
            tree()
            print tree.string
            self.mySimulation.evaluate_fitness(tree)
        
if __name__ == '__main__':
    unittest.main()
