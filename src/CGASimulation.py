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
        database = Utilities.readDatabase(databaseFile)
        self.indices = database['indices']
        self.singleFrequencies = database['singleFrequencies']
        self.jointFrequencies = database['jointFrequencies']
        if not os.path.exists(pdbFile):
            raise IOError, 'something is wrong with your pdb file; check yourself . . .'
        self.distances = Utilities.calculateAtomicDistances(pdbFile)
        self.cgag = CGAGenerator.CGAGenerator()
        # parameters specific to the simulation
        self.population = []
        self.time = 0
        self.forestSize = forestSize
        
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
    
    def evaluate(self):
        """Evaluates the fitness of each species in the population."""
        pass

    def calculateAccuracy(self):
        """Calculates the accuracy from an input"""
        pass
#        accuracy = []
#        weights = 0.0
#        proteinDiameter = mean(self.distances.values()) - min(self.distances.values())
#        proteinMinimum = min(self.distances.values())
#        for i,j in something:
#            accuracy.append(weightedValue*((self.distances[(i,j)]-proteinMinimum)/proteinDiameter))
#            weights += weightedValue
#        finalAccuracy = 1 - sum(accuracy)/weights        
    

class CGASimulationTests(unittest.TestCase):
    def setUp(self):
        self.mySimulation = CGASimulation('../tests/pdz_test.db', '../tests/1iu0.pdb',forestSize=10)
        self.mySimulation.populate(treetype='fixed',treeSize=5)

    def testPopulation(self):
        print "\n\n----- testing population generation -----"
        for x in self.mySimulation.population:
            x()
            print 'String rep : %s' %(x.string)

    def testSimulation(self):
        print "\n\n----- testing simulation -----"
        print self.mySimulation.indices
        print 'Type of P_i  : ', type(MATH.asarray(self.mySimulation.singleFrequencies[20]))
        print 'Type of P_ij : ', type(MATH.asarray(self.mySimulation.jointFrequencies[(7,20)]))
        print (MATH.asarray(self.mySimulation.singleFrequencies[20]))*log(MATH.asarray(self.mySimulation.singleFrequencies[20]))
        #print self.mySimulation.jointFrequencies[(7,20)]*log(self.mySimulation.jointFrequencies[(7,20)])
        #print log(self.mySimulation.jointFrequencies[(7,20)])
        print self.mySimulation.distances[(7, 20)]
    
        
if __name__ == '__main__':
    unittest.main()