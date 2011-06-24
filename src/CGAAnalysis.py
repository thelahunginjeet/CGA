#!/usr/bin/env/python

import unittest,os
import sqlite3
import numpy as np

class CGAAnalysis(object):
    """Class to perform post-hoc analysis on database files created by CGASimulation.  Reads simulation parameters from the
    cgapar table, and performs calculations/extracts information in the cgarun table."""
    def __init__(self,path,dbFile):
        self.path = path
        self.dbFile = dbFile
        self.connection = sqlite3.connect(os.path.join(self.path,self.dbFile))
        self.cursor = self.connection.cursor()
        # fetch run parameters and fields; load them into a dictionary
        self.pardict = {}
        fields = [x[1] for x in self.cursor.execute("""PRAGMA table_info(cgapar)""").fetchall()]
        # below might fail on an empty database; insert error checking
        values = [x for x in self.cursor.execute("""SELECT * FROM cgapar""").fetchall()[0]]
        for i in range(0,len(fields)):
            self.pardict[fields[i]] = values[i]
        # these hold data
        self.generations = None
        self.fitness = None
        self.parsimony = None
        self.finitewts = None
        self.maxfitness = None
        self.meanfitness = None
        self.offline = None
        self.online = None
        
        
    def __repr__(self):
        """String representation of the run, with values for parameters and flags."""
        repstr = 'Run database : ' + os.path.join(self.path,self.dbFile) + '\n'
        repstr += 'Parameters :\n'
        for k in sorted(self.pardict):
            repstr += '\t%s : %s\n' % (k,self.pardict[k])
        return repstr
    
    
    def get_generations(self):
        """Determines the generations at which data was taken."""
        if self.generations is None:
            self.generations = np.unique([x[0] for x in self.cursor.execute("""SELECT generation FROM cgarun""").fetchall()])
        return self.generations
    
    def get_fitness(self):
        """Fetches and returns all fitness values recorded during the run."""
        if self.fitness is None:
            self.fitness = [x[0] for x in self.cursor.execute("""SELECT fitness FROM cgarun""").fetchall()]
        return self.fitness
    
    def get_parsimony(self):
        """Fetches and returns all parsimony values recorded during the run."""
        if self.parsimony is None:
            self.parsimony = [x[0] for x in self.cursor.execute("""SELECT parsimony FROM cgarun""").fetchall()]
        return self.parsimony
    
    def get_finite_weights(self):
        """Fetches and returns all fraction-of-weights-which-are-finite terms."""
        if self.finitewts is None:
            self.finitewts = [x[0] for x in self.cursor.execute("""SELECT finitewts FROM cgarun""").fetchall()]
        return self.finitewts
    
    def get_mean_parsimony(self):
        if self.generations is None:
            _g = self.get_generations()
        meanpars = []
        _pars = self.get_parsimony()
        if not self.pardict.has_key('forestSize'):
            print 'ERROR.  Cannot determine size of forest.'
            return None
        else:
            npop = self.pardict['forestSize']
        for i in range(0,len(self.generations)):
            chunk = [x for x in _pars[npop*i:npop*i+npop]]
            meanpars.append(np.mean(chunk))
        return np.asarray(meanpars)
    
    
    def get_mean_finite_weights(self):
        if self.generations is None:
            _g = self.get_generations()
        meanfw = []
        _fw = [x[0] for x in self.cursor.execute("""SELECT finitewts FROM cgarun""").fetchall()]
        if not self.pardict.has_key('forestSize'):
            print 'ERROR.  Cannot determine size of forest.'
            return None
        else:
            npop = self.pardict['forestSize']
        for i in range(0,len(self.generations)):
            chunk = [x for x in _fw[npop*i:npop*i+npop]]
            meanfw.append(np.mean(chunk))
        return np.asarray(meanfw)
            
    
    
    def get_mean_fitness(self):
        """Calculates and returns the average fitness (over the forest), per sampling step, ignoring trees with negative fitness."""
        # need generations to determine sampling rate
        if self.generations is None:
            _g = self.get_generations()
        if self.meanfitness is None:
            self.meanfitness = []
            _fit = self.get_fitness()
            if not self.pardict.has_key('forestSize'):
                print 'ERROR.  Cannot determine size of forest.'
                return None
            else:
                npop = self.pardict['forestSize']
            for i in range(0,len(self.generations)):
                chunk = [x for x in _fit[npop*i:npop*i+npop] if x >= 0] # drops negatives
                self.meanfitness.append(np.mean(chunk))
            self.meanfitness = np.asarray(self.meanfitness)
        return self.meanfitness
    
    
    def get_max_fitness(self):
        """Calculates and returns the maximum fitness (over the forest), per sampling step."""
        if self.generations is None:
            _g = self.get_generations()
        if self.maxfitness is None:
            self.maxfitness = []
            _fit = self.get_fitness()
            if not self.pardict.has_key('forestSize'):
                print 'ERROR.  Cannot determine size of forest.'
                return None
            else:
                npop = self.pardict['forestSize']
            for i in range(0,len(self.generations)):
                chunk = [x for x in _fit[npop*i:npop*i+npop]]
                self.maxfitness.append(np.max(chunk))
            self.maxfitness = np.asarray(self.maxfitness)
        return self.maxfitness
    
    
    def get_offline_fitness(self):
        """Offline performance at time t is defined as the average value, over t, of the best fitness
        that has been seen up to t."""
        if self.offline is None:
            self.offline = []
            mf = self.get_max_fitness()
            for i in range(1,len(mf)):
                self.offline.append(np.max(mf[0:i])) 
            self.offline = (np.asarray(self.offline).cumsum())/xrange(1,len(self.offline)+1)
        return self.offline
        
    
    def get_online_fitness(self):
        """Online performance at time t is defined by the average fitness of all individuals that have
        been evaluated up to t steps.  This implementation assumes all members of the forest at each 
        generation are used in the average fitness, which is not strictly true."""
        if self.online is None:
            mf = self.get_mean_fitness()
            self.online = mf.cumsum()/xrange(1,len(mf)+1)
        return self.online
    

class CGAAnalysisTests(unittest.TestCase):
    def setUp(self):
        # database file to test - need something in ../tests directory (just rename one 
        #    created by a run)
        self.analysis = CGAAnalysis(path='../tests',dbFile='testdb.sqldb')
    
    def testRepr(self):
        print "\n\n----- testing string representation of class -----"
        print self.analysis
    
    def testGetGenerations(self):
        print "\n\n----- testing fetch and return of generations -----"
        g = self.analysis.get_generations()
        print 'Sampling rate: ', g[1]-g[0]
        print 'Time range: %d . . . %d' %(g[0],g[-1])
        
    def testGetMeanFitness(self):
        print "\n\n----- testing calculation of mean fitness -----"
        meanf = self.analysis.get_mean_fitness()
        
    def testGetMaxFitness(self):
        print "\n\n----- testing calculation of max fitness -----"
        maxf = self.analysis.get_max_fitness()
        
    def testGetOnlineFitness(self):
        print "\n\n----- testing calculation of online fitness -----"
        online = self.analysis.get_online_fitness()
        
    def testGetOfflineFitness(self):
        print "\n\n----- testing calculation of offline fitness -----"
        offline = self.analysis.get_offline_fitness()
        
    def testPlotFitness(self):
        print "\n\n----- plotting calculated fitness values -----"
        import pylab
        g = self.analysis.get_generations()
        meanf = self.analysis.get_mean_fitness()
        maxf = self.analysis.get_max_fitness()
        online = self.analysis.get_online_fitness()
        offline = self.analysis.get_offline_fitness()
        pylab.subplot(221)
        pylab.plot(g,meanf,'k-')
        pylab.ylabel('Mean Fitness')
        pylab.subplot(222)
        pylab.plot(g,maxf,'k-')
        pylab.ylabel('Max Fitness')
        pylab.subplot(223)
        pylab.plot(g,online,'k-')
        pylab.ylabel('Online Fitness')
        pylab.subplot(224)
        pylab.plot(g[1:],offline,'k-')
        pylab.ylabel('Offline Fitness')
        pylab.show()
        

if __name__ == 'main':
    unittest.main()
        