#!/usr/bin/env/python

import unittest,os
import sqlite3
import numpy as np

"""Exceptions for the CGAAnalysis class."""
class CGADBFieldException(Exception):
    def __init__(self,dbField):
        print "Field %s does not exist in run database. Check table_info(cgarun)." % dbField

class CGADBException(Exception):
    def __init__(self,dbFileName):
        print "Database %s does not exist." % dbFileName

class CGAAnalysis(object):
    """Class to perform post-hoc analysis on database files created by CGASimulation.  Parameter
    values are read from the run database."""
    def __init__(self,path,dbFileName):
        self.path = path
        self.dbFileName = dbFileName
        if not os.path.exists(os.path.join(self.path,self.dbFileName)):
            raise CGADBException(os.path.join(self.path,self.dbFileName))
        self.connection = sqlite3.connect(os.path.join(self.path,self.dbFileName))
        self.cursor = self.connection.cursor()
        self.pardict = {}
        fields = [x[1] for x in self.cursor.execute("""PRAGMA table_info(cgapar)""").fetchall()]
        # below might fail on an empty database; insert error checking
        values = [x for x in self.cursor.execute("""SELECT * FROM cgapar""").fetchall()[0]]
        for i in range(0,len(fields)):
            self.pardict[fields[i]] = values[i]
        # these hold data so re-fetching from the database (slow) is not necessary
        self.dataDict = {}
        
        
    def __repr__(self):
        """String representation of the run, with values for parameters and flags."""
        repstr = 'Run database : ' + os.path.join(self.path,self.dbFileName) + '\n'
        repstr += 'Parameters :\n'
        for k in sorted(self.pardict):
            repstr += '\t%s : %s\n' % (k,self.pardict[k])
        return repstr
    
    
    def db_fetch(self,dbField):
        """Does a simple fetch of the contents of a field in the run database; fields not present in the table will
        fail, and the contents are dumped into the data dictionary, then returned.  The data dictionary is keyed on 
        the attribute names: if the key is already present, the stored value is simply returned, rather than
        re-fetched."""
        if self.dataDict.has_key(dbField):
            pass
        else:
            dbSelect = """SELECT %s FROM cgarun""" % dbField
            try:
                _temp = [x[0] for x in self.cursor.execute(dbSelect).fetchall()]
            except:
                raise CGADBFieldException(dbField)
            if dbField is 'generation':
                _temp = np.unique(_temp)
            self.dataDict[dbField] = _temp
        return self.dataDict[dbField]
                
      
    def db_calc(self,dbField,rfunc=np.mean,negOnly=True):
        """Computes functions of desired fields in the database; these should be 'reducing' functions 
        (mean, min, max) which operate on all members of the population at a given step.  Results are 
        stored keyed on (rfunc.func_name)_dbField' in the dataDictionary, and returned by this function.
        (So any custom function you use needs to have the attribute func_name.) Any necessary db_fetches 
        are performed first.  If you want to return only the negative values - for example in fitness
        functions where these are remapped NaN, Inf, and -Infs - set negOnly to True."""
        key = rfunc.func_name+dbField
        if self.dataDict.has_key(key):
            pass
        else:
            _g = self.db_fetch('generation')
            _temp = self.db_fetch(dbField)
            calcvals = []
            if not self.pardict.has_key('forestSize'):
                print 'ERROR.  Cannot determine size of forest.'
                return None
            else:
                npop = self.pardict['forestSize']
            for i in range(0,len(_g)):
                if negOnly:
                    chunk = [x for x in _temp[npop*i:npop*i+npop] if x < 0]
                else:
                    chunk = [x for x in _temp[npop*i:npop*i+npop]]
                try:
                    f = rfunc(chunk)
                except ValueError:
                    f = np.nan
                calcvals.append(f)
            self.dataDict[key] = np.asarray(calcvals)
        return self.dataDict[key]
       
    '''
    def offline_fitness_calc(self):
        """Offline performance at time t is defined as the average value, over t, of the best fitness
        that has been seen up to t."""
        if not self.dataDict.has_key('offline_fitness'):
            offline = []
            mf = self.db_calc('fitness',np.max)
            for i in range(1,len(mf)):
                offline.append(np.max(mf[0:i])) 
            offline = (np.asarray(offline).cumsum())/xrange(1,len(offline)+1)
            self.dataDict['offline_fitness'] = offline
        return self.dataDict['offline_fitness']
        
    
    def online_fitness_calc(self):
        """Online performance at time t is defined by the average fitness of all individuals that have
        been evaluated up to t steps.  This implementation assumes all members of the forest at each 
        generation are used in the average fitness, which is not strictly true."""
        if not self.dataDict.has_key('online_fitness'):
            mf = self.db_calc('fitness',np.mean)
            online = mf.cumsum()/xrange(1,len(mf)+1)
            self.dataDict['online_fitness'] = online
        return self.dataDict['online_fitness']
    '''

class CGAAnalysisTests(unittest.TestCase):
    def setUp(self):
        # database file to test - need something in ../tests directory (just rename one 
        #    created by a run)
        self.analysis = CGAAnalysis(path='../tests',dbFileName='testdb.sqldb')
        self.dbFields = ['generation','weighted_accuracy','parsimony','finitewts']
        self.funcs = [np.mean,np.max]
    
    def testRepr(self):
        print "\n\n----- testing string representation of class -----"
        print self.analysis
    
    def testDBFetch(self):
        print "\n\n----- testing fetch and return of cgarun info -----"
        for d in self.dbFields:
            temp = self.analysis.db_fetch(d)
            print 'Contents of %s (length %d): %f, . . . , %f' % (d,len(temp),temp[0],temp[-1])
    
    def testDBCalc(self):
        print "\n\n----- testing calculations using db fields -----"
        for d in self.dbFields:
            if d is not 'generation':
                for f in self.funcs:
                    if d is 'finitewts':
                        temp = self.analysis.db_calc(d,f,False)
                    else:
                        temp = self.analysis.db_calc(d,f,True)
                    print 'Contents of %s(%s) : %f, . . ., %f' %(f.func_name,d,temp[0],temp[-1])
    
    """
    def testOnlineFitnessCalc(self):
        print "\n\n----- testing calculation of online fitness -----"
        online = self.analysis.online_fitness_calc()
        print 'Online fitness : %f, . . ., %f' %(online[0],online[-1])
        
    def testOfflineFitnessCalc(self):
        print "\n\n----- testing calculation of offline fitness -----"
        offline = self.analysis.offline_fitness_calc()
        print 'Offline fitness : %f, . . ., %f' %(offline[0],offline[-1])   
    """
    
    def testPlotFitness(self):
        print "\n\n----- plotting calculated fitness values -----"
        import pylab
        g = self.analysis.db_fetch('generation')
        pylab.subplot(321)
        pylab.plot(g,self.analysis.db_calc('weighted_accuracy',np.mean,True),'k-')
        pylab.ylabel('Mean Fitness')
        pylab.xlim([g[0],g[-1]])
        pylab.subplot(322)
        pylab.plot(g,self.analysis.db_calc('weighted_accuracy',np.max,True),'k-')
        pylab.ylabel('Max Fitness')
        pylab.xlim([g[0],g[-1]])
        pylab.subplot(323)
        pylab.plot(g,self.analysis.db_calc('parsimony',np.mean,True),'k-')
        pylab.ylabel('Mean Parsimony')
        pylab.xlim([g[0],g[-1]])
        pylab.subplot(324)
        pylab.plot(g,self.analysis.db_calc('parsimony',np.max,True),'k-')
        pylab.ylabel('Max Parsimony')
        pylab.xlim([g[0],g[-1]])
        pylab.subplot(325)
        pylab.plot(g,self.analysis.db_calc('finitewts',np.mean,False),'k-')
        pylab.ylabel('Mean FW')
        pylab.xlim([g[0],g[-1]])
        pylab.ylim([0,1.1])
        pylab.subplot(326)
        pylab.plot(g,self.analysis.db_calc('finitewts',np.max,False),'k-')
        pylab.ylabel('Max FW')
        pylab.xlim([g[0],g[-1]])
        pylab.ylim([0,1.1])
        pylab.show()
        

if __name__ == 'main':
    unittest.main()
        