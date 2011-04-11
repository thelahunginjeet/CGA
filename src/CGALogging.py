#!/usr/bin/env python

'''Simple implementation of Subject/Observer Design Pattern in python.  I have modified Subject/Observer to optionally
accept and pass through kwargs dicts so that the observers do not need to know the subject's attributes and you can 
update dictionaries of data based on key/value pairs.'''


import types, unittest, numpy
from CGAFunctions import DataMethodFactory
from scipy import mean
from numpy import isnan, isinf, int

class Subject(object):
    '''Standard observer pattern, but I have modified notify to accept a kwargs dict; the update() function doesn't
    need to look at them (and they could be skipped).
    '''
    def __init__(self):
        self._observers = []
    
    def attach(self, observer):
        if not observer in self._observers:
            self._observers.append(observer)
            
    def detach(self, observer):
        try:
            self._observers.remove(observer)
        except ValueError:
            pass
    
    def notify(self, modifier=None, **kwargs):
        for observer in self._observers:
            if modifier != observer:
                observer.update(self,**kwargs)


class Observer(object):
    '''Putting a pass in update() allows anything derived from Observer to be attached
    and updated without failure, though potentially without any results either, if
    update() is not overloaded.
    '''
    def __init__(self):
        pass

    def update(self, subject, **kwargs):
        pass

    
class DataLogger(Observer):
    '''The data logger keeps a dict of lists, for example:
            dataLogger.data['time'] = [0,1,2,3,...]
            dataLogger.data['val_1'] = [0.1,0.2,0.3,...]
            dataLogger.data['val_2'] = [1.1,2.3,4.5,...]
        the lists can be lists of anything (even other lists).  The notification/update cycle
        is called thusly:
            subject.notify('time'=0,'val_1'=0.1,'val_2'=1.1)
    '''
    def __init__(self):
        super(DataLogger,self).__init__()
        self.data = dict()
    
    def update(self,subject,**kwargs):
        for k in kwargs:
            if self.data.has_key(k):
                self.data[k].append(kwargs[k])
            else:
                self.data[k] = []
                self.data[k].append(kwargs[k])


class SqliteLogger(Observer):
    """This is a logger that uses sqlite3 to log to a database file two tables.  The first
    is called 'cgafunctions' and has the format:
            TEXT function
            TEXT latex
            INTEGER generation
            REAL fitness
            REAL minFitness
            REAL meanFitness
            REAL maxFitness 
                + number of function/data nodes per type (e.g. INTEGER f_log)
    'cgafunctions' is updated/added to during the run. The second table, which stores 
    simulation parameters (only logged at the beginning, never updated) is called
    'cgaruns' and has format:
            REAL probGrow
            REAL probPrune
            REAL probMutate
            REAL probCross
            REAL treePar (interpretation depends on treeType)
            INTEGER forestSize
            TEXT treeType
            TEXT selectionMethod
            DATE stamp
    """
    
    
    def __init__(self, dbFile):
        assert type(dbFile) is str
        super(SqliteLogger, self).__init__()
        # prepare the variable function columns by pointing to correct getter
        self.funcs = {}
        dataFactory = DataMethodFactory()
        for data in dataFactory.data:
            if data in ('p_i', 'p_j', 'p_ij'):
                self.funcs[dataFactory.data[data][0]] = data
        for unary in dataFactory.unary:
            self.funcs[dataFactory.unary[unary][0]] = unary
        for binary in dataFactory.binary:
            self.funcs[dataFactory.binary[binary][0]] = binary
        self.forder = sorted(self.funcs.values())
        fcolumns = ""
        self.COLUMNS = "(function, latex, generation, fitness, min_fitness, mean_fitness, max_fitness, "
        for func in self.forder:
            fcolumns += "f_%s INTEGER,"%(func)
            self.COLUMNS += "f_%s, "%(func)
        self.COLUMNS = self.COLUMNS[:-2] + ")" 
        self.QUESTIONS = "(%s)"%(((len(self.forder) + 7)*'?,')[:-1])
        # this is for the metadata 
        self.RUNCOLS = "(generation,probGrow,probPrune,probMutate,probCross,treePar_p,treePar_r,forestSize,treeType,selectionMethod)"
        # fire up the sqlite
        import sqlite3
        self.connection = sqlite3.connect(dbFile)
        # this is the function table
        try:
            with self.connection:
                self.connection.execute("""CREATE TABLE IF NOT EXISTS cgafunctions (
                                                function TEXT UNIQUE PRIMARY KEY, 
                                                latex TEXT, 
                                                generation INTEGER,
                                                fitness REAL,
                                                min_fitness REAL,
                                                mean_fitness REAL,
                                                max_fitness REAL,
                                                %s);"""%(fcolumns[:-1]))
        except sqlite3.IntegrityError:
            print "there was a problem initializing your database . . ."
        # metadata table
        try:
            with self.connection:
                # table creation
                self.connection.execute("""CREATE TABLE IF NOT EXISTS cgaruns (
                                                generation INTEGER UNIQUE PRIMARY KEY,
                                                probGrow REAL,
                                                probPrune REAL,
                                                probMutate REAL,
                                                probCross REAL,
                                                treePar_p REAL,
                                                treePar_r REAL,
                                                forestSize INTEGER,
                                                treeType TEXT,
                                                selectionMethod TEXT);""")
        except sqlite3.IntegrityError:
            print "there was a problem initializing your database . . ."
            
                
    def update(self, subject, **kwargs):
        try:
            minFit = min([x.fitness for x in subject.population if not isnan(x.fitness) and not isinf(x.fitness)])
        except ValueError:
            minFit = 0.0
        try:
            meanFit = mean([x.fitness for x in subject.population if not isnan(x.fitness) and not isinf(x.fitness)])
        except ValueError:
            meanFit = 0.0
        try:
            maxFit = max([x.fitness for x in subject.population if not isnan(x.fitness) and not isinf(x.fitness)])        
        except ValueError:
            maxFit = 0.0
        funcs = {}.fromkeys(self.forder)
        for chromosome in subject.population:
            # reset number dictionary
            for f in funcs:
                funcs[f] = 0
            tree = chromosome.tree
            function, latex, generation, fitness = tree.getString(), tree.getLatex(), kwargs['time'], chromosome.fitness
            nodes = tree.getNodes()
            # determine number of each node
            for node in [n for n in nodes if n.string in self.funcs]:
                funcs[self.funcs[node.string]] += 1
            values = [function, latex, generation, fitness, minFit, meanFit, maxFit] + [funcs[x] for x in self.forder]            
            self.connection.execute("""INSERT OR REPLACE INTO cgafunctions %s VALUES %s"""%(self.COLUMNS, self.QUESTIONS), values)
        self.connection.commit()
        # update the cgaruns table - only write the majority of the metadata for generation zero
        if int(kwargs['time']) == 0:
            runvals = (kwargs['time'],subject.pG,subject.pP,subject.pM,subject.pC,subject.treeGenDict['p'],subject.treeGenDict['r'],subject.forestSize,subject.treeGenDict['treetype'],subject.selectionMethod)
            self.connection.execute("""INSERT OR REPLACE INTO cgaruns %s VALUES %s"""%(self.RUNCOLS,runvals))
        else:
            self.connection.execute("""INSERT OR REPLACE INTO cgaruns (generation) VALUES (%d)""" % kwargs['time'])
        self.connection.commit()
        

class CGALoggingTests(unittest.TestCase):
    def setUp(self):
        # a simple class to log
        class LogMe(Subject):
            def __init__(self):
                super(LogMe,self).__init__()
                self.ascalar = 1.0
                self.astring = 'mystring'
                self.alist = [1.0,2.0,3.0]
                self.afunction = lambda x: x**2
        self.logMe = LogMe()
                
    def testNotification(self):
        obs = DataLogger()
        self.logMe.attach(obs)
        self.logMe.notify(value=self.logMe.ascalar,text=self.logMe.astring,thelist=self.logMe.alist,thefunc=self.logMe.afunction)
        self.logMe.notify(value=self.logMe.ascalar+1)
        self.logMe.notify(thefunc=lambda x : x**3)
        print 'Logged data : '
        for k in obs.data.keys():
            print k,obs.data[k]
        self.logMe.detach(obs)
        
    def testSqliteLogger(self):
        obs = SqliteLogger('../tests/test.db')
        self.logMe.attach(obs)
        # can't notify because .logMe isn't a simulation object
#        self.logMe.notify()
        
        
if __name__ == 'main':
    pass