#!/usr/bin/env python

'''Simple implementation of Subject/Observer Design Pattern in python.  I have modified Subject/Observer to optionally
accept and pass through kwargs dicts so that the observers do not need to know the subject's attributes and you can 
update dictionaries of data based on key/value pairs.'''


import types, unittest, numpy, os
from CGAFunctions import DataMethodFactory
from scipy import mean
from numpy import isnan, isinf, int
import datetime

def digitize_value(value):
    '''Unbound method to convert nan, inf, and -inf to floating point values (or whatever you
    want, I suppose).  Change the convert dictionary to change the conversion.'''
    convert = {'nan':-111,'inf':-11,'-inf':-1}
    if convert.has_key(str(value)):
        return convert[str(value)]
    else:
        return value


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
    """This is a logger that uses sqlite3 to log to two separate databases.  The first is called 
    'run_####.sqldb', in which #### is a unique number (timestamp) created when the observer is
    created, designed to index different runs.  The run database has two tables:
        TABLE cgarun: records the whole population at each notify()
            INTEGER generation
            REAL fitness
            TEXT function
            TEXT latex
        
        TABLE cgapar: records simulation parameters (logged only once)
            REAL probGrow 
            REAL probPrune 
            REAL probMutate 
            REAL probCross 
            TEXT tree_type 
            REAL tree_p 
            REAL tree_r 
            INTEGER forestSize
            TEXT selectionMethod
            TEXT elitism
            INTEGER eliteN
            TEXT fitnessMethod
        
    The second database is called 'cgafunctions.sqldb'.  It has a single table cgafunctions which
    counts all the functions we have ever found, over multiple runs.  
        TABLE cgafunctions: logged every time there is an update
            TEXT UNIQUE PRIMARY KEY function
            TEXT latex
            REAL fitness
            INTEGER count
    count is appropriately incremented (using a query followed by an update).
    """
    
    
    def __init__(self,path):
        assert type(path) is str
        super(SqliteLogger, self).__init__()
        self.FUNCOLS = "(function, latex, fitness, count)"
        self.PARCOLS = "(probGrow, probPrune, probMutate, probCross, tree_type, tree_p, tree_r, forestSize, selectionMethod, elitism, eliteN, fitnessMethod)"
        self.RUNCOLS = "(generation, fitness, function, latex)"
        # fire up sqlite3 and access/create the databases
        import sqlite3
        tStamp = '%s'%(datetime.datetime.now())
        tSP = tStamp.split(' ')
        runDBFileName = 'run_'+tSP[0]+'_'+tSP[1]+'.sqldb'
        funDBFileName = 'cgafunctions.sqldb'
        
        self.runconnection = sqlite3.connect(os.path.join(path,runDBFileName))
        self.runcursor = self.runconnection.cursor()
        self.funconnection = sqlite3.connect(os.path.join(path,funDBFileName))
        self.funcursor = self.funconnection.cursor()
        # this is the function table
        try:
            with self.funconnection:
                self.funconnection.execute("""CREATE TABLE IF NOT EXISTS cgafunctions (
                                                function TEXT UNIQUE PRIMARY KEY, 
                                                latex TEXT, 
                                                fitness REAL,
                                                count INTEGER);""")
        except sqlite3.IntegrityError:
            print "there was a problem initializing your table . . ."
        # metadata table
        try:
            with self.runconnection:
                # table creation
                self.runconnection.execute("""CREATE TABLE IF NOT EXISTS cgapar (
                                                probGrow REAL,
                                                probPrune REAL,
                                                probMutate REAL,
                                                probCross REAL,
                                                tree_type TEXT,
                                                tree_p REAL,
                                                tree_r REAL,
                                                forestSize INTEGER,
                                                selectionMethod TEXT,
                                                elitism TEXT,
                                                eliteN INTEGER,
                                                fitnessMethod TEXT);""")
        except sqlite3.IntegrityError:
            print "there was a problem initializing your table . . ."
        # run table
        try:
            with self.runconnection:
                # table creation
                self.runconnection.execute("""CREATE TABLE IF NOT EXISTS cgarun (
                                                generation INTEGER,
                                                fitness REAL,
                                                function TEXT,
                                                latex TEXT);""")
        except sqlite3.IntegrityError:
            print 'there was a problem initializing your table . . .'
            
                
    def update(self, subject, **kwargs):
        for chromosome in subject.population:
            # insert everything from population
            tree = chromosome.tree
            function, latex, generation, fitness = tree.getString(), tree.getLatex(), kwargs['time'], digitize_value(chromosome.fitness)
            runvals = tuple([generation,fitness,function,latex])
            # run table update
            self.runconnection.execute("""INSERT OR REPLACE INTO cgarun %s VALUES %s"""%(self.RUNCOLS,runvals))
            self.runconnection.commit()
            # inserts/updates into the uniques table
            self.funcursor.execute("""SELECT count FROM cgafunctions WHERE function='%s'"""%function)
            cTup = self.funcursor.fetchall()
            if len(cTup) == 0:
                # first appearance - count should be 1
                funvals = tuple([function,latex,fitness,1])
                self.funcursor.execute("""INSERT INTO cgafunctions %s VALUES %s"""%(self.FUNCOLS,funvals))
            else:
                # has appeared before - update counter
                self.funcursor.execute("""UPDATE cgafunctions SET count=%d+1 WHERE function='%s'"""%(cTup[0][0],function))
            self.funconnection.commit()
        # update the cgaruns table - only occurs during first notify()
        if int(kwargs['time']) == 0:
            parvals = (subject.pG,subject.pP,subject.pM,subject.pC,subject.treeGenDict['treetype'],subject.treeGenDict['p'],subject.treeGenDict['r'],subject.forestSize,subject.selectionMethod,subject.elitism.__repr__(),subject.eliteN,subject.fitnessMethod)
            self.runconnection.execute("""INSERT OR REPLACE INTO cgapar %s VALUES %s"""%(self.PARCOLS,parvals))
            self.runconnection.commit()
        

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
        obs = SqliteLogger('../tests')
        self.logMe.attach(obs)
        
        
if __name__ == 'main':
    pass