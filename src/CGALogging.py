#!/usr/bin/env python

'''Simple implementation of Subject/Observer Design Pattern in python.  I have modified Subject/Observer to optionally
accept and pass through kwargs dicts so that the observers do not need to know the subject's attributes and you can 
update dictionaries of data based on key/value pairs.'''


import types, unittest, numpy
from hashlib import md5 
from CGAFunctions import DataMethodFactory


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
    """This is a logger that uses sqlite3 to log to a database file with the format:
            INTEGER md5-hash;
            TEXT function
            TEXT latex
            INTEGER generation
            INTEGER (one for each) f_log, f_+, etc."""
    
    def __init__(self, dbFile):
        assert type(dbFile) is str
        import sqlite3
        super(SqliteLogger, self).__init__()
        self.connection = sqlite3.connect(dbFile)
        # prepare the variable function columns
#        self.funcs = sorted(DataMethodFactory().data.keys())
#        self.funcs += sorted(DataMethodFactory().unary.keys())
#        self.funcs += sorted(DataMethodFactory().binary.keys())
#        self.funcs += sorted(DataMethodFacctor().scalars.keys())
        try:
            with self.connection:
                self.connection.execute("""CREATE TABLE IF NOT EXISTS cgafunctions (
                                                md5 TEXT UNIQUE PRIMARY KEY,
                                                function TEXT, 
                                                latex TEXT, 
                                                generation INTEGER,
                                                fitness REAL);""")
        except sqlite3.IntegrityError:
            print "there was a problem initializing your database . . ."
                
                
    def update(self, subject, **kwargs):
        print "number of trees : ", len(subject.population)
        for chromosome in subject.population:
            tree = chromosome.tree
            function = "'" + tree.getString() + "'"
            hash = md5(function).hexdigest()
            latex = "'" + tree.getLatex() + "'"
            generation = kwargs['time']
            fitness = chromosome.fitness
            values = (hash, function, latex, generation, fitness)
            self.connection.execute("""INSERT OR REPLACE INTO cgafunctions (md5, function, latex, generation, fitness) VALUES 
                                        (?, ?, ?, ?, ?)""", values)
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
        self.logMe.notify()
        
        
if __name__ == 'main':
    pass