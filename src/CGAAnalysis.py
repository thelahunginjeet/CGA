#!/usr/bin/env/python

import unittest,os
import sqlite3
import numpy as np
from pyparsing import *

"""Utility functions used in the parser."""
def flatten(l):
    """Generator that flattens a list."""
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el,basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def unique(l, idfun=repr):
    """There are lots of ways to return a list of unique items in a sequence; here is one."""
    seen = {}
    return [seen.setdefault(idfun(e),e) for e in l if idfun(e) not in seen]

def strip_brackets(s,brack='{}'):
    """Replaces open/close brackets with whitespace.  brack should be '{}','[]','()' or
    somesuch.  Useful to remove nested brackets in strings."""
    return s.replace(brack[0],'').replace(brack[1],'')

def fix_latex(ls,fakeCC='$'):
    """To avoid lots of unnecessary brackets, we use a different command character in the db
    latex strings and then replace it with the proper character afterwards."""
    return ls.strip().replace(fakeCC,'\\')


"""Exceptions for the CGAAnalysis classes."""
class CGADBFieldException(Exception):
    def __init__(self,dbField):
        print "Field %s does not exist in database. Check table_info(dbName)." % dbField

class CGADBException(Exception):
    def __init__(self,dbFileName):
        print "Database %s does not exist." % dbFileName

class CGADBCmpException(Exception):
    def __init__(self,cmp):
        print "Comparator %s not allowed for database fetches." % cmp
        

class CGAFunctionParser(object):
    """Post-hoc parsing and tokenizing of strings written to the database.  Current operator
    and variable list (not at all in BNF, just a shorthand - [] = optional):
        variable :: p_ij, p_i, p_j, or N
        constant :: E or PI
        fnumber  :: [+/-]+digits+[.digits]
        operand  :: variable, constant, or fnumber
        prefix   :: tr, det, s_m, exp, sum_ij, log, tanh, sinh, cosh
        postfix  :: **2 or ^T
        plusop   :: + or -
        multop   :: * or /
    """
    def __init__(self):
        self._variable = Literal("p_ij") | Literal("p_i") | Literal("p_j") | Literal("N")
        self._constant = Literal("PI") | Literal("E")
        self._fnumber = Combine(Optional("-")+Word(nums)+Optional("."+Word(nums)))
        self._operand = self._variable | self._constant | self._fnumber
        
        self._prefix = Literal("tr") | Literal("exp") | Literal("sum_ij") | Literal("log") | Literal("tanh") | Literal("sinh") | Literal("cosh") | Literal("det") | Literal("s_m")
        self._postfix = Literal("**2") | Literal("^T")
        self._plusop = oneOf('+ -')
        self._multop = oneOf('* /')
        
        self._expr = operatorPrecedence(self._operand,
            [(self._postfix, 1, opAssoc.LEFT),\
             (self._prefix, 1, opAssoc.RIGHT), \
             (self._plusop, 2, opAssoc.LEFT), \
             (self._multop, 2, opAssoc.LEFT),])
    
    def __repr__(self):
        """The pieces of the grammar in pyparsing have representations, so we can pass those along to
        generate some kind of readable representation of the parser's grammar."""
        repStr = 'Grammar for CGA Parser :\n'
        repStr += '   variables  :: %s\n' % strip_brackets(self._variable.__repr__())
        repStr += '   constants  :: %s\n' % strip_brackets(self._constant.__repr__())
        repStr += '   numbers    :: %s\n' % strip_brackets(self._fnumber.__repr__())
        repStr += '   operands   :: variables | constants | numbers\n'
        repStr += '   prefix op  :: %s\n' % strip_brackets(self._prefix.__repr__())
        repStr += '   postfix op :: %s\n' % strip_brackets(self._postfix.__repr__())
        repStr += '   plus op    :: %s\n' % strip_brackets(self._plusop.__repr__())
        repStr += '   mult op    :: %s\n' % strip_brackets(self._multop.__repr__())
        return repStr
    
    def parse(self,t):
        return self._expr.parseString(t)
        
    
    def tokenize(self,t):
        tokenList = [x for x in flatten(self._expr.parseString(t))]
        # count the tokens and return a dictionary
        tokens = {}.fromkeys(unique(tokenList))
        for k in tokens:
            tokens[k] = tokenList.count(k)
        return tokens


class CGAFunctionAnalysis(object):
    """Performs post-hoc analysis on the cgafunctions database, which logs everything ever seen
    throughout multiple runs.  The filename will typically be cgafunctions.sqldb, but this is 
    not required - the SQL table name (cgafunctions) is, however."""
    def __init__(self,path,dbFileName='cgafunctions.sqldb'):
        self.path = path
        self.dbFileName = dbFileName
        fName = os.path.join(self.path,self.dbFileName)
        if not os.path.exists(fName):
            raise CGADBException(fName)
        self.connection = sqlite3.connect(fName)
        self.cursor = self.connection.cursor()
        self.fields = [x[1] for x in self.cursor.execute("""PRAGMA table_info(cgafunctions)""").fetchall()]
        self.fieldMap = {}
        count = 0
        # gives the fieldname to position mapping
        for x in self.fields:
            self.fieldMap[x] = count
            count += 1

        
    def __repr__(self):
        """String representation of contents of the database."""
        repstr = 'Function database: ' + os.path.join(self.path,self.dbFileName) + '\n'
        repstr += 'Fields logged:\n'
        for k in self.fields:
            repstr += '\t%s\n' % k
        return repstr
    
    def select_records(self,dbField,cmp,cutoff,negOnly=True):
        """This database can be large, so unlike CGAAnalysis(), the fetching function here uses
        SQL to pull only records matching the desired criteria, rather than prepulling everything
        and using python to sort it.  Usage example:
            dbField = parsimony
            cmp = gt
            cutoff = -10
            negOnly = True
        will pull records from the function database whose parsimony value is greater than -10, 
        dropping any records with positive parsimony.  All values in the table for each record
        are returned, in a dictionary keyed on fields. 
        
        dbField here acts as a field to select on, not the only field which is returned (unlike 
        in CGAAnalysis).  'negOnly' is mostly useful for ignoring fitness values mapped to special 
        values : 1,11,111 for nan,-inf,inf, for example.
        
        Allowed values for kwargs:
            dbField : any allowed db field listed in self.fields
            cmp     : 'gt','lt','gteq', 'lteq', 'eq'
            cutoff  : any floating point number or integer
            negOnly : True or False
            
        """
        if dbField not in self.fields:
            raise CGADBFieldException(dbField)
        if cmp not in ('gt','lt','lteq','gteq','eq'):
            raise CGADBCmpException(cmp)
        op = None
        if cmp == 'gt':
            op = '>'
        elif cmp == 'lt':
            op = '<'
        elif cmp == 'gteq':
            op = '>='
        elif cmp == 'lteq':
            op = '<='
        elif cmp == 'eq':
            op = '=='
        # now do the pull
        dbSelect = """SELECT * FROM cgafunctions WHERE %s %s %f""" %(dbField,op,cutoff)
        _temp = self.cursor.execute(dbSelect).fetchall()
        if negOnly:
            records = [x for x in _temp if x[self.fieldMap[dbField]] <= 0]
        else:
            records = _temp
        recordDict = dict().fromkeys(self.fields)
        for k in recordDict:
            recordDict[k] = list(zip(*records)[self.fieldMap[k]])
            if k == 'latex':
                recordDict[k] = map(fix_latex,recordDict[k])
        return recordDict
        
            
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
        self.functions = CGAFunctionAnalysis(path='../tests',dbFileName='cgafunctions.sqldb')
        self.dbFields = ['generation','weighted_accuracy','parsimony','finitewts']
        self.funcs = [np.mean,np.max]
        self.parser = CGAFunctionParser()
        
    def testParsing(self):
        print "\n\n----- testing string representation of classes -----"
        tests = [ "p_i**2",\
                       "cosh(p_i)^T + E", \
                       "(p_i*(p_i)^T)**2",\
                       "exp((p_i+2.10987719))",\
                       "s_m(tanh(-1.0))",\
                       "det(PI)", \
                       "p_ij*(1/N)", \
                       "(p_i+p_j)",\
                       "log(p_i*(p_ij+E))",\
                       "tr(p_ij/(p_i+(p_j)^T))",\
                       "sum_ij((p_ij*(p_i+PI)+tr(PI))*(p_i*(p_j)^T))"]
        for t in tests:
            try:
                print t, "->", self.parser.parse(t)
            except:
                print "FAILED"
                
    
    def testRepr(self):
        print "\n\n----- testing string representation of classes -----"
        print 'Repr. of CGAAnalysis:'
        print self.analysis
        print 'Repr. of CGAFunctionAnalysis:'
        print self.functions
        print 'Repr. of CGAParser:'
        print self.parser
    
    def testDBFetch(self):
        print "\n\n----- testing fetch and return of cgarun info -----"
        for d in self.dbFields:
            temp = self.analysis.db_fetch(d)
            print 'Contents of %s (length %d): %f, . . . , %f' % (d,len(temp),temp[0],temp[-1])
    
    def testFunctionFetch(self):
        print "\n\n----- testing fetch of cgafunction records -----"
        recordDict = self.functions.select_records('parsimony',cmp='gt',cutoff=-10,negOnly=False)
        print 'Number of fields returned : %d' %len(recordDict)
        print 'Number of records : %d' % len(recordDict['parsimony'])
        print 'Min, max of selected field in records : %f, %f' % (np.min(recordDict['parsimony']),np.max(recordDict['parsimony']))
        
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
        