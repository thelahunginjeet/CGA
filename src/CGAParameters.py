#!/usr/bin/env python

import unittest

# TODO:
#  - better parameter checking of values in CGAParameters; things are currently pretty unsafe
#  - should all parameters be assumed to be dictionaries (CGAParameters is basically a dict of dicts)?

class CGAParameters(object):
    """A dict of dicts for simulation parameters; has default values necessary for CGASimulations.  These
    can be overwritten before the simulation begins, and CGASimulation does some checking on the 
    parameters."""
    def __init__(self):
        self.treep = {'forestSize':30,'treetype':'fixed','p':5,'r':0.6}
        self.mutation = {'pG':0.025, 'pP':0.025, 'pHC':0.1, 'pM':0.05, 'pC':0.7}
        self.timing = {'timeSteps':100,'sampGen':10}
        self.selection = {'method':'pareto_tournament'}
        self.fitness = {'weighted_accuracy':True,'parsimony':True,'finitewts':False}
        self.parFields = ['treep','mutation','selection','fitness','timing']

    def set(self,name=None,**kwargs):
        """ set syntax is as follows:
                set(string={...}) sets CGAParameters.string = {...}
                set('string',x=0.1,y=0.2,...) sets individual key/value pairs in parameters['string']
        """
        for k in kwargs:
            if name is None:
                setattr(self,k,kwargs[k])
            else:
                if hasattr(self,name):
                    getattr(self,name)[k] = kwargs[k]
                else:
                    print 'ERROR: Parameter field \'%s\' does not exist.' % name
            
        
    def __repr__(self):
        repstr = '-Parameter values (%d fields set)-\n' % len(self.parFields)
        for k in self.parFields:
            repstr += '   -%s-\n' % k
            attToPrint = getattr(self,k)
            for l in attToPrint:
                repstr += '     %s : %s\n' % (l,attToPrint[l])
        return repstr
    

class CGAParametersTests(unittest.TestCase):
    def setUp(self):
        self.cgap = CGAParameters()
    
    def testDictSet(self):
        print"\n\n----- testing update via full dictionary -----"
        print 'BEFORE DICT UPDATE:'
        print self.cgap
        self.cgap.set(treep={'forestSize':30, 'treetype':'fixed','p':5,'r':0.6})
        self.cgap.set(mutation={'pG':0.025, 'pP':0.025, 'pHC':0.1, 'pM':0.05, 'pC':0.7})
        self.cgap.set(fitness={'weighted_accuracy':True,'parsimony':True})
        print 'AFTER DICT UPDATE:'
        print self.cgap
    
    def testKeyValueSet(self):
        print"\n\n----- testing update via key/value pairs -----"
        self.cgap.set(treep={'forestSize':30, 'treetype':'fixed','p':5,'r':0.6})
        self.cgap.set(mutation={'pG':0.025, 'pP':0.025, 'pHC':0.1, 'pM':0.05, 'pC':0.7})
        print 'BEFORE KEY VALUE UPDATES:'
        print self.cgap
        self.cgap.set('treep',forestSize=10000)
        self.cgap.set('mutation',pG=0.1)
        self.cgap.set('selection',method='pareto_tournament')
        print 'AFTER KEY/VALUE UPDATES:'
        print self.cgap


if __name__ == '__main__':
    unittest.main()