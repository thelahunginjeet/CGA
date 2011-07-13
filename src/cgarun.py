#!/usr/bin/env python

"""
scr_cgarun.py

Example script to do offline runs of the CGA.
"""

# if you want to collect data, import CGALogging
import CGASimulation, CGALogging, CGAParameters
from numpy import max

def run_cga():
    # which db and pdb file to use?
    proteinDBFileName = '../tests/pdz_test.db'
    pdbFileName = '../tests/1iu0.pdb'
    # usage of new CGAParameters object; any necessary missing parameters are filled in by 
    #    default.  You can set single parameters in a field, or multiple parameters as a
    #    dict.
    cgap = CGAParameters.CGAParameters()
    cgap.set('timing',timeSteps=1500)
    cgap.set(fitness={'weighted_accuracy':True,'parsimony':True,'finitewts':False})
    mySim = CGASimulation.CGASimulation(databaseFile=proteinDBFileName, pdbFile=pdbFileName,cgap=cgap)
     # create and attach a DataLogger
    dataLogger = CGALogging.SqliteLogger('../tests')
    mySim.attach(dataLogger)
    # generate an initial population
    mySim.populate()
    # now start running and logging data (this bit is awkward; the simulation
    #    knows how long it is going to run, so we could just modify advance()
    #    to move one or more steps ahead, or add another wrapper).
    for n in range(0, mySim.cgap.timing['timeSteps']):
        print 'working on generation %d' % n
        mySim.advance()
    mySim.detach(dataLogger)
    
    
if __name__ == '__main__':
    rundata = run_cga()
