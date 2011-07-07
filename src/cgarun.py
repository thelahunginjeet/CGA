#!/usr/bin/env python

"""
scr_cgarun.py

Example script to do offline runs of the CGA.
"""

# if you want to collect data, import CGALogging
import CGASimulation, CGALogging
from numpy import max

def run_cga():
    # which db and pdb file to use?
    proteinDBFileName = '../tests/pdz_test.db'
    pdbFileName = '../tests/1iu0.pdb'
    # how many generations to run for?
    nGen = 1500
    # set up a simulation
    #    treeType:
    #        fixed = fixed number of nodes
    #        p = tree size (#nonterminal nodes)
    #        r = bias parameter (prob. of Binary vs Unary node)
    mySim = CGASimulation.CGASimulation(databaseFile=proteinDBFileName, pdbFile=pdbFileName, sampGen=10, selectionMethod='pareto_tournament', fitnessMethod='topN_weighted_accuracy', forestSize=30, treeGenDict={'treetype':'fixed','p':15,'r':0.5})
     # create and attach a DataLogger
    dataLogger = CGALogging.SqliteLogger('../tests')
    mySim.attach(dataLogger)
    # generate an initial population
    mySim.populate()
    # now start running and logging data
    for n in range(0, nGen):
        print 'working on generation %d' % n
        mySim.advance()
    mySim.detach(dataLogger)
    
    
if __name__ == '__main__':
    rundata = run_cga()
