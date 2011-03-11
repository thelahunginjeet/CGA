#!/usr/bin/env python

"""
scr_cgarun.py

Example script to do offline runs of the CGA.
"""

# if you want to collect data, import CGALogging
import CGASimulation,CGALogging
# for plotting at the end
import pylab

def run_cga():
    # which db and pdb file to use?
    dataBaseFileName = '../tests/pdz_test.db'
    pdbFileName = '../tests/1iu0.pdb'
    # how many generations to run for?
    nGen = 200
    # set up a simulation
    mySim = CGASimulation.CGASimulation(databaseFile=dataBaseFileName,pdbFile=pdbFileName,forestSize=50)
     # create and attach a DataLogger
    dataLogger = CGALogging.DataLogger()
    mySim.attach(dataLogger)
    # generate an initial population (default will be exponentially distributed tree sizes
    mySim.populate()
    # now start running and logging data
    for n in range(0,nGen):
        mySim.advance()
    # plot some stuff
    pylab.subplot(221)
    pylab.plot(dataLogger.data['time'],dataLogger.data['maxSize'],'r-')
    pylab.plot(dataLogger.data['time'],dataLogger.data['minSize'],'b-',hold=True)
    pylab.title('Tree Size Bounds (Min/Max)')
    pylab.subplot(222)
    pylab.plot(dataLogger.data['time'],dataLogger.data['maxFit'],'k-')
    pylab.title('Maximum Fitness')
    pylab.subplot(223)
    pylab.plot(dataLogger.data['time'],dataLogger.data['meanFit'],'k-')
    pylab.title('Mean Fitness')
    pylab.xlabel('Generation Number')
    pylab.subplot(224)
    pylab.plot(dataLogger.data['time'],dataLogger.data['wellFormed'],'k-')
    pylab.title('Well-formed Fraction')
    pylab.xlabel('Generation Number')
    pylab.show()
    # detach the observer, but return the data contents
    mySim.detach(dataLogger)
    return dataLogger.data
    
if __name__ == 'main':
    rundata = run_cga()