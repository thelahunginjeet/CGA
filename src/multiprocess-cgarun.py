#!/usr/bin/env python

"""
Multiprocessing script to search parameter space for the CGA 
"""

import CGASimulation, CGALogging, CGAParameters
from multiprocessing import Pool
import time, numpy

# get around stupid warnings
numpy.seterr(all='ignore')


proteinDBFileName = '../tests/pdz_test.db'
pdbFileName = '../tests/1iu0.pdb'
databaseDirectory = '../../parameterdb'


def worker(parameters):
	"""Individual CGA worker for a given set of parameters
		parameters -> (pC, pHC)"""
	global proteinDBFileName, pdbFileName
	assert type(parameters) is tuple and len(parameters) == 2
	# set the parameters for the simulation
	cgap = CGAParameters.CGAParameters()
	cgap.set('timing', timeSteps=300, sampGen=5)
	cgap.set(fitness={'weighted_accuracy':True, 'parsimony':True, 'finitewts':False})
	cgap.set('mutation', pC=parameters[0], pHC=parameters[1])

	# create simulation and logging
	t = time.clock()
	print 'Starting worker (pC=%.2f, pHC=%.2f) at %f'%(parameters[0], parameters[1], t)
	sim = CGASimulation.CGASimulation(databaseFile=proteinDBFileName, pdbFile=pdbFileName, cgap=cgap)
	logger = CGALogging.SqliteLogger(databaseDirectory)
	sim.attach(logger)

	# generate an initial population
	sim.populate()
	for n in range(sim.cgap.timing['timeSteps']):
		sim.advance()
		print "worker (pC=%.2f, pHC=%.2f) @ generation %d" % (parameters[0], parameters[1], n)
	sim.detach(logger)
	print 'Finished worker (pC=%.2f, pHC=%.2f) at %f' % (parameters[0], parameters[1], time.clock()-t)


def test(parameters):
	print "Starting with parameters : ", parameters


def spawn():
	"""Spawn multiple processes to search parameter space for pC and pHC
		pC -> default @ 0.7
		pHC -> default @ 0.1"""
	parameters = [] 
	for pC in numpy.arange(0.65, 0.85, 0.04):
		for pHC in numpy.arange(0.05, 0.15, 0.02):
			parameters.append((float(pC), float(pHC)))
	# launch the workers
	pool = Pool(5)
	pool.map(worker, parameters)

    
if __name__ == '__main__':
	spawn()
