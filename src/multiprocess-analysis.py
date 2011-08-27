#!/usr/bin/env python

"""
Script to analyze the results of a multiprocess run.
"""

from CGAAnalysis import CGAAnalysis
import pylab, glob, os, numpy


databases = "../../parameterdb/run*.sqldb"


def order():
	files = glob.glob(databases)
	parameters = []
	for f in files:
		path, dbfile = os.path.split(f)
		analysis = CGAAnalysis(path, dbfile)
		parameters.append((analysis.pardict['pC'], analysis.pardict['pHC']))
	return sorted(zip(parameters, files))


def main():
	pairs = order()
	files, parameters = [x[1] for x in pairs], [x[0] for x in pairs]
	figure = pylab.figure(figsize=(12, 12))
	for i in range(5):
		for j in range(5):
			path, dbfile = os.path.split(files[5*i + j])
			analysis = CGAAnalysis(path, dbfile)
			pC, pHC = parameters[5*i + j]
			subplot = figure.add_subplot(5, 5, 5*i + j + 1)
			g = analysis.db_fetch('generation')
			subplot.plot(g, analysis.db_calc('parsimony', numpy.max, True), 'k-')
			subplot.plot(g, analysis.db_calc('parsimony', numpy.mean, True), 'r-')
			#--------------------------------------------------
			# subplot.plot(g, analysis.db_calc('weighted_accuracy', numpy.max, True), 'k-')
			# subplot.plot(g, analysis.db_calc('weighted_accuracy', numpy.mean, True), 'r-')
			#-------------------------------------------------- 
			if j == 0:
				subplot.set_ylabel("pC=%.2f"%(float(pC)))
				subplot.set_yticks([-2000, -1500, -1000, -500, 0])
				# subplot.set_yticks([-1.5, -1.0, -0.5, 0.0])
			else:
				subplot.set_yticks([])
			if i == 4:
				subplot.set_xlabel("pHC=%.2f"%(float(pHC)))
			subplot.set_ylim([-2000, 0])
			subplot.text(45, -1900, "%.1fMb"%(os.path.getsize(files[5*i+j])/1000000.0,))
			# subplot.set_ylim([-1.5, 0.0])
			subplot.set_xticks([])
	pylab.show()
	

if __name__ == '__main__':
	# main()
	order()
