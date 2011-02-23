#!/usr/bin/env python

import unittest, operator, random
import numpy as MATH
from numpy.random import randint, uniform
from CGAPreprocessing import Utilities

class Data(object):
	"""Simple container for keeping information about a data source"""
	def __init__(self, string, latex, data):
		self.string = string
		self.latex = latex
		# note : self.function just return data object for bottom of recursion for consistency
		self.function = data
		

class Function(object):
	"""Simple container for keeping information about a function"""
	def __init__(self, string, latex, method):
		self.string = string
		self.latex = latex # use raw strings for string with escape characters
		self.function = method
		

class Functions(dict):
	"""General function class to subclass"""
	def __init__(self):
		pass
		
	def returnRandom(self):
		return self[self.keys()[randint(0,len(self))]]
		

class ImmutableData(Functions):
	"""Simple class for immutable data items; this is where we would store problem-specific 
	loaded-in data, or constants like e, pi, etc.  This class can't deal with adding a random value; 
	that still has to be done elsewhere."""	
	def __init__(self):
		super(ImmutableData, self).__init__()
		self['e'] = Data("e", r'e', MATH.e)
		self['pi'] = Data("pi", r'\pi', MATH.pi)
		self['1'] = Data("1", r'1', 1.0)
		self['1/2'] = Data("1/2", r'\frac{1}{2}', 0.5)
		self['1/N'] = Data("1/N", r'\frac{1}{N}', 1./20.)

		
class ProteinData(Functions):
	"""Simple class to hold the protein frequency data; the initial (2x2) matrix associated with the
	data elements ensures initial evaluteability but is just a dummy.  For real applications, external
	functions update the node value to the current column from the MSA upon evaluation/looping.
		  note : the indices indicate the ith column, jth column (>= i), and the joint i,j"""
	def __init__(self, databaseFile):
		super(ProteinData, self).__init__()
		self['p_i'] = Data("p_i", r'\rho_{i}', MATH.array([[0.1,0.1],[0.9,0.9]]))
		self['p_j'] = Data("p_j", r'\rho_{j}', MATH.array([[0.3,0.3],[0.7,0.7]]))
		self['p_ij'] = Data("p_ij", r'\rho_{ij}', MATH.array([[0.2,0.4],[0.8,0.6]]))


class ScalarizingFunctions(Functions):
	"""This is a special set of functions that convert arrays into scalars; the final step in almost all
	scoring methods will be double summation over amino acids, sine the input data is in the form
	of matrices and all other operators preserve this shape.  Therefore, for the coevolving
	residue prediction problem, the root of the tree *must* be a scalarizing operator.  However, the 
	scalarizing functions are written to work on scalars as well (just returning the scalar), in order
	to run simple tests on non-protein data. """
	def __init__(self):	
		super(ScalarizingFunctions, self).__init__()
		self['tr'] = Function("tr(%s)", r'{\mathrm Tr}\left\{%s\right\}', self.trace)
		self['sum_ij'] = Function("sum_ij(%s)", r'\Sigma_{ij}\left(%s\right)', self.dsum)
		
	def trace(self,x):
		try:
			y = x.trace()
		except ValueError:
			y = x
		return y
	
	def dsum(self,x):
		try:
			y = sum(sum(x))
		except TypeError:
			y = x
		return y
	

class UnaryFunctions(Functions):
	"""Simple class for unary functions; these all operate array-wise, returning an array if used on
	an array."""	
	def __init__(self):
		super(UnaryFunctions, self).__init__()
		self['sin'] = Function("sin(%s)", r'\sin\left(%s\right)', MATH.sin)
		self['exp'] = Function("exp(%s)", r'\exp\left(%s\right)', MATH.exp)
		self['log'] = Function("log(%s)", r'\log\left(%s\right)', MATH.log)
		self['tanh'] = Function("tanh(%s)",r'\tanh\left(%s\right)', MATH.tanh)

						
class BinaryFunctions(Functions):
	"""Simple class for binary functions; these also operate element-wise on arrays."""
	def __init__(self):
		super(BinaryFunctions, self).__init__()
		self['+'] = Function("(%s+%s)", r'\left(%s+%s\right)', MATH.add)
		self['-'] = Function("(%s-%s)", r'\left(%s-%s\right)', MATH.subtract)
		self['*'] = Function("(%s*%s)", r'%s\dot%s', MATH.multiply)
		self['/'] = Function("(%s/%s)", r'\frac{%s}{%s}', MATH.divide)
		

class CGAFunctionsTests(unittest.TestCase):
	def setUp(self):
		self.binaryFunctions = BinaryFunctions()
		self.unaryFunctions = UnaryFunctions()
		self.data = ImmutableData()
		self.scalFunctions = ScalarizingFunctions()
		self.proteinData = ProteinData('../tests/pdz_test.db')
	
	def testName(self):
		print "\n----- testing function names -----"
		self.assertEquals(self.binaryFunctions['+'].string, '(%s+%s)')
		self.assertEquals(self.unaryFunctions['log'].string,'log(%s)')
		self.assertEquals(self.scalFunctions['tr'].string, 'tr(%s)')
		self.assertEquals(self.data['pi'].string,'pi')
		self.assertEquals(self.proteinData['p_ij'].string,'p_ij')
		
	def testData(self):
		print "\n----- testing data containers -----"
		print 'e = ',self.data['e'].function
		print 'pi = ',self.data['pi'].function
		print '1 = ',self.data['1'].function
		print '1/2 = ',self.data['1/2'].function
		print '1/N = ',self.data['1/N'].function
		print 'p_i = ',self.proteinData['p_i'].function
		print 'p_j = ',self.proteinData['p_j'].function
		print 'p_ij = ',self.proteinData['p_ij'].function

if __name__ == '__main__':
	unittest.main()