#!/usr/bin/env python

import unittest, operator, random
import numpy as MATH
from numpy.random import randint, uniform
from CGAPreprocessing import Utilities


class DataFunction(object):
	"""Simple container for keeping track of data or function (unbound method)
		- for simplicity, self.function contains the data or function
		- e.g. self.function = numpy.pi -or- self.function = tanh"""
	def __init__(self, string, latex, dataOrFunction):
		self.string, self.latex, self.function = string, latex, dataOrFunction


class DataMethodFactory(dict):
	"""Main class with static methods for getting data and functions"""
	def __init__(self):
		# the data functions
		self.data = {}
		self.data['e'] = ("e", r'e', MATH.e)
		self.data['pi'] = ("pi", r'\pi', MATH.pi)
		self.data['1'] = ("1", r'1', 1.0)
		self.data['1/2'] = ("1/2", r'\frac{1}{2}', 0.5)
		self.data['1/N'] = ("1/N", r'\frac{1}{N}', 1./20.)
		self.data['p_i'] = ("p_i", r'\rho_{i}', -1.0)
		self.data['p_j'] = ("p_j", r'\rho_{j}', -1.0)
		self.data['p_ij'] = ("p_ij", r'\rho_{ij}', -1.0)
		self.DATA = len(self.data)
		
		# the unary functions
		self.unary = {}
		self.unary['sin'] = ("sin(%s)", r'\sin\left(%s\right)', MATH.sin)
		self.unary['exp'] = ("exp(%s)", r'\exp\left(%s\right)', MATH.exp)
		self.unary['log'] = ("log(%s)", r'\log\left(%s\right)', MATH.log)
		self.unary['tanh'] = ("tanh(%s)",r'\tanh\left(%s\right)', MATH.tanh)
		self.unary['transpose'] = ("(%s)^T",r'%s^T',MATH.transpose)
		self.UNARY = len(self.unary)

		# the binary functions
		self.binary = {}
		self.binary['+'] = ("(%s+%s)", r'\left(%s+%s\right)', MATH.add)
		self.binary['-'] = ("(%s-%s)", r'\left(%s-%s\right)', MATH.subtract)
		self.binary['*'] = ("(%s*%s)", r'%s\dot%s', MATH.multiply)
		self.binary['/'] = ("(%s/%s)", r'\frac{%s}{%s}', MATH.divide)
		self.BINARY = len(self.binary)
		
		# the scalarizing functions
		self.scalars = {}
		self.scalars['tr'] = ("tr(%s)", r'{\mathrm Tr}\left\{%s\right\}', nantrace)
		self.scalars['sum_ij'] = ("sum_ij(%s)", r'\Sigma_{ij}\left(%s\right)', dsum)
		self.SCALARS = len(self.scalars)

		# the scalarizing functions 
		# MOVE THESE TO A UTILITY CLASS
		def nantrace(self, x):
			# nan compatible trace
			try:
				y = MATH.nansum(x.diagonal())
			except:
				y = x
			return y

		def dsum(self, x):
			# nan compatible sub of the elements of a matrix
			try:
				y = MATH.nansum(x)
			except IndexError:
				y = x
			return y

		@staticmethod
		def random(dictionary, size):
			return dictionary[dictionary.keys()[randint(0, size)]]

		@staticmethod
		def randomDF(name, dictionary, size):
			assert name is None or name in dictionary
			if name is None:
				a, b, c = dictionary[name]
				return DataFunction(a, b, c)
			a, b, c = DataMethodFactory.random(dictionary, size)
			return DataFunction(a, b, c)
					
		def getData(self, name=None):
			"""Method to return a named data element or if None a random data element"""
			return randomDF(name, self.data, self.DATA)
			
		def getUnary(self, name=None):
			"""Method to return a named unary operator or if None a random unary operator"""
			return randomDF(name, self.unary, self.UNARY)
			
		def getBinary(self, name=None):
			"""Method to return a named binary operator or if None a random binary operator"""
			return randomDF(name, self.binary, self.BINARY)
		 
		def getScalar(self, name=None):
			"""Method to return a named scalarizing operator or if None a random scalarizing operator"""
			return randomDF(name, self.scalars, self.SCALARS)


#class Data(object):
#	"""Simple container for keeping information about a data source"""
#	def __init__(self, string, latex, data):
#		self.string = string
#		self.latex = latex
#		# note : self.function just return data object for bottom of recursion for consistency
#		self.function = data
		

#class Function(object):
#	"""Simple container for keeping information about a function"""
#	def __init__(self, string, latex, method):
#		self.string = string
#		self.latex = latex # use raw strings for string with escape characters
#		self.function = method
		

#class Functions(dict):
#	"""General function class to subclass"""
#	def returnRandom(self):
#		return self[self.keys()[randint(0,len(self))]]
		

#class ImmutableData(Functions):
#	"""Simple class for immutable data items; this is where we would store problem-specific 
#	loaded-in data, or constants like e, pi, etc.  This class can't deal with adding a random value; 
#	that still has to be done elsewhere.  Protein data is also stored here, but not (initially)
#	with real values (just dummies).  An adaptor (CGASimulation) has to plug in real values from a 
#	database and compute the N(N-1)/2 values needed to determine accuracy."""	
#	def __init__(self):
#		self['e'] = Data("e", r'e', MATH.e)
#		self['pi'] = Data("pi", r'\pi', MATH.pi)
#		self['1'] = Data("1", r'1', 1.0)
#		self['1/2'] = Data("1/2", r'\frac{1}{2}', 0.5)
#		self['1/N'] = Data("1/N", r'\frac{1}{N}', 1./20.)
#		self['p_i'] = Data("p_i", r'\rho_{i}', -1.0)
#		self['p_j'] = Data("p_j", r'\rho_{j}', -1.0)
#		self['p_ij'] = Data("p_ij", r'\rho_{ij}', -1.0)
		

#class ScalarizingFunctions(Functions):
#	"""This is a special set of functions that convert arrays into scalars; the final step in almost all
#	scoring methods will be double summation over amino acids, sine the input data is in the form
#	of matrices and all other operators preserve this shape.  Therefore, for the coevolving
#	residue prediction problem, the root of the tree *must* be a scalarizing operator.  However, the 
#	scalarizing functions are written to work on scalars as well (just returning the scalar), in order
#	to run simple tests on non-matrix input data. """
#	def __init__(self):	
#		self['tr'] = Function("tr(%s)", r'{\mathrm Tr}\left\{%s\right\}', self.nantrace)
#		self['sum_ij'] = Function("sum_ij(%s)", r'\Sigma_{ij}\left(%s\right)', self.dsum)
#
#		
#	def nantrace(self,x):
#		# nan compatible trace
#		try:
#			y = MATH.nansum(x.diagonal())
#		except:
#			y = x
#		return y
#
#	
#	def dsum(self,x):
#		# nan compatible sub of the elements of a matrix
#		try:
#			y = MATH.nansum(x)
#		except IndexError:
#			y = x
#		return y
	

#class UnaryFunctions(Functions):
#	"""Simple class for unary functions; these all operate array-wise, returning an array if used on
#	an array."""	
#	def __init__(self):
#		self['sin'] = Function("sin(%s)", r'\sin\left(%s\right)', MATH.sin)
#		self['exp'] = Function("exp(%s)", r'\exp\left(%s\right)', MATH.exp)
#		self['log'] = Function("log(%s)", r'\log\left(%s\right)', MATH.log)
#		self['tanh'] = Function("tanh(%s)",r'\tanh\left(%s\right)', MATH.tanh)
#		self['transpose'] = Function("(%s)^T",r'%s^T',MATH.transpose)

						
#class BinaryFunctions(Functions):
#	"""Simple class for binary functions; these also operate element-wise on arrays."""
#	def __init__(self):
#		self['+'] = Function("(%s+%s)", r'\left(%s+%s\right)', MATH.add)
#		self['-'] = Function("(%s-%s)", r'\left(%s-%s\right)', MATH.subtract)
#		self['*'] = Function("(%s*%s)", r'%s\dot%s', MATH.multiply)
#		self['/'] = Function("(%s/%s)", r'\frac{%s}{%s}', MATH.divide)
#		

class CGAFunctionsTests(unittest.TestCase):
	def setUp(self):
		self.methodFactory = DataMethodFactory()
#		self.binaryFunctions = BinaryFunctions()
#		self.unaryFunctions = UnaryFunctions()
#		self.data = ImmutableData()
#		self.scalFunctions = ScalarizingFunctions()
	
	def testName(self):
		print "\n----- testing function names -----"
		self.assertEquals(self.methodFactory.getData('pi').string, 'pi')
#		self.assertEquals(self.binaryFunctions['+'].string, '(%s+%s)')
#		self.assertEquals(self.unaryFunctions['log'].string,'log(%s)')
#		self.assertEquals(self.scalFunctions['tr'].string, 'tr(%s)')
#		self.assertEquals(self.data['pi'].string,'pi')
#		self.assertEquals(self.data['p_ij'].string,'p_ij')
		
		
if __name__ == '__main__':
	unittest.main()