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
		self.scalars['tr'] = ("tr(%s)", r'{\mathrm Tr}\left\{%s\right\}', DataMethodFactory.nantrace)
		self.scalars['sum_ij'] = ("sum_ij(%s)", r'\Sigma_{ij}\left(%s\right)', DataMethodFactory.dsum)
		self.SCALARS = len(self.scalars)
		
		# reverse dictionaries - it's one-to-one, so no problems; the reverse dictionaries allow
		#	you to use node.string to access members.  This is important for copying nodes.
		# TODO - these don't work for mutable data
		self.rev_data = dict((self.data[k][0],k) for k in self.data)
		self.rev_unary = dict((self.unary[k][0],k) for k in self.unary)
		self.rev_binary = dict((self.binary[k][0],k) for k in self.binary)
		self.rev_scalars = dict((self.scalars[k][0],k) for k in self.scalars)

	@staticmethod
	def random(dictionary, size):
		return dictionary[dictionary.keys()[randint(0, size)]]

	@staticmethod
	def randomDF(name, dictionary, size):
		assert name is None or name in dictionary
		if name in dictionary:
			a, b, c = dictionary[name]
			return DataFunction(a, b, c)
		a, b, c = DataMethodFactory.random(dictionary, size)
		return DataFunction(a, b, c)
		
	def getData(self, name=None):
		"""Method to return a named data element or if None a random data element"""
		if self.data.has_key(name) or name is None:
			# regular forward access
			return DataMethodFactory.randomDF(name, self.data, self.DATA)
		else:
			# reverse access
			return DataMethodFactory.randomDF(self.rev_data[name], self.data, self.DATA)
	
	def getUnary(self, name=None):
		"""Method to return a named unary operator or if None a random unary operator"""
		if self.unary.has_key(name) or name is None:
			return DataMethodFactory.randomDF(name, self.unary, self.UNARY)
		else:
			return DataMethodFactory.randomDF(self.rev_unary[name], self.unary, self.UNARY)
	
	def getBinary(self, name=None):
		"""Method to return a named binary operator or if None a random binary operator"""
		if self.binary.has_key(name) or name is None:
			return DataMethodFactory.randomDF(name, self.binary, self.BINARY)
		else:
			return DataMethodFactory.randomDF(self.rev_binary[name], self.binary, self.BINARY)

	def getScalar(self, name=None):
		"""Method to return a named scalarizing operator or if None a random scalarizing operator"""
		if self.scalars.has_key(name) or name is None:
			return DataMethodFactory.randomDF(name, self.scalars, self.SCALARS)
		else:
			return DataMethodFactory.randomDF(self.rev_scalars[name], self.scalars, self.SCALARS)

	@staticmethod
	def nantrace(self, x):
		# nan compatible trace
		try:
			y = MATH.nansum(x.diagonal())
		except:
			y = x
		return y

	@staticmethod
	def dsum(self, x):
		# nan compatible sub of the elements of a matrix
		try:
			y = MATH.nansum(x)
		except IndexError:
			y = x
		return y


class CGAFunctionsTests(unittest.TestCase):
	def setUp(self):
		self.methodFactory = DataMethodFactory()
	
	def testName(self):
		print "\n----- testing function names -----"
		self.assertEquals(self.methodFactory.getData('pi').string, 'pi')
		self.assertEquals(self.methodFactory.getScalar('tr').string, 'tr(%s)')
		self.assertEquals(self.methodFactory.getUnary('log').string, 'log(%s)')
		self.assertEquals(self.methodFactory.getBinary('+').string, '(%s+%s)')
	
	def testReverseAccess(self):
		print "\n----- testing reverse dictionary access -----"
		scalar = self.methodFactory.getScalar('tr')
		unary = self.methodFactory.getUnary('log')
		binary = self.methodFactory.getBinary('+')
		self.assertEquals(scalar.string, self.methodFactory.getScalar(scalar.string).string)
		self.assertEquals(unary.string, self.methodFactory.getUnary(unary.string).string)
		self.assertEquals(binary.string, self.methodFactory.getBinary(binary.string).string)
	
if __name__ == '__main__':
	unittest.main()