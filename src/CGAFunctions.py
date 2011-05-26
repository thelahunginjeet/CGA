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
		# the data functions (with random number generator)
		self.data = {}
		self.data['e'] = ("e", r' e', MATH.e)
		self.data['pi'] = ("pi", r'$pi', MATH.pi)
		# ephemeral random number
		rnum = randint(-1,2)*uniform() + randint(-3, 4) 
		self.data[str(rnum)] = (str(rnum), r' %s'%str(rnum), rnum)
		self.data['1/N'] = ("(1/N)", r'$frac{1}{N}', 1./20.)
		self.data['p_i'] = ("p_i", r'$rho_{i}', -1.0)
		self.data['p_j'] = ("p_j", r'$rho_{j}', -1.0)
		self.data['p_ij'] = ("p_ij", r'$rho_{ij}', -1.0)
		self.DATA = len(self.data)
		
		# the unary functions
		self.unary = {}
		self.unary['exp'] = ("exp(%s)", r'$exp$left( %s$right) ', MATH.exp)
		self.unary['log'] = ("log(%s)", r'$log$left( %s$right) ', MATH.log)
		self.unary['tanh'] = ("tanh(%s)", r'$tanh$left( %s$right) ', MATH.tanh)
		self.unary['sinh'] = ("sinh(%s)", r'$sinh$left( %s$right) ', MATH.sinh)
		self.unary['cosh'] = ("cosh(%s)", r'$cosh$left( %s$right) ', MATH.cosh)
		self.unary['transpose'] = ("(%s)^T",r'$left( %s$right)^{T} ',MATH.transpose)
		self.unary['square'] = ("(%s)**2", r'$left( %s$right)^{2} ', self.sqr)
		self.UNARY = len(self.unary)

		# the binary functions
		self.binary = {}
		self.binary['add'] = ("(%s+%s)", r'$left( %s + %s$right) ', MATH.add)
		self.binary['subtract'] = ("(%s-%s)", r'$left( %s- %s$right) ', MATH.subtract)
		self.binary['multiply'] = ("(%s*%s)", r'$left(%s$cdot %s$right) ', MATH.multiply)
		self.binary['divide'] = ("(%s/%s)", r'$frac{%s}{%s}', MATH.divide)
		self.BINARY = len(self.binary)
		
		# the scalarizing functions
		self.scalars = {}
		self.scalars['tr'] = ("tr(%s)", r' {$mathrm{Tr}}$left( %s$right) ', self.nantrace)
		self.scalars['sum_ij'] = ("sum_ij(%s)", r'$Sigma_{ij}$left( %s$right) ', self.dsum)
		self.SCALARS = len(self.scalars)
		
		# reverse dicts (wait for it . . . ah) for accessing by .string
		self.atad = dict((self.data[k][0],k) for k in self.data)
		self.yranu = dict((self.unary[k][0],k) for k in self.unary)
		self.yranib = dict((self.binary[k][0],k) for k in self.binary)
		self.sralacs = dict((self.scalars[k][0],k) for k in self.scalars)

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
		"""Method to return a named data element, if None a random data element, or a copy of an 
			ephemeral random number (either a float or int)"""
		if name in self.data or name is None:
			return DataMethodFactory.randomDF(name, self.data, self.DATA)
		elif name in self.atad:
			return DataMethodFactory.randomDF(self.atad[name], self.data, self.DATA)
		else:
			# asking for ephemeral random number (either float or int)			
			if '.' in name:
				return DataFunction(name, name, float(name))
			return DataFunction(name, name, int(name))
	
	def getUnary(self, name=None):
		"""Method to return a named unary operator or if None a random unary operator"""
		if name in self.unary or name is None:
			return DataMethodFactory.randomDF(name, self.unary, self.UNARY)
		return DataMethodFactory.randomDF(self.yranu[name], self.unary, self.UNARY)
	
	def getBinary(self, name=None):
		"""Method to return a named binary operator or if None a random binary operator"""
		if name in self.binary or name is None:
			return DataMethodFactory.randomDF(name, self.binary, self.BINARY)
		return DataMethodFactory.randomDF(self.yranib[name], self.binary, self.BINARY)

	def getScalar(self, name=None):
		"""Method to return a named scalarizing operator or if None a random scalarizing operator"""
		if self.scalars.has_key(name) or name is None:
			return DataMethodFactory.randomDF(name, self.scalars, self.SCALARS)
		return DataMethodFactory.randomDF(self.sralacs[name], self.scalars, self.SCALARS)
	

	@staticmethod
	def nantrace(x):
		# nan compatible trace
		try:
			y = MATH.nansum(x.diagonal())
		except:
			y = x
		return y

	@staticmethod
	def dsum(x):
		# nan compatible sum of the elements of a matrix
		try:
			y = MATH.nansum(x)
		except IndexError:
			y = x
		return y
	
	@ staticmethod
	def sqr(x):
		return x**2


class CGAFunctionsTests(unittest.TestCase):
	def setUp(self):
		self.methodFactory = DataMethodFactory()
	
	def testName(self):
		print "\n----- testing function names -----"
		self.assertEquals(self.methodFactory.getData('pi').string, 'pi')
		self.assertEquals(self.methodFactory.getScalar('tr').string, 'tr(%s)')
		self.assertEquals(self.methodFactory.getUnary('log').string, 'log(%s)')
		self.assertEquals(self.methodFactory.getBinary('add').string, '(%s+%s)')
		
	def testFunctions(self):
		print "\n----- testing function evaluation -----"
		mynansum = self.methodFactory.getScalar('sum_ij')
		mynantrace = self.methodFactory.getScalar('tr')
		mylog = self.methodFactory.getUnary('log')
		mysum = self.methodFactory.getBinary('add')
		self.assertAlmostEquals(0.0,mylog.function(1.0))
		self.assertAlmostEquals(2.0,mysum.function(1.0,1.0))
		self.assertAlmostEquals(2.0,mynansum.function([1.0,MATH.nan,1.0,MATH.nan]))
		self.assertAlmostEquals(2.0,mynansum.function(MATH.eye(2)))
	
	def testReverseAccess(self):
		print "\n----- testing reverse dictionary access -----"
		N = 10
		for i in range(N):
			binary = self.methodFactory.getBinary()
			binary2 = self.methodFactory.getBinary(binary.string)
			if binary.string == binary2.string:
				self.assertEquals((binary.string, binary.latex, binary.function), (binary.string, binary.latex, binary.function))
			else:	
				self.assertNotEquals((binary.string, binary.latex, binary.function), (binary.string, binary.latex, binary.function))
			self.assertNotEquals(binary, binary2)
		for i in range(N):
			unary = self.methodFactory.getUnary()
			unary2 = self.methodFactory.getUnary(unary.string)
			if unary.string == unary2.string:
				self.assertEquals((unary.string, unary.latex, unary.function), (unary.string, unary.latex, unary.function))
			else:
				self.assertNotEquals((unary.string, unary.latex, unary.function), (unary.string, unary.latex, unary.function))
			self.assertNotEquals(unary, unary2)
		for i in range(N):
			data = self.methodFactory.getData()
			data2 = self.methodFactory.getData(data.string)
			if data.string == data2.string:
				self.assertEquals((data.string, data.latex, data.function), (data.string, data.latex, data.function))
			else:				
				self.assertNotEquals((data.string, data.latex, data.function), (data.string, data.latex, data.function))	
			self.assertNotEquals(data, data2)
		for i in range(N):
			scalar = self.methodFactory.getScalar()
			scalar2 = self.methodFactory.getScalar(scalar.string)
			if scalar.string == scalar2.string:
				self.assertEquals((scalar.string, scalar.latex, scalar.function), (scalar.string, scalar.latex, scalar.function))
			else:
				self.assertNotEquals((scalar.string, scalar.latex, scalar.function), (scalar.string, scalar.latex, scalar.function))
			self.assertNotEquals(scalar, scalar2)

	
if __name__ == '__main__':
	unittest.main()
