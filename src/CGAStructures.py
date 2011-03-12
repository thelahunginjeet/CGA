#!/usr/bin/env python

# TODO : 
#	- use the __call__ method for memoization for all of the evaluation method calls
#	- use abstract base classes so that methods are implemented across subclasses (abc module)
#	- ZeroDivisionErrors need to be checked during the function evaluations


# NOTES :
#	- garbage collection might be a problem with the trees; I'm not sure there's any good way around
#	  	this (after reading a bit about del and such on forums).  Our best bet is probably shown in 
#	  	CGAGenerator, where we can cut external references by setting the parent of a subtree we 
#		want to delete to None.  This way we still have circular references, but they are contained.

import unittest
import CGAFunctions
from networkx import Graph, draw
import pylab

class Node(object):
	"""General node object to subclass."""
	def __init__(self,function):
		self.latex = function.latex
		self.string = function.string
		self.function = function.function
		self.nxstring = function.string+'_'+id(self).__repr__()
		# a way to deal with whether I'm the root or not; all nodes are not root by default
		self.header = False
		
	def __repr__(self):
		"""String representation of a node used for debugging."""
		return str(type(self)).split('.')[1].split("'")[0] + " : " + self.string
		
	def _evalNodes(self):
		return [self]
		
	def _evalString(self):
		return self.string
	
	def _evalFunction(self):
		return None
		
	def _evalLatex(self):
		return self.latex
		
	def _evalEdges(self):
		return []
		
	def getChildren(self):
		"""Return the left and right children of a node"""
		return (None, None)
	
	def setChildren(self,left=None,right=None):
		"""This needs to be overloaded if your node can actually have children."""
		pass
	
	def getIdentity(self):
		"""Return whether or not the node is a left (0) or right (1) node"""
		if hasattr(self, 'identity'):
			return self.identity
		else:
			return None
	
	def getHeader(self):
		"""Return whether or not the node is the root (True) or not (False)"""
		if hasattr(self,'header'):
			return self.header
		else:
			return None
		
	def setIdentity(self, integer):
		"""Set the node as being a left node (0) or a right node (1)"""
		if integer == 0:
			self.identity = 0
		elif integer == 1:
			self.identity = 1
		else:
			raise TypeError, "you have attempted to set a Node identity as an object different than 0 or 1; think again . . ."
	
	def setHeader(self, bool):
		"""Set the node as being the root (True) or not (False)"""
		assert bool in (True, False)
		self.header = bool

	
class EmptyNode(Node):
	"""Node containing nothing that is subclassed to Data, Unary and Binary Nodes."""
	def __init__(self):		
		self.latex = 'empty'
		self.string = 'empty'
		self.function = None
		self.nxstring = None
		

class DataNode(Node):
	"""General data or constant node that is terminal"""
	def __init__(self, data):
		super(DataNode, self).__init__(data)
		
	def _evalFunction(self):
		"""Return the data object so that it can be evaluated"""
		return self.function
	
	def replaceData(self, realData):
		"""Replaces the default data with real bonafide data"""
		self.function = realData
	
class ScalarNode(Node):
	"""General node for scalarizing function operations; this class is entirely redundant (unfortunately) with UnaryNode,
	but we need it to ensure we don't swap the root node for a non-scalarizing node and lose evaluability (so we can
	type check and only replace a ScalarNode with another one)."""
	def __init__(self, function):
		super(ScalarNode, self).__init__(function)
		self.left = EmptyNode()
		self.left.parent = self
		self.left.setIdentity(0)
		
	def _evalNodes(self):
		return [self] + self.left._evalNodes()
		
	def _evalFunction(self):
		leval = self.left._evalFunction()
		if leval is None:
			return None
		else:
			return self.function(leval)
		
	def _evalString(self):
		return self.string % (self.left._evalString())
		
	def _evalLatex(self):
		return self.latex % (self.left._evalLatex())
		
	def _evalEdges(self):
		return [(self.nxstring,self.left.nxstring)] + self.left._evalEdges()
		
	def setChildren(self, left=None, right=None):
		if left is not None:
			if isinstance(left, Node):
				self.left = left
				left.parent = self
				left.setIdentity(0)
			else:
				raise TypeError, "you have foolishly tried to set a Unary left child as a non-Node object; think again . . ."

	def getChildren(self):
		return (self.left, None)
	

class UnaryNode(Node):
	"""General node for unary function operations"""
	def __init__(self, function):
		super(UnaryNode, self).__init__(function)
		self.left = EmptyNode()
		self.left.parent = self
		self.left.setIdentity(0)
		
	def _evalNodes(self):
		return [self] + self.left._evalNodes()
		
	def _evalFunction(self):
		leval = self.left._evalFunction()
		if leval is None:
			return None
		else:
			return self.function(leval)
		
	def _evalString(self):
		return self.string % (self.left._evalString())
		
	def _evalLatex(self):
		return self.latex % (self.left._evalLatex())
		
	def _evalEdges(self):
		return [(self.nxstring,self.left.nxstring)] + self.left._evalEdges()
		
	def setChildren(self, left=None, right=None):
		if left is not None:
			if isinstance(left, Node):
				self.left = left
				left.parent = self
				left.setIdentity(0)
			else:
				raise TypeError, "you have foolishly tried to set a Unary left child as a non-Node object; think again . . ."

	def getChildren(self):
		return (self.left, None)
				

class BinaryNode(Node):
	"""General node for binary function operations"""
	def __init__(self, function):
		super(BinaryNode, self).__init__(function)
		self.left = EmptyNode()
		self.left.parent = self
		self.left.setIdentity(0)
		self.right = EmptyNode()
		self.right.parent = self
		self.right.setIdentity(1)
		
	def _evalNodes(self):
		return [self] + self.left._evalNodes() + self.right._evalNodes()

	def _evalFunction(self):
		leval = self.left._evalFunction()
		reval = self.right._evalFunction()
		if leval is None or reval is None:
			return None
		else:
			return self.function(leval, reval)

	def _evalString(self):
		return self.string % (self.left._evalString(), self.right._evalString())

	def _evalLatex(self):
		return self.latex % (self.left._evalLatex(), self.right._evalLatex())
		
	def _evalEdges(self):
		return [(self.nxstring,self.left.nxstring),(self.nxstring,self.right.nxstring)] + self.left._evalEdges() + self.right._evalEdges()
				
	def setChildren(self, left=None, right=None):
		if left is not None:
			if isinstance(left, Node):				
				self.left = left
				left.parent = self
				left.setIdentity(0)
			else:
				raise TypeError, "you have foolishly tried to set a Binary left child as a non-Node object; think again . . ."		
		if right is not None:
			if isinstance(right, Node):
				self.right = right
				right.parent = self
				right.setIdentity(1)
			else:
				raise TypeError, "you have foolishly tried to set a Binary right child as a non-Node object; think again . . ."

	def getChildren(self):
		return (self.left, self.right)
			
	
class AlgorithmTree(object):
	"""General tree structure for recursion"""
	def __init__(self, root):
		self.root = root
		self.root.setHeader(True)
		self.root.setIdentity(0)
		self.graph = Graph()
		self.update()
	
	def __call__(self):
		"""Make the tree callable to evaluate expressions"""
		self.update()
		self.evaluateFunction()
		self.evaluateString()
		self.evaluateLatex()
		self.evaluateGraph()
	
	def __repr__(self):
		"""String representation a tree."""
		self()	# make sure the tree is evaluated first
		output = "function eval : %s\nstring eval : %s\nLaTeX eval : %s\nEdges eval : %s" \
			%(self.function, self.string, self.latex, self.graph.edges())
		return output
	
	def update(self):
		"""General function to update any properties"""
		self.evaluateNodes()
	
	def evaluateNodes(self):
		"""Recurse the tree and keep track of all nodes"""
		self.nodes = self.root._evalNodes()
		self.termini = [x for x in self.nodes if isinstance(x, DataNode) or isinstance(x, EmptyNode)]
		
	def evaluateFunction(self):
		"""Recurse the tree and evaluate the function; if there are still empty nodes in
		the tree, this will be meaningless."""
		self.function = self.root._evalFunction()			
		
	def evaluateString(self):
		"""Recurse the tree and evaluate the string expression"""
		self.string = self.root._evalString()
	
	def evaluateLatex(self):
		"""Recurse the tree and evaluate the LaTeX expression"""
		self.latex = self.root._evalLatex()
		
	def evaluateGraph(self):
		"""Clears the current graph, recurses the tree to get the edgelist, and then 
		returns the graph composed of those edges."""
		self.graph.clear()
		self.graph.add_edges_from(self.root._evalEdges())
						
	
class AlgorithmTreeTests(unittest.TestCase):
	"""Test suite for making AlgorithmTree operations"""	
	def setUp(self):
		self.binaryFunctions = CGAFunctions.BinaryFunctions()
		self.unaryFunctions = CGAFunctions.UnaryFunctions()
		self.data = CGAFunctions.ImmutableData()
	
	def testNxNodeNaming(self):
		print "\n----- testing unique node names for nx.graph -----"
		root = BinaryNode(self.binaryFunctions['/'])
		print 'Root node namestring : %s' % root.nxstring
		node1 = UnaryNode(self.unaryFunctions['log'])
		print 'Node1 namestring : %s' % node1.nxstring
		node2 = UnaryNode(self.unaryFunctions['log'])
		print 'Node2 namestring : %s' % node2.nxstring
		self.assertNotEquals(node1.nxstring,node2.nxstring)
				
	def testRecursion(self): 
		print "\n----- testing tree recursion -----"
		root = BinaryNode(self.binaryFunctions['/'])
		node1 = UnaryNode(self.unaryFunctions['exp'])
		node2 = UnaryNode(self.unaryFunctions['log'])
		node3 = UnaryNode(self.unaryFunctions['sin'])
		node4 = UnaryNode(self.unaryFunctions['log'])
		# immutable data - stored
		constant1 = DataNode(self.data['pi'])
		# mutable data - might require one function call then the value is fixed
		value = CGAFunctions.uniform()
		constant2 = DataNode(CGAFunctions.Data(str(value),str(value),value))
		tree = AlgorithmTree(root)
		root.setChildren(node1, node2)
		node1.setChildren(node3)
		node3.setChildren(node4)
		node2.setChildren(constant1)
		node4.setChildren(constant2)
		tree()
		print tree
		
	def testGraph(self):
		print "\n----- testing graph drawing -----"
		root = BinaryNode(self.binaryFunctions['/'])
		node1 = UnaryNode(self.unaryFunctions['exp'])
		node2 = UnaryNode(self.unaryFunctions['log'])
		node3 = UnaryNode(self.unaryFunctions['sin'])
		node4 = UnaryNode(self.unaryFunctions['log'])
		# immutable data - stored
		constant1 = DataNode(self.data['pi'])
		# mutable data - might require one function call then the value is fixed
		value = CGAFunctions.uniform()
		constant2 = DataNode(CGAFunctions.Data(str(value),str(value),value))
		tree = AlgorithmTree(root)
		root.setChildren(node1, node2)
		node1.setChildren(node3)
		node3.setChildren(node4)
		node2.setChildren(constant1)
		node4.setChildren(constant2)
		tree()
		draw(tree.graph)
		pylab.show()
		
if __name__ == '__main__':
	unittest.main()