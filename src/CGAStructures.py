#!/usr/bin/env python

import unittest
import CGAFunctions
from networkx import Graph, draw
import pylab


class Node(object):
	"""General node object to subclass that acts almost as an abstract class / interface"""
	def __init__(self, function):
		assert type(function) in (CGAFunctions.Data, CGAFunctions.Function)
		self.latex, self.string, self.function = function.latex, function.string, function.function
		self.header = False
		self.nxstring = function.string+'_'+id(self).__repr__()
		
	def __repr__(self):
		"""String representation of a node used for debugging."""
		return str(type(self)).split('.')[1].split("'")[0] + " : " + self.string
		
	def _evalNodes(self):
		raise Exception, "base class _evalNodes() should never be called"
		
	def _evalString(self):
		return self.string
	
	def _evalFunction(self):
		return self.function
		
	def _evalLatex(self):
		return self.latex
		
	def _evalEdges(self):
		return []
	
	def clean(self):
		self = None
		return self
		
	def getChildren(self):
		"""Return the left and right children of a node"""
		return (None, None)
	
	def setChildren(self,left=None,right=None):
		"""This needs to be overloaded if your node can actually have children."""
		pass
	
	def getIdentity(self):
		"""Return whether or not the node is a left (0) or right (1) node"""
		assert hasattr(self, 'identity')
		return self.identity

	def setIdentity(self, integer):
		"""Set the node as being a left node (0) or a right node (1)"""
		assert integer in (0, 1)
		self.identity = integer
		
	def getHeader(self):
		"""Return whether or not the node is the root (True) or not (False)"""
		assert hasattr(self, 'header')
		return self.header
	
	def setHeader(self, bool):
		"""Set the node as being the root (True) or not (False)"""
		assert bool in (True, False)
		self.header = bool
		

class DataNode(Node):
	"""General data or constant node that is terminal"""
	def __init__(self, data=None):
		assert data is None or type(data) is CGAFunctions.Data
		if data is None:
			super(DataNode, self).__init__(CGAFunctions.ImmutableData().returnRandom())
		else:
			super(DataNode, self).__init__(data)
			
	def _evalNodes(self):
		return [self]
	
	def replaceData(self, realData):
		"""Replaces the default data with real bonafide data"""
		assert realData is not None
		self.function = realData

	
class ScalarNode(Node):
	"""General node for scalarizing function operations; this class is entirely redundant (unfortunately) with UnaryNode,
	but we need it to ensure we don't swap the root node for a non-scalarizing node and lose evaluability (so we can
	type check and only replace a ScalarNode with another one)."""
	def __init__(self, function):
		super(ScalarNode, self).__init__(function)
		self.left = DataNode()
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
		self.left = DataNode()
		self.left.parent = self
		self.left.setIdentity(0)
	
	def clean(self):
		self.left = self.left.clean()
		return self.left
	
	def _evalNodes(self):
		return [self] + self.left._evalNodes()
		
	def _evalFunction(self):
		leval = self.left._evalFunction()
		assert leval is not None
		return self.function(leval)
#		if leval is None:
#			return None
#		else:
#			return self.function(leval)
		
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
		self.left = DataNode()
#		self.left = EmptyNode()
		self.left.parent = self
		self.left.setIdentity(0)
		self.right = DataNode()
#		self.right = EmptyNode()
		self.right.parent = self
		self.right.setIdentity(1)
		
	def clean(self):
		self.left = self.left.clean()
		self.right = self.right.clean()
		return self.right  # just need to return None, so either left or right
		
	def _evalNodes(self):
		return [self] + self.left._evalNodes() + self.right._evalNodes()

	def _evalFunction(self):
		leval = self.left._evalFunction()
		reval = self.right._evalFunction()
		assert leval is not None and reval is not None
		return self.function(leval, reval)
#		if leval is None or reval is None:
#			return None
#		else:
#			return self.function(leval, reval)

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
	
	# might need to simply have getters to return things like the terminii so that update isn't always called
	def evaluateNodes(self):
		"""Recurse the tree and return a list of all of the nodes in the tree"""
		self.nodes = self.root._evalNodes()
		self.termini = [x for x in self.nodes if isinstance(x, DataNode)]
		
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
		
	def testClean(self):
		print "\n----- testing delete on graph -----"
		root = BinaryNode(self.binaryFunctions['/'])
		node1 = UnaryNode(self.unaryFunctions['exp'])
		node2 = UnaryNode(self.unaryFunctions['log'])
		node3 = UnaryNode(self.unaryFunctions['sin'])
		node4 = UnaryNode(self.unaryFunctions['log'])
		constant1 = DataNode(self.data['pi'])
		value = CGAFunctions.uniform()
		constant2 = DataNode(CGAFunctions.Data(str(value),str(value),value))
		tree = AlgorithmTree(root)
		root.setChildren(node1, node2)
		node1.setChildren(node3)
		node3.setChildren(node4)
		node2.setChildren(constant1)
		node4.setChildren(constant2)
		tree()
		root.setChildren(DataNode(), None)
		node1.clean()
		tree()
		draw(tree.graph)
		pylab.show()

		
if __name__ == '__main__':
	unittest.main()