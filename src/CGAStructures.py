#!/usr/bin/env python

import unittest
import CGAFunctions
from networkx import Graph, draw
import pylab


class Node(object):
	"""General node object to subclass that acts almost as an abstract class / interface"""
	def __init__(self, function):
		assert type(function) is CGAFunctions.DataFunction
		self.latex, self.string, self.function = function.latex, function.string, function.function
		self.header = False
		self.nxstring = function.string + '_' + str(id(self))
		
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
	
	def copy(self):
		"""Returns a new instance of a node identical to the calling node; must be
		overridden in base classes."""
		raise Exception, "base class copy() should never be called"
	
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
		if data is None:
			super(DataNode, self).__init__(CGAFunctions.DataMethodFactory().getData())
		else:
			super(DataNode, self).__init__(data)
			
	def _evalNodes(self):
		return [self]
	
	def replaceData(self, realData):
		"""Replaces the default data with real bonafide data"""
		assert realData is not None
		self.function = realData
		
	def copy(self):
		"""Copy method for a data node is simple, as data node are always terminal."""
		return DataNode(CGAFunctions.DataMethodFactory().getData(self.string))

	
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
	
	def copy(self):
		"""Returns a copy (new instance) of the node.  Also copies the node's children."""
		newScalar = ScalarNode(CGAFunctions.DataMethodFactory().getScalar(self.string))
		newScalar.setChildren(self.left.copy())
		return newScalar
		
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
	
	def copy(self):
		"""Returns a copy (new instance) of the node.  Also copies the node's children."""
		newUnary = UnaryNode(CGAFunctions.DataMethodFactory().getUnary(self.string))
		newUnary.setChildren(self.left.copy())
		return newUnary
	
	def _evalNodes(self):
		return [self] + self.left._evalNodes()
		
	def _evalFunction(self):
		leval = self.left._evalFunction()
		assert leval is not None
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
		self.left = DataNode()
		self.left.parent = self
		self.left.setIdentity(0)
		self.right = DataNode()
		self.right.parent = self
		self.right.setIdentity(1)
		
	def clean(self):
		self.left = self.left.clean()
		self.right = self.right.clean()
		return self.right  # just need to return None, so either left or right
	
	def copy(self):
		"""Returns a copy (new instance) of the node. Also copies the node's children."""
		newBinary = BinaryNode(CGAFunctions.DataMethodFactory().getBinary(self.string))
		newBinary.setChildren(self.left.copy(),self.right.copy())
		return newBinary
	
	def _evalNodes(self):
		return [self] + self.left._evalNodes() + self.right._evalNodes()

	def _evalFunction(self):
		leval = self.left._evalFunction()
		reval = self.right._evalFunction()
		assert leval is not None and reval is not None
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
	"""General tree structure for recursion.  Tree members (function, graph, string, latex) that are
	assembled from the constituent nodes are queried using accessor methods, which are called
	when the tree has changed. """
	def __init__(self, root):
		self.root = root
		self.root.setHeader(True)
		self.root.setIdentity(0)
		self._isCurrent = {'nodes':False,'string':False,'function':False,'latex':False,'graph':False}
		self._graph = Graph()
		self._function = None
		self._string = None
		self._latex = None

	#def _update(self):
	#	"""Updates all the tree values, and sets the current flag to true."""
	#	self._nodes = self.root._evalNodes()
	#	self._termini = [x for x in self._nodes if isinstance(x, DataNode)] 
	#	self._function = self.root._evalFunction()
	#	self._string = self.root._evalString()
	#	self._latex = self.root._evalLatex()
	#	self._graph.clear()
	#	self._graph.add_edges_from(self.root._evalEdges())
	#	self._isCurrent = True
	
	def __repr__(self):
		"""String representation of a tree."""
		output = "function eval : %s\nstring eval : %s\nLaTeX eval : %s\nEdges eval : %s" \
			%(self.getFunction(), self.getString(), self.getLatex(), self.getGraph().edges())
		return output
	
	def _evaluateNodes(self):
		"""Recurse the tree and update the nodes and the tree termini"""
		self._nodes = self.root._evalNodes()
		self._termini = [x for x in self._nodes if isinstance(x, DataNode)]
	
	def copy(self):
		"""Recursive copy of the entire tree; returns a tree."""
		newTree = AlgorithmTree(self.root.copy())
		return newTree
	
	def getGraph(self):
		"""Either returns the stored (correct) graph, or if the nodes need updating the graph is cleared,
		the tree recursed to get the edgelist, and then the graph composed of those edges is returned."""
		if self._isCurrent['graph']:
			return self._graph
		self._evaluateNodes()
		self._graph.clear()
		self._graph.add_edges_from(self.root._evalEdges())
		self._isCurrent['graph'] = True
		return self._graph
	
	def getFunction(self):
		"""If the nodes are not current, recurse the tree and evaluate the tree's function.
		Otherwise, just return the correct, stored function."""
		if self._isCurrent['function']:
			return self._function
		self._evaluateNodes()
		self._function = self.root._evalFunction()
		self._isCurrent['function'] = True
		return self._function			
	
	def getString(self):
		"""If the nodes are not current, recurse the tree and evaluate the string expression for the tree.
		Otherwise, just return the correct, stored string expression."""
		if self._isCurrent['string']:
			return self._string
		self._evaluateNodes()
		self._string = self.root._evalString()
		self._isCurrent['string'] = True
		return self._string
	
	def getLatex(self):
		"""If the nodes are not current, recurse the tree and evaluate the latex expression for the tree.
		Otherwise, just return the correct, stored latex expression."""
		if self._isCurrent['latex']:
			return self._latex
		self._evaluateNodes()
		self._latex = self.root._evalLatex()
		self._isCurrent['latex'] = True
		return self._latex
	
	def getNodes(self):
		"""If the nodes are not current, recurse the tree and return a list of all of the nodes in the tree.
		Otherwise, just return a list of the nodes."""
		if self._isCurrent['nodes']:
			return self._nodes
		self._evaluateNodes()
		self._isCurrent['nodes'] = True
		return self._nodes
	
	def getTermini(self):
		"""If the nodes are not current, recurse the tree to update the nodes are return the terminal nodes.
		Otherwise, just return the currently stored list."""
		if self._isCurrent['nodes']:
			return self._termini
		self._evaluateNodes()
		return self._termini
	
	def nodesHaveChanged(self):
		"""You are responsible for calling this if you modify the tree.  (Not sure if this can
		be made automatic or whether it is worth it)."""
		# ensure that accessors re-evaluate at next request
		for k in self._isCurrent:
			self._isCurrent[k] = False		
		
		
	#ef evaluateFunction(self):
	#	"""Recurse the tree and evaluate the function; if there are still empty nodes in
	#	the tree, this will be meaningless."""
	#	self._function = self.root._evalFunction()			
		
	#def evaluateString(self):
	#	"""Recurse the tree and evaluate the string expression"""
	#	self._string = self.root._evalString()
	
	#def evaluateLatex(self):
	#	"""Recurse the tree and evaluate the LaTeX expression"""
	#	self._latex = self.root._evalLatex()
		
	#def evaluateGraph(self):
	#	"""Clears the current graph, recurses the tree to get the edgelist, and then 
	#	returns the graph composed of those edges."""
	#	self._graph.clear()
	#	self._graph.add_edges_from(self.root._evalEdges())
						
	
class AlgorithmTreeTests(unittest.TestCase):
	"""Test suite for making AlgorithmTree operations"""	
	def setUp(self):
		self.methodFactory = CGAFunctions.DataMethodFactory()
	
	def testNxNodeNaming(self):
		print "\n----- testing unique node names for nx.graph -----"
		root = BinaryNode(self.methodFactory.getBinary('/'))
		print 'Root node namestring : %s' % root.nxstring
		node1 = UnaryNode(self.methodFactory.getUnary('log'))
		print 'Node1 namestring : %s' % node1.nxstring
		node2 = UnaryNode(self.methodFactory.getUnary('log'))
		print 'Node2 namestring : %s' % node2.nxstring
		self.assertNotEquals(node1.nxstring,node2.nxstring)
				
	def testRecursion(self): 
		print "\n----- testing tree recursion -----"
		root = BinaryNode(self.methodFactory.getBinary('/'))
		node1 = UnaryNode(self.methodFactory.getUnary('exp'))
		node2 = UnaryNode(self.methodFactory.getUnary('log'))
		node3 = UnaryNode(self.methodFactory.getUnary('sin'))
		node4 = UnaryNode(self.methodFactory.getUnary('log'))
		# immutable data - stored
		constant1 = DataNode(self.methodFactory.getData('pi'))
		# mutable data - might require one function call then the value is fixed
		value = CGAFunctions.uniform()
		constant2 = DataNode(CGAFunctions.DataFunction(str(value),str(value),value))
		tree = AlgorithmTree(root)
		root.setChildren(node1, node2)
		node1.setChildren(node3)
		node3.setChildren(node4)
		node2.setChildren(constant1)
		node4.setChildren(constant2)
		print tree
		
	def testGraph(self):
		print "\n----- testing graph drawing -----"
		root = BinaryNode(self.methodFactory.getBinary('/'))
		node1 = UnaryNode(self.methodFactory.getUnary('exp'))
		node2 = UnaryNode(self.methodFactory.getUnary('log'))
		node3 = UnaryNode(self.methodFactory.getUnary('sin'))
		node4 = UnaryNode(self.methodFactory.getUnary('log'))
		# immutable data - stored
		constant1 = DataNode(self.methodFactory.getData('pi'))
		# mutable data - might require one function call then the value is fixed
		value = CGAFunctions.uniform()
		constant2 = DataNode(CGAFunctions.DataFunction(str(value),str(value),value))
		tree = AlgorithmTree(root)
		root.setChildren(node1, node2)
		node1.setChildren(node3)
		node3.setChildren(node4)
		node2.setChildren(constant1)
		node4.setChildren(constant2)
		draw(tree.getGraph())
		pylab.show()
		
	def testClean(self):
		print "\n----- testing delete on graph -----"
		root = BinaryNode(self.methodFactory.getBinary('/'))
		node1 = UnaryNode(self.methodFactory.getUnary('exp'))
		node2 = UnaryNode(self.methodFactory.getUnary('log'))
		node3 = UnaryNode(self.methodFactory.getUnary('sin'))
		node4 = UnaryNode(self.methodFactory.getUnary('log'))
		# immutable data - stored
		constant1 = DataNode(self.methodFactory.getData('pi'))
		# mutable data - might require one function call then the value is fixed
		value = CGAFunctions.uniform()
		constant2 = DataNode(CGAFunctions.DataFunction(str(value),str(value),value))
		tree = AlgorithmTree(root)
		root.setChildren(node1, node2)
		node1.setChildren(node3)
		node3.setChildren(node4)
		node2.setChildren(constant1)
		node4.setChildren(constant2)
		print 'Tree before node cleaning:'
		print tree
		root.setChildren(DataNode(), None)
		tree.nodesHaveChanged()
		node1.clean()
		print 'Tree after node cleaning:'
		print tree
		draw(tree.getGraph())
		pylab.show()
	
	def testCopyDataNode(self):
		print "\n----- testing copy() of data node -----"
		dNode = DataNode(self.methodFactory.getData('pi'))
		# do a copy
		newDataNode = dNode.copy()
		# compare addresses
		self.assertNotEquals(id(dNode),id(newDataNode))
		# set the original to None
		dNode = None
		# original should be None but copy should not
		self.assertEquals(dNode,None)
		self.assertNotEquals(newDataNode,None)
		
	def testCopyUnaryNode(self):
		print "\n----- testing copy() of unary node -----"
		# make a unary node with a data child
		uNode = UnaryNode(self.methodFactory.getUnary('exp'))
		dNode = DataNode(self.methodFactory.getData('pi'))
		uNode.setChildren(dNode)
		# copy the unary node
		newUNode = uNode.copy()
		# check addresses
		self.assertNotEquals(id(uNode),id(newUNode))
		self.assertNotEquals(id(uNode.getChildren()[0]),id(newUNode.getChildren()[0]))
		
	def testCopyBinaryNode(self):
		print "\n----- testing copy() of binary node -----"
		# make a binary node with data children
		bNode = BinaryNode(self.methodFactory.getBinary('+'))
		dNode1 = DataNode(self.methodFactory.getData('pi'))
		dNode2 = DataNode(self.methodFactory.getData('1'))
		bNode.setChildren(dNode1,dNode2)
		# copy the node
		newBNode = bNode.copy()
		self.assertNotEquals(id(newBNode),id(bNode))
		self.assertNotEquals(id(newBNode.getChildren()[0]),id(bNode.getChildren()[0]))
		self.assertNotEquals(id(newBNode.getChildren()[1]),id(bNode.getChildren()[1]))
		
	def testRecursiveCopy(self):
		print "\n----- testing recursive copy() -----"
		# make a small subtree (no root, just linked nodes)
		uNode1 = UnaryNode(self.methodFactory.getUnary('exp'))
		uNode2 = UnaryNode(self.methodFactory.getUnary('tanh'))
		dNode = DataNode(self.methodFactory.getData('pi'))
		uNode1.setChildren(uNode2)
		uNode2.setChildren(dNode)
		# copy uNode1
		newUNode = uNode1.copy()
		print 'Original subtree: ', uNode1,uNode1.getChildren()[0],uNode1.getChildren()[0].getChildren()[0]
		print 'Copied subtree: ', newUNode,newUNode.getChildren()[0],newUNode.getChildren()[0].getChildren()[0]
		# check addresses
		self.assertNotEquals(id(uNode1),id(newUNode))
		self.assertNotEquals(id(newUNode.getChildren()[0]),id(uNode2))
		self.assertNotEquals(id(newUNode.getChildren()[0].getChildren()[0]),id(dNode))
		
	def testTreeCopy(self):
		print "\n----- testing recursive copy() of algorithm tree -----"
		root = BinaryNode(self.methodFactory.getBinary('/'))
		node1 = UnaryNode(self.methodFactory.getUnary('exp'))
		node2 = UnaryNode(self.methodFactory.getUnary('log'))
		node3 = UnaryNode(self.methodFactory.getUnary('sin'))
		node4 = UnaryNode(self.methodFactory.getUnary('log'))
		# immutable data - stored
		constant1 = DataNode(self.methodFactory.getData('pi'))
		constant2 = DataNode(self.methodFactory.getData('1'))
		tree = AlgorithmTree(root)
		root.setChildren(node1, node2)
		node1.setChildren(node3)
		node3.setChildren(node4)
		node2.setChildren(constant1)
		node4.setChildren(constant2)
		#tree()
		print "Original tree (node,id): "
		print [(n.string,id(n)) for n in tree.getNodes()]
		print "Root, root header value:", tree.root,tree.root.getHeader()
		newTree = tree.copy()
		root.clean()
		#newTree()
		print "Copy tree (node,id): "
		print [(n.string,id(n)) for n in newTree.getNodes()]
		print "Root, root header value:", newTree.root,newTree.root.getHeader()
		
	def testStateVariable(self):
		print "\n----- testing state variable for node updates  -----"
		root = ScalarNode(self.methodFactory.getScalar('tr'))
		node1 = BinaryNode(self.methodFactory.getBinary('-'))
		# immutable data - stored
		constant1 = DataNode(self.methodFactory.getData('e'))
		constant2 = DataNode(self.methodFactory.getData('pi'))
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(constant1,constant2)
		# save some old tree data
		oldString = tree.getString()
		oldFunction = tree.getFunction()
		oldLatex = tree.getLatex()
		# stick a new node in there
		newNode = BinaryNode(self.methodFactory.getBinary('*'))
		root.setChildren(newNode,None)
		newNode.setChildren(constant1,constant2)
		node1.clean()
		# status has not changed, so string/etc. are not current!
		self.assertEquals(tree.getString(),oldString)
		self.assertEquals(tree.getLatex(),oldLatex)
		self.assertEquals(tree.getFunction(),oldFunction)
		print 'Assertion passed, get() values equal old values.'
		# now register the update
		tree.nodesHaveChanged()
		# now the old information should NOT be the same as the new stuff
		self.assertNotEquals(tree.getString(),oldString)
		self.assertNotEquals(tree.getLatex(),oldLatex)
		self.assertNotEquals(tree.getFunction(),oldFunction)
		print 'Assertion passed, change properly registered everywhere.'

		
if __name__ == '__main__':
	unittest.main()