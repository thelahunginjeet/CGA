#!/usr/bin/env python

# TODO :
#	- check the root being selected in all of the function calls
#	- remove deprecated versions of mutation operators

# NOTES :
#   - as of 1/4/11, all "atomic" and GA-type functions pass basic tests, behaving as expected

import unittest
from numpy.random import uniform
from numpy.random import randint
import random, copy
from CGAFunctions import BinaryFunctions, UnaryFunctions, ScalarizingFunctions, ImmutableData
from CGAStructures import EmptyNode, DataNode, BinaryNode, UnaryNode, ScalarNode, AlgorithmTree


class CGAGenerator(object):
	"""General class to make generate random trees, make mutations, crossovers, etc.  
	The four "atomic" methods are:
		1. _extend()  : replace a terminal node with a functional node
		2. _delete()  : remove an entire subtree, terminating the tree at that node
		3. _replace() : in situ replacement of nodes, preserving identity
		4. _swap()    : swap two subtrees
	These base methods are used to make more GA-type operations, which are designed to ensure
	return of an evaluateable tree (the atomic methods do not guarantee this)."""
	def __init__(self):
		pass
	
	@staticmethod
	def _extend(tNode,newNode):
		"""Grows a tree by replacing a terminal node (any one in tree.termini) with the 
		node newNode.  
			Properties:
				-This will never involve the root, so no root checking is required. 
				-newNode can be a dataNode, in which case it will be a terminal leaf (replacing an empty node)
				-The new terminal nodes are not filled in.
				-Should not cause a GC problem; tNode has no children."""
		# DEBUG - maybe want to check that terminalNode is either EmptyNode or DataNode, and pass if not?
		identity = tNode.getIdentity()
		if identity == 0:
			tNode.parent.setChildren(newNode,None)
			tNode.parent = None
			del tNode
		elif identity == 1:
			tNode.parent.setChildren(None,newNode)
			tNode.parent = None
			del tNode
		else:
			raise TypeError, "The identity of your node (%s) isn't equal to 0 or 1; try again . . ." % (identity)
		
	@staticmethod
	def _delete(fNode):
		"""Truncates a tree by replacing the given functional Node with an empty node.  This effectively
		removes an entire subtree.
			Properties:
				-If functionalNode == root, return without deletion, as this would just delete the
					entire tree.
				-Might cause a GC problem; None-ing fNode's .parent variable is an attempt to 
					remove external references and hope GC does its job."""
		header = fNode.getHeader()
		if header == True:
			pass
		else:
			identity = fNode.getIdentity()
			if identity == 0:
				fNode.parent.setChildren(EmptyNode(),None)
				fNode.parent = None
				fNode.setChildren(None,None)
				del fNode
			elif identity == 1:
				fNode.parent.setChildren(None,EmptyNode())
				fNode.parent = None
				fNode.setChildren(None,None)
				del fNode
			else:
				raise TypeError, "The identity of your node (%s) isn't equal to 0 or 1; try again . . ." % (identity)
	
	@staticmethod
	def _replace(tree,oldNode,newNode):
		"""Replaces one node with another of the same type, in situ.  So Binary->Binary and Unary->Unary, and
		Data->Data.
			Properties:
				-Checking that the types of the new and old nodes are the same is necessary.
				-The root can be the oldNode; we need to ensure the new node is known as the root.
				-Unfortunately, to ensure proper behavior when oldNode == root, you have to pass
					the tree into this function as well.
				-GC should be OK here;  parents and children for the oldNode are set to None before
					refcount decrement, so no external references should linger."""
		if type(oldNode) != type(newNode):
			raise TypeError, "You cannot do in-position replacement with a node of a different type."
		header = oldNode.getHeader()
		identity = oldNode.getIdentity()
		nodeChildren = oldNode.getChildren()
		if header == True:
			newNode.setHeader(True)
			newNode.setIdentity(identity)
			newNode.setChildren(nodeChildren[0],nodeChildren[1])
			tree.root = newNode
			oldNode.setChildren(None,None)
			del oldNode
		else:
			if identity == 0:
				oldNode.parent.setChildren(newNode,None)
				newNode.setChildren(nodeChildren[0],nodeChildren[1])
				oldNode.parent = None
				oldNode.setChildren(None,None)
				del oldNode
			elif identity == 1:
				oldNode.parent.setChildren(None, newNode)
				newNode.setChildren(nodeChildren[0],nodeChildren[1])
				oldNode.parent = None
				oldNode.setChildren(None,None)
				del oldNode
			else:
				raise TypeError, "the identity of your node (%s) isn't equal to 0 or 1; try again . . ."%(identity)
			
	@staticmethod
	def _swap(nodeOne,nodeTwo):
		"""Swaps two subtrees.  The nodes may be any node except the root; swapping roots (or one root)
		involves a lot of painful checking, and isn't that useful.
			Properties:
				-If either nodeOne or nodeTwo is the root, the function just returns with no crossover.
				-GC should be fine in this case; we don't actually want to delete anything, just
					update the references."""
		oneHeader = nodeOne.getHeader()
		twoHeader = nodeTwo.getHeader()
		if oneHeader == True or twoHeader == True:
			# don't allow this, just return
			return
		oneIdentity = nodeOne.getIdentity()
		twoIdentity = nodeTwo.getIdentity()
		tmpOne = copy.copy(nodeOne)
		tmpTwo = copy.copy(nodeTwo)
		tmpOne.parent, tmpTwo.parent = None, None
		if (oneIdentity, twoIdentity) == (0, 0):
			nodeOne.parent.setChildren(tmpTwo)
			nodeTwo.parent.setChildren(tmpOne)
		elif (oneIdentity, twoIdentity) == (0, 1):
			nodeOne.parent.setChildren(tmpTwo)
			nodeTwo.parent.setChildren(None, tmpOne)
		elif (oneIdentity, twoIdentity) == (1, 0):
			nodeOne.parent.setChildren(None, tmpTwo)
			nodeTwo.parent.setChildren(tmpOne)
		elif (oneIdentity, twoIdentity) == (1, 1):
			nodeOne.parent.setChildren(None, tmpTwo)
			nodeTwo.parent.setChildren(None, tmpOne)
		else:
			raise TypeError, "You have some node identity issues (%s, %s); try again . . ."%(oneIdentity, twoIdentity)		
				
							
	@staticmethod
	def _createRandomFunctionalNode():	
		if uniform() < 0.5:
			node = BinaryNode(BinaryFunctions().returnRandom())
		else:
			node = UnaryNode(UnaryFunctions().returnRandom())
		return node
	
	@staticmethod
	def _createRandomScalarNode():
		return ScalarNode(ScalarizingFunctions().returnRandom())

	@staticmethod
	def _createRandomDataNode():
		return DataNode(ImmutableData().returnRandom())
	
	@staticmethod
	def _fillTreeWithRandomData(tree):
		"""Replaces all terminal data nodes with stored, random data nodes."""
		tree.update()
		for node in tree.termini:
			if type(node) == EmptyNode:
				dNode = CGAGenerator._createRandomDataNode()
				CGAGenerator._extend(node,dNode)
		tree.update()
			
	@staticmethod
	def generate(number=10):
		"""Function to generate and return a random tree with fixed number of internal nodes added.  The root has to 
		be a scalarizing node."""
		root = CGAGenerator._createRandomScalarNode()
		tree = AlgorithmTree(root)
		while number > 0:
			tree.update()
			tNode = random.choice(tree.termini)
			fNode = CGAGenerator._createRandomFunctionalNode()
			CGAGenerator._extend(tNode,fNode)
			number -= 1
		CGAGenerator._fillTreeWithRandomData(tree)
		return tree

	@staticmethod
	def expgenerate(p):
		"""Generates and returns a random tree; nonterminal nodes are extended with probability p and 
		terminated with probability 1-p."""
		root = CGAGenerator._createRandomScalarNode()
		tree = AlgorithmTree(root)
		ntries = 1
		while True:
			tree.update()
			freenodes = [x for x in tree.termini if type(x) is EmptyNode]
			if len(freenodes) == 0:
				break # all nodes are terminal
			tNode = random.choice(freenodes)
			r = uniform()
			if r < p:
				# extend
				fNode = CGAGenerator._createRandomFunctionalNode()
				CGAGenerator._extend(tNode,fNode)
			else:
				# terminate at a data node
				dNode = CGAGenerator._createRandomDataNode()
				CGAGenerator._extend(tNode,dNode)
			ntries += 1
		return tree

	@staticmethod
	def point_mutate(tree,node):
		"""In-place switch of the node in tree, with a random node of the same type."""
		tree.update()		
		if isinstance(node, UnaryNode):
			newNode = UnaryNode(UnaryFunctions().returnRandom())
		elif isinstance(node, DataNode):
			newNode = DataNode(ImmutableData().returnRandom())
		elif isinstance(node, EmptyNode):
			newNode = EmptyNode()
		elif isinstance(node, BinaryNode):
			newNode = BinaryNode(BinaryFunctions().returnRandom())
		elif isinstance(node, ScalarNode):
			newNode = ScalarNode(ScalarizingFunctions().returnRandom())
		else:
			raise TypeError, "there seems to be a mystery node; try again . . ."
		CGAGenerator._replace(tree,node,newNode)
		tree.update()
		
	@staticmethod
	def prune(tree,node):
		"""Accepts one input tree and removes the subtree rooted at node; node is replaced with a terminal 
		data node, so the tree remains evaluateable.  If the node chosen is a data node or the root,
		nothing will be done - the tree will be unmodified."""
		tree.update()
		if type(node) is DataNode or type(node) is EmptyNode:
			pass
		else:
			CGAGenerator._delete(node)
			CGAGenerator._fillTreeWithRandomData(tree)
		tree.update()
		
	@staticmethod
	def grow(tree,tNode):
		"""Takes an input tree and extends it by changing node into a new (random) Unary or 
		Binary node.  New random data nodes are added as terminal leaves."""
		# TODO - check that the input node is terminal data; right now this can do things unsafe
		tree.update()
		fNode = CGAGenerator._createRandomFunctionalNode()
		CGAGenerator._extend(tNode, fNode)
		dNode1 = CGAGenerator._createRandomDataNode()
		dNode2 = CGAGenerator._createRandomDataNode()
		fNode.setChildren(dNode1,dNode2)
		tree.update()
		
	@staticmethod
	def single_crossover(treeOne,nodeOne,treeTwo,nodeTwo):
		"""Takes two input trees and swaps the subtrees rooted in nodeOne, nodeTwo respectively.
		So after this operation:
			treeOne's nodeOne subtree = treeTwo's nodeTwo subtree
			treeTwo's nodeTwo subtree = treeOne's nodeOne subtree.
		If either nodeOne or nodeTwo are their respective tree's root, the crossover will not be performed."""
		treeOne.update()
		treeTwo.update()
		# check for roots; don't do anything if one is the root
		if nodeOne.getHeader() == True or nodeTwo.getHeader() == True:
			pass
		else:	
			CGAGenerator._swap(nodeOne, nodeTwo)
			treeOne.update()
			treeTwo.update()

	'''
	DEPRECATED METHOD
	@staticmethod
	def mutate(tree):
		"""Takes an input tree and does an in-place switch of a random node.
		DEPRECATED."""
		tree.update()		
		treeNode = random.choice(tree.nodes)
		if isinstance(treeNode, UnaryNode):
			newNode = UnaryNode(UnaryFunctions().returnRandom())
		elif isinstance(treeNode, DataNode):
			newNode = DataNode(ImmutableData().returnRandom())
		elif isinstance(treeNode, EmptyNode):
			newNode = EmptyNode()
		elif isinstance(treeNode, BinaryNode):
			newNode = BinaryNode(BinaryFunctions().returnRandom())
		elif isinstance(treeNode, ScalarNode):
			newNode = ScalarNode(ScalarizingFunctions().returnRandom())
		else:
			raise TypeError, "there seems to be a mystery node; try again . . ."
		CGAGenerator._replace(tree,treeNode,newNode)
		tree.update()
	'''	
	'''
	DEPRECATED METHOD	
	@staticmethod
	def grow(tree):
		"""Takes an input tree and extends it by changing a random terminal data node into a new Unary or 
		Binary node.  New random data nodes are added as terminal leaves."""
		tree.update()
		tNode = random.choice(tree.termini)
		fNode = CGAGenerator._createRandomFunctionalNode()
		CGAGenerator._extend(tNode, fNode)
		dNode1 = CGAGenerator._createRandomDataNode()
		dNode2 = CGAGenerator._createRandomDataNode()
		fNode.setChildren(dNode1,dNode2)
		tree.update()
		
	
	@staticmethod
	def crossover(treeOne, treeTwo):
		"""Takes two input trees and swaps a random non-root node in treeOne with a random 
		non-root node in treeTwo."""
		treeOne.update()
		treeTwo.update()
		while True:
			nodeOne = random.choice(treeOne.nodes)
			if nodeOne.getHeader() == False:
				break
		while True:
			nodeTwo = random.choice(treeTwo.nodes)
			if nodeTwo.getHeader() == False:
				break
		CGAGenerator._swap(nodeOne, nodeTwo)
		treeOne.update()
		treeTwo.update()
		
	'''
	'''
	DEPRECATED METHOD
	@staticmethod
	def prune(tree):
		"""Accepts one input tree and deletes a randomly chosen subtree; the node chosen is replaced with a 
		terminal data node, so the tree remains evaluateable.  If the node chosen is a data node or the root,
		nothing will be done - the tree will be unmodified."""
		tree.update()
		ranNode = random.choice(tree.nodes)
		print ranNode
		if type(ranNode) is DataNode or type(ranNode) is EmptyNode:
			pass
		else:
			CGAGenerator._delete(ranNode)
			CGAGenerator._fillTreeWithRandomData(tree)
		tree.update()
	'''

	@staticmethod
	def duplicate():
		pass
		

class CGAGeneratorTests(unittest.TestCase):
	def setUp(self):
		self.scalFunctions = ScalarizingFunctions()
		self.binaryFunctions = BinaryFunctions()
		self.unaryFunctions = UnaryFunctions()
		self.data = ImmutableData()
		# set up the tree
		self.root = ScalarNode(self.scalFunctions['sum_ij'])
		self.node1 = BinaryNode(self.binaryFunctions['+'])
		self.node2 = BinaryNode(self.binaryFunctions['-'])
		self.node3 = UnaryNode(self.unaryFunctions['log'])
		# immutable data - stored
		self.constant1 = DataNode(self.data['e'])
		self.constant2 = DataNode(self.data['pi'])
		self.constant3 = DataNode(self.data['1'])
		# set up the tree
		self.testTree = AlgorithmTree(self.root)
		self.root.setChildren(self.node1)
		self.node1.setChildren(self.node2,self.node3)
		self.node2.setChildren(self.constant1,self.constant2)
		self.node3.setChildren(self.constant3)
		
	def testExtend(self):
		print "\n\n----- testing extension -----"
		random.choice(self.testTree.termini)
		print 'Tree before extension: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		newNode = BinaryNode(self.binaryFunctions['*'])
		CGAGenerator._extend(self.constant1,newNode)
		print 'Tree after extension: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		
	def testDelete(self):
		print "\n\n----- testing deletion -----"
		print 'Tree before deletion: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		CGAGenerator._delete(self.node1)
		print 'Tree after deletion: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		
	def testNonRootReplace(self):
		print "\n\n----- testing replacement -----"
		print 'Tree before replacement: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		newNode = BinaryNode(self.binaryFunctions['*'])
		CGAGenerator._replace(self.testTree,self.node1,newNode)
		print 'Tree after replacement: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		
	def testRootReplace(self):
		print "\n\n----- testing replacement of root node -----"
		print 'Tree before root replacement: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		newNode = ScalarNode(self.scalFunctions['tr'])
		CGAGenerator._replace(self.testTree,self.root,newNode)
		print 'Tree after root replacement: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		
	def testSwap(self):
		print "\n\n----- testing subtree swap -----"
		# need another tree to swap with the test tree
		root = UnaryNode(self.unaryFunctions['log'])
		node1 = BinaryNode(self.binaryFunctions['*'])
		constant1 = DataNode(self.data['1/2'])
		constant2 = DataNode(self.data['e'])
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(constant1,constant2)
		print 'Trees before the swap: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		tree.update()
		tree()
		print tree
		CGAGenerator._swap(self.node1,node1)
		print 'Trees after the swap: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		tree.update()
		tree()
		print tree
		
	def testGenerate(self):
		print "\n\n----- testing generate(tree) -----"
		tree = CGAGenerator.generate(3)
		tree()		
		print tree
	
	def testExpGenerate(self):
		print "\n\n----- testing expgenerate(tree) -----"
		tree = CGAGenerator.expgenerate(0.6)
		tree()
		print tree
		
	def testPointMutate(self):
		print "\n\n----- testing point_mutate(tree,node) -----"
		print "Tree before mutation: "
		self.testTree.update()
		self.testTree()
		print self.testTree
		CGAGenerator.point_mutate(self.testTree,random.choice(self.testTree.nodes))
		print "Tree after mutation: "
		self.testTree.update()
		self.testTree()
		print self.testTree
		
	def testRepeatedPointMutate(self):
		print "\n\n----- testing repeated runs of point_mutate(tree,node) -----"
		print "Tree before mutation: "
		self.testTree.update()
		self.testTree()
		print self.testTree
		for i in xrange(0,10000):
			CGAGenerator.point_mutate(self.testTree,random.choice(self.testTree.nodes))
		print "Tree after mutation: "
		self.testTree.update()
		self.testTree()
		print self.testTree
		
	""" DEPRECATED TEST
	def testMutate(self):
		print "\n\n----- testing mutate(tree) -----"
		print "Tree before mutation: "
		self.testTree.update()
		self.testTree()
		print self.testTree
		CGAGenerator.mutate(self.testTree)
		print "Tree after mutation: "
		self.testTree.update()
		self.testTree()
		print self.testTree
	"""
	
	def testGrow(self):
		print "\n\n----- testing grow(tree,tNode) -----"
		print "Tree before extension: "
		self.testTree.update()
		self.testTree()
		print self.testTree
		CGAGenerator.grow(self.testTree,self.constant2)
		print "Tree after extension: "
		self.testTree.update()
		self.testTree()
		print self.testTree
	
	
	""" DEPRECATED TEST
	def testGrow(self):
		print "\n\n----- testing grow(tree) -----"
		print "Tree before extension: "
		self.testTree.update()
		self.testTree()
		print self.testTree
		CGAGenerator.grow(self.testTree)
		print "Tree after extension: "
		self.testTree.update()
		self.testTree()
		print self.testTree
	"""
	
	def testCrossover(self):
		print "\n\n----- testing crossover(tree1,node1,tree2,node2) -----"
		# need another tree to crossover with the test tree
		root = UnaryNode(self.unaryFunctions['log'])
		node1 = BinaryNode(self.binaryFunctions['*'])
		constant1 = DataNode(self.data['1/2'])
		constant2 = DataNode(self.data['e'])
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(constant1,constant2)
		print 'Trees before the crossover: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		tree.update()
		tree()
		print tree
		CGAGenerator.single_crossover(self.testTree,self.node2,tree,node1)
		print 'Trees after the crossover: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		tree.update()
		tree()
		print tree
		
	"""
	DEPRECATED TEST
	def testCrossover(self):
		print "\n\n----- testing crossover(tree1,tree2) -----"
		# need another tree to crossover with the test tree
		root = UnaryNode(self.unaryFunctions['log'])
		node1 = BinaryNode(self.binaryFunctions['*'])
		constant1 = DataNode(self.data['1/2'])
		constant2 = DataNode(self.data['e'])
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(constant1,constant2)
		print 'Trees before the crossover: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		tree.update()
		tree()
		print tree
		CGAGenerator.crossover(self.testTree,tree)
		print 'Trees after the crossover: '
		self.testTree.update()
		self.testTree()
		print self.testTree
		tree.update()
		tree()
		print tree
	"""
	
	def testPrune(self):
		print "\n\n----- testing prune(tree,node) -----"
		# need to make something that isn't just roots and data nodes so
		# 	we are likely to pick a truncateable node
		root = UnaryNode(self.unaryFunctions['log'])
		node1 = UnaryNode(self.unaryFunctions['exp'])
		node2 = UnaryNode(self.unaryFunctions['tanh'])
		node3 = UnaryNode(self.unaryFunctions['sin'])
		node4 = BinaryNode(self.binaryFunctions['*'])
		constant1 = DataNode(self.data['1/2'])
		constant2 = DataNode(self.data['e'])
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(node2)
		node2.setChildren(node3)
		node3.setChildren(node4)
		node4.setChildren(constant1,constant2)
		print "Tree before pruning: "
		tree.update()
		tree()
		print tree
		CGAGenerator.prune(tree,node3)
		print "Tree after pruning: "
		tree.update()
		tree()
		print tree
	
	"""DEPRECATED TEST
	def testPrune(self):
		print "\n\n----- testing prune(tree) -----"
		# need to make something that isn't just roots and data nodes so
		# 	we are likely to pick a truncateable node
		root = UnaryNode(self.unaryFunctions['log'])
		node1 = UnaryNode(self.unaryFunctions['exp'])
		node2 = UnaryNode(self.unaryFunctions['tanh'])
		node3 = UnaryNode(self.unaryFunctions['sin'])
		node4 = BinaryNode(self.binaryFunctions['*'])
		constant1 = DataNode(self.data['1/2'])
		constant2 = DataNode(self.data['e'])
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(node2)
		node2.setChildren(node3)
		node3.setChildren(node4)
		node4.setChildren(constant1,constant2)
		print "Tree before pruning: "
		tree.update()
		tree()
		print tree
		CGAGenerator.prune(tree)
		print "Tree after pruning: "
		tree.update()
		tree()
		print tree
	"""	
		
if __name__ == '__main__':
	unittest.main()