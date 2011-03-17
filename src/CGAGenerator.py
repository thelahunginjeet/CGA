#!/usr/bin/env python

# TODO :
#	- check the root being selected in all of the function calls
#	- remove deprecated versions of mutation operators

import unittest
from numpy.random import uniform
from numpy.random import randint
import random, copy
import CGAFunctions
from CGAStructures import DataNode, BinaryNode, UnaryNode, ScalarNode, AlgorithmTree


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
		raise TypeError, "this is a utility class with static methods; don't initialize me"
	
	@staticmethod
	def _extend(tNode, newNode):
		"""Grows a tree by replacing a terminal node (any one in tree.termini) with the 
		node newNode.  
			Properties:
				-This will never involve the root, so no root checking is required. 
				-newNode can be a dataNode, in which case it will be a terminal leaf (replacing an empty node)
				-The new terminal nodes are not filled in.
				-Should not cause a GC problem; tNode has no children."""
		assert tNode.getHeader() is False
		identity = tNode.getIdentity()
		if identity:  # right node
			tNode.parent.setChildren(None, newNode)
		else:  # left node
			tNode.parent.setChildren(newNode, None)
		tNode.parent = None
		tNode.clean()
		
	@staticmethod
	def _delete(fNode):
		"""Truncates a tree by replacing the given functional Node with a newly created, random
		 data node.  This effectively removes an entire subtree.
			Properties:
				-If functionalNode == root, return without deletion, as this would just delete the
					entire tree.
				-Might cause a GC problem; None-ing fNode's .parent variable is an attempt to 
					remove external references and hope GC does its job."""
		assert type(fNode) is BinaryNode or UnaryNode or ScalarNode
		if fNode.getHeader(): # don't do this to the root
			pass
		else:
			dNode = DataNode()
			if fNode.getIdentity():  # right node
				fNode.parent.setChildren(None, dNode)
			else:  # left node
				fNode.parent.setChildren(dNode, None)
			fNode.parent = None
			fNode.clean()
	
	@staticmethod
	def _replace(tree, oldNode, newNode):
		"""Replaces one node with another of the same type.  So Binary->Binary, Unary->Unary, 
		Data->Data, and Scalar->Scalar
			Properties:
				-Checking that the types of the new and old nodes are the same is necessary.
				-The root can be the oldNode; we need to ensure the new node is known as the root.
				-Unfortunately, to ensure proper behavior when oldNode == root, you have to pass
					the tree into this function as well.
				-GC should be OK here;  parents and children for the oldNode are set to None before
					refcount decrement, so no external references should linger."""
		assert type(oldNode) is type(newNode)
		childLeft, childRight = oldNode.getChildren()
		if oldNode.getHeader():
			newNode.setHeader(True)
			tree.root = newNode
		else:
			newNode.setHeader(False)
			if oldNode.getIdentity():
				oldNode.parent.setChildren(None, newNode)
			else:
				oldNode.parent.setChildren(newNode, None)
		newNode.setIdentity(oldNode.getIdentity())
		newNode.setChildren(childLeft, childRight)
		# DO NOT USE .clean() HERE!  IT WILL None NODES TOO FAR DOWN IN THE TREE
		oldNode.parent = None
		oldNode.setChildren(None,None)

							
	@staticmethod
	def _swap(nodeOne, nodeTwo):
		"""Swaps two subtrees.  The nodes may be any node except the root; swapping roots (or one root)
		involves a lot of painful checking, and isn't that useful.
			Properties:
				-If either nodeOne or nodeTwo is the root, the function just returns with no crossover.
				-GC should be fine in this case; we don't actually want to delete anything, just
					update the references."""
		if nodeOne.getHeader() or nodeTwo.getHeader():
			return None
		oneIdentity, twoIdentity = nodeOne.getIdentity(), nodeTwo.getIdentity()
		# make copies of the subtrees to cut in
		tmpOne, tmpTwo = nodeOne.copy(),nodeTwo.copy()
		tmpOne.parent, tmpTwo.parent = None, None
		if oneIdentity:  # one right
			nodeOne.parent.setChildren(None, tmpTwo)
			if twoIdentity:  # two right
				nodeTwo.parent.setChildren(None, tmpOne)
			else:  # two left
				nodeTwo.parent.setChildren(tmpOne, None)
		else:  # one left
			nodeOne.parent.setChildren(tmpTwo)
			if twoIdentity:  # two right
				nodeTwo.parent.setChildren(None, tmpOne)
			else:  # two left
				nodeTwo.parent.setChildren(tmpOne, None)
		nodeOne.parent = None
		nodeTwo.parent = None
		nodeOne.clean()
		nodeTwo.clean()
											
	@staticmethod
	def _getRandomFunctionalNode():
		if uniform < 0.5:
			return BinaryNode(CGAFunctions.DataMethodFactory().getBinary())
		return UnaryNode(CGAFunctions.DataMethodFactory().getUnary())
				
	@staticmethod
	def generate(number=10):
		"""Function to generate and return a random tree with fixed number of internal nodes added.  The root has to 
		be a scalarizing node."""
		root = ScalarNode(CGAFunctions.DataMethodFactory().getScalar())
		tree = AlgorithmTree(root)
		while number > 0:
			tree.update()
			tNode = random.choice(tree.termini)
			fNode = CGAGenerator._getRandomFunctionalNode()
			CGAGenerator._extend(tNode,fNode)
			number -= 1
		return tree

	@staticmethod
	def expgenerate(p):
		"""Generates and returns a random tree.  Current tree is extended at a *single* node with
		probability p, and process terminates with probability 1-p."""
		root = ScalarNode(CGAFunctions.DataMethodFactory().getScalar())
		tree = AlgorithmTree(root)
		extend = True
		while extend:
			tree.update()
			tNode = random.choice([x for x in tree.termini])
			fNode = CGAGenerator._getRandomFunctionalNode()
			CGAGenerator._extend(tNode,fNode)
			if uniform() < p:
				extend = True
			else:
				extend = False
		return tree

	@staticmethod
	def point_mutate(tree,node):
		"""Replaces node with a random node of the same type."""
		tree.update()		
		if isinstance(node, UnaryNode):
			newNode = UnaryNode(CGAFunctions.DataMethodFactory().getUnary())
		elif isinstance(node, DataNode):
			newNode = DataNode(CGAFunctions.DataMethodFactory().getData())
		elif isinstance(node, BinaryNode):
			newNode = BinaryNode(CGAFunctions.DataMethodFactory().getBinary())
		elif isinstance(node, ScalarNode):
			newNode = ScalarNode(CGAFunctions.DataMethodFactory().getScalar())
		else:
			raise TypeError, "there seems to be a mystery node; try again . . ."
		CGAGenerator._replace(tree, node, newNode)
		tree.update()
		
	@staticmethod
	def prune(tree,node):
		"""Accepts one input tree and removes the subtree rooted at node; node is replaced with a terminal 
		data node, so the tree remains evaluateable.  If the node chosen is a data node or the root,
		nothing will be done - the tree will be unmodified."""
		tree.update()
		if type(node) is DataNode:
			pass
		else:
			CGAGenerator._delete(node)
		tree.update()
		
	@staticmethod
	def grow(tree,tNode):
		"""Takes an input tree and extends it by changing node into a new (random) Unary or 
		Binary node.  New random data nodes are added as terminal leaves."""
		assert type(tNode) is DataNode
		tree.update()
		fNode = CGAGenerator._getRandomFunctionalNode()
		CGAGenerator._extend(tNode, fNode)
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
		if nodeOne.getHeader() or nodeTwo.getHeader():
			pass
		else:	
			CGAGenerator._swap(nodeOne, nodeTwo)
			treeOne.update()
			treeTwo.update()


class CGAGeneratorTests(unittest.TestCase):
	def setUp(self):
		self.methodFactory = CGAFunctions.DataMethodFactory()
		# set up the tree
		self.root = ScalarNode(self.methodFactory.getScalar('sum_ij'))
		self.node1 = BinaryNode(self.methodFactory.getBinary('+'))
		self.node2 = BinaryNode(self.methodFactory.getBinary('-'))
		self.node3 = UnaryNode(self.methodFactory.getUnary('log'))
		# immutable data - stored
		self.constant1 = DataNode(self.methodFactory.getData('e'))
		self.constant2 = DataNode(self.methodFactory.getData('pi'))
		self.constant3 = DataNode(self.methodFactory.getData('1'))
		# set up the tree
		self.testTree = AlgorithmTree(self.root)
		self.root.setChildren(self.node1)
		self.node1.setChildren(self.node2,self.node3)
		self.node2.setChildren(self.constant1,self.constant2)
		self.node3.setChildren(self.constant3)
		
	def testExtend(self):
		print "\n\n----- testing extension -----"
		self.testTree.update()
		self.testTree()
		print 'Tree before extension (node,id): '
		print [(n.string,id(n)) for n in self.testTree.nodes]
		self.testTree()
		newNode = BinaryNode(self.methodFactory.getBinary('*'))
		CGAGenerator._extend(self.constant1,newNode)
		print 'Tree after extension (node,id): '
		self.testTree.update()
		self.testTree()
		print [(n.string,id(n)) for n in self.testTree.nodes]
		
	def testGrow(self):
		print "\n\n----- testing grow(tree,tNode) -----"
		self.testTree.update()
		self.testTree()
		print "Tree before extension (node,id): "
		print [(n.string,id(n)) for n in self.testTree.nodes]
		CGAGenerator.grow(self.testTree,self.constant2)
		self.testTree.update()
		self.testTree()
		print "Tree after extension (node,id): "
		print [(n.string,id(n)) for n in self.testTree.nodes]
		
	def testDelete(self):
		print "\n\n----- testing deletion -----"
		self.testTree.update()
		self.testTree()
		print 'Tree before deletion (node, id): '
		print [(n.string,id(n)) for n in self.testTree.nodes]
		CGAGenerator._delete(self.node2)
		self.testTree.update()
		self.testTree()
		print 'Tree after deletion (node, id): '
		print [(n.string,id(n)) for n in self.testTree.nodes]
	
	def testPrune(self):
		print "\n\n----- testing prune(tree,node) -----"
		self.testTree.update()
		self.testTree()
		print "Tree before pruning (node,id): "
		print [(n.string,id(n)) for n in self.testTree.nodes]
		CGAGenerator.prune(self.testTree,self.node3)
		self.testTree.update()
		self.testTree()
		print "Tree after pruning (node,id): "
		print [(n.string,id(n)) for n in self.testTree.nodes]
		
	def testNonRootReplace(self):
		print "\n\n----- testing non-root replacement -----"
		self.testTree.update()
		self.testTree()
		print 'Tree before replacement (node,id): '
		print [(n.string,id(n)) for n in self.testTree.nodes]
		newNode = BinaryNode(self.methodFactory.getBinary('*'))
		CGAGenerator._replace(self.testTree, self.node1, newNode)
		self.testTree.update()
		self.testTree()
		print 'Tree after replacement (node,id): '
		print [(n.string,id(n)) for n in self.testTree.nodes]
	
	def testRootReplace(self):
		print "\n\n----- testing replacement of root node -----"
		self.testTree.update()
		self.testTree()
		print 'Tree before root replacement (node,id): '
		print [(n.string,id(n)) for n in self.testTree.nodes]
		print 'Root : ', self.testTree.root
		print 'Root header : ', self.testTree.root.getHeader()
		newNode = ScalarNode(self.methodFactory.getScalar('tr'))
		CGAGenerator._replace(self.testTree,self.testTree.root,newNode)
		self.testTree.update()
		self.testTree()
		print 'Tree after root replacement (node, id): '
		print [(n.string,id(n)) for n in self.testTree.nodes]
		print 'Root : ', self.testTree.root
		print 'Root header : ', self.testTree.root.getHeader()
		
	def testSwap(self):
		print "\n\n----- testing subtree swap -----"
		root = ScalarNode(self.methodFactory.getScalar('tr'))
		node1 = UnaryNode(self.methodFactory.getUnary('log'))
		node2 = BinaryNode(self.methodFactory.getBinary('*'))
		constant1 = DataNode(self.methodFactory.getData('1/2'))
		constant2 = DataNode(self.methodFactory.getData('e'))
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(node2)
		node2.setChildren(constant1,constant2)
		self.testTree.update()
		self.testTree()
		tree.update()
		tree()
		print 'Trees before the swap (node,id): '
		print 'Tree One:'
		print [(n.string,id(n)) for n in self.testTree.nodes]
		print 'Tree Two:'
		print [(n.string,id(n)) for n in tree.nodes]
		CGAGenerator._swap(self.node2,node1)
		self.testTree.update()
		self.testTree()
		tree.update()
		tree()
		print 'Trees after the swap (node,id): '
		print 'Tree One:'
		print [(n.string,id(n)) for n in self.testTree.nodes]
		print 'Tree Two:'
		print [(n.string,id(n)) for n in tree.nodes]
	
	def testCrossover(self):
		print "\n\n----- testing single_crossover(tree1,tree1Node,tree2,tree2Node) -----"
		root = ScalarNode(self.methodFactory.getScalar('tr'))
		node1 = UnaryNode(self.methodFactory.getUnary('log'))
		node2 = BinaryNode(self.methodFactory.getBinary('*'))
		constant1 = DataNode(self.methodFactory.getData('1/2'))
		constant2 = DataNode(self.methodFactory.getData('e'))
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(node2)
		node2.setChildren(constant1,constant2)
		self.testTree.update()
		self.testTree()
		tree.update()
		tree()
		print 'Trees before the swap (node,id): '
		print 'Tree One:'
		print [(n.string,id(n)) for n in self.testTree.nodes]
		print 'Tree Two:'
		print [(n.string,id(n)) for n in tree.nodes]
		CGAGenerator.single_crossover(self.testTree,self.node2,tree,node1)
		self.testTree.update()
		self.testTree()
		tree.update()
		tree()
		print 'Trees after the swap (node,id): '
		print 'Tree One:'
		print [(n.string,id(n)) for n in self.testTree.nodes]
		print 'Tree Two:'
		print [(n.string,id(n)) for n in tree.nodes]
		
	def testGenerate(self):
		print "\n\n----- testing generate(tree) -----"
		tree = CGAGenerator.generate(3)
		tree()		
		print tree
	
	def testExpGenerate(self):
		print "\n\n----- testing expgenerate(tree) -----"
		tree = CGAGenerator.expgenerate(0.85)
		tree()
		print tree
	
	def testPointMutate(self):
		print "\n\n----- testing point_mutate(tree,node) -----"
		self.testTree.update()
		self.testTree()
		print "Tree before mutation (node,id): "
		print [(n.string,id(n)) for n in self.testTree.nodes]
		CGAGenerator.point_mutate(self.testTree,random.choice(self.testTree.nodes))
		self.testTree.update()
		self.testTree()
		print "Tree after mutation (node,id): "
		print [(n.string,id(n)) for n in self.testTree.nodes]

	def testRepeatedPointMutate(self):
		print "\n\n----- testing repeated runs of point_mutate(tree,node) -----"
		self.testTree.update()
		self.testTree()
		print "Tree before mutation (node,id): "
		print [(n.string,id(n)) for n in self.testTree.nodes]
		for i in xrange(0,10000):
			CGAGenerator.point_mutate(self.testTree,random.choice(self.testTree.nodes))
			self.testTree.update()
		self.testTree.update()
		self.testTree()
		print "Tree after mutation (node,id): "
		print [(n.string,id(n)) for n in self.testTree.nodes]
	
		
		
if __name__ == '__main__':
	unittest.main()