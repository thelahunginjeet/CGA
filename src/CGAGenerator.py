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
	return of an evaluateable tree (the atomic methods do not guarantee this).  There are also
	methods to generate special trees (mutual information, OMES), to check fitness functions
	and such for them."""
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
		# MIGHT NEED TO USE A .COPY() METHOD HERE TO INSURE POINTER INDEPENDENCE
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
#		oldNode.clean()
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
		tmpOne, tmpTwo = nodeOne.copy(), nodeTwo.copy()
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
		nodeOne.parent, nodeTwo.parent = None, None
		nodeOne.clean()
		nodeTwo.clean()
											
	@staticmethod
	def _getRandomFunctionalNode(r=0.5):
		flip = uniform()
		if flip < r:
			return BinaryNode(CGAFunctions.DataMethodFactory().getBinary())
		else:
			if flip > 0.9:
				return ScalarNode(CGAFunctions.DataMethodFactory().getScalar())
			return UnaryNode(CGAFunctions.DataMethodFactory().getUnary())
				
	@staticmethod
	def generate(number=10,r=0.5):
		"""Function to generate and return a random tree with fixed number of internal nodes added.  
		The root has to be a scalarizing node.  Biasing of the tree towards (or away from) binary
		nodes can be controlled with the parameter r."""
		root = ScalarNode(CGAFunctions.DataMethodFactory().getScalar())
		tree = AlgorithmTree(root)
		while number > 0:
			tNode = random.choice(tree.getTermini())
			fNode = CGAGenerator._getRandomFunctionalNode(r)
			CGAGenerator._extend(tNode,fNode)
			number -= 1
		return tree	


	@staticmethod
	def expgenerate(p,r=0.5):
		"""Generates and returns a random tree.  Current tree is extended at a *single* node with
		probability p, and process terminates with probability 1-p."""
		root = ScalarNode(CGAFunctions.DataMethodFactory().getScalar())
		tree = AlgorithmTree(root)
		extend = True
		while extend:
			tNode = random.choice(tree.getTermini())
			fNode = CGAGenerator._getRandomFunctionalNode(r)
			CGAGenerator._extend(tNode,fNode)
			if uniform() < p:
				extend = True
			else:
				extend = False
		return tree

	
	@staticmethod
	def generate_special_tree(treename):
		"""Generates and returns a "special" tree.  Allowed treenames are:
				['MI','OMES']."""
		if treename == 'MI':
			root = ScalarNode(CGAFunctions.DataMethodFactory().getScalar('sum_ij'))
			n1 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('subtract'))
			n2 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('multiply'))
			n3 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('add'))
			d1 = DataNode(CGAFunctions.DataMethodFactory().getData('p_ij'))
			n4 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('log'))
			n5 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('multiply'))
			n6 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('multiply'))
			d2 = DataNode(CGAFunctions.DataMethodFactory().getData('p_ij'))
			d3 = DataNode(CGAFunctions.DataMethodFactory().getData('1/N'))
			n7 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('multiply'))
			d4 = DataNode(CGAFunctions.DataMethodFactory().getData('1/N'))
			n8 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('multiply'))
			d5 = DataNode(CGAFunctions.DataMethodFactory().getData('p_i'))
			n9 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('log'))
			n10 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('transpose'))
			n11 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('log'))
			d6 = DataNode(CGAFunctions.DataMethodFactory().getData('p_i'))
			d7 = DataNode(CGAFunctions.DataMethodFactory().getData('p_j'))
			n12 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('transpose'))
			d8 = DataNode(CGAFunctions.DataMethodFactory().getData('p_j'))
			# set up the tree
			specialTree = AlgorithmTree(root)
			root.setChildren(n1)
			n1.setChildren(n2,n3)
			n2.setChildren(d1,n4)
			n4.setChildren(d2)
			n3.setChildren(n5,n6)
			n5.setChildren(d3,n7)
			n6.setChildren(d4,n8)
			n7.setChildren(d5,n9)
			n8.setChildren(n10,n11)
			n9.setChildren(d6)
			n10.setChildren(d7)
			n11.setChildren(n12)
			n12.setChildren(d8)			
		elif treename == 'OMES':
			root = ScalarNode(CGAFunctions.DataMethodFactory().getScalar('sum_ij'))
			n1 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('divide'))
			n2 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('square'))
			n3 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('multiply'))
			n4 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('subtract'))
			d1 = DataNode(CGAFunctions.DataMethodFactory().getData('p_i'))
			n5 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('transpose'))
			d2 = DataNode(CGAFunctions.DataMethodFactory().getData('p_ij'))
			n6 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('multiply'))
			d3 = DataNode(CGAFunctions.DataMethodFactory().getData('p_j'))
			d4 = DataNode(CGAFunctions.DataMethodFactory().getData('p_i'))
			n7 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('transpose'))
			d5 = DataNode(CGAFunctions.DataMethodFactory().getData('p_j'))
			# set up the tree
			specialTree = AlgorithmTree(root)
			root.setChildren(n1)
			n1.setChildren(n2,n3)
			n2.setChildren(n4)
			n3.setChildren(d1,n5)
			n4.setChildren(d2,n6)
			n5.setChildren(d3)
			n6.setChildren(d4,n7)
			n7.setChildren(d5)
		elif treename == 'BADTR':
			root = ScalarNode(CGAFunctions.DataMethodFactory().getScalar('sum_ij'))
			n1 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('divide'))
			d1 = DataNode(CGAFunctions.DataMethodFactory().getData('p_i'))
			n2 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('subtract'))
			n3 = BinaryNode(CGAFunctions.DataMethodFactory().getBinary('divide'))
			n4 = ScalarNode(CGAFunctions.DataMethodFactory().getScalar('tr'))
			d2 = DataNode(CGAFunctions.DataMethodFactory().getData('e'))
			n5 = UnaryNode(CGAFunctions.DataMethodFactory().getUnary('transpose'))
			d3 = DataNode(CGAFunctions.DataMethodFactory().getData('p_i'))
			d4 = d2 = DataNode(CGAFunctions.DataMethodFactory().getData('e'))
			# set up the tree
			specialTree = AlgorithmTree(root)
			root.setChildren(n1)
			n1.setChildren(d1,n2)
			n2.setChildren(n3,n4)
			n3.setChildren(d2,n5)
			n4.setChildren(d3)
			n5.setChildren(d4)
		else:
			specialTree = None
		return specialTree

	
	@staticmethod
	def point_mutate(tree,node):
		"""Replaces node with a random node of the same type."""
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
		
	@staticmethod
	def prune(tree, node):
		"""Accepts one input tree and removes the subtree rooted at node; node is replaced with a terminal 
		data node, so the tree remains evaluateable.  If the node chosen is a data node or the root,
		nothing will be done - the tree will be unmodified."""
		if type(node) is DataNode:
			pass
		else:
			CGAGenerator._delete(node)
		
	@staticmethod
	def grow(tree, tNode):
		"""Takes an input tree and extends it by changing node into a new (random) Unary or 
		Binary node.  New random data nodes are added as terminal leaves."""
		assert type(tNode) is DataNode
		fNode = CGAGenerator._getRandomFunctionalNode(0.5)
		CGAGenerator._extend(tNode, fNode)
		
	@staticmethod
	def single_crossover(nodeOne, nodeTwo):
		"""Takes two input trees and swaps the subtrees rooted in nodeOne, nodeTwo respectively.
		So after this operation:
			treeOne's nodeOne subtree = treeTwo's nodeTwo subtree
			treeTwo's nodeTwo subtree = treeOne's nodeOne subtree.
		If either nodeOne or nodeTwo are their respective tree's root, the crossover will not be performed."""
		if nodeOne.getHeader() or nodeTwo.getHeader():
			pass
		else:	
			CGAGenerator._swap(nodeOne, nodeTwo)


class CGAGeneratorTests(unittest.TestCase):
	def setUp(self):
		self.methodFactory = CGAFunctions.DataMethodFactory()
		# set up the tree
		self.root = ScalarNode(self.methodFactory.getScalar('sum_ij'))
		self.node1 = BinaryNode(self.methodFactory.getBinary('add'))
		self.node2 = BinaryNode(self.methodFactory.getBinary('subtract'))
		self.node3 = UnaryNode(self.methodFactory.getUnary('log'))
		# immutable data - stored
		self.constant1 = DataNode(self.methodFactory.getData('e'))
		self.constant2 = DataNode(self.methodFactory.getData('pi'))
		self.constant3 = DataNode(self.methodFactory.getData('1/N'))
		# set up the tree
		self.testTree = AlgorithmTree(self.root)
		self.root.setChildren(self.node1)
		self.node1.setChildren(self.node2,self.node3)
		self.node2.setChildren(self.constant1,self.constant2)
		self.node3.setChildren(self.constant3)
		
	def testExtend(self):
		print "\n\n----- testing extension -----"
		print 'Tree before extension (node,id): '
		self.testTree()
		newNode = BinaryNode(self.methodFactory.getBinary('multiply'))
		CGAGenerator._extend(self.constant1,newNode)
		print 'Tree after extension (node,id): '
		self.testTree()
		
	def testGrow(self):
		print "\n\n----- testing grow(tree,tNode) -----"
		print "Tree before extension (node,id): "
		self.testTree()
		CGAGenerator.grow(self.testTree,self.constant2)
		print "Tree after extension (node,id): "
		self.testTree()
		
	def testDelete(self):
		print "\n\n----- testing deletion -----"
		print 'Tree before deletion (node, id): '
		self.testTree()
		CGAGenerator._delete(self.node2)
		print 'Tree after deletion (node, id): '
		self.testTree()
	
	def testPrune(self):
		print "\n\n----- testing prune(tree,node) -----"
		print "Tree before pruning (node,id): "
		self.testTree()
		CGAGenerator.prune(self.testTree,self.node3)
		print "Tree after pruning (node,id): "
		self.testTree()
		
	def testNonRootReplace(self):
		print "\n\n----- testing non-root replacement -----"
		print 'Tree before replacement (node,id): '
		self.testTree()
		newNode = BinaryNode(self.methodFactory.getBinary('multiply'))
		CGAGenerator._replace(self.testTree, self.node1, newNode)
		print 'Tree after replacement (node,id): '
		self.testTree()
	
	def testRootReplace(self):
		print "\n\n----- testing replacement of root node -----"
		print 'Tree before root replacement (node,id): '
		self.testTree()
		print 'Root : ', self.testTree.root
		print 'Root header : ', self.testTree.root.getHeader()
		newNode = ScalarNode(self.methodFactory.getScalar('tr'))
		CGAGenerator._replace(self.testTree,self.testTree.root,newNode)
		print 'Tree after root replacement (node, id): '
		self.testTree()
		print 'Root : ', self.testTree.root
		print 'Root header : ', self.testTree.root.getHeader()
		
	def testSwap(self):
		print "\n\n----- testing subtree swap -----"
		root = ScalarNode(self.methodFactory.getScalar('tr'))
		node1 = UnaryNode(self.methodFactory.getUnary('log'))
		node2 = BinaryNode(self.methodFactory.getBinary('multiply'))
		constant1 = DataNode(self.methodFactory.getData('1/N'))
		constant2 = DataNode(self.methodFactory.getData('e'))
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(node2)
		node2.setChildren(constant1,constant2)
		print 'Trees before the swap (node,id): '
		print 'Tree One:'
		self.testTree()		
		print 'Tree Two:'
		tree()
		CGAGenerator._swap(self.node2, node1)
		print 'Trees after the swap (node,id): '
		print 'Tree One:'
		self.testTree()
		print 'Tree Two:'
		tree()
	
	def testCrossover(self):
		print "\n\n----- testing single_crossover(tree1Node, tree2Node) -----"
		root = ScalarNode(self.methodFactory.getScalar('tr'))
		node1 = UnaryNode(self.methodFactory.getUnary('log'))
		node2 = BinaryNode(self.methodFactory.getBinary('multiply'))
		constant1 = DataNode(self.methodFactory.getData('1/N'))
		constant2 = DataNode(self.methodFactory.getData('e'))
		tree = AlgorithmTree(root)
		root.setChildren(node1)
		node1.setChildren(node2)
		node2.setChildren(constant1,constant2)
		print 'Trees before the swap (node,id): '
		print 'Tree One:'
		self.testTree()
		print 'Tree Two:'
		tree()
		CGAGenerator.single_crossover(self.node2, node1)
		print 'Trees after the swap (node,id): '
		print 'Tree One:'
		self.testTree()
		print 'Tree Two:'
		tree()
		
	def testGenerate(self):
		print "\n\n----- testing generate(tree) -----"
		tree = CGAGenerator.generate(3)
		tree()		
	
	def testExpGenerate(self):
		print "\n\n----- testing expgenerate(tree) -----"
		tree = CGAGenerator.expgenerate(0.85)
		tree()
	
	def testPointMutate(self):
		print "\n\n----- testing point_mutate(tree,node) -----"
		print "Tree before mutation (node,id): "
		self.testTree()
		CGAGenerator.point_mutate(self.testTree,random.choice(self.testTree.getNodes()))
		print "Tree after mutation (node,id): "
		self.testTree()

	def testRepeatedPointMutate(self):
		print "\n\n----- testing repeated runs of point_mutate(tree,node) -----"
		print "Tree before mutation (node,id): "
		self.testTree()
		for i in xrange(0,10000):
			CGAGenerator.point_mutate(self.testTree, random.choice(self.testTree.getNodes()))
		print "Tree after mutation (node,id): "
		self.testTree()
	
	def testSpecialTreeGeneration(self):
		print "\n\n----- testing generation of 'special' trees -----"
		treeset = ['MI','OMES','BADTR','DOES NOT EXIST']
		for t in treeset:
			tree = CGAGenerator.generate_special_tree(t)
			if tree is not None:
				print '%s : %s', (t,tree.getString(),tree.getLatex())
			else:
				print 'Tree %s does not exist.' % t
		
if __name__ == '__main__':
	unittest.main()