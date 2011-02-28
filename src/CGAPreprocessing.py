#!/usr/bin/env python

import unittest, re, cPickle, os
import Bio.PDB as PDB
from numpy import matrix, zeros
from itertools import izip

class Utilities(object):
    """Simple class of utility functions"""
    def __init__(self):
        pass
    
    @staticmethod
    def readFastaSequences(seqFile):
        """Read in sequences from fasta file and return dictionary"""
        sequences = dict()
        seqFile = open(seqFile,'r')
        listInput = seqFile.read().split('>')[1:]
        seqFile.close()
        for seqRecord in listInput:
            reSeq = re.findall('.+',seqRecord)
            seq = str()
            for seqPart in reSeq[1:]:
                seq += seqPart
            sequences[reSeq[0].strip()]=seq
        return sequences
    
    @staticmethod
    def loadMSA(msaFile):
        """Function to load a multiple sequence alignment and return """
        sequences = Utilities.readFastaSequences(msaFile)
        # make sure all of the sequences are uppercase for algorithms
        for s in [x for x in sequences if x != '#=GC RF']:
            sequences[s] = sequences[s].upper()
        # check to make sure that every sequence is the same length (aka an alignment)
        if len({}.fromkeys([len(x) for x in sequences.values()])) > 1:
            raise ValueError, "your MSA doesn't have matrix dimension; perhaps you meant another file . . ."
        else:
            dimensions = (len(sequences),len(sequences.values()[0]))
            # make a fixed list of sequences to iterate through excluding HMMER reference
            sequenceList = [sequences[x] for x in sequences if x != '#=GC RF']
            columns = {}
            for i in xrange(dimensions[1]):
                columns[i] = [x[i] for x in sequenceList]
        return sequences, columns
    
    @staticmethod            
    def translatePositions(sequences, columns):
        """Function to make a mapping (dict) between the MSA position and individual sequence positions"""
        mapping = dict()
        for column in columns:
            position = dict() # Uses physical numbering of positions and NOT Python numbering
            for sequence in sequences:
                if sequences[sequence][column] is '-':
                    position[sequence] = None
                else:
                    position[sequence] = len(sequences[sequence][:column+1])-sequences[sequence][:column+1].count('-')
            mapping[column] = position
        return mapping
      
    @staticmethod      
    def calculateGapFrequency(columns):
        """Function to calculate the gap frequency as a function of column position"""
        gaps = dict()
        for column in columns:
            gaps[column] = columns[column].count('-')/float(len(columns[column]))
        return gaps
    
    @staticmethod
    def writeDatabase(databaseFile, dictionary):
        """Writes a database file of the single and joint frequencies
            - dict['singleFrequencies'] = dict of matrices (20x1)
            - dict['jointFrequencies'] = dict of matrices (20x20)"""
        databaseFile = open(databaseFile, 'wb')
        cPickle.dump(dictionary, databaseFile, -1)
        databaseFile.close()
        
    @staticmethod
    def readDatabase(databaseFile):
        """Reads in a database file containing indices, single frequencies, and joint frequencies"""
        if os.path.exists(databaseFile):
            databaseFile = open(databaseFile,'rb')
            dictionary = cPickle.load(databaseFile)
            return dictionary
        else:
            raise IOError, 'your database file does not exist; please check your file name . . .'

    @staticmethod
    def calculateAtomicDistances(pdbFile):
        """Calculates the residue-residue distances from a PDB file"""
        """From an input PDB file, compute the residue-residue distance and return a dictionary on the edges that
         gives the minimal distance between two residues in Angstroms, calculated as the minimum distance between 
         Cb atoms (Ca if residue == Gly).  Ignores waters or any other funky ligands."""
        strucParser = PDB.PDBParser()
        pdbStruc = strucParser.get_structure('pdb', pdbFile)
        modelNumber = pdbStruc.child_dict.keys()[0]
        chain = pdbStruc[modelNumber].child_dict.keys()[0]
        residues = pdbStruc[modelNumber][chain].child_dict
        distances = {}
        aminoAcids =  ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        for residue1 in [x for x in residues if residues[x].resname in aminoAcids]:
            for residue2 in [x for x in residues if x != residue1 and residues[x].resname in aminoAcids]:
                if residues[residue1].resname == 'GLY' and residues[residue2].resname != 'GLY':
                    distances[(residue1[1], residue2[1])] = residues[residue1]['CA'] - residues[residue2]['CB']
                elif residues[residue1].resname != 'GLY' and residues[residue2].resname == 'GLY':
                    distances[(residue1[1], residue2[1])] =  residues[residue1]['CB'] - residues[residue2]['CA']
                elif residues[residue1].resname == 'GLY' and residues[residue2].resname == 'GLY':
                    distances[(residue1[1], residue2[1])] = residues[residue1]['CA'] - residues[residue2]['CA']
                else:
                    distances[(residue1[1], residue2[1])] = residues[residue1]['CB'] - residues[residue2]['CB']
        return distances


class MSA(object):
    """This is a simple multiple sequence alignment class to be used with MSAAlgorithms"""
    def __init__(self, alnFile, canonSequence, databaseFile=None, gapThreshold=0.20):
        """Simple constructor that loads the sequences, columns, and dimensions and translates the positions"""
        sequences, columns = Utilities.loadMSA(alnFile)
        if canonSequence not in sequences:
            raise KeyError, "your canonical sequence isnt' in the input MSA; please check your canonical sequence header . . ."
        mapping = Utilities.translatePositions(sequences, columns)
        gaps = Utilities.calculateGapFrequency(columns)
        # check to see if a reference HMMER sequence is present
        if sequences.has_key('#=GC RF'):
            hmmer = True
        else:
            hmmer = False
        # calculate position frequencies using the canonical sequence for mapping
        aminoAcids = tuple('ACDEFGHIKLMNPQRSTVWY')
        self.indices = []
        self.singleFrequencies = {}
        for column in columns:
            pmap = mapping[column][canonSequence]
            if pmap is not None and gaps[column] < gapThreshold:
                self.indices.append(pmap)
                self.singleFrequencies[pmap] = matrix(zeros([20, 1], dtype='float64'))
                for i in range(len(aminoAcids)):
                    observations = float(len([x for x in columns[column] if x is not '-']))
                    self.singleFrequencies[pmap][i, 0] = columns[column].count(aminoAcids[i])/observations
        self.jointFrequencies = {}
        for column1 in columns:
            for column2 in [x for x in columns if x >= column1]:
                pmap1 = mapping[column1][canonSequence]
                pmap2 = mapping[column2][canonSequence]
                if pmap1 is not None and pmap2 is not None and gaps[column1] < gapThreshold and gaps[column2] < gapThreshold:
                    self.jointFrequencies[(pmap1, pmap2)] = matrix(zeros([20, 20], dtype='float64'))
                    for i in range(len(aminoAcids)):
                        for j in range(len(aminoAcids)):
                            pairs = [p for p in izip(columns[column1], columns[column2]) if '-' not in p]
                            observations = float(len(pairs))
                            self.jointFrequencies[(pmap1, pmap2)][i, j] = pairs.count((aminoAcids[i], aminoAcids[j]))/observations
#        self.indices = [mapping[x][canonSequence] for x in columns if mapping[x][canonSequence] is not None]
        if databaseFile is not None:
            dictionary = {'indices':self.indices, 'singleFrequencies':self.singleFrequencies, 'jointFrequencies':self.jointFrequencies}
            Utilities.writeDatabase(databaseFile, dictionary)


class CGAFunctionsTests(unittest.TestCase):
    def setUp(self):
        pass
       
    def testMSA(self):
        # this will write a test database file
        MSA('../tests/pdz_test.aln', '1IU0', '../tests/pdz_test.db')
    
        
    
if __name__ == '__main__':
    unittest.main()
