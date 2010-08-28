#!/usr/bin/env python
"""Tests for Genetic Algorithm Repair code.

This tests classes which are designed for repairing organisms after
mutation and crossover.
"""
# standard library
import unittest

# biopython
from Bio.Alphabet import Alphabet
from Bio.Seq import MutableSeq

# local stuff
from Bio.NeuralNetwork.Gene.Schema import Schema
from Bio.GA.Organism import Organism
from Bio.GA.Repair.Stabilizing import AmbiguousRepair



class TestAlphabet(Alphabet):
    """Simple test alphabet.
    """
    alphabet_matches = {"1": "1",
                        "2": "2",
                        "3": "3",
                        "*": "123"}
                        
    letters = ["1", "2", "3", "*"]

def test_fitness(genome):
    """Simple class for calculating fitnesses.
    """
    return 1

class AmbiguousRepairTest(unittest.TestCase):
    """Test for the ability to repair too many ambiguous genes in a genome.
    """
    def setUp(self):
        alphabet = TestAlphabet()
        test_genome = MutableSeq("11*22*33*", alphabet)
        self.organism = Organism(test_genome, test_fitness)
        
        self.ambig_info = Schema(alphabet.alphabet_matches)

    def test_single_repair(self):
        """Test repair of a single ambiguous position in a genome.
        """
        repairer = AmbiguousRepair(self.ambig_info, 2)

        for repair_attempt in range(5):
            new_org = repairer.repair(self.organism)
            new_genome_seq = new_org.genome.toseq()

            assert new_genome_seq.count("*") == 2, \
                   "Did not repair genome, got %s" % new_genome_seq.tostring()

    def test_multiple_repair(self):
        """Test repair of multiple ambiguous positions in a genome.
        """
        repairer = AmbiguousRepair(self.ambig_info, 0)

        for repair_attempt in range(5):
            new_org = repairer.repair(self.organism)
            new_genome_seq = new_org.genome.toseq()

            assert new_genome_seq.count("*") == 0, \
                   "Did not repair genome, got %s" % new_genome_seq.tostring()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
