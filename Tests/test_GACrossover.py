#!/usr/bin/env python
"""Tests different Genetic Algorithm crossover classes.
"""
# standard library
import sys
import string

# biopython
from Bio.Seq import MutableSeq

# local stuff
from Bio.GA.Organism               import Organism
from Bio.GA.Crossover.General      import SafeFitnessCrossover
from Bio.GA.Crossover.GeneralPoint import GeneralPointCrossover
from Bio.GA.Crossover.GeneralPoint import InterleaveCrossover
from Bio.GA.Crossover.TwoPoint     import TwoPointCrossover
from Bio.GA.Crossover.Point        import SinglePointCrossover
from Bio.GA.Crossover.Uniform      import UniformCrossover


# PyUnit
import unittest

def run_tests(argv):
    ALL_TESTS = [SinglePointTest, FourPointTest, InterleaveTest,
                 TwoPointTest, UniformTest, SafeFitnessTest]
    
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    
    for test in ALL_TESTS:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        runner.run(cur_suite)

class TestAlphabet:
    """Simple test alphabet.
    """
    letters = ["1", "2", "3"]

    def contains( self, oalpha ): return 1
    
def test_fitness(genome):
    """Simple class for calculating fitnesses.
    """
    seq_genome = genome.toseq()
    return int(seq_genome.data)

class SinglePointTest(unittest.TestCase):
    """Test simple point crossovers.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        genome_1 = MutableSeq("11111", self.alphabet)
        self.org_1 = Organism(genome_1, test_fitness)

        genome_2 = MutableSeq("22222", self.alphabet)
        self.org_2 = Organism(genome_2, test_fitness)
        
        self.crossover  = SinglePointCrossover(1.0)
	
    def t_basic_crossover(self):
        """Test basic point crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.crossover.do_crossover(self.org_1,
			                                   self.org_2)

        assert new_org_1.genome != start_genome_1 and \
               new_org_2.genome != start_genome_2, \
               "Did not perform a crossover when expected."

        assert new_org_1 != self.org_1 and new_org_2 != self.org_2, \
               "Returned an exact copy of the original organism."

class UniformTest(unittest.TestCase):
    """Test simple point crossovers.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        genome_1 = MutableSeq("11111111", self.alphabet)
        self.org_1 = Organism(genome_1, test_fitness)

        genome_2 = MutableSeq("22222222", self.alphabet)
        self.org_2 = Organism(genome_2, test_fitness)
        
	genome_3 = MutableSeq("333", self.alphabet)
	self.org_3 = Organism(genome_3, test_fitness)
	
        self.crossover = UniformCrossover(1.0, 0.8)

    def t_basic_crossover(self):
        """Test basic uniform crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.crossover.do_crossover(self.org_1,
                                                           self.org_2)

        assert new_org_1.genome != start_genome_1 and \
               new_org_2.genome != start_genome_2, \
               "Did not perform a crossover when expected."

        assert new_org_1 != self.org_1 and new_org_2 != self.org_2, \
               "Returned an exact copy of the original organism."
	return

    def t_ds_prop_uniform_crossover(self):
	"""Test properties of differing genome length, uniform crossovers.
	"""
        new_org_1, new_org_2 = self.crossover.do_crossover(self.org_1,
			                                   self.org_3)


	assert len(new_org_1.genome) > len(new_org_2.genome), \
	       "Strings are of wrong sizes after uniform crossover."

        assert string.count( new_org_2.genome.tostring(), "1" ) ==  \
	       string.count( new_org_1.genome.tostring(), "3" ),    \
	       "There should be equal distributions of the smaller string"

	assert self.org_1.genome[len(new_org_2.genome):] == \
	       new_org_1.genome[len(new_org_2.genome):], \
	       "Uniform should not touch non-overlapping elements of genome"
	return
    def t_ss_prop_uniform_crossover(self):
	"""Test properties of equal genome length, uniform crossovers.
	"""
        new_org_1, new_org_2 = self.crossover.do_crossover(self.org_1,
			                                   self.org_2)


	assert len(new_org_1.genome) == len(new_org_2.genome), \
	       "Strings are of different sizes after symmetric crossover."

        assert string.count( new_org_1.genome.tostring(), "1" ) ==  \
	       string.count( new_org_2.genome.tostring(), "2" ) and \
	       string.count( new_org_1.genome.tostring(), "2" ) ==  \
 	       string.count( new_org_2.genome.tostring(), "1" ),    \
	       "There should be equal, inverse distributions"
	return

class InterleaveTest(unittest.TestCase):
    """Test 'simple' 4-point crossovers.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        genome_1 = MutableSeq("11111", self.alphabet)
        self.org_1 = Organism(genome_1, test_fitness)

        genome_2 = MutableSeq("22222", self.alphabet)
        self.org_2 = Organism(genome_2, test_fitness)
        
        genome_3 = MutableSeq("333333333", self.alphabet)
        self.org_3 = Organism(genome_3, test_fitness)
        
	self._crossover  = InterleaveCrossover(1.0)

    def t_basic_crossover(self):
        """Test basic interleave crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self._crossover.do_crossover(self.org_1,
	                                                    self.org_2)

	assert new_org_1.genome != start_genome_1 and \
               new_org_2.genome != start_genome_2, \
               "Did not perform a crossover when expected."

        assert new_org_1 != self.org_1 and new_org_2 != self.org_2, \
               "Returned an exact copy of the original organism."
	return
	       
    def t_prop_sym_crossover(self):
	"""Test properties of interleave point crossover.
	"""
        new_org_1, new_org_2 = self._crossover.do_crossover(self.org_1,
				                            self.org_2)

	assert len(new_org_1.genome) == len(new_org_2.genome), \
	       "Strings are of different sizes after symmetric crossover."

        assert string.count( new_org_1.genome.tostring(), "1" ) ==  \
	       string.count( new_org_2.genome.tostring(), "2" ) and \
	       string.count( new_org_1.genome.tostring(), "2" ) ==  \
 	       string.count( new_org_2.genome.tostring(), "1" ),    \
	       "There should be equal, inverse distributions"
	       
	assert new_org_1.genome.tostring() == "12121" and  \
	       new_org_2.genome.tostring() == "21212",     \
	       "Did not interleave."
	       
	return

    def t_prop_asym_crossover(self):
	"""Test basic interleave crossover with asymmetric genomes.
	"""
        start_genome_1 = self.org_1.genome[:]
        start_genome_3 = self.org_3.genome[:]

	new_org_1, new_org_3 = self._crossover.do_crossover(self.org_1,
	                                                    self.org_3)
	
	assert new_org_1.genome != start_genome_1 and \
               new_org_3.genome != start_genome_3, \
               "Did not perform a crossover when expected."

        assert new_org_1 != self.org_1 and new_org_3 != self.org_3, \
               "Returned an exact copy of the original organism."
	       
	assert new_org_1.genome.tostring() == "13131" and  \
	       new_org_3.genome.tostring() == "31313333",  \
	       "Did not interleave with growth."
	return
    
class FourPointTest(unittest.TestCase):
    """Test 'simple' 4-point crossovers.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        genome_1 = MutableSeq("11111", self.alphabet)
        self.org_1 = Organism(genome_1, test_fitness)

        genome_2 = MutableSeq("22222", self.alphabet)
        self.org_2 = Organism(genome_2, test_fitness)
        
	self.sym_crossover  = GeneralPointCrossover(3,1.0)
	self.asym_crossover = GeneralPointCrossover(4,1.0)
	
    def t_basic_crossover(self):
        """Test basic 4-point crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.sym_crossover.do_crossover(self.org_1,
	                                                       self.org_2)

	assert new_org_1.genome != start_genome_1 and \
               new_org_2.genome != start_genome_2, \
               "Did not perform a crossover when expected."

        assert new_org_1 != self.org_1 and new_org_2 != self.org_2, \
               "Returned an exact copy of the original organism."
	return
	       
    def t_prop_sym_crossover(self):
	"""Test properties of symmetric 4-point crossover.
	"""
        new_org_1, new_org_2 = self.sym_crossover.do_crossover(self.org_1,
	                                                       self.org_2)

	assert len(new_org_1.genome) == len(new_org_2.genome), \
	       "Strings are of different sizes after symmetric crossover."

	
        assert string.count( new_org_1.genome.tostring(), "1" ) ==  \
	       string.count( new_org_2.genome.tostring(), "2" ) and \
	       string.count( new_org_1.genome.tostring(), "2" ) ==  \
 	       string.count( new_org_2.genome.tostring(), "1" ),    \
	       "There should be equal, inverse distributions"
	       
	return
    
    def t_basic_asym_crossover(self):
        """Test basic asymmetric 2-point crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.asym_crossover.do_crossover(self.org_1,
	                                                        self.org_2)

        assert new_org_1.genome != start_genome_1 and \
               new_org_2.genome != start_genome_2, \
               "Did not perform a crossover when expected."

        assert new_org_1 != self.org_1 and new_org_2 != self.org_2, \
               "Returned an exact copy of the original organism."

	return
    
    
class TwoPointTest(unittest.TestCase):
    """Test simple 2-point crossovers.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        genome_1 = MutableSeq("11111111", self.alphabet)
        self.org_1 = Organism(genome_1, test_fitness)

        genome_2 = MutableSeq("22222222", self.alphabet)
        self.org_2 = Organism(genome_2, test_fitness)
        
	self.asym_crossover = TwoPointCrossover(1.0)

    def t_basic_asym_crossover(self):
        """Test basic asymmetric 2-point crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.asym_crossover.do_crossover(self.org_1,
	                                                        self.org_2)

        assert new_org_1.genome != start_genome_1 and \
               new_org_2.genome != start_genome_2, \
               "Did not perform a crossover when expected."

        assert new_org_1 != self.org_1 and new_org_2 != self.org_2, \
               "Returned an exact copy of the original organism."

	return
    
class TestCrossover:
    """Provide basic crossover functionality for testing SafeFitness.
    """
    def __init__(self):
        # whether or not to produce new organisms with lower fitness
        # higher fitness, or the same organism
        self.type = "lower"

    def do_crossover(self, org_1, org_2):
        seq_org1 = org_1.genome.toseq()
        seq_org2 = org_2.genome.toseq()
        org1_genome = seq_org1.data
        org2_genome = seq_org2.data

        new_org_1 = org_1.copy()
        new_org_2 = org_2.copy()
        
        if self.type == "same":
            return new_org_1, new_org_2
        elif self.type == "lower":
            new_org1_genome = str(int(org1_genome) - 1)
            new_org2_genome = str(int(org2_genome) - 1)

            new_org_1.genome = MutableSeq(new_org1_genome,
                                          org_1.genome.alphabet)
            new_org_2.genome = MutableSeq(new_org2_genome,
                                          org_2.genome.alphabet)
        elif self.type == "higher":
            new_org1_genome = str(int(org1_genome) + 1)
            new_org2_genome = str(int(org2_genome) + 1)
        else:
            raise ValueError("Got type %s" % self.type)

        new_org_1.genome = MutableSeq(new_org1_genome,
                                      org_1.genome.alphabet)
        new_org_2.genome = MutableSeq(new_org2_genome,
                                      org_2.genome.alphabet)

        return new_org_1, new_org_2
                
class SafeFitnessTest(unittest.TestCase):
    """Tests for crossovers which do not reduce fitness.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        genome_1 = MutableSeq("2", self.alphabet)
        self.org_1 = Organism(genome_1, test_fitness)

        genome_2 = MutableSeq("2", self.alphabet)
        self.org_2 = Organism(genome_2, test_fitness)

        self.test_crossover = TestCrossover()

    def t_keep_higher(self):
        """Make sure we always keep higher fitness when specified.
        """
        crossover = SafeFitnessCrossover(self.test_crossover)

        self.test_crossover.type = "same"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        assert (new_org_1 == self.org_1 and new_org_2 == self.org_2), \
               "Did not retain organism for same fitness."

        self.test_crossover.type = "lower"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        assert (new_org_1 == self.org_1 and new_org_2 == self.org_2), \
               "Did not retain organism when crossover had lower fitness."

        self.test_crossover.type = "higher"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        assert (new_org_1.fitness > self.org_1.fitness and
                new_org_2.fitness > self.org_2.fitness), \
                "Did not get new organism when it had higher fitness."

    def t_keep_lower(self):
        """Make sure we do normal crossover functionality when specified.
        """
        crossover = SafeFitnessCrossover(self.test_crossover, 1.0)

        self.test_crossover.type = "same"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        assert (new_org_1 == self.org_1 and new_org_2 == self.org_2), \
               "Did not retain organism for same fitness."

        self.test_crossover.type = "lower"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        assert (new_org_1 != self.org_1 and new_org_2 != self.org_2), \
               "Did not retain lower fitness organism in crossover."

        self.test_crossover.type = "higher"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        assert (new_org_1.fitness > self.org_1.fitness and
                new_org_2.fitness > self.org_2.fitness), \
                "Did not get new organism under higher fitness conditions."

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
