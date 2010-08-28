#!/usr/bin/env python
"""Tests different Genetic Algorithm crossover classes.
"""
# standard library
import unittest

# biopython
from Bio.Seq import MutableSeq
from Bio.Alphabet import SingleLetterAlphabet

# local stuff
from Bio.GA.Organism import Organism
from Bio.GA.Crossover.General import SafeFitnessCrossover
from Bio.GA.Crossover.GeneralPoint import GeneralPointCrossover
from Bio.GA.Crossover.GeneralPoint import InterleaveCrossover
from Bio.GA.Crossover.TwoPoint import TwoPointCrossover
from Bio.GA.Crossover.Point import SinglePointCrossover
from Bio.GA.Crossover.Uniform import UniformCrossover


class TestAlphabet(SingleLetterAlphabet):
    """Simple test alphabet.
    """
    letters = ["1", "2", "3"]

    def contains(self, oalpha):
        return True

def test_fitness(genome):
    """Simple class for calculating fitnesses.
    """
    seq_genome = genome.toseq()
    return int(seq_genome.tostring())

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

    def test_basic_crossover(self):
        """Test basic point crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.crossover.do_crossover(self.org_1,
                                                           self.org_2)
        self.assertNotEqual(str(new_org_1.genome), str(start_genome_1),
                            "Did not perform a crossover when expected.")
        self.assertNotEqual(str(new_org_2.genome), str(start_genome_2),
                            "Did not perform a crossover when expected.")
        
        self.assertNotEqual(str(new_org_1), str(self.org_1),
                            "Returned an exact copy of the original organism.")
        self.assertNotEqual(str(new_org_2), str(self.org_2),
                            "Returned an exact copy of the original organism.")

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

    def test_basic_crossover(self):
        """Test basic uniform crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.crossover.do_crossover(self.org_1,
                                                           self.org_2)

        self.assertNotEqual(str(new_org_1.genome), str(start_genome_1),
                            "Did not perform a crossover when expected.")
        self.assertNotEqual(str(new_org_2.genome), str(start_genome_2),
                            "Did not perform a crossover when expected.")
    
        self.assertNotEqual(str(new_org_1), str(self.org_1),
                            "Returned an exact copy of the original organism.")
        self.assertNotEqual(str(new_org_2), str(self.org_2),
                            "Returned an exact copy of the original organism.")

    def test_ds_prop_uniform_crossover(self):
        """Test properties of differing genome length, uniform crossovers.
        """
        new_org_1, new_org_2 = self.crossover.do_crossover(self.org_1,
                                                           self.org_3)


        self.assertTrue(len(new_org_1.genome) > len(new_org_2.genome),
                     "Strings are of wrong sizes after uniform crossover.")

        self.assertEqual(new_org_2.genome.tostring().count("1"),
                         new_org_1.genome.tostring().count("3"),
                         "There should be equal distributions of the smaller string")

        self.assertEqual(str(self.org_1.genome[len(new_org_2.genome):]),
                         str(new_org_1.genome[len(new_org_2.genome):]),
                         "Uniform should not touch non-overlapping elements of genome")
    
    def test_ss_prop_uniform_crossover(self):
        """Test properties of equal genome length, uniform crossovers.
        """
        new_org_1, new_org_2 = self.crossover.do_crossover(self.org_1,
                                                           self.org_2)

        self.assertEqual(len(new_org_1.genome), len(new_org_2.genome),
                         "Strings are of different sizes after uniform crossover.")

        self.assertEqual(new_org_1.genome.tostring().count("1"),
                         new_org_2.genome.tostring().count("2"),
                         "There should be equal, inverse distributions")
        self.assertEqual(new_org_1.genome.tostring().count("2") ,
                         new_org_2.genome.tostring().count("1"),
                         "There should be equal, inverse distributions")


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

    def test_basic_crossover(self):
        """Test basic interleave crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self._crossover.do_crossover(self.org_1,
                                                            self.org_2)

        self.assertNotEqual(str(new_org_1.genome), str(start_genome_1),
                            "Did not perform a crossover when expected.")
        self.assertNotEqual(str(new_org_2.genome), str(start_genome_2),
                            "Did not perform a crossover when expected.")

        self.assertNotEqual(str(new_org_1), str(self.org_1),
                            "Returned an exact copy of the original organism.")
        self.assertNotEqual(str(new_org_2), str(self.org_2),
                            "Returned an exact copy of the original organism.")
       
    def test_prop_sym_crossover(self):
        """Test properties of interleave point crossover."""
        new_org_1, new_org_2 = self._crossover.do_crossover(self.org_1,
                                                            self.org_2)

        self.assertEqual(len(new_org_1.genome), len(new_org_2.genome),
            "Strings are of different sizes after interleave point crossover.")

        self.assertEqual(new_org_1.genome.tostring().count("1"),
                         new_org_2.genome.tostring().count("2"),
                         "There should be equal, inverse distributions")
        self.assertEqual(new_org_1.genome.tostring().count("2") ,
                         new_org_2.genome.tostring().count("1"),
                         "There should be equal, inverse distributions")
       
        self.assertEqual(new_org_1.genome.tostring(), "12121",
                         "Did not interleave.")
        self.assertEqual(new_org_2.genome.tostring(), "21212",
                         "Did not interleave.")

    def test_prop_asym_crossover(self):
        """Test basic interleave crossover with asymmetric genomes."""
        start_genome_1 = self.org_1.genome[:]
        start_genome_3 = self.org_3.genome[:]

        new_org_1, new_org_3 = self._crossover.do_crossover(self.org_1,
                                                            self.org_3)

        self.assertNotEqual(str(new_org_1.genome), str(start_genome_1),
                            "Did not perform a crossover when expected.")
        self.assertNotEqual(str(new_org_3.genome), str(start_genome_3),
                            "Did not perform a crossover when expected.")

        self.assertNotEqual(str(new_org_1), str(self.org_1),
                            "Returned an exact copy of the original organism.")
        self.assertNotEqual(str(new_org_3), str(self.org_3),
                            "Returned an exact copy of the original organism.")

        self.assertEqual(new_org_1.genome.tostring(), "13131",
                         "Did not interleave with growth.")
        self.assertEqual(new_org_3.genome.tostring(), "31313333",
                         "Did not interleave with growth.")
    
class FourPointTest(unittest.TestCase):
    """Test 'simple' 4-point crossovers.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        genome_1 = MutableSeq("11111", self.alphabet)
        self.org_1 = Organism(genome_1, test_fitness)

        genome_2 = MutableSeq("22222", self.alphabet)
        self.org_2 = Organism(genome_2, test_fitness)
        
        self.sym_crossover = GeneralPointCrossover(3,1.0)
        self.asym_crossover = GeneralPointCrossover(4,1.0)

    def test_basic_crossover(self):
        """Test basic 4-point crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.sym_crossover.do_crossover(self.org_1,
                                                               self.org_2)

        self.assertNotEqual(str(new_org_1.genome), str(start_genome_1),
                            "Did not perform a crossover when expected.")
        self.assertNotEqual(str(new_org_2.genome), str(start_genome_2),
                            "Did not perform a crossover when expected.")

        self.assertNotEqual(str(new_org_1), str(self.org_1),
                            "Returned an exact copy of the original organism.")
        self.assertNotEqual(str(new_org_2), str(self.org_2),
                            "Returned an exact copy of the original organism.")
               
    def test_prop_sym_crossover(self):
        """Test properties of symmetric 4-point crossover.
        """
        new_org_1, new_org_2 = self.sym_crossover.do_crossover(self.org_1,
                                                               self.org_2)

        self.assertEqual(len(new_org_1.genome), len(new_org_2.genome),
                         "Strings are of different sizes after symmetric crossover.")

        self.assertEqual(new_org_1.genome.tostring().count("1"),
                         new_org_2.genome.tostring().count("2"),
                         "There should be equal, inverse distributions")
        self.assertEqual(new_org_1.genome.tostring().count("2") ,
                         new_org_2.genome.tostring().count("1"),
                         "There should be equal, inverse distributions")
    
    def test_basic_asym_crossover(self):
        """Test basic asymmetric 2-point crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.asym_crossover.do_crossover(self.org_1,
                                                                self.org_2)

        self.assertNotEqual(str(new_org_1.genome), str(start_genome_1),
                            "Did not perform a crossover when expected.")
        self.assertNotEqual(str(new_org_2.genome), str(start_genome_2),
                            "Did not perform a crossover when expected.")

        self.assertNotEqual(str(new_org_1), str(self.org_1),
                            "Returned an exact copy of the original organism.")
        self.assertNotEqual(str(new_org_2), str(self.org_2),
                            "Returned an exact copy of the original organism.")
    
    
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

    def test_basic_asym_crossover(self):
        """Test basic asymmetric 2-point crossover functionality.
        """
        start_genome_1 = self.org_1.genome[:]
        start_genome_2 = self.org_2.genome[:]
        
        new_org_1, new_org_2 = self.asym_crossover.do_crossover(self.org_1,
                                                                self.org_2)

        self.assertNotEqual(str(new_org_1.genome), str(start_genome_1),
                            "Did not perform a crossover when expected.")
        self.assertNotEqual(str(new_org_2.genome), str(start_genome_2),
                            "Did not perform a crossover when expected.")

        self.assertNotEqual(str(new_org_1), str(self.org_1),
                            "Returned an exact copy of the original organism.")
        self.assertNotEqual(str(new_org_2), str(self.org_2),
                            "Returned an exact copy of the original organism.")

    
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
        org1_genome = seq_org1.tostring()
        org2_genome = seq_org2.tostring()

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

    def test_keep_higher(self):
        """Make sure we always keep higher fitness when specified.
        """
        crossover = SafeFitnessCrossover(self.test_crossover)

        self.test_crossover.type = "same"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        self.assertEqual(str(new_org_1), str(self.org_1),
                         "Did not retain organism for same fitness.")
        self.assertEqual(str(new_org_2), str(self.org_2),
                         "Did not retain organism for same fitness.")

        self.test_crossover.type = "lower"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        self.assertEqual(str(new_org_1), str(self.org_1),
                         "Did not retain organism when crossover had lower fitness.")
        self.assertEqual(str(new_org_2), str(self.org_2),
                         "Did not retain organism when crossover had lower fitness.")

        self.test_crossover.type = "higher"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        self.assertTrue(new_org_1.fitness > self.org_1.fitness and \
                     new_org_2.fitness > self.org_2.fitness,
                     "Did not get new organism when it had higher fitness.")

    def test_keep_lower(self):
        """Make sure we do normal crossover functionality when specified.
        """
        crossover = SafeFitnessCrossover(self.test_crossover, 1.0)

        self.test_crossover.type = "same"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        self.assertEqual(str(new_org_1), str(self.org_1),
                         "Did not retain organism for same fitness.")
        self.assertEqual(str(new_org_2), str(self.org_2),
                         "Did not retain organism for same fitness.")

        self.test_crossover.type = "lower"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        self.assertNotEqual(str(new_org_1), str(self.org_1),
                         "Retained lower fitness organism in crossover.")
        self.assertNotEqual(str(new_org_2), str(self.org_2),
                         "Retained lower fitness organism in crossover.")
        
        self.test_crossover.type = "higher"
        new_org_1, new_org_2 = crossover.do_crossover(self.org_1, self.org_2)

        self.assertTrue(new_org_1.fitness > self.org_1.fitness and \
                     new_org_2.fitness > self.org_2.fitness,
                     "Did not get new organism under higher fitness conditions.")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
