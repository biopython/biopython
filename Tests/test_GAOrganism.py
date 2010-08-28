#!/usr/bin/env python
"""Tests for an Organism in a Genetic Algorithm population.
"""
# standard library
import unittest

# Biopython
from Bio import Alphabet
from Bio.Seq import MutableSeq


# local stuff
from Bio.GA import Organism

    

# -- utility functions
class TestAlphabet(Alphabet.Alphabet):
    """Simple alphabet for test purposes.
    """
    letters = ["1", "2", "3", "4"]

def genome_generator():
    """Generate a genome for testing purposes.
    """
    return MutableSeq("1234", TestAlphabet())

def fitness_calculator(genome):
    """Calculate fitness for testing purposes.
    """
    assert isinstance(genome, MutableSeq), "Expected MutableSeq for a genome."

    regular_seq = genome.toseq()
    return int(regular_seq.tostring())

class CreatePopulationTest(unittest.TestCase):
    """Tests for utility functions for creating populations.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()

    def test_function_population(self):
        """Create a population using a function to generate genomes.
        """
        num_orgs = 10
        new_pop = Organism.function_population(genome_generator,
                                               num_orgs, fitness_calculator)

        assert len(new_pop) == num_orgs, "Expected %s organisms, got %s" \
               % (num_orgs, len(new_pops))

        for org in new_pop:
            assert isinstance(org, Organism.Organism), \
                   "Expected to get an organism, got %r" % org

            exp_fit = fitness_calculator(org.genome)
            assert org.fitness == exp_fit, \
                   "Expected fitness of %s, got %s" % (org.fitness, exp_fit)

    def test_random_population(self):
        """Create a population randomly from a alphabet.
        """
        num_orgs = 10
        genome_size = 5
        new_pop = Organism.random_population(self.alphabet, genome_size,
                                             num_orgs, fitness_calculator)

        assert len(new_pop) == num_orgs, "Expected %s organisms, got %s" \
               % (num_orgs, len(new_pops))

        for org in new_pop:
            assert isinstance(org, Organism.Organism), \
                   "Expected to get an organism, got %r" % org

            exp_fit = fitness_calculator(org.genome)
            assert org.fitness == exp_fit, \
                   "Expected fitness of %s, got %s" % (org.fitness, exp_fit)

            assert len(org.genome) == genome_size, \
                   "Expected genome size of %s, got %s" % (len(org.genome),
                                                           genome_size)

    def test_random_population_types(self):
        """Creating a random population with different types of alphabets.
        """
        class DoubleAlphabet:
            letters = [1.0, 2.0]

        class CharacterAlphabet:
            letters = ["a", "b"]

        class IntegerAlphabet:
            letters = [1, 2]

        def test_fitness(genome):
            return 2

        all_alphabets = [DoubleAlphabet(), CharacterAlphabet(),
                         IntegerAlphabet()]

        for alphabet in all_alphabets:
            new_pop = Organism.random_population(alphabet, 5, 10,
                                                 test_fitness)

class OrganismTest(unittest.TestCase):
    """Tests for an organism in a GA population.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        self.genome = MutableSeq("1234", self.alphabet)
        self.organism = Organism.Organism(self.genome, fitness_calculator)

    def test_organism_basic(self):
        """Exercise basic organism functionality.
        """
        same_genome = MutableSeq("1234", self.alphabet)
        same_organism = Organism.Organism(same_genome, fitness_calculator)

        dif_genome = MutableSeq("1111", self.alphabet)
        dif_organism = Organism.Organism(dif_genome, fitness_calculator)

        assert str(self.organism) == str(same_organism), \
               "Comparison doesn't work for identical organisms."

        assert str(self.organism) != str(dif_organism), \
               "Comparison doesn't work for different organism."

    def test_organism_fitness(self):
        """Test the ability to deal with the fitness of the genome.
        """
        assert self.organism.fitness == 1234, \
               "Unexpected fitness %s" % self.organism.fitness
        
        new_genome = MutableSeq("1111", self.alphabet)
        self.organism.genome = new_genome
        self.organism.recalculate_fitness()

        assert self.organism.fitness == 1111, \
               "Unexpected fitness %s" % self.organism.fitness

    def test_organism_copy(self):
        """Test copying of organisms.
        """
        new_organism = self.organism.copy()

        new_organism.genome.append("1")

        assert str(new_organism.genome) != str(self.organism.genome), \
               "Did not provide a copy of the organism."

    def test_provide_fitness(self):
        """Test that providing a pre-calculated fitness works.
        """
        def fitness_calc(genome):
            raise ValueError("Should not have been executed.")

        genome = self.organism.genome

        # make sure not supplying fitness works
        try:
            new_org = Organism.Organism(genome, fitness_calc)
            raise AssertionError("Did not calculate fitness when expected.")
        except ValueError:
            pass

        # make sure supplying fitness works
        new_org = Organism.Organism(genome, fitness_calc, 50)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
