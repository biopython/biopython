#!/usr/bin/env python
"""Tests for Genetic Algorithm mutation functionality.
"""
# standard library
import unittest

# biopython
from Bio.Seq import MutableSeq
from Bio.Alphabet import SingleLetterAlphabet

# local stuff
from Bio.GA.Organism import Organism
from Bio.GA.Mutation.General import SafeFitnessMutation
from Bio.GA.Mutation.Simple import ConversionMutation
from Bio.GA.Mutation.Simple import SinglePositionMutation


class TestAlphabet(SingleLetterAlphabet):
    """Simple test alphabet.
    """
    letters = ["1", "2", "3"]

def test_fitness(genome):
    """Simple class for calculating fitnesses.
    """
    seq_genome = genome.toseq()
    return int(seq_genome.tostring())

class MutationHelper:
    """Mixin class which provides useful functions for testing mutations.
    """
    num_trials = 500
    
    def _always_mutate(self, mutator, expected_percent):
        """Test the ability of a mutator to always mutate.

        Arguments:

        o mutator - The mutation class we're testing.

        o expected_percent - The minimum percent of mutations we expect
        to see under 'always mutate.' This will depend on how many letters
        are in the alphabet and other factors.
        """
        num_mutations = 0
        for trial in range(self.num_trials):
            new_org = mutator.mutate(self.organism)

            # if we see a visilble mutation, mark it down
            if new_org != self.organism:
                num_mutations += 1

        percent_mutants = float(num_mutations) / float(self.num_trials)
        assert percent_mutants > expected_percent, \
               "Did not recieve an acceptable number of mutations."

    def _never_mutate(self, mutator):
        """Test that a mutator does not cause unexpected mutations.
        """
        for trial in range(self.num_trials):
            new_org = mutator.mutate(self.organism)
            assert new_org == self.organism, "Unexpected mutation found"

class ConversionTest(unittest.TestCase, MutationHelper):
    """Test mutation which just converts one gene in the chromosome.
    """
    def setUp(self):
        genome = MutableSeq("1111", TestAlphabet())
        self.organism = Organism(genome, test_fitness)

    def test_always_mutate(self):
        """Test ability to cause mutations.
        """
        mutator = ConversionMutation(mutation_rate = 1.0)

        # when we mutate randomly by chance, we expect to get 2/3
        # visible mutations (there are three letters in the alphabet and
        # one change cannot be observed since it is a mutation back to itself)
        # For a four letter genome, the chance of being exactly the same
        # after mutations is about .01, so being better than 90% different
        # if a reasonable expectation.
        expected_percent = .9

        self._always_mutate(mutator, expected_percent)

    def test_never_mutate(self):
        """Make sure we do not mutate at unexpected times.
        """
        mutator = ConversionMutation(mutation_rate = 0.0)
        self._never_mutate(mutator)

class SinglePositionTest(unittest.TestCase, MutationHelper):
    """Test mutations at a single position in a genome.
    """
    def setUp(self):
        genome = MutableSeq("1111", TestAlphabet())
        self.organism = Organism(genome, test_fitness)

    def test_always_mutate(self):
        """Test ability to cause mutations.
        """
        mutator = SinglePositionMutation(mutation_rate = 1.0)

        # when we mutate randomly by chance, we expect to get 2/3
        # visible mutations (there are three letters in the alphabet and
        # one change cannot be observed since it is a mutation back to itself)
        expected_percent = .6

        self._always_mutate(mutator, expected_percent)

    def test_never_mutate(self):
        """Make sure we do not mutate at unexpected times.
        """
        mutator = SinglePositionMutation(mutation_rate = 0.0)
        self._never_mutate(mutator)

class TestMutator:
    """Provide basic mutator ability.
    """
    def __init__(self):
        self.type = "lower"

    def mutate(self, org):
        org_genome_seq = org.genome.toseq()
        old_org_genome = org_genome_seq.tostring()
        
        new_org = org.copy()
        
        if self.type == "same":
            return new_org
        elif self.type == "lower":
            new_org.genome = MutableSeq(str(int(old_org_genome) - 1),
                                        org_genome_seq.alphabet)
            return new_org
        elif self.type == "higher":
            new_org.genome = MutableSeq(str(int(old_org_genome) + 1),
                                        org_genome_seq.alphabet)
            return new_org
        else:
            raise ValueError("Got type %s" % self.type)

class SafeFitnessTest(unittest.TestCase):
    """Test mutation which does not allow decreases in fitness.
    """
    def setUp(self):
        self.alphabet = TestAlphabet()
        genome = MutableSeq("2", self.alphabet)
        self.org = Organism(genome, test_fitness)

        self.test_mutator = TestMutator()

    def test_keep_higher(self): 
        """Make sure we always keep the higher fitness.
        """
        mutator = SafeFitnessMutation(self.test_mutator)

        self.test_mutator.type = "same"
        new_org = mutator.mutate(self.org)
        assert (new_org == self.org), \
               "Did not retain organism for same fitness."

        self.test_mutator.type = "lower"
        new_org = mutator.mutate(self.org)
        assert (new_org == self.org), \
               "Did not retain organism when crossover had lower fitness."

        self.test_mutator.type = "higher"
        new_org = mutator.mutate(self.org)
        assert (new_org.fitness > self.org.fitness), \
                "Did not get new organism when it had higher fitness."

    def test_keep_new(self):
        """Make sure we always keep the new organism when specified.
        """
        mutator = SafeFitnessMutation(self.test_mutator, 1.0)

        self.test_mutator.type = "same"
        new_org = mutator.mutate(self.org)
        assert (new_org == self.org), \
               "Did not retain organism for same fitness."

        self.test_mutator.type = "lower"
        new_org = mutator.mutate(self.org)
        assert (new_org.fitness < self.org.fitness), \
               "Did not get new organism when it had lower fitness."

        self.test_mutator.type = "higher"
        new_org = mutator.mutate(self.org)
        assert (new_org.fitness > self.org.fitness), \
                "Did not get new organism under higher fitness conditions."

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
