#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Genetic Algorithm classes that provide selection capabilities.
"""
# standard library
import random
import unittest

# biopython
from Bio.Seq import MutableSeq
from Bio.Alphabet import SingleLetterAlphabet

# local stuff
from Bio.GA.Organism import Organism
from Bio.GA.Selection.Diversity import DiversitySelection
from Bio.GA.Selection.Tournament import TournamentSelection
from Bio.GA.Selection.RouletteWheel import RouletteWheelSelection


# --- helper classes and functions

class TestAlphabet(SingleLetterAlphabet):
    """Simple test alphabet.
    """                        
    letters = ["0", "1", "2", "3"]

def test_fitness(genome):
    """Simple class for calculating fitnesses.
    """
    genome_seq = genome.toseq()
    return int(genome_seq.tostring())

class NoSelection:
    """A simple 'selection' class that just returns the generated population.
    """
    def select(self, population):
        return population

class NoMutation:
    """Simple 'mutation' class that doesn't do anything.
    """
    def mutate(self, org):
        return org.copy()

class NoCrossover:
    """Simple 'crossover' class that doesn't do anything.
    """
    def do_crossover(self, org_1, org_2):
        return org_1.copy(), org_2.copy()

class NoRepair:
    """Simple 'repair' class that doesn't do anything.
    """
    def repair(self, org):
        return org.copy()

def random_genome():
    """Return a random genome string.
    """
    alphabet = TestAlphabet()

    new_genome = ""
    for letter in range(3):
        new_genome += random.choice(alphabet.letters)

    return MutableSeq(new_genome, alphabet)

def random_organism():
    """Generate a random organism.
    """
    genome = random_genome()
    return Organism(genome, test_fitness)

# --- the actual test classes

class DiversitySelectionTest(unittest.TestCase):
    """Test selection trying to maximize diversity.
    """
    def setUp(self):
        self.selector = DiversitySelection(NoSelection(), random_genome)

    def test_get_new_organism(self):
        """Getting a new organism not in the new population.
        """
        org = random_organism()
        old_pop = [org]
        new_pop = []

        new_org = self.selector._get_new_organism(new_pop, old_pop)
        self.assertEqual(new_org, org,
                         "Got an unexpected organism %s" % new_org)

    def test_no_retrieve_organism(self):
        """Test not getting an organism already in the new population.
        """
        org = random_organism()
        old_pop = [org]
        new_pop = [org]

        new_org = self.selector._get_new_organism(new_pop, old_pop)
        #assert new_org != org, "Got organism already in the new population."
        #TODO - Why was the above commented out?

    def test_selection(self):
        """Test basic selection on a small population.
        """
        pop = [random_organism() for org_num in range(50)]

        new_pop = self.selector.select(pop)

        self.assertEqual(len(new_pop), len(pop),
                         "Did not maintain population size.")

class TournamentSelectionTest(unittest.TestCase):
    """Test selection based on a tournament style scheme.
    """
    def setUp(self):
        self.selector = TournamentSelection(NoMutation(), NoCrossover(),
                                            NoRepair(), 2)

    def test_select_best(self):
        """Ensure selection of the best organism in a population of 2.
        """
        #Create any two non equal organisms
        org_1 = random_organism()
        while True:
            org_2 = random_organism()
            if org_2.fitness != org_1.fitness:
                break
        #Sort them so org_1 is most fit
        if org_2.fitness > org_1.fitness:
            org_1, org_2 = org_2, org_1
        self.assertTrue(org_1.fitness > org_2.fitness)
        
        pop = [org_1, org_2]
        new_pop = self.selector.select(pop)
        for org in new_pop:
            self.assertEqual(org, org_1,
                             "Got a worse organism selected.")

        #Just to make sure the selector isn't doing something
        #silly with the order, try this with the input reserved:
        pop = [org_2, org_1]
        new_pop = self.selector.select(pop)
        for org in new_pop:
            self.assertEqual(org, org_1,
                             "Got a worse organism selected.")

    def test_selection(self):
        """Test basic selection on a small population.
        """
        pop = [random_organism() for org_num in range(50)]
        new_pop = self.selector.select(pop)

        self.assertEqual(len(new_pop), len(pop),
                         "Did not maintain population size.")


class RouletteWheelSelectionTest(unittest.TestCase):
    """Test selection using a roulette wheel selection scheme.
    """
    def setUp(self):
        self.selector = RouletteWheelSelection(NoMutation(), NoCrossover(),
                                               NoRepair())

    def test_select_best(self):
        """Ensure selection of a best organism in a population of 2.
        """
        worst_genome = MutableSeq("0", TestAlphabet())
        worst_org = Organism(worst_genome, test_fitness)

        better_genome = MutableSeq("1", TestAlphabet())
        better_org = Organism(better_genome, test_fitness)

        new_pop = self.selector.select([worst_org, better_org])
        for org in new_pop:
            self.assertEqual(org, better_org,
                             "Worse organism unexpectly selected.")

    def test_selection(self):
        """Test basic selection on a small population.
        """
        pop = [random_organism() for org_num in range(50)]
        new_pop = self.selector.select(pop)

        self.assertEqual(len(new_pop), len(pop),
                         "Did not maintain population size.")

        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
