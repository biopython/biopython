# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Select individuals into a new population trying to maintain diversity.

This selection mechanism seeks to try and get both high fitness levels
and high diversity levels in the population.
"""
# standard modules
import random
import math

# biopython
from Bio.Seq import MutableSeq

# local modules
from .Abstract import AbstractSelection
from .Tournament import TournamentSelection


class DiversitySelection(AbstractSelection):
    """Implement diversity selection.

    Diversity selection is performed by trying to select individuals
    from the population that aren't already in the new_population. A group
    of selected individuals is then subjected to selection using
    a passed selection routine.

    If new individuals can not be selected, new individuals will be
    randomly generated and inserted into the population.
    """
    def __init__(self, internal_selector, genome_generator):
        """Initialize a diversity selector.

        Arguments:

        o internal_selector - A selection object that will be used to select
        individuals based on fitness, perform crossover, mutation and repair.

        o genome_generator - A function that, when called, will return a
        genome to be used for a new organism. The genome returned must
        be a MutableSeq() object.
        """
        self._internal_selector = internal_selector
        self._genome_generator = genome_generator

        self.sub_pop_percent = .1
        self.random_tries = 10

    def _get_new_organism(self, new_pop, old_pop):
        """Get a new organism from old_pop that isn't in new_pop.

        This attempts to select an organism from old_pop that isn't in
        new_pop. If we can't do this in the number of tries specified
        by the class attribute random_tries, we generate a new random
        organism and return that.
        """
        # try to pick an organism that isn't in the population
        new_org = None
        num_tries = 0
        while new_org is None and num_tries < self.random_tries:
            chance_org = random.choice(old_pop)

            if chance_org not in new_pop:
                new_org = chance_org

            num_tries += 1

        # if we don't get an organism, generate a random one
        if new_org is None:
            new_org = old_pop[0].copy()
            random_genome = self._genome_generator()
            new_org.genome = random_genome
            new_org.recalculate_fitness()

        return new_org

    def select(self, population):
        """Perform selection on the current population, encouraging diversity.
        """
        new_population = []

        while len(new_population) < len(population):
            # generate a sub population
            sub_pop_size = int(math.ceil(len(population) *
                                         self.sub_pop_percent))
            sub_pop = []
            for individual in range(sub_pop_size):
                new_org = self._get_new_organism(new_population, population)
                sub_pop.append(new_org)

            # put the new sub population through selection, mutation
            # and all of that good stuff
            new_sub_pop = self._internal_selector.select(sub_pop)

            new_population.extend(new_sub_pop)

        # return the new population, which should have the same number
        # of individuals we started with.
        return new_population[:len(population)]
