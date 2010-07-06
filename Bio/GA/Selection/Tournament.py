"""Provide Tournament style selection.

This implements selection based on a tournament style. In this model of
selection, two individuals are randomly chosen from the population, and
the organism with the higher fitness is considered the 'winner' and moves
to the next generation.
"""
# standard modules
import random

# local modules
from Abstract import AbstractSelection

class TournamentSelection(AbstractSelection):
    """Implement tournament style selection.
    """
    def __init__(self, mutator, crossover, repairer, num_competitors = 2):
        """Initialize the tournament selector.

        Arguments:

        o num_competitors-- The number of individiuals that should be
        involved in a selection round. By default we just have two
        individuals (head to head!).

        See AbstractSelection for a description of the arguments to
        the initializer.
        """
        AbstractSelection.__init__(self, mutator, crossover, repairer)

        if num_competitors < 2:
            raise ValueError("Must have at least 2 competitors!")
        
        self._num_competitors = num_competitors

    def select(self, population):
        """Perform selection on the population using the Tournament model.

        Arguments:

        o population -- A population of organisms on which we will perform
        selection. The individuals are assumed to have fitness values which
        are due to their current genome (ie. the fitness is up to date).
        """
        # we want to create a new population of the same size as the original
        new_population = []

        while len(new_population) < len(population):
            # select two individuals using tournament selection
            new_orgs = []
            # select for two individuals
            for round_num in range(2):
                competitors = []
                while len(competitors) < self._num_competitors:
                    new_org = random.choice(population)
                    if new_org not in competitors:
                        competitors.append(new_org)
                                    
                # sort the competitors by fitness, this will put them
                # from lowest to highest
                competitors.sort(key = lambda org: org.fitness)

                # get the best organism
                new_orgs.append(competitors[-1])

            assert len(new_orgs) == 2, "Expected two organisms to be selected"

            # do mutation and crossover to get the new organisms
            new_org_1, new_org_2 = self.mutate_and_crossover(new_orgs[0],
                                                             new_orgs[1])

            new_population.extend([new_org_1, new_org_2])

        return new_population
                
