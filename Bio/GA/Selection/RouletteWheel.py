"""Implement Roulette Wheel selection on a population.

This implements Roulette Wheel selection in which individuals are
selected from a population randomly, with their proportion of selection
based on their relative fitness in the population.
"""
# standard modules
import random
import copy

# local modules
from Abstract import AbstractSelection

class RouletteWheelSelection(AbstractSelection):
    """Roulette wheel selection proportional to individuals fitness.

    The implements a roulette wheel selector that selects individuals
    from the population, and performs mutation and crossover on
    the selected individuals.
    """
    def __init__(self, mutator, crossover, repairer = None):
        """Initialize the selector.

        Arguments:

        o mutator -- A Mutation object which will perform mutation
        on an individual.

        o crossover -- A Crossover object which will take two
        individuals and produce two new individuals which may
        have had crossover occur.

        o repairer -- A class which can do repair on rearranged genomes
        to eliminate infeasible individuals. If set at None, so repair
        will be done.
        """
        AbstractSelection.__init__(self, mutator, crossover, repairer)

    def select(self, population):
        """Perform selection on the population based using a Roulette model.

        Arguments:

        o population -- A population of organisms on which we will perform
        selection. The individuals are assumed to have fitness values which
        are due to their current genome.
        """
        # set up the current probabilities for selecting organisms
        # from the population
        prob_wheel = self._set_up_wheel(population)
        probs = prob_wheel.keys()
        probs.sort()
        
        # now create the new population with the same size as the original
        new_population = []

        for pair_spin in range(len(population) // 2):
            # select two individuals using roulette wheel selection
            choice_num_1 = random.random()
            choice_num_2 = random.random()

            # now grab the two organisms from the probabilities
            chosen_org_1 = None
            chosen_org_2 = None
            prev_prob = 0
            for cur_prob in probs:
                if choice_num_1 > prev_prob and choice_num_1 <= cur_prob:
                    chosen_org_1 = prob_wheel[cur_prob]
                if choice_num_2 > prev_prob and choice_num_2 <= cur_prob:
                    chosen_org_2 = prob_wheel[cur_prob]

                prev_prob = cur_prob

            assert chosen_org_1 is not None, "Didn't select organism one"
            assert chosen_org_2 is not None, "Didn't select organism two"

            # do mutation and crossover to get the new organisms
            new_org_1, new_org_2 = self.mutate_and_crossover(chosen_org_1,
                                                             chosen_org_2)
            
            new_population.extend([new_org_1, new_org_2])

        return new_population

    def _set_up_wheel(self, population):
        """Set up the roulette wheel based on the fitnesses.

        This creates a fitness proportional 'wheel' that will be used for
        selecting based on random numbers.

        Returns:
        
        o A dictionary where the keys are the 'high' value that an
        individual will be selected. The low value is determined by
        the previous key in a sorted list of keys. For instance, if we
        have a sorted list of keys like:

        [.1, .3, .7, 1]

        Then the individual whose key is .1 will be selected if a number
        between 0 and .1 is chosen, the individual whose key is .3 will
        be selected if the number is between .1 and .3, and so on.

        The values of the dictionary are the organism instances.
        """
        # first sum up the total fitness in the population
        total_fitness = 0
        for org in population:
            total_fitness += org.fitness

        # now create the wheel dictionary for all of the individuals
        wheel_dict = {}
        total_percentage = 0
        for org in population:
            org_percentage = float(org.fitness) / float(total_fitness)

            # the organisms chance of being picked goes from the previous
            # percentage (total_percentage) to the previous percentage
            # plus the organisms specific fitness percentage
            wheel_dict[total_percentage + org_percentage] = copy.copy(org)

            # keep a running total of where we are at in the percentages
            total_percentage += org_percentage

        return wheel_dict
