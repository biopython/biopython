"""Perform point crossovers between the genomes of two organisms.

These modules perform simple crossover between two genomes like:

genome 1 --       A B C D
                     |
                     |
genome 2 --       a b c d

After crossover:

new genome 1 --  A B c d
new genome 2 --  a b C D

This is just a demonstration of single point crossover, but more complicated
cases are possible.
"""
# standard modules
import random

class SinglePointCrossover:
    """Perform single point crossover between genomes at some defined rates.

    This performs a single crossover between two genomes at some
    defined frequency. The location of the crossover is chosen randomly
    if the crossover meets the probability to occur.
    """
    def __init__(self, crossover_prob = .1):
        """Initialize to do crossovers at the specified probability.
        """
        self._crossover_prob = crossover_prob
        # a random number generator to test if we have a crossover
        self._crossover_rand = random.Random()
        # a random number generator to pick a crossover location
        self._location_rand = random.Random()

    def do_crossover(self, org_1, org_2):
        """Potentially do a crossover between the two organisms.
        """
        new_org_1 = org_1.copy()
        new_org_2 = org_2.copy()
        
        # determine if we have a crossover
        crossover_chance = self._crossover_rand.random()
        if crossover_chance <= self._crossover_prob:
            # choose a place to do the crossover
            crossover_loc = \
              self._location_rand.choice(range(len(new_org_1.genome)))

            # temporarily get the part of the first genome, we'll
            # need to give it to genome 2 in a second
            temp_genome = new_org_1.genome[crossover_loc:]

            # transfer the part of the genome from 2 to 1
            new_org_1.genome[crossover_loc:] = new_org_2.genome[crossover_loc:]

            # now give genome 2 the stuff from genome 1
            new_org_2.genome[crossover_loc:] = temp_genome

        return new_org_1, new_org_2


        
        

         
