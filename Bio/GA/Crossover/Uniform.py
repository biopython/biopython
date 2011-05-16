"""Perform uniform crossovers between the genomes of two organisms.


genome 1 --       A B C D E F G
                  .   . .   .
genome 2 --       a b c d e f g

After crossover:

new genome 1 --  a B c d E f G
new genome 2 --  A b C D e F g

Uniform Crossover is a standard crossover technique for
rapid mutation-behavior.  
"""
# standard modules
import random

class UniformCrossover(object):
    """Perform single point crossover between genomes at some defined rates.

    This performs a single crossover between two genomes at some
    defined frequency. The location of the crossover is chosen randomly
    if the crossover meets the probability to occur.
    """
    def __init__(self, crossover_prob = .1, uniform_prob = 0.7):
        """Initialize to do uniform crossover at the specified probability and frequency.
        """
        self._crossover_prob = crossover_prob
        self._uniform_prob   = uniform_prob
        return
        
    def do_crossover(self, org_1, org_2):
        """Potentially do a crossover between the two organisms.
        """
        new_org_1 = org_1.copy()
        new_org_2 = org_2.copy()
        
        # determine if we have a crossover
        crossover_chance = random.random()
        if crossover_chance <= self._crossover_prob:
            minlen = min(len(new_org_1.genome),len(new_org_2.genome))
            for i in range( minlen ):
                uniform_chance = random.random()
                if uniform_chance <= self._uniform_prob:
                    # cycle element
                    temp = new_org_1.genome[i]
                    new_org_1.genome[i] = new_org_2.genome[i]
                    new_org_2.genome[i] = temp
            
        return new_org_1, new_org_2
