"""Perform two-point crossovers between the genomes of two organisms.

This module performs single-point crossover between two genomes.

SinglePointCrossover:
genome 1 --       A B C*D E F
genome 2 --       a b c*d e f

new genome 1 --   A B C d e f
new genome 2 --   a b c D E F

"""
# standard modules
from GeneralPoint import TwoCrossover

class SinglePointCrossover(TwoCrossover):
    """Perform point crossover between genomes at some defined rate.

    This performs a crossover between two genomes at some defined 
    frequency.  Length of genome is preserved, as the crossover 
    point is the same for either genome.
    """
    def __init__(self, crossover_prob = .1):
        """Initialize to do crossovers at the specified probability.
        """
        TwoCrossover.__init__(self, 1, crossover_prob)
