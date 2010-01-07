"""Perform two-point crossovers between the genomes of two organisms.

This module performs two-point crossover between two genomes.
There are two flavors: OnePointCrossover (Point) and TwoPointCrossover.

TwoPointCrossover is the minimal crossover technique that
facilitates diverse genome length.  Do not use this if you need to
maintain consistent genome length.

TwoPointCrossover:
genome 1 --       A B*C D E F
genome 2 --       a b c*d e f

new genome 1 --   A B d e f
new genome 2 --   a b c C D E F

"""
# standard modules
from GeneralPoint import TwoCrossover

class TwoPointCrossover(TwoCrossover):
    """Perform two point crossover between genomes at some defined rate.

    This performs a crossover between two genomes at some defined frequency.
    The location of the points of crossover are chosen randomly if the
    crossover meets the probability to occur.  
    """
    def __init__(self, crossover_prob = .1):
        """Initialize to do crossovers at the specified probability.
        """
        TwoCrossover.__init__(self, 2, crossover_prob)
