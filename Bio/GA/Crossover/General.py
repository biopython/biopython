"""General functionality for crossover that doesn't apply.

This collects Crossover stuff that doesn't deal with any specific
type of crossover.
"""
# standard library
import random

# local stuff
from Bio.GA.Organism import Organism

class SafeFitnessCrossover(object):
    """Perform crossovers, but do not allow decreases in organism fitness.

    This doesn't actually do any crossover work, but instead relies on
    another class to do the crossover and just checks that newly created
    organisms do not have less fitness. This is useful for cases where
    crossovers can 
    """
    def __init__(self, actual_crossover, accept_less = 0.0):
        """Initialize to do safe crossovers.

        Arguments:

        o actual_crossover - A Crossover class which actually implements
        crossover functionality.

        o accept_less - A probability to accept crossovers which
        generate less fitness. This allows you to accept some
        crossovers which reduce fitness, but not all of them.
        """
        self._crossover = actual_crossover
        self._accept_less_percent = accept_less
        self._accept_less_rand = random.Random()

    def do_crossover(self, org_1, org_2):
        """Perform a safe crossover between the two organism.
        """
        new_org_1, new_org_2 = self._crossover.do_crossover(org_1, org_2)

        return_orgs = []

        for start_org, new_org in ((org_1, new_org_1),
                                   (org_2, new_org_2)):
            new_org.recalculate_fitness()

            # if the starting organism has a better fitness,
            # keep it, minding the acceptance of less favorable change policy
            if start_org.fitness > new_org.fitness:
                accept_change = self._accept_less_rand.random()
                if accept_change <= self._accept_less_percent:
                    return_orgs.append(new_org)
                else:
                    return_orgs.append(start_org)
            else:
                return_orgs.append(new_org)

        assert len(return_orgs) == 2, "Should have two organisms to return."

        return return_orgs
