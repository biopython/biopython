"""General functionality for mutations.
"""
# standard library
import random

# local stuff
from Bio.GA.Organism import Organism

class SafeFitnessMutation(object):
    """Perform mutations, but do not allow decreases in organism fitness.

    This doesn't actually do any mutation work, but just checks that
    newly create organisms do not have lower fitnesses. 
    """
    def __init__(self, actual_mutation, accept_less = 0.0):
        """Initialize to do safe mutations

        Arguments:

        o actual_mutation - A Mutation class which actually implements
        mutation. functionality.

        o accept_less - A probability to accept mutations which
        generate lower fitness. This allows you to accept some
        crossovers which reduce fitness, but not all of them.
        """
        self._mutation = actual_mutation
        self._accept_less_percent = accept_less
        self._accept_less_rand = random.Random()

    def mutate(self, org):
        """Perform safe mutation of the specified organism.
        """
        new_org = self._mutation.mutate(org)
        new_org.recalculate_fitness()

        if org.fitness > new_org.fitness:
            accept_less_chance = self._accept_less_rand.random()
            if accept_less_chance <= self._accept_less_percent:
                return new_org
            else:
                return org
        else:
            return new_org
