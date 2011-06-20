"""Base selection class from which all Selectors should derive.
"""

class AbstractSelection(object):
    """Base class for Selector classes.

    This classes provides useful functions for different selector classes
    and also defines the functions that all selector classes must
    implement.

    This class should not be used directly, but rather should be subclassed.
    """
    def __init__(self, mutator, crossover, repairer = None):
        """Initialize a selector.

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
        self._mutator = mutator
        self._crossover = crossover
        self._repairer = repairer

    def mutate_and_crossover(self, org_1, org_2):
        """Perform mutation and crossover on the two organisms.

        This uses the classes mutator and crossover functions to
        perform the manipulations.

        If a repair class is available, then the rearranged genomes will
        be repaired to make them feasible.
        
        The newly created individuals are returned.
        """
        # first crossover the two organisms
        cross_org_1, cross_org_2 = self._crossover.do_crossover(org_1, org_2)

        # now perform mutation on the two organisms
        final_org_1 = self._mutator.mutate(cross_org_1)
        final_org_2 = self._mutator.mutate(cross_org_2)

        # if we have a repair class, perform repair
        if self._repairer is not None:
            final_org_1 = self._repairer.repair(final_org_1)
            final_org_2 = self._repairer.repair(final_org_2)

        return final_org_1, final_org_2

    def select(self, population):
        raise NotImplementedError("Derived classes must implement.")
