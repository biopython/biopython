# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Methods for performing repairs that will Stabilize genomes.

These methods perform repair to keep chromosomes from drifting too far in
any direction -- ie. bring them back to a stabilizing center. This may be
useful in cases where fitness functions alone won't keep chromosomes in
check.
"""
# standard library
import random


class AmbiguousRepair(object):
    """Perform repair to reduce the number of Ambiguous genes in a genome.

    In cases where ambiguous genes are allowed in a genome (for example,
    where you have a wild card character like '*' that will match
    anything), these can come to dominate a genome since, really, the
    best fitness is someting like '*******'. This repair protects against
    that by changing ambiguous characters into some non-ambiguous gene.
    """
    def __init__(self, ambig_finder, num_ambiguous):
        """Initialize the repair class.

        Arguments:

        o ambig_finder - A class implementing the function find_ambiguous
        which will return a list of all ambiguous positions in a sequence.
        It also must have the function all_unambiguous, which will return
        all allowed unambiguous letters.

        o num_ambiguous - The minimum number of ambiguous items that are
        allowed in a genome. If there are more than this present, repair
        will be performed.
        """
        self._ambig_finder = ambig_finder
        self._num_ambiguous = num_ambiguous
        self._alphabet_letters = ambig_finder.all_unambiguous()

    def repair(self, organism):
        """Perform a repair to remove excess ambiguous genes.
        """
        new_org = organism.copy()

        # start getting rid of ambiguous items
        while True:
            # first find all of the ambigous items
            seq_genome = new_org.genome.toseq()
            all_ambiguous = self._ambig_finder.find_ambiguous(str(seq_genome))

            # if we have less then the number of ambiguous allowed, stop
            if len(all_ambiguous) <= self._num_ambiguous:
                break

            # remove an ambiguous item and replace it with a non-ambiguous
            to_change = random.choice(all_ambiguous)
            new_gene = random.choice(self._alphabet_letters)
            new_org.genome[to_change] = new_gene

        return new_org
