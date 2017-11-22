# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Perform Simple mutations on an organism's genome."""
# standard modules
import random

from Bio._py3k import range


class SinglePositionMutation(object):
    """Perform a conversion mutation, but only at a single point in the genome.

    This does not randomize the genome as much as ConversionMutation, since
    only one change is allowed per genome at the specified mutation rate.
    """

    def __init__(self, mutation_rate=0.001):
        """Initialize a mutator.

        Arguments:
         - mutation_rate - The chance of a mutation happening once in the
           genome.

        """
        self._mutation_rate = mutation_rate
        # Originally there was an RNG for each case. Python 3.6 was
        # released with a bug under Windows where all RNGs created
        # within a 15ms window have the same initial seed.
        # To reduce possible correlations, re-use a single RNG.
        rng = random.Random()
        # a random number generator to test if we have a mutation
        self._mutation_rand = rng
        # a random number generator to switch to a new alphabet letter
        self._switch_rand = rng
        # a random number generator to find the mutation position
        self._pos_rand = rng

    def mutate(self, organism):
        """Mutate the organism's genome."""
        mutated_org = organism.copy()
        gene_choices = mutated_org.genome.alphabet.letters

        mutation_chance = self._mutation_rand.random()
        if mutation_chance <= self._mutation_rate:
            # pick a gene position to mutate at
            mutation_pos = \
                         self._pos_rand.choice(list(range(len(mutated_org.genome))))

            # get a new letter to replace the position at
            new_letter = self._switch_rand.choice(gene_choices)

            mutated_org.genome[mutation_pos] = new_letter

        return mutated_org


class ConversionMutation(object):
    """Potentially mutate any item to another in the alphabet.

    This just performs switching mutation -- changing one gene of a genome
    to any other potential gene, at some defined frequency. If the organism
    is determined to mutate, then the alphabet item it is equally likely
    to switch to any other letter in the alphabet.
    """

    def __init__(self, mutation_rate=0.001):
        """Inititialize a mutator.

        Arguments:
         - mutation_rate -- The chance of a mutation happening at any
           position in the genome.

        """
        self._mutation_rate = mutation_rate
        # a random number generator to test if we have a mutation
        self._mutation_rand = random.Random()
        # a random number generator to switch to a new alphabet letter
        self._switch_rand = random.Random()

    def mutate(self, organism):
        """Mutate the organism's genome."""
        mutated_org = organism.copy()

        gene_choices = mutated_org.genome.alphabet.letters

        # potentially mutate any gene in the genome
        for gene_index in range(len(mutated_org.genome)):
            mutation_chance = self._mutation_rand.random()
            # if we have a mutation
            if mutation_chance <= self._mutation_rate:
                # get a new letter
                new_letter = self._switch_rand.choice(gene_choices)
                mutated_org.genome[gene_index] = new_letter

        return mutated_org
