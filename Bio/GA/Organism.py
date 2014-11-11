# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Deal with an Organism in a Genetic Algorithm population.
"""
# standard modules
import sys  # for Python 3 hack
import random
import array

# Sequence objects from Biopython
from Bio.Seq import MutableSeq


def function_population(new_genome, num_organisms, fitness_calculator):
    """Generate a population given a function to create genomes

    Arguments:

    o new_genome - A function or callable object that will return
    a genome that can be used for a new organism. This new genome
    should be a MutableSeq object with a specified alphabet.

    o num_organisms - The number of individuals we want in the population.

    o fitness_calculator -- A function that will calculate the fitness
    of the organism when given the organisms genome.
    """
    all_orgs = []

    for org_num in range(num_organisms):
        cur_genome = new_genome()
        all_orgs.append(Organism(cur_genome, fitness_calculator))

    return all_orgs


def random_population(genome_alphabet, genome_size, num_organisms,
                      fitness_calculator):
    """Generate a population of individuals with randomly set genomes.

    Arguments:

    o genome_alphabet -- An Alphabet object describing all of the
    possible letters that could potentially be in the genome of an
    organism.

    o genome_size -- The size of each organisms genome.

    o num_organism -- The number of organisms we want in the population.

    o fitness_calculator -- A function that will calculate the fitness
    of the organism when given the organisms genome.
    """
    all_orgs = []

    # a random number generator to get letters for the genome
    letter_rand = random.Random()

    # figure out what type of characters are in the alphabet
    if isinstance(genome_alphabet.letters[0], str):
        if sys.version_info[0] == 3:
            alphabet_type = "u"  # Use unicode string on Python 3
        else:
            alphabet_type = "c"  # Use byte string on Python 2
    elif isinstance(genome_alphabet.letters[0], int):
        alphabet_type = "i"
    elif isinstance(genome_alphabet.letters[0], float):
        alphabet_type = "d"
    else:
        raise ValueError(
            "Alphabet type is unsupported: %s" % genome_alphabet.letters)

    for org_num in range(num_organisms):
        new_genome = MutableSeq(array.array(alphabet_type), genome_alphabet)

        # generate the genome randomly
        for gene_num in range(genome_size):
            new_gene = letter_rand.choice(genome_alphabet.letters)
            new_genome.append(new_gene)

        # add the new organism with this genome
        all_orgs.append(Organism(new_genome, fitness_calculator))

    return all_orgs


class Organism(object):
    """Represent a single individual in a population.

    Attributes:

    o genome -- The genome of the organism. This is a Bio.MutableSeq
    object that has the sequence of the genome, and the alphabet
    describing all elements that can be a part of the genome.

    o fitness -- The calculate fitness of the organism. This fitness is
    based on the last time it was calculated using the fitness_calculator.
    So... the fitness could potentially be out of date with the real genome
    if you are not careful to recalculate it after changes with
    recalculate_fitness()
    """
    def __init__(self, genome, fitness_calculator, start_fitness=None):
        """Initialize an organism

        Arguments:

        o genome -- A MutableSeq object representing the sequence of the
        genome.

        o fitness_calculator -- A function that will calculate the fitness
        of the organism when given the organisms genome.

        o start_fitness - the starting fitness corresponding with the
        given genome. If not supplied, the fitness will be calculated
        using fitness_calculator.
        """
        assert isinstance(genome, MutableSeq), "Genome must be a MutableSeq"

        self.genome = genome
        self._fitness_calc = fitness_calculator

        # calculate the fitness of the genome
        if start_fitness is None:
            self.fitness = self._fitness_calc(self.genome)
        else:
            self.fitness = start_fitness

    def __str__(self):
        """Provide a string output for debugging.
        """
        return "Genome: %s; Fitness %s" % (str(self.genome), self.fitness)

    def __eq__(self, other):
        """Compare organisms by their genomes (as strings of letters).
        """
        # See Bio/Seq.py and the comments there about shifting to
        # using simple string equality. Previously Seq objects used
        # object equality, while MutableSeq objects used alphabet
        # aware string equality.
        return str(self.genome) == str(other.genome)

    def __ne__(self, other):
        """Compare organisms by their genomes (as strings of letters).
        """
        return str(self.genome) != str(other.genome)

    def __lt__(self, other):
        """Compare organisms by their genomes (as strings of letters).
        """
        return str(self.genome) < str(other.genome)

    def __le__(self, other):
        """Compare organisms by their genomes (as strings of letters).
        """
        return str(self.genome) <= str(other.genome)

    def __gt__(self, other):
        """Compare organisms by their genomes (as strings of letters).
        """
        return str(self.genome) > str(other.genome)

    def __ge__(self, other):
        """Compare organisms by their genomes (as strings of letters).
        """
        return str(self.genome) >= str(other.genome)

    def copy(self):
        """Return a copy of the organism.

        This makes it easy to duplicate an organism before changing it.
        """
        copy_genome = self.genome[:]
        return Organism(copy_genome, self._fitness_calc, self.fitness)

    def recalculate_fitness(self):
        """Calculate and reset the fitness of the current genome

        This should be called after the genome is updated to ensure that
        fitness always stays in sync with the current genome.
        """
        self.fitness = self._fitness_calc(self.genome)
