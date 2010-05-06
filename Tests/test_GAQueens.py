#!/usr/bin/env python
"""Try out the N-queens problem for an arbitrary number of queens.

This program uses Genetic Algorithms to try to solve the N queens
problem, in which you place N different queens on an N by N chess
board such that no two queens will be attacking each other.

We represent queens on the board as a tuple like (1, 2, 3, 4, 5)
which would be 5 queens diaganol across the board.

Usage:
python test_GAQueens.py <Number of Queens to place>

where <Number of Queens to place> is just a number specifying how many
queens you want to try to calculate this for.

When called as part of the Biopython unit test suite, 5 queens are used.
"""
# standard library
import sys
import math
import random
import copy
import time

# Biopython
from Bio import Alphabet

# Genetic Algorithm stuff
from Bio.GA.Evolver import GenerationEvolver
from Bio.GA import Organism
from Bio.GA.Mutation.Simple import ConversionMutation
from Bio.GA.Crossover.Point import SinglePointCrossover
from Bio.GA.Selection.RouletteWheel import RouletteWheelSelection
from Bio.GA.Selection.Tournament import TournamentSelection

VERBOSE = 0


def main(num_queens):

    print "Calculating for %s queens..." % num_queens

    num_orgs = 1000
    print "Generating an initial population of %s organisms..." % num_orgs
    queen_alphabet = QueensAlphabet(num_queens)

    start_population = Organism.random_population(queen_alphabet, num_queens,
                                                  num_orgs, queens_fitness)

    print "Evolving the population and searching for a solution..."

    mutator = QueensMutation(mutation_rate = 0.05)
    crossover = QueensCrossover(queens_fitness, crossover_prob = .2,
                                max_crossover_size = 4)
    repair = QueensRepair()
    # rw_selector = RouletteWheelSelection(mutator, crossover, repair)
    t_selector = TournamentSelection(mutator, crossover, repair, 5)

    start_time = time.ctime(time.time())
    evolver = GenerationEvolver(start_population, t_selector)
    evolved_pop = evolver.evolve(queens_solved)
    end_time = time.ctime(time.time())

    unique_solutions = []
    for organism in evolved_pop:
        if organism.fitness == num_queens:
            if organism not in unique_solutions:
                unique_solutions.append(organism)

    if VERBOSE:
        print "Search started at %s and ended at %s" % (start_time, end_time)
        for organism in unique_solutions:
            print "We did it!", organism
            display_board(organism.genome)
        

def display_board(genome):
    """Display a genome in the N-queens problem.

    Inspired by the display function in the queens.py solution to the N-queens
    problem in the Python demo scripts.
    """
    print '+-' + '--'*len(genome) + '+'

    for row in range(len(genome)):
        print '|',
        for genome_item in genome:
            if genome_item == row:
                print 'Q',
            else:
                print '.',
        print '|'

    print '+-' + '--'*len(genome) + '+'

def queens_solved(organisms):
    """Determine if we have solved the problem.

    We just search through the population for an organism that has a
    fitness that is equal to the number of queens in the population.
    If so, we have a solution, otherwise we need to keep looking.
    """
    for org in organisms:
        if org.fitness == len(org.genome):
            return 1

    # if we got here we didn't do it
    return 0
     
def queens_fitness(genome):
    """Calculate the fitness of an organization of queens on the chessboard.

    Arguments:

    o genome -- A MutableSeq object specifying an organism genome.

    The number returned is the number of unattacked queens on the board.
    """
    fitness = 0

    # check each queen on the board
    for check_queen_col in range(len(genome)):
        is_attacked = 0
        # check against all other queens on the board
        for other_queen_col in range(len(genome)):
            # only check a queen if it isn't exactly the same queen
            if check_queen_col != other_queen_col:
                # get the row for the two queens we are comparing
                check_queen_row = int(genome[check_queen_col])
                other_queen_row = int(genome[other_queen_col])
                
                # a queen is attacked if it is in a row with another queen
                if check_queen_row == other_queen_row:
                    is_attacked = 1
                    break
                # or it is attacked if it is diaganol to another queen
                elif (abs(check_queen_row - other_queen_row) ==
                      abs(check_queen_col - other_queen_col)):
                    is_attacked = 1
                    break

        if not(is_attacked):
            fitness += 1

    return fitness

class QueensAlphabet(Alphabet.Alphabet):
    def __init__(self, num_queens):
        """Initialize with the number of queens we are calculating for.
        """
        # set up the letters for the alphabet
        assert 0 <= num_queens <= 9
        self.letters = "".join(str(i) for i in range(num_queens))

# --- Problem specific crossover, mutation and repair operations
class QueensRepair:
    """A repair function to help create correct N-Queens solutions.

    This attempts to help generate correct solutions by offering some
    amount of repair to remove queens that are located in the same rows.
    After repair, a sequence should have no queens in the same row.

    So, if you start with something infeasible like (1, 2, 2, 3, 3, 4),
    after running it through repair you'll get a feasible individual
    like (1, 2, 5, 3, 6, 4). This should greatly reduce the number of
    individuals that need to be searched through in a population.
    """
    def __init__(self, repair_prob = 1):
        """Initialize the repairer.

        Arguments:

        o repair_prob -- The probability that we'll repair a genome.
        By default, we always repair.
        """
        self._repair_prob = repair_prob

    def _get_duplicates(self, genome):
        """Return all of the letters in the genome that are duplicated.

        This checks every letter in the genome (which are the rows of
        the chessboard, in this case), and adds them to a list of duplicated
        items if there is more than one of them, and then returns this list.
        """
        duplicates = []
        for item in genome.alphabet.letters:
            if genome.count(str(item)) > 1:
                duplicates.append(item)

        return duplicates

    def _get_unused(self, genome):
        """Return all of the letters in the genome which are unused.

        This checks the letters in the genome (which are th rows on the
        chessboard) and returns all items which are not used.
        """
        unused = []
        for item in genome.alphabet.letters:
            if genome.count(str(item)) == 0:
                unused.append(item)

        return unused

    def repair(self, organism):
        """Repair the specified genome to make it feasible.

        Arguments:

        o organism -- The Organism object we are going to perform the
        repair on.
        """
        # check if we should repair or not
        repair_chance = random.random()
        if repair_chance <= self._repair_prob:
            while 1:
                # get the duplicated items we need to work on
                duplicated_items = self._get_duplicates(organism.genome)

                if len(duplicated_items) == 0:
                    break

                # take the first duplicated element, and convert it to
                # a row that is not already taken
                duplicated_pos = organism.genome.index(duplicated_items[0])

                free_rows = self._get_unused(organism.genome)
                assert len(free_rows) > 0, "Unexpected lack of empty rows"

                new_item = random.choice(free_rows)
                organism.genome[duplicated_pos] = new_item

        return organism
        
class QueensCrossover:
    """Crossover operation to help in solving the N-Queens problem.

    This tries to perform smarter crossovers by picking out regions of
    the genome that have high fitness.

    It scans through both genomes in the crossover with a window half the
    size of the genome, and finds the region with the highest fitness in
    both genomes. It then recombines these high fitness windows to form
    the new genome that is returned.
    """
    def __init__(self, fitness_func, crossover_prob = .1,
                 max_crossover_size = 4):
        """Initialize to do N-Queens optimized crossover.

        Arguments:

        o fitness_func -- A function that can calculate the fitness of
        a genome.

        o crossover_prob -- The probability of having a crossover
        between two passed in organisms.

        o max_crossover_size -- The maximum crossover size of the 'best' region
        to search for.
        """
        self._crossover_prob = crossover_prob
        self._fitness_calc = fitness_func
        self._max_crossover_size = max_crossover_size

    def do_crossover(self, org_1, org_2):
        """Perform a crossover between two organisms.
        """
        new_org_1 = org_1.copy()
        new_org_2 = org_2.copy()
        
        # determine if we have a crossover
        crossover_chance = random.random()
        if crossover_chance <= self._crossover_prob:
            # find the region of highest probability in both orgs
            best_1, rest_1 = self._find_best_region(new_org_1.genome,
                                                    make_best_larger = 1)
            best_2, rest_2 = self._find_best_region(new_org_2.genome,
                                                    make_best_larger = 0)

            assert len(best_1) + len(best_2) == len(rest_1) + len(rest_2), \
                   "Did not preserve genome length!"
            
            new_org_1.genome = best_1 + best_2
            new_org_2.genome = rest_1 + rest_2

        return new_org_1, new_org_2

    def _find_best_region(self, genome, make_best_larger = 1):
        """Find the best region in the given genome.

        Arguments:

        o genome -- A MutableSeq object specifying the genome of an organism

        o make_best_larger -- A flag to determine whether the best region
        we should search for should be the larger region of the split
        caused by crossover or the smaller region. This makes it easy
        to split two genomes, recombine them, and get a solution that
        makes sense.

        Returns:
        o Two MutableSeq objects. They are both half of the size of the passed
        genome. The first is the highest fitness region of the genome and the
        second is the rest of the genome.
        """
        first_region = max(len(genome) / 2, self._max_crossover_size)
        second_region = len(genome) - first_region
        
        if make_best_larger:
            region_size = max(first_region, second_region)
        else:
            region_size = min(first_region, second_region)

        # loop through all of the segments and find the best fitness segment
        
        # represent best_fitness as a three tuple with the coordinates of
        # the start and end as the first two elements, and the fitness of
        # the region as the last element. Start with a value that
        # will overridden right away
        best_fitness = [0, 0, -1]
        for start_index in range(len(genome) - region_size):
            region_fitness = \
             self._fitness_calc(genome[start_index: start_index + region_size])

            if region_fitness > best_fitness[2]:
                best_fitness = [start_index, start_index + region_size,
                                region_fitness]

        # get the the two regions and return 'em
        best_region = genome[best_fitness[0]:best_fitness[1]]
        rest_region = genome[0:best_fitness[0]] + genome[best_fitness[1]:]

        return best_region, rest_region
            

class QueensMutation:
    """Mutation operation to help in the N-Queens problem.

    This performs mutation, but instead of randomly mutating a single
    item to any other, it tries to mutate it to a row that is not already
    taken at some other position in the genome. This thus tries to
    generate more 'correct' mutations that will help achieve the solution.
    """
    def __init__(self, mutation_rate = 0.001):
        """Inititialize a mutator.

        Arguments:

        o mutation_rate -- The change of a mutation happening at any
        position in the genome.
        """
        self._mutation_rate = mutation_rate

    def mutate(self, organism):
        """Mutate the genome trying to put in 'helpful' mutations.
        """
        new_org = organism.copy()
        gene_choices = list(new_org.genome.alphabet.letters)

        # potentially mutate any gene in the genome
        for gene_index in range(len(new_org.genome)):
            mutation_chance = random.random()
            # if we have a mutation
            if mutation_chance <= self._mutation_rate:
                # find only choices that are not already taken elsewhere
                # in the genome
                gene_choices = list(new_org.genome.alphabet.letters)

                for gene in new_org.genome:
                    if gene in gene_choices:
                        gene_choices.remove(gene)

                # if there are no choices left, we are stuck going for random
                if len(gene_choices) == 0:
                    gene_choices = list(new_org.genome.alphabet.letters)
                
                # get a new letter with the left-over choices
                new_letter = random.choice(gene_choices)
                new_org.genome[gene_index] = new_letter

        return new_org
 
num_queens = 5

if __name__ == "__main__":
    if len(sys.argv) == 2:
        num_queens = int(sys.argv[1])
    elif len(sys.argv) > 2:
        print "Usage:"
        print "python test_GAQueens.py <Number of Queens to place>\n"
        print "where <Number of Queens to place> is an optional parameter"
        print "specifying how many queens you want to try to calculate"
        print "this for. The default number of queens to place is 5."
        sys.exit(1)

main(num_queens)
