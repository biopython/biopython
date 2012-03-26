"""Deal with Motifs or Signatures allowing ambiguity in the sequences.

This class contains Schema which deal with Motifs and Signatures at
a higher level, by introducing `don't care` (ambiguity) symbols into
the sequences. For instance, you could combine the following Motifs:

'GATC', 'GATG', 'GATG', 'GATT'

as all falling under a schema like 'GAT*', where the star indicates a
character can be anything. This helps us condense a whole ton of
motifs or signatures.
"""
# standard modules
import random
import re

# biopython
from Bio import Alphabet
from Bio.Seq import MutableSeq

# neural network libraries
from Pattern import PatternRepository

# genetic algorithm libraries
from Bio.GA import Organism
from Bio.GA.Evolver import GenerationEvolver
from Bio.GA.Mutation.Simple import SinglePositionMutation
from Bio.GA.Crossover.Point import SinglePointCrossover
from Bio.GA.Repair.Stabilizing import AmbiguousRepair
from Bio.GA.Selection.Tournament import TournamentSelection
from Bio.GA.Selection.Diversity import DiversitySelection

class Schema(object):
    """Deal with motifs that have ambiguity characters in it.

    This motif class allows specific ambiguity characters and tries to
    speed up finding motifs using regular expressions.

    This is likely to be a replacement for the Schema representation,
    since it allows multiple ambiguity characters to be used.
    """
    def __init__(self, ambiguity_info):
        """Initialize with ambiguity information.

        Arguments:
        
        o ambiguity_info - A dictionary which maps letters in the motifs to
        the ambiguous characters which they might represent. For example,
        {'R' : 'AG'} specifies that Rs in the motif can match a A or a G.
        All letters in the motif must be represented in the ambiguity_info
        dictionary.
        """
        self._ambiguity_info = ambiguity_info

        # a cache of all encoded motifs
        self._motif_cache = {}

    def encode_motif(self, motif):
        """Encode the passed motif as a regular expression pattern object.
        
        Arguments:

        o motif - The motif we want to encode. This should be a string.
        
        Returns:
        A compiled regular expression pattern object that can be used
        for searching strings.
        """
        regexp_string = ""

        for motif_letter in motif:
            try:
                letter_matches = self._ambiguity_info[motif_letter]
            except KeyError:
                raise KeyError("No match information for letter %s"
                               % motif_letter)

            if len(letter_matches) > 1:
                regexp_match = "[" + letter_matches + "]"
            elif len(letter_matches) == 1:
                regexp_match = letter_matches
            else:
                raise ValueError("Unexpected match information %s"
                                 % letter_matches)

            regexp_string += regexp_match

        return re.compile(regexp_string)

    def find_ambiguous(self, motif):
        """Return the location of ambiguous items in the motif.

        This just checks through the motif and compares each letter
        against the ambiguity information. If a letter stands for multiple
        items, it is ambiguous.
        """
        ambig_positions = []
        for motif_letter_pos in range(len(motif)):
            motif_letter = motif[motif_letter_pos]
            try:
                letter_matches = self._ambiguity_info[motif_letter]
            except KeyError:
                raise KeyError("No match information for letter %s"
                               % motif_letter)

            if len(letter_matches) > 1:
                ambig_positions.append(motif_letter_pos)

        return ambig_positions

    def num_ambiguous(self, motif):
        """Return the number of ambiguous letters in a given motif.
        """
        ambig_positions = self.find_ambiguous(motif)
        return len(ambig_positions)

    def find_matches(self, motif, query):
        """Return all non-overlapping motif matches in the query string.

        This utilizes the regular expression findall function, and will
        return a list of all non-overlapping occurances in query that
        match the ambiguous motif.
        """
        try:
            motif_pattern = self._motif_cache[motif]
        except KeyError:
            motif_pattern = self.encode_motif(motif)
            self._motif_cache[motif] = motif_pattern

        return motif_pattern.findall(query)

    def num_matches(self, motif, query):
        """Find the number of non-overlapping times motif occurs in query.
        """
        all_matches = self.find_matches(motif, query)
        return len(all_matches)

    def all_unambiguous(self):
        """Return a listing of all unambiguous letters allowed in motifs.
        """
        all_letters = sorted(self._ambiguity_info)
        unambig_letters = []

        for letter in all_letters:
            possible_matches = self._ambiguity_info[letter]
            if len(possible_matches) == 1:
                unambig_letters.append(letter)

        return unambig_letters

# --- helper classes and functions for the default SchemaFinder

# -- Alphabets

class SchemaDNAAlphabet(Alphabet.Alphabet):
    """Alphabet of a simple Schema for DNA sequences.

    This defines a simple alphabet for DNA sequences that has a single
    character which can match any other character.

    o G,A,T,C - The standard unambiguous DNA alphabet.

    o * - Any letter
    """
    letters = ["G", "A", "T", "C", "*"]
    
    alphabet_matches = {"G" : "G",
                        "A" : "A",
                        "T" : "T",
                        "C" : "C",
                        "*" : "GATC"}

# -- GA schema finder

class GeneticAlgorithmFinder(object):
    """Find schemas using a genetic algorithm approach.

    This approach to finding schema uses Genetic Algorithms to evolve
    a set of schema and find the best schema for a specific set of
    records.

    The 'default' finder searches for ambiguous DNA elements. This
    can be overridden easily by creating a GeneticAlgorithmFinder
    with a different alphabet.
    """
    def __init__(self, alphabet = SchemaDNAAlphabet()):
        """Initialize a finder to get schemas using Genetic Algorithms.

        Arguments:

        o alphabet -- The alphabet which specifies the contents of the
        schemas we'll be generating. This alphabet must contain the
        attribute 'alphabet_matches', which is a dictionary specifying
        the potential ambiguities of each letter in the alphabet. These
        ambiguities will be used in building up the schema.
        """
        self.alphabet = alphabet

        self.initial_population = 500
        self.min_generations = 10

        self._set_up_genetic_algorithm()

    def _set_up_genetic_algorithm(self):
        """Overrideable function to set up the genetic algorithm parameters.

        This functions sole job is to set up the different genetic
        algorithm functionality. Since this can be quite complicated, this
        allows cusotmizablity of all of the parameters. If you want to
        customize specially, you can inherit from this class and override
        this function.
        """
        self.motif_generator = RandomMotifGenerator(self.alphabet)
        
        self.mutator = SinglePositionMutation(mutation_rate = 0.1)
        self.crossover = SinglePointCrossover(crossover_prob = 0.25)
        self.repair = AmbiguousRepair(Schema(self.alphabet.alphabet_matches),
                                      4)
        self.base_selector = TournamentSelection(self.mutator, self.crossover,
                                                 self.repair, 2)
        self.selector = DiversitySelection(self.base_selector,
                                           self.motif_generator.random_motif)

    def find_schemas(self, fitness, num_schemas):
        """Find the given number of unique schemas using a genetic algorithm

        Arguments:

        o fitness - A callable object (ie. function) which will evaluate
        the fitness of a motif.

        o num_schemas - The number of unique schemas with good fitness
        that we want to generate.
        """
        start_population = \
           Organism.function_population(self.motif_generator.random_motif,
                                        self.initial_population,
                                        fitness)
        finisher = SimpleFinisher(num_schemas, self.min_generations)

        # set up the evolver and do the evolution
        evolver = GenerationEvolver(start_population, self.selector)
        evolved_pop = evolver.evolve(finisher.is_finished)

        # convert the evolved population into a PatternRepository
        schema_info = {}
        for org in evolved_pop:
            # convert the Genome from a MutableSeq to a Seq so that
            # the schemas are just strings (and not array("c")s)
            seq_genome = org.genome.toseq()
            schema_info[seq_genome.tostring()] = org.fitness

        return PatternRepository(schema_info)

# -- fitness classes

class DifferentialSchemaFitness(object):
    """Calculate fitness for schemas that differentiate between sequences.
    """
    def __init__(self, positive_seqs, negative_seqs, schema_evaluator):
        """Initialize with different sequences to evaluate

        Arguments:
        
        o positive_seq - A list of SeqRecord objects which are the 'positive'
        sequences -- the ones we want to select for.

        o negative_seq - A list of SeqRecord objects which are the 'negative'
        sequences that we want to avoid selecting.

        o schema_evaluator - An Schema class which can be used to
        evaluate find motif matches in sequences.
        """
        self._pos_seqs = positive_seqs
        self._neg_seqs = negative_seqs
        self._schema_eval = schema_evaluator

    def calculate_fitness(self, genome):
        """Calculate the fitness for a given schema.

        Fitness is specified by the number of occurances of the schema in
        the positive sequences minus the number of occurances in the
        negative examples.

        This fitness is then modified by multiplying by the length of the
        schema and then dividing by the number of ambiguous characters in
        the schema. This helps select for schema which are longer and have
        less redundancy.
        """
        # convert the genome into a string
        seq_motif = genome.toseq()
        motif = seq_motif.tostring()
        
        # get the counts in the positive examples
        num_pos = 0
        for seq_record in self._pos_seqs:
            cur_counts = self._schema_eval.num_matches(motif,
                                                      seq_record.seq.tostring())
            num_pos += cur_counts

        # get the counts in the negative examples
        num_neg = 0
        for seq_record in self._neg_seqs:
            cur_counts = self._schema_eval.num_matches(motif,
                                                      seq_record.seq.tostring())

            num_neg += cur_counts

        num_ambiguous = self._schema_eval.num_ambiguous(motif)
        # weight the ambiguous stuff more highly
        num_ambiguous = pow(2.0, num_ambiguous)
        # increment num ambiguous to prevent division by zero errors.
        num_ambiguous += 1

        motif_size = len(motif)
        motif_size = motif_size * 4.0

        discerning_power = num_pos - num_neg
        
        diff = (discerning_power * motif_size) / float(num_ambiguous)
        return diff

class MostCountSchemaFitness(object):
    """Calculate a fitness giving weight to schemas that match many times.

    This fitness function tries to maximize schemas which are found many
    times in a group of sequences.
    """
    def __init__(self, seq_records, schema_evaluator):
        """Initialize with sequences to evaluate.

        Arguments:
        
        o seq_records -- A set of SeqRecord objects which we use to
        calculate the fitness.

        o schema_evaluator - An Schema class which can be used to
        evaluate find motif matches in sequences.
        """
        self._records = seq_records
        self._evaluator = schema_evaluator

    def calculate_fitness(self, genome):
        """Calculate the fitness of a genome based on schema matches.

        This bases the fitness of a genome completely on the number of times
        it matches in the set of seq_records. Matching more times gives a
        better fitness
        """
        # convert the genome into a string
        seq_motif = genome.toseq()
        motif = seq_motif.tostring()
        
        # find the number of times the genome matches
        num_times = 0
        for seq_record in self._records:
            cur_counts = self._evaluator.num_matches(motif,
                                                     seq_record.seq.tostring())
            num_times += cur_counts

        return num_times

# -- Helper classes
class RandomMotifGenerator(object):
    """Generate a random motif within given parameters.
    """
    def __init__(self, alphabet, min_size = 12, max_size = 17):
        """Initialize with the motif parameters.

        Arguments:

        o alphabet - An alphabet specifying what letters can be inserted in
        a motif.

        o min_size, max_size - Specify the range of sizes for motifs.
        """
        self._alphabet = alphabet
        self._min_size = min_size
        self._max_size = max_size

    def random_motif(self):
        """Create a random motif within the given parameters.
        
        This returns a single motif string with letters from the given
        alphabet. The size of the motif will be randomly chosen between
        max_size and min_size.
        """
        motif_size = random.randrange(self._min_size, self._max_size)

        motif = ""
        for letter_num in range(motif_size):
            cur_letter = random.choice(self._alphabet.letters)
            motif += cur_letter

        return MutableSeq(motif, self._alphabet)

class SimpleFinisher(object):
    """Determine when we are done evolving motifs.

    This takes the very simple approach of halting evolution when the
    GA has proceeded for a specified number of generations and has
    a given number of unique schema with positive fitness.
    """
    def __init__(self, num_schemas, min_generations = 100):
        """Initialize the finisher with its parameters.

        Arguments:

        o num_schemas -- the number of useful (positive fitness) schemas
        we want to generation

        o min_generations -- The minimum number of generations to allow
        the GA to proceed.
        """
        self.num_generations = 0

        self.num_schemas = num_schemas
        self.min_generations = min_generations

    def is_finished(self, organisms):
        """Determine when we can stop evolving the population.
        """
        self.num_generations += 1
        # print "generation %s" % self.num_generations

        if self.num_generations >= self.min_generations:
            all_seqs = []
            for org in organisms:
                if org.fitness > 0:
                    if org.genome not in all_seqs:
                        all_seqs.append(org.genome)

            if len(all_seqs) >= self.num_schemas:
                return 1

        return 0
# ---

class SchemaFinder(object):
    """Find schema in a set of sequences using a genetic algorithm approach.

    Finding good schemas is very difficult because it takes forever to
    enumerate all of the potential schemas. This finder using a genetic
    algorithm approach to evolve good schema which match many times in
    a set of sequences.

    The default implementation of the finder is ready to find schemas
    in a set of DNA sequences, but the finder can be customized to deal
    with any type of data.
    """
    def __init__(self, num_schemas = 100,
                 schema_finder = GeneticAlgorithmFinder()):
        self.num_schemas = num_schemas
        self._finder = schema_finder

        self.evaluator = Schema(self._finder.alphabet.alphabet_matches)

    def find(self, seq_records):
        """Find well-represented schemas in the given set of SeqRecords.
        """
        fitness_evaluator = MostCountSchemaFitness(seq_records,
                                                   self.evaluator)

        return self._finder.find_schemas(fitness_evaluator.calculate_fitness,
                                         self.num_schemas)

    def find_differences(self, first_records, second_records):
        """Find schemas which differentiate between the two sets of SeqRecords.
        """
        fitness_evaluator = DifferentialSchemaFitness(first_records,
                                                      second_records,
                                                      self.evaluator)

        return self._finder.find_schemas(fitness_evaluator.calculate_fitness,
                                         self.num_schemas)

class SchemaCoder(object):
    """Convert a sequence into a representation of ambiguous motifs (schemas).

    This takes a sequence, and returns the number of times specified
    motifs are found in the sequence. This lets you represent a sequence
    as just a count of (possibly ambiguous) motifs.
    """
    def __init__(self, schemas, ambiguous_converter):
        """Initialize the coder to convert sequences

        Arguments:

        o schema - A list of all of the schemas we want to search for
        in input sequences.

        o ambiguous_converter - An Schema class which can be
        used to convert motifs into regular expressions for searching.
        """
        self._schemas = schemas
        self._converter = ambiguous_converter

    def representation(self, sequence):
        """Represent the given input sequence as a bunch of motif counts.

        Arguments:

        o sequence - A Bio.Seq object we are going to represent as schemas.

        This takes the sequence, searches for the motifs within it, and then
        returns counts specifying the relative number of times each motifs
        was found. The frequencies are in the order the original motifs were
        passed into the initializer.
        """
        schema_counts = []

        for schema in self._schemas:
            num_counts = self._converter.num_matches(schema, sequence.tostring())
            schema_counts.append(num_counts)

        # normalize the counts to go between zero and one
        min_count = 0
        max_count = max(schema_counts)

        # only normalize if we've actually found something, otherwise
        # we'll just return 0 for everything
        if max_count > 0:
            for count_num in range(len(schema_counts)):
                schema_counts[count_num] = (float(schema_counts[count_num]) -
                                           float(min_count)) / float(max_count)

        return schema_counts

def matches_schema(pattern, schema, ambiguity_character = '*'):
    """Determine whether or not the given pattern matches the schema.

    Arguments:

    o pattern - A string representing the pattern we want to check for
    matching. This pattern can contain ambiguity characters (which are
    assumed to be the same as those in the schema).

    o schema - A string schema with ambiguity characters.

    o ambiguity_character - The character used for ambiguity in the schema.
    """
    if len(pattern) != len(schema):
        return 0

    # check each position, and return a non match if the schema and pattern
    # are non ambiguous and don't match
    for pos in range(len(pattern)):
        if (schema[pos] != ambiguity_character and
            pattern[pos] != ambiguity_character and
            pattern[pos] != schema[pos]):
            
            return 0

    return 1

class SchemaFactory(object):
    """Generate Schema from inputs of Motifs or Signatures.
    """
    def __init__(self, ambiguity_symbol = '*'):
        """Initialize the SchemaFactory

        Arguments:

        o ambiguity_symbol -- The symbol to use when specifying that
        a position is arbitrary.
        """
        self._ambiguity_symbol = ambiguity_symbol

    def from_motifs(self, motif_repository, motif_percent, num_ambiguous):
        """Generate schema from a list of motifs.

        Arguments:

        o motif_repository - A MotifRepository class that has all of the
        motifs we want to convert to Schema.

        o motif_percent - The percentage of motifs in the motif bank which
        should be matches. We'll try to create schema that match this
        percentage of motifs.

        o num_ambiguous - The number of ambiguous characters to include
        in each schema. The positions of these ambiguous characters will
        be randomly selected.
        """
        # get all of the motifs we can deal with
        all_motifs = motif_repository.get_top_percentage(motif_percent)

        # start building up schemas
        schema_info = {}
        # continue until we've built schema matching the desired percentage
        # of motifs
        total_count = self._get_num_motifs(motif_repository, all_motifs)
        matched_count = 0
        assert total_count > 0, "Expected to have motifs to match"
        while (float(matched_count) / float(total_count)) < motif_percent:
            
            new_schema, matching_motifs = \
                        self._get_unique_schema(schema_info.keys(),
                                                all_motifs, num_ambiguous)

            # get the number of counts for the new schema and clean up
            # the motif list
            schema_counts = 0
            for motif in matching_motifs:
                # get the counts for the motif
                schema_counts += motif_repository.count(motif)

                # remove the motif from the motif list since it is already
                # represented by this schema
                all_motifs.remove(motif)


            # all the schema info
            schema_info[new_schema] = schema_counts

            matched_count += schema_counts

            # print "percentage:", float(matched_count) / float(total_count)

        return PatternRepository(schema_info)

    def _get_num_motifs(self, repository, motif_list):
        """Return the number of motif counts for the list of motifs.
        """
        motif_count = 0
        for motif in motif_list:
            motif_count += repository.count(motif)

        return motif_count

    def _get_unique_schema(self, cur_schemas, motif_list, num_ambiguous):
        """Retrieve a unique schema from a motif.

        We don't want to end up with schema that match the same thing,
        since this could lead to ambiguous results, and be messy. This
        tries to create schema, and checks that they do not match any
        currently existing schema.
        """
        # create a schema starting with a random motif
        # we'll keep doing this until we get a completely new schema that
        # doesn't match any old schema
        num_tries = 0
        
        while 1:
            # pick a motif to work from and make a schema from it
            cur_motif = random.choice(motif_list)
            
            num_tries += 1
                
            new_schema, matching_motifs = \
                        self._schema_from_motif(cur_motif, motif_list,
                                                num_ambiguous)

            has_match = 0
            for old_schema in cur_schemas:
                if matches_schema(new_schema, old_schema,
                                  self._ambiguity_symbol):
                    has_match = 1

            # if the schema doesn't match any other schema we've got
            # a good one
            if not(has_match):
                break

            # check for big loops in which we can't find a new schema
            assert num_tries < 150, \
                   "Could not generate schema in %s tries from %s with %s" \
                   % (num_tries, motif_list, cur_schemas)

        return new_schema, matching_motifs

    def _schema_from_motif(self, motif, motif_list, num_ambiguous):
        """Create a schema from a given starting motif.

        Arguments:

        o motif - A motif with the pattern we will start from.

        o motif_list - The total motifs we have.to match to.

        o num_ambiguous - The number of ambiguous characters that should
        be present in the schema.

        Returns:

        o A string representing the newly generated schema.

        o A list of all of the motifs in motif_list that match the schema.
        """
        assert motif in motif_list, \
               "Expected starting motif present in remaining motifs."

        # convert random positions in the motif to ambiguous characters
        # convert the motif into a list of characters so we can manipulate it
        new_schema_list = list(motif)
        for add_ambiguous in range(num_ambiguous):
            # add an ambiguous position in a new place in the motif
            while 1:
                ambig_pos = random.choice(range(len(new_schema_list)))

                # only add a position if it isn't already ambiguous
                # otherwise, we'll try again
                if new_schema_list[ambig_pos] != self._ambiguity_symbol:
                    new_schema_list[ambig_pos] = self._ambiguity_symbol
                    break

        # convert the schema back to a string
        new_schema = ''.join(new_schema_list)

        # get the motifs that the schema matches
        matched_motifs = []
        for motif in motif_list:
            if matches_schema(motif, new_schema, self._ambiguity_symbol):
                matched_motifs.append(motif)

        return new_schema, matched_motifs
            
    def from_signatures(self, signature_repository, num_ambiguous):
        raise NotImplementedError("Still need to code this.")
