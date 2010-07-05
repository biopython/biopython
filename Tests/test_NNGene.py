#!/usr/bin/env python

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Test the different representations of Genes.

This exercises the Motif, Schema and Signature methods of representing
genes, as well as generic Pattern methods.
"""
# standard library
import os
import unittest

# Biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# stuff we are testing
from Bio.NeuralNetwork.Gene import Schema
from Bio.NeuralNetwork.Gene import Motif
from Bio.NeuralNetwork.Gene import Signature
from Bio.NeuralNetwork.Gene import Pattern

VERBOSE = 0


# --- Tests for Pattern

class PatternIOTest(unittest.TestCase):
    """Tests for reading and writing patterns to a file.
    """
    def setUp(self):
        self.alphabet = IUPAC.ambiguous_dna
        self.test_file = os.path.join("NeuralNetwork", "patternio.txt")
        #Remove any existing copy of the output file,
        if os.path.isfile(self.test_file):
            os.remove(self.test_file)
        self.pattern_io = Pattern.PatternIO(self.alphabet)

    def tearDown(self):
        #Clean up by removing our output file,
        if os.path.isfile(self.test_file):
            os.remove(self.test_file)

    def test_motif(self):
        """Reading and writing motifs to a file
        """
        # write to a file
        motifs = ["GAC", "AAA", "TTT", "GGG"]
        output_handle = open(self.test_file, "w")
        self.pattern_io.write(motifs, output_handle)
        output_handle.close()

        # read 'em back
        input_handle = open(self.test_file, "r")
        read_motifs = self.pattern_io.read(input_handle)
        input_handle.close()
        assert read_motifs == motifs, \
               "Failed to get back expected motifs %s, got %s" \
               % (motifs, read_motifs)

        # write seqs
        seq_motifs = []
        for motif in motifs:
            seq_motifs.append(Seq(motif, self.alphabet))
        output_handle = open(self.test_file, "w")
        self.pattern_io.write_seq(seq_motifs, output_handle)
        output_handle.close()

        # read the seqs back
        input_handle = open(self.test_file, "r")
        read_motifs = self.pattern_io.read(input_handle)
        input_handle.close()
        assert read_motifs == motifs, \
               "Failed to get back expected motifs %s from seqs, got %s" \
               % (motifs, read_motifs)

    def test_schema(self):
        """Reading and writing schemas to a file.
        """
        schemas = ["GTR", "GAC"]
        # write out the schemas
        output_handle = open(self.test_file, "w")
        self.pattern_io.write(schemas, output_handle)
        output_handle.close()

        # read back the schemas
        input_handle = open(self.test_file, "r")
        read_schemas = self.pattern_io.read(input_handle)
        input_handle.close()
        assert schemas == read_schemas, \
               "Read incorrect schemas %s, expected %s." \
               % (read_schemas, schemas)

        # --- make sure inappropriate alphabets are reported
        schemas = ["GTR", "G*C"] # '*' not in the unambigous alphabet
        output_handle = open(self.test_file, "w")
        self.pattern_io.write(schemas, output_handle)
        output_handle.close()

        input_handle = open(self.test_file, "r")
        try:
            read_schemas = self.pattern_io.read(input_handle)
            raise AssertionError("Did not report error on bad alphabet.")
        except ValueError:
            pass # expected behavior
        except:
            raise AssertionError("Got unexpected error while reading.")

        input_handle.close()

    def test_signature(self):
        """Reading and writing signatures to a file.
        """
        signatures = [("GAC", "GAC"), ("AAA", "TTT")]
        output_handle = open(self.test_file, "w")
        self.pattern_io.write(signatures, output_handle)
        output_handle.close()

        input_handle = open(self.test_file, "r")
        read_sigs = self.pattern_io.read(input_handle)
        input_handle.close()
        assert read_sigs == signatures, \
               "Got back unexpected signatures %s, wanted %s" \
               % (read_sigs, signatures)

class PatternRepositoryTest(unittest.TestCase):
    """Tests for retrieving info from a repository of patterns.
    """ 
    def setUp(self):
        self.motifs = {"GATC" : 30,
                       "GGGG" : 10,
                       "GTAG" : 0,
                       "AAAA" : -10,
                       "ATAT" : -20}

        self.repository = Pattern.PatternRepository(self.motifs)

    def test_get_all(self):
        """Retrieve all patterns from a repository.
        """
        all_motifs = self.repository.get_all()

        assert all_motifs == ["GATC", "GGGG", "GTAG", "AAAA", "ATAT"], \
               "Unexpected motifs returned %s" % all_motifs

    def test_get_random(self):
        """Retrieve random patterns from the repository.
        """
        for num_patterns in range(5):
            patterns = self.repository.get_random(num_patterns)
            assert len(patterns) == num_patterns, \
                   "Got unexpected number of patterns %s, expected %s" \
                   % (len(patterns), num_patterns)

            for pattern in patterns:
                assert pattern in self.motifs.keys(), \
                       "Got unexpected pattern %s" % pattern

    def test_get_top_percentage(self):
        """Retrieve the top percentge of patterns from the repository.
        """
        for num_patterns, percentage in ((1, 0.2), (2, .4), (5, 1.0)):
            patterns = self.repository.get_top_percentage(percentage)
            assert len(patterns) == num_patterns, \
                   "Got unexpected number of patterns %s, expected %s" \
                   % (len(patterns), num_patterns)

            for pattern in patterns:
                assert pattern in self.motifs.keys(), \
                       "Got unexpected pattern %s" % pattern      

    def test_get_top(self):
        """Retrieve a certain number of the top patterns.
        """
        for num_patterns in range(5):
            patterns = self.repository.get_top(num_patterns)
            assert len(patterns) == num_patterns, \
                   "Got unexpected number of patterns %s, expected %s" \
                   % (len(patterns), num_patterns)

            for pattern in patterns:
                assert pattern in self.motifs.keys(), \
                       "Got unexpected pattern %s" % pattern       

    def test_get_differing(self):
        """Retrieve patterns from both sides of the list (top and bottom).
        """
        patterns = self.repository.get_differing(2, 2)
        assert patterns == ["GATC", "GGGG", "AAAA", "ATAT"], \
               "Got unexpected patterns %s" % patterns

    def test_remove_polyA(self):
        """Test the ability to remove A rich patterns from the repository.
        """
        patterns = self.repository.get_all()
        assert len(patterns) == 5, "Unexpected starting: %s" % patterns

        self.repository.remove_polyA()
        
        patterns = self.repository.get_all()
        assert len(patterns) == 3, "Unexpected ending: %s" % patterns
        assert patterns == ["GATC", "GGGG", "GTAG"], \
               "Unexpected patterns: %s" % patterns

    def test_count(self):
        """Retrieve counts for particular patterns in the repository.
        """
        num_times = self.repository.count("GGGG")
        assert num_times == 10, \
               "Did not count item in the respository: %s" % num_times

        num_times = self.repository.count("NOT_IN_THERE")
        assert num_times == 0, \
               "Counted items not in repository: %s" % num_times

# --- Tests for motifs

class MotifFinderTest(unittest.TestCase):
    """Tests for finding motifs from sequences.
    """
    def setUp(self):
        test_file = os.path.join('NeuralNetwork', 'enolase.fasta')
        diff_file = os.path.join('NeuralNetwork', 'repeat.fasta')

        self.test_records = []
        self.diff_records = []

        # load the records
        for file, records in ((test_file, self.test_records),
                              (diff_file, self.diff_records)):

            handle = open(file, 'r')

            iterator = SeqIO.parse(handle, "fasta",
                                   alphabet=IUPAC.unambiguous_dna)
            while 1:
                try:
                    seq_record = iterator.next()
                except StopIteration:
                    break
                if seq_record is None:
                    break

                records.append(seq_record)

            handle.close()

        self.motif_finder = Motif.MotifFinder()

    def test_find(self):
        """Find all motifs in a set of sequences.
        """
        motif_repository = self.motif_finder.find(self.test_records, 8)
        top_motif = motif_repository.get_top(1)

        assert top_motif[0] == 'TTGGAAAG', \
               "Got unexpected motif %s" % top_motif[0]

    def test_find_differences(self):
        """Find the difference in motif counts between two sets of sequences.
        """
        motif_repository = \
               self.motif_finder.find_differences(self.test_records,
                                                  self.diff_records, 8)

        top, bottom = motif_repository.get_differing(1, 1)

        assert top == "TTGGAAAG", "Got unexpected top motif %s" % top
        assert bottom == "AATGGCAT", "Got unexpected bottom motif %s" % bottom

class MotifCoderTest(unittest.TestCase):
    """Test the ability to encode sequences as a set of motifs.
    """
    def setUp(self):
        motifs = ["GAG", "GAT", "GCC", "ATA"]

        self.match_strings = (("GATCGCC", [0.0, 1.0, 1.0, 0.0]),
                              ("GATGATCGAGCC", [.5, 1.0, .5, 0.0]))
        
        self.coder = Motif.MotifCoder(motifs)

        
    def test_representation(self):
        """Convert a sequence into its motif representation.
        """
        for match_string, expected in self.match_strings:
            seq_to_code = Seq(match_string, IUPAC.unambiguous_dna)
            matches = self.coder.representation(seq_to_code)

            assert matches == expected, \
                   "Did not match representation, expected %s, got %s" \
                   % (expected, matches)

# --- Tests for schemas

class SchemaTest(unittest.TestCase):
    """Matching ambiguous motifs with multiple ambiguity characters.
    """
    def setUp(self):
        ambiguity_chars = {"G" : "G",
                           "A" : "A",
                           "T" : "T",
                           "C" : "C",
                           "R" : "AG",
                           "*" : "AGTC"}

        self.motif_coder = Schema.Schema(ambiguity_chars)

        self.match_string = "GATAG"
        self.match_info = [("GA", ["GA"]),
                           ("GATAG", ["GATAG"]),
                           ("GA*AG", ["GATAG"]),
                           ("GATRG", ["GATAG"]),
                           ("*A", ["GA", "TA"])]

    def test_find_matches(self):
        """Find all matches in a sequence.
        """
        for motif, expected in self.match_info:
            found_matches = self.motif_coder.find_matches(motif,
                                                          self.match_string)
            assert found_matches == expected, "Expected %s, got %s" \
                   % (expected, found_matches)

    def test_num_matches(self):
        """Find how many matches are present in a sequence.
        """
        for motif, expected in self.match_info:
            num_matches = self.motif_coder.num_matches(motif,
                                                       self.match_string)
            assert num_matches == len(expected), \
                   "Expected %s, got %s" % (num_matches, len(expected))

    def test_find_ambiguous(self):
        """Find the positions of ambiguous items in a sequence.
        """
        ambig_info = (("GATC", []),
                      ("G***", [1, 2, 3]),
                      ("GART", [2]),
                      ("*R*R", [0, 1, 2, 3]))

        for motif, expected in ambig_info:
            found_positions = self.motif_coder.find_ambiguous(motif)
            assert found_positions == expected, \
                   "Expected %s, got %s for %s" % (expected, found_positions,
                                                   motif)
        
    def test_num_ambiguous(self):
        """Find the number of ambiguous items in a sequence.
        """
        ambig_info = (("GATC", 0),
                      ("G***", 3),
                      ("GART", 1),
                      ("*R*R", 4))

        for motif, expected in ambig_info:
            found_num = self.motif_coder.num_ambiguous(motif)
            assert found_num == expected, \
                   "Expected %s, got %s for %s" % (expected, found_num, motif)

    def test_motif_cache(self):
        """Make sure motif compiled regular expressions are cached properly.
        """
        test_motif = "GATC"

        self.motif_coder.find_matches(test_motif, "GATCGATC")

        self.assertTrue(test_motif in self.motif_coder._motif_cache,
                     "Did not find motif cached properly.")

        # make sure we don't bomb out if we use the same motif twice
        self.motif_coder.find_matches(test_motif, "GATCGATC")

    def test_all_unambiguous(self):
        """Return all unambiguous characters that can be in a motif.
        """
        found_unambig = self.motif_coder.all_unambiguous()

        expected = ["A", "C", "G", "T"]
        assert found_unambig == expected, \
               "Got %s, expected %s" % (found_unambig, expected)

class SchemaFinderTest(unittest.TestCase):
    """Test finding schemas from a set of sequences.
    """
    def setUp(self):
        test_file = os.path.join('NeuralNetwork', 'enolase.fasta')
        diff_file = os.path.join('NeuralNetwork', 'repeat.fasta')

        self.test_records = []
        self.diff_records = []

        # load the records
        for file, records in ((test_file, self.test_records),
                              (diff_file, self.diff_records)):

            handle = open(file, 'r')
            records.extend(SeqIO.parse(handle, "fasta",
                                       alphabet=IUPAC.unambiguous_dna))
            handle.close()

        self.num_schemas = 2
        schema_ga = Schema.GeneticAlgorithmFinder()
        schema_ga.min_generations = 1
        self.finder = Schema.SchemaFinder(num_schemas = self.num_schemas,
                                          schema_finder = schema_ga)

    def test_find(self):
        """Find schemas from sequence inputs.
        """
        # this test takes too long
        if VERBOSE:
            repository = self.finder.find(self.test_records + self.diff_records)
            schemas = repository.get_all()

            assert len(schemas) >= self.num_schemas, "Got too few schemas."

    def test_find_differences(self):
        """Find schemas that differentiate between two sets of sequences.
        """
        # this test takes too long
        if VERBOSE:
            repository = self.finder.find_differences(self.test_records,
                                                      self.diff_records)
            schemas = repository.get_all()

            assert len(schemas) >= self.num_schemas, "Got too few schemas."
        
class SchemaCoderTest(unittest.TestCase):
    """Test encoding sequences as a grouping of motifs.
    """
    def setUp(self):
        ambiguity_chars = {"G" : "G",
                           "A" : "A",
                           "T" : "T",
                           "C" : "C",
                           "R" : "AG",
                           "*" : "AGTC"}

        motif_representation = Schema.Schema(ambiguity_chars)
        motifs = ("GA", "GATAG", "GA*AG", "GATRG", "*A")
        self.motif_coder = Schema.SchemaCoder(motifs,
                                              motif_representation)

        self.match_strings = [("GATAG", [.5, .5, .5, .5, 1.0]),
                              ("GAGAGATA", [float(3) / float(4), 0,
                                            float(1) / float(4), 0,
                                            1])]

    def test_representation(self):
        """Convert a string into a representation of motifs.
        """
        for match_string, expected in self.match_strings:
            match_seq = Seq(match_string, IUPAC.unambiguous_dna)
            found_rep = self.motif_coder.representation(match_seq)
            assert found_rep == expected, "Got %s, expected %s" % \
                   (found_rep, expected)
    
class SchemaMatchingTest(unittest.TestCase):
    """Matching schema to strings works correctly.
    """
    def shortDescription(self):
        return "%s:%s" % (self.__class__.__name__, self.__doc__)
    
    def runTest(self):
        match = Schema.matches_schema("GATC", "AAAAA")
        assert match == 0, "Expected no match because of length differences"

        match = Schema.matches_schema("GATC", "GAT*")
        assert match == 1, "Expected match"

        match = Schema.matches_schema("GATC", "GATC")
        assert match == 1, "Expected match"

        match = Schema.matches_schema("GATC", "C*TC")
        assert match == 0, "Expected no match because of char mismatch."

        match = Schema.matches_schema("G*TC", "*TTC")
        assert match == 1, "Expected match because of ambiguity."

class SchemaFactoryTest(unittest.TestCase):
    """Test the SchemaFactory for generating Schemas.
    """
    def __init__(self, method):
        unittest.TestCase.__init__(self, method)

        # a cached schema bank, so we don't have to load it multiple times
        self.schema_bank = None
    
    def setUp(self):
        self.factory = Schema.SchemaFactory()

        self.test_file = os.path.join(os.getcwd(), "NeuralNetwork", "enolase.fasta")

        ambiguity_chars = {"G" : "G",
                           "A" : "A",
                           "T" : "T",
                           "C" : "C",
                           "R" : "AG",
                           "*" : "AGTC"}

        self.schema = Schema.Schema(ambiguity_chars)

    def test_easy_from_motifs(self):
        """Generating schema from a simple list of motifs.
        """
        motifs = {"GATCGAA" : 20,
                  "GATCGAT" : 15,
                  "GATTGAC" : 25,
                  "TTTTTTT" : 10}

        motif_bank = Pattern.PatternRepository(motifs)

        schema_bank = self.factory.from_motifs(motif_bank, .5, 2)
        if VERBOSE:
            print "\nSchemas:"
            for schema in schema_bank.get_all():
                print "%s: %s" % (schema, schema_bank.count(schema))

    def test_hard_from_motifs(self):
        """Generating schema from a real life set of motifs.
        """
        schema_bank = self._load_schema_repository()

        if VERBOSE:
            print "\nSchemas:"
            for schema in schema_bank.get_top(5):
                print "%s: %s" % (schema, schema_bank.count(schema))

    def _load_schema_repository(self):
        """Helper function to load a schema repository from a file.

        This also caches a schema bank, to prevent having to do this
        time consuming operation multiple times.
        """
        # if we already have a cached repository, return it
        if self.schema_bank is not None:
            return self.schema_bank
        
        # otherwise, we'll read in a new schema bank

        # read in the all of the motif records
        motif_handle = open(self.test_file, 'r')
        seq_records = list(SeqIO.parse(motif_handle, "fasta",
                                       alphabet=IUPAC.unambiguous_dna))
        motif_handle.close()

        # find motifs from the file
        motif_finder = Motif.MotifFinder()
        motif_size = 9

        motif_bank = motif_finder.find(seq_records, motif_size)
        
        schema_bank = self.factory.from_motifs(motif_bank, .1, 2)

        # cache the repository
        self.schema_bank = schema_bank

        return schema_bank

    def test_schema_representation(self):
        """Convert sequences into schema representations.
        """
        # get a set of schemas we want to code the sequence in
        schema_bank = self._load_schema_repository()
        top_schemas = schema_bank.get_top(25)
        schema_coder = Schema.SchemaCoder(top_schemas, self.schema)

        # get the sequences one at a time, and encode them
        fasta_handle = open(self.test_file, 'r')
        for seq_record in SeqIO.parse(fasta_handle, "fasta",
                                      alphabet=IUPAC.unambiguous_dna):
            schema_values = schema_coder.representation(seq_record.seq)
            if VERBOSE:
                print "Schema values:", schema_values
        fasta_handle.close()

# --- Tests for Signatures
class SignatureFinderTest(unittest.TestCase):
    """Test the ability to find signatures in a set of sequences.
    """
    def setUp(self):
        test_file = os.path.join('NeuralNetwork', 'enolase.fasta')

        self.test_records = []

        # load the records
        handle = open(test_file, 'r')
        self.test_records = list(SeqIO.parse(handle, "fasta",
                                             alphabet=IUPAC.unambiguous_dna))
        handle.close()

        self.sig_finder = Signature.SignatureFinder()

    def test_find(self):
        """Find signatures from sequence inputs.
        """
        repository = self.sig_finder.find(self.test_records, 6, 9)
        top_sig = repository.get_top(1)

        assert top_sig[0] == ('TTGGAA', 'TGGAAA'), \
               "Unexpected signature %s" % top_sig[0]

class SignatureCoderTest(unittest.TestCase):
    """Test the ability to encode sequences as a set of signatures.
    """
    def setUp(self):
        signatures = [("GAC", "GAC"), ("AAA", "TTT"), ("CAA", "TTG")]

        self.coder = Signature.SignatureCoder(signatures, 9)

        self.test_seqs = [("GACAAAGACTTT", [1.0, 1.0, 0.0]),
                          ("CAAAGACGACTTTAAATTT", [0.5, 1.0, 0.0]),
                          ("AAATTTAAAGACTTTGAC", [1.0 / 3.0, 1.0, 0.0]),
                          ("GACGAC", [1.0, 0.0, 0.0]),
                          ("GACAAAAAAAAAGAC", [1.0, 0.0, 0.0]),
                          ("GACAAAAAAAAAAGAC", [0.0, 0.0, 0.0])]

    def test_representation(self):
        """Convert a sequence into its signature representation.
        """
        for seq_string, expected in self.test_seqs:
            test_seq = Seq(seq_string, IUPAC.unambiguous_dna)
            predicted = self.coder.representation(test_seq)

            assert predicted == expected, \
                   "Non-expected representation %s for %s, wanted %s" \
                   % (predicted, seq_string, expected)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
