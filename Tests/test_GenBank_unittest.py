# Copyright 2013 by Kai Blin.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import warnings
from os import path, getcwd
from collections import namedtuple

from Bio import SeqIO
from Bio import BiopythonParserWarning

class GenBankTests(unittest.TestCase):
    def test_invalid_product_line_raises_value_error(self):
        "Test GenBank parsing invalid product line raises ValueError"
        def parse_invalid_product_line():
            rec = SeqIO.read(path.join('GenBank', 'invalid_product.gb'),
                             'genbank')
        self.assertRaises(ValueError, parse_invalid_product_line)

#define a named tuple to make tests more explicit
File = namedtuple('File', ['path', 'name'])

#set base directory
gb_file_dir = path.join(getcwd(), 'GenBank')

#well behaved files
test_files = ['noref.gb', 'cor6_6.gb', 'iro.gb', 'pri1.gb', 'arab1.gb',
              'protein_refseq.gb', 'extra_keywords.gb', 'one_of.gb',
              'NT_019265.gb', 'origin_line.gb', 'blank_seq.gb',
              'dbsource_wrap.gb', 'gbvrl1_start.seq', 'NC_005816.gb',
              'empty_feature_qualifier.gb']
test_files = [File(path.join(gb_file_dir, f), f) for f in test_files]

#files that induce warnings
warn_files = ['no_end_marker.gb', 'wrong_sequence_indent.gb',
              'invalid_locus_line_spacing.gb', 'invalid_misc_feature.gb', 
              '1MRR_A.gp']
warn_files = [File(path.join(gb_file_dir, f), f) for f in warn_files]


class GenBankTestsManyFiles(unittest.TestCase):
    
    def generic_parse_one(self,filetuple):
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            ft = next(SeqIO.parse(filetuple.path, 'genbank'))
            # Verify some things
            assert len(w) == 0

    def generic_parse_seq_id(self,filetuple):
        rec = next(SeqIO.parse(filetuple.path, 'genbank'))
        seq = rec.seq
        id = rec.id

    def warn_inducing_parse(self,filetuple):
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            ft = next(SeqIO.parse(filetuple.path, 'genbank'))
            # Verify some things

#set up tests that expect no warnings
for filetuple in test_files:
    name = filetuple.name.split(".")[0]

    def funct(fn):
        f = lambda x : x.generic_parse_one(fn)
        f.__doc__ = "Checking nucleotide file %s" % fn.name
        return f

    def funct2(fn):
        f = lambda x : x.generic_parse_seq_id(fn)
        f.__doc__ = "Checking nucleotide file %s" % fn.name
        return f

    setattr(GenBankTestsManyFiles, "test_nuc_%s"%name, funct(filetuple))
    setattr(GenBankTestsManyFiles, "test_nuc_seq_and_id%s"%name, \
            funct2(filetuple))
    del funct

#set up tests that expect a warning
for filetuple in warn_files:
    name = filetuple.name.split(".")[0]

    def funct(fn):
        f = lambda x : x.warn_inducing_parse(fn)
        f.__doc__ = "Checking nucleotide file %s" % fn.name
        return f
    
    setattr(GenBankTestsManyFiles, "test_warnings_from_%s"%name, \
            funct(filetuple))
    del funct


#tests that check granular consumer behavior
import re
from Bio.GenBank import _FeatureConsumer
from Bio.GenBank.utils import FeatureValueCleaner
from Bio.GenBank.Scanner import GenBankScanner
from Bio._py3k import StringIO

class ConsumerBehaviorTest(unittest.TestCase):
    
    recordfile = "brca_FJ940752.gb" 

    def setUp(self):
        self.handle.seek(0)
        self.scanner = GenBankScanner(debug=0)
        self.consumer = _FeatureConsumer(use_fuzziness=1,
                            feature_cleaner=FeatureValueCleaner())

    @classmethod
    def setUpClass(cls):
        cls.handle = open(path.join('GenBank', cls.recordfile), 'r')
    
    @classmethod
    def tearDownClass(cls):
        cls.handle.close()
                
    def test_handle_assignment(self):
        self.scanner.set_handle(self.handle)
        self.assertEqual(self.scanner.line, "")

    def test_complete_(self):
        scanner = self.scanner
        consumer = self.consumer
        #step 0: manually setting the handle
        # same as >>> scanner.set_handle(self.handle)
        scanner.handle = self.handle
        scanner.line = ""
        #step 1: find_start
        scanner.find_start()
        self.assertEqual(scanner.line[0:20], "LOCUS       FJ940752")
        self.assertTrue(bool(re.match(r"^LOCUS(\s)*",scanner.line)))
        #step 2: feed_first_line (expect no advancement)
        scanner._feed_first_line(self.consumer, line=scanner.line)
        self.assertEqual(scanner.line[0:20], "LOCUS       FJ940752")
        self.assertTrue(bool(re.match(r"^LOCUS(\s)*",scanner.line)))
        #step 3: feed all header data
        scanner._feed_header_lines(consumer, scanner.parse_header())
        self.assertEqual(scanner.line[0:40], "FEATURES             Location/Qualifiers")
        self.assertTrue(bool(re.match(r"^FEATURES(\s)*Location/Qualifiers",scanner.line)))
        #step 4: feed the features
        featuretuple = scanner.parse_features(skip=False)
        testline = scanner.line
        scanner._feed_feature_table(consumer, featuretuple)
        self.assertEqual(testline[0:6], "ORIGIN")
        self.assertTrue(bool(re.match(r"^ORIGIN",testline)))
        #step 5: feed the footer
        misc_lines, sequence_string = scanner.parse_footer()
        self.assertEqual("//", scanner.line)
        self.assertTrue(bool(re.match(r"^//",scanner.line)))
        #step 6: finish the consumer
        scanner._feed_misc_lines(consumer, misc_lines)
        consumer.sequence(sequence_string)

    def test_get_start_by_cut(self):
        #get positions
        position = 0
        while True:
            line = self.handle.readline()
            if re.match(r"^FEATURES(\s)*Location/Qualifiers",line):
                position = self.handle.tell()
                break
            if not line:
                raise ValueError("record has bad header")
        self.handle.seek(0)
        top = StringIO(self.handle.read(position-0))
        #use the fake thing
        scanner = self.scanner
        consumer = self.consumer
        scanner.set_handle(top)
        scanner.find_start()
        scanner._feed_first_line(self.consumer, line=scanner.line)
        scanner._feed_header_lines(consumer, scanner.parse_header())
        self.assertEqual(consumer.data.id, "FJ940752")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
