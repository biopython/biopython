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
        parse = lambda ft: next(SeqIO.parse(ft.path, 'genbank'))
        self.assertNoWarn(BiopythonParserWarning, parse, filetuple)

    def generic_parse_seq_id(self,filetuple):
        rec = next(SeqIO.parse(filetuple.path, 'genbank'))
        seq = rec.seq
        id = rec.id

    def warn_inducing_parse(self,filetuple):
        parse = lambda ft: next(SeqIO.parse(ft.path, 'genbank'))
        self.assertWarns(BiopythonParserWarning, parse, filetuple)

    def assertWarns(self, warning, callable, *args, **kwds):
        with warnings.catch_warnings(record=True) as warning_list:
            warnings.simplefilter('always')
            result = callable(*args, **kwds)
            self.assertTrue(any(item.category == warning for item in warning_list))

    def assertNoWarn(self, warning, callable, *args, **kwds):
        with warnings.catch_warnings(record=True) as warning_list:
            warnings.simplefilter('always')
            result = callable(*args, **kwds)
            self.assertTrue(not any(item.category == warning for item in warning_list))

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

for filetuple in warn_files:
    name = filetuple.name.split(".")[0]

    def funct(fn):
        f = lambda x : x.warn_inducing_parse(fn)
        f.__doc__ = "Checking nucleotide file %s" % fn.name
        return f

    setattr(GenBankTestsManyFiles, "test_warnings_from_%s"%name, \
            funct(filetuple))
    del funct


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
