# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.SeqIO.convert(...) function."""
import unittest
from StringIO import StringIO

from Bio import AlignIO
from Bio.Alphabet import generic_protein, generic_nucleotide, generic_dna


#Top level function as this makes it easier to use for debugging:
def check_convert(in_filename, in_format, out_format, alphabet=None):
    #Write it out using parse/write
    handle = StringIO()
    aligns = list(AlignIO.parse(open(in_filename), in_format, None, alphabet))
    try:
        count = AlignIO.write(aligns, handle, out_format)
    except ValueError:
        count = 0
    #Write it out using convert passing filename and handle
    handle2 = StringIO()
    try:
        count2 = AlignIO.convert(in_filename, in_format, handle2, out_format, alphabet)
    except ValueError:
        count2 = 0
    assert count == count2
    assert handle.getvalue() == handle2.getvalue()
    #Write it out using convert passing handle and handle
    handle2 = StringIO()
    try:
        count2 = AlignIO.convert(open(in_filename), in_format, handle2, out_format, alphabet)
    except ValueError:
        count2 = 0
    assert count == count2
    assert handle.getvalue() == handle2.getvalue()
    #TODO - convert passing an output filename?

class ConvertTests(unittest.TestCase):
    """Cunning unit test where methods are added at run time."""
    def simple_check(self, filename, in_format, out_format, alphabet):
        check_convert(filename, in_format, out_format, alphabet)

tests = [
    ('Clustalw/hedgehog.aln', "clustal", None),
    ('Nexus/test_Nexus_input.nex', "nexus", None),
    ('Stockholm/simple.sth', "stockholm", None),
    ('GFF/multi.fna', "fasta", generic_nucleotide),
    ("Quality/example.fastq", "fastq", None),
    ("Quality/example.fastq", "fastq-sanger", generic_dna),
    ('Fasta/output001.m10', "fasta-m10", None),
    ('IntelliGenetics/VIF_mase-pro.txt', "ig", generic_protein),
    ('NBRF/clustalw.pir', "pir", None),
    ]
output_formats = ["fasta"] + sorted(AlignIO._FormatToWriter)

for filename, in_format, alphabet in tests:
    for out_format in output_formats:
        def funct(fn,fmt1, fmt2, alpha):
            f = lambda x : x.simple_check(fn, fmt1, fmt2, alpha)
            f.__doc__ = "Convert %s from %s to %s" % (fn, fmt1, fmt2)
            return f
        setattr(ConvertTests, "test_%s_%s_to_%s" \
                % (filename.replace("/","_").replace(".","_"), in_format, out_format),
                funct(filename, in_format, out_format, alphabet))
    del funct

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
