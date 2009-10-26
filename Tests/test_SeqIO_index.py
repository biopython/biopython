# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Additional unit tests for Bio.SeqIO.convert(...) function."""
import os
import unittest
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO._index import _FormatToIndexedDict
from Bio.Alphabet import generic_protein, generic_nucleotide, generic_dna

BINARY_FORMATS = ["sff", "sff-trim"]

class IndexDictTests(unittest.TestCase) :
    """Cunning unit test where methods are added at run time."""
    def simple_check(self, filename, format, alphabet) :
        if format in BINARY_FORMATS :
            mode = "rb"
        else :
            mode = "r"
        id_list = [rec.id for rec in \
                   SeqIO.parse(open(filename, mode), format, alphabet)]
        rec_dict = SeqIO.index(filename, format, alphabet)
        self.assertEqual(set(id_list), set(rec_dict.keys()))
        #This is redundant, I just want to make sure len works:
        self.assertEqual(len(id_list), len(rec_dict))
        #Make sure boolean evaluation works
        self.assertEqual(bool(id_list), bool(rec_dict))
        for key in id_list :
            self.assert_(key in rec_dict)
            self.assertEqual(key, rec_dict[key].id)
            self.assertEqual(key, rec_dict.get(key).id)
        #Check non-existant keys,
        try :
            rec = rec_dict[chr(0)]
            raise ValueError("Accessing a non-existant key should fail")
        except KeyError :
            pass
        self.assertEqual(rec_dict.get(chr(0)), None)
        self.assertEqual(rec_dict.get(chr(0), chr(1)), chr(1))
        #Now check iteritems...
        for key, rec in rec_dict.iteritems() :
            self.assert_(key in id_list)
            self.assert_(isinstance(rec, SeqRecord))
            self.assertEqual(rec.id, key)
        #Now check non-defined methods...
        self.assertRaises(NotImplementedError, rec_dict.values)
        self.assertRaises(NotImplementedError, rec_dict.popitem)
        self.assertRaises(NotImplementedError, rec_dict.pop, chr(0))
        self.assertRaises(NotImplementedError, rec_dict.pop, chr(0), chr(1))
        self.assertRaises(NotImplementedError, rec_dict.clear)
        self.assertRaises(NotImplementedError, rec_dict.__setitem__, "X", None)
        self.assertRaises(NotImplementedError, rec_dict.copy)
        self.assertRaises(NotImplementedError, rec_dict.fromkeys, [])
        #Done
            
tests = [
    ("Ace/contig1.ace", "ace", generic_dna),
    ("Ace/consed_sample.ace", "ace", None),
    ("Ace/seq.cap.ace", "ace", generic_dna),
    ("Quality/wrapping_original_sanger.fastq", "fastq", None),
    ("Quality/example.fastq", "fastq", None),
    ("Quality/example.fastq", "fastq-sanger", generic_dna),
    ("Quality/tricky.fastq", "fastq", generic_nucleotide),
    ("Quality/sanger_faked.fastq", "fastq-sanger", generic_dna),
    ("Quality/solexa_faked.fastq", "fastq-solexa", generic_dna),
    ("Quality/illumina_faked.fastq", "fastq-illumina", generic_dna),
    ("EMBL/U87107.embl", "embl", None),
    ("EMBL/TRBG361.embl", "embl", None),
    ("GenBank/NC_000932.faa", "fasta", generic_protein),
    ("GenBank/NC_005816.faa", "fasta", generic_protein),
    ("GenBank/NC_005816.ffn", "fasta", generic_dna),
    ("GenBank/NC_005816.fna", "fasta", generic_dna),
    ("GenBank/NC_005816.gb", "gb", None),
    ("GenBank/cor6_6.gb", "genbank", None),
    ("IntelliGenetics/vpu_nucaligned.txt", "ig", generic_nucleotide),
    ("IntelliGenetics/TAT_mase_nuc.txt", "ig", None),
    ("IntelliGenetics/VIF_mase-pro.txt", "ig", generic_protein),
    ("Phd/phd1", "phd", generic_dna),
    ("Phd/phd2", "phd", None),
    ("Phd/phd_solexa", "phd", generic_dna),
    ("Phd/phd_454", "phd", generic_dna),
    ("NBRF/B_nuc.pir", "pir", generic_nucleotide),
    ("NBRF/Cw_prot.pir", "pir", generic_protein),
    ("NBRF/clustalw.pir", "pir", None),
    ("SwissProt/sp001", "swiss", None),
    ("SwissProt/sp010", "swiss", None),
    ("SwissProt/sp016", "swiss", None),
    ("Roche/E3MFGYR02_random_10_reads.sff", "sff", generic_dna),
    ("Roche/E3MFGYR02_index_at_start.sff", "sff", generic_dna),
    ("Roche/E3MFGYR02_index_in_middle.sff", "sff", generic_dna),
    ("Roche/greek.sff", "sff", generic_nucleotide),
    ("Roche/paired.sff", "sff", None),
    ]
for filename, format, alphabet in tests :
    assert format in _FormatToIndexedDict
    def funct(fn,fmt,alpha) :
        f = lambda x : x.simple_check(fn, fmt, alpha)
        f.__doc__ = "Index %s file %s" % (fmt, fn)
        return f
    setattr(IndexDictTests, "test_%s_%s" \
            % (filename.replace("/","_").replace(".","_"), format),
            funct(filename, format, alphabet))
    del funct

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
