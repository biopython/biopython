# Copyright 2009-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Additional unit tests for Bio.SeqIO.index(...) function."""
import sys
if sys.version_info[0] >= 3:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "Skipping since currently this is very slow on Python 3.")
    
import os
import unittest
from StringIO import StringIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO._index import _FormatToIndexedDict
from Bio.Alphabet import generic_protein, generic_nucleotide, generic_dna

from seq_tests_common import compare_record

def add_prefix(key):
    """Dummy key_function for testing index code."""
    return "id_" + key

class IndexDictTests(unittest.TestCase):
    """Cunning unit test where methods are added at run time."""
    def simple_check(self, filename, format, alphabet):
        if format in SeqIO._BinaryFormats:
            mode = "rb"
        else :
            mode = "r"
        id_list = [rec.id for rec in \
                   SeqIO.parse(open(filename, mode), format, alphabet)]
        #Without key_function
        rec_dict = SeqIO.index(filename, format, alphabet)
        self.check_dict_methods(rec_dict, id_list, id_list)
        #Check with key_function
        key_list = [add_prefix(id) for id in id_list]
        rec_dict = SeqIO.index(filename, format, alphabet, add_prefix)
        self.check_dict_methods(rec_dict, key_list, id_list)
        #Done
    
    def check_dict_methods(self, rec_dict, keys, ids):
        self.assertEqual(set(keys), set(rec_dict.keys()))
        #This is redundant, I just want to make sure len works:
        self.assertEqual(len(keys), len(rec_dict))
        #Make sure boolean evaluation works
        self.assertEqual(bool(keys), bool(rec_dict))
        for key,id in zip(keys, ids):
            self.assertTrue(key in rec_dict)
            self.assertEqual(id, rec_dict[key].id)
            self.assertEqual(id, rec_dict.get(key).id)
        #Check non-existant keys,
        assert chr(0) not in keys, "Bad example in test"
        try:
            rec = rec_dict[chr(0)]
            raise ValueError("Accessing a non-existent key should fail")
        except KeyError:
            pass
        self.assertEqual(rec_dict.get(chr(0)), None)
        self.assertEqual(rec_dict.get(chr(0), chr(1)), chr(1))
        if hasattr(dict, "iteritems"):
            #Python 2.x
            for key, rec in rec_dict.iteritems():
                self.assertTrue(key in keys)
                self.assertTrue(isinstance(rec, SeqRecord))
                self.assertTrue(rec.id in ids)
            #Now check non-defined methods...
            self.assertRaises(NotImplementedError, rec_dict.items)
            self.assertRaises(NotImplementedError, rec_dict.values)
        else:
            #Python 3
            assert not hasattr(rec_dict, "iteritems")
            for key, rec in rec_dict.iteritems():
                self.assertTrue(key in id_list)
                self.assertTrue(isinstance(rec, SeqRecord))
                self.assertTrue(rec.id in ids)
            for rec in rec_dict.itervalues():
                self.assertTrue(key in id_list)
                self.assertTrue(isinstance(rec, SeqRecord))
                self.assertTrue(rec.id in ids)
        #Check the following fail
        self.assertRaises(NotImplementedError, rec_dict.popitem)
        self.assertRaises(NotImplementedError, rec_dict.pop, chr(0))
        self.assertRaises(NotImplementedError, rec_dict.pop, chr(0), chr(1))
        self.assertRaises(NotImplementedError, rec_dict.clear)
        self.assertRaises(NotImplementedError, rec_dict.__setitem__, "X", None)
        self.assertRaises(NotImplementedError, rec_dict.copy)
        self.assertRaises(NotImplementedError, rec_dict.fromkeys, [])

    def get_raw_check(self, filename, format, alphabet):
        if format in SeqIO._BinaryFormats:
            #This means SFF at the moment, which does not get
            #implement the get_raw method
            return
        handle = open(filename, "rU")
        raw_file = handle.read()
        handle.close()
        #Also checking the key_function here
        id_list = [rec.id.lower() for rec in \
                   SeqIO.parse(filename, format, alphabet)]
        rec_dict = SeqIO.index(filename, format, alphabet,
                               key_function = lambda x : x.lower())
        self.assertEqual(set(id_list), set(rec_dict.keys()))
        self.assertEqual(len(id_list), len(rec_dict))
        for key in id_list:
            self.assertTrue(key in rec_dict)
            self.assertEqual(key, rec_dict[key].id.lower())
            self.assertEqual(key, rec_dict.get(key).id.lower())
            raw = rec_dict.get_raw(key)
            self.assertTrue(raw.strip())
            self.assertTrue(raw in raw_file)
            if format in ["ig", "uniprot-xml"]:
               #These have a header structure and can't be parsed
               #individually (at least, not right now).
               continue
            rec1 = rec_dict[key]
            rec2 = SeqIO.read(StringIO(raw), format, alphabet)
            self.assertEqual(True, compare_record(rec1, rec2))

    def test_duplicates_index(self):
        """Index file with duplicate identifers with Bio.SeqIO.index()"""
        self.assertRaises(ValueError, SeqIO.index, "Fasta/dups.fasta", "fasta")

    def test_duplicates_to_dict(self):
        """Index file with duplicate identifers with Bio.SeqIO.to_dict()"""
        handle = open("Fasta/dups.fasta", "rU")
        iterator = SeqIO.parse(handle, "fasta")
        self.assertRaises(ValueError, SeqIO.to_dict, iterator)
        handle.close()

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
    ("EMBL/epo_prt_selection.embl", "embl", None),
    ("EMBL/U87107.embl", "embl", None),
    ("EMBL/TRBG361.embl", "embl", None),
    ("EMBL/A04195.imgt", "embl", None), #Not a proper EMBL file, an IMGT file
    ("EMBL/A04195.imgt", "imgt", None),
    ("GenBank/NC_000932.faa", "fasta", generic_protein),
    ("GenBank/NC_005816.faa", "fasta", generic_protein),
    ("GenBank/NC_005816.tsv", "tab", generic_protein),
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
    ("SwissProt/multi_ex.txt", "swiss", None),
    ("SwissProt/multi_ex.xml", "uniprot-xml", None),
    ("SwissProt/multi_ex.fasta", "fasta", None),
    ("Roche/E3MFGYR02_random_10_reads.sff", "sff", generic_dna),
    ("Roche/E3MFGYR02_random_10_reads.sff", "sff-trim", generic_dna),
    ("Roche/E3MFGYR02_index_at_start.sff", "sff", generic_dna),
    ("Roche/E3MFGYR02_index_in_middle.sff", "sff", generic_dna),
    ("Roche/E3MFGYR02_alt_index_at_start.sff", "sff", generic_dna),
    ("Roche/E3MFGYR02_alt_index_in_middle.sff", "sff", generic_dna),
    ("Roche/E3MFGYR02_alt_index_at_end.sff", "sff", generic_dna),
    ("Roche/E3MFGYR02_no_manifest.sff", "sff", generic_dna),
    ("Roche/greek.sff", "sff", generic_nucleotide),
    ("Roche/greek.sff", "sff-trim", generic_nucleotide),
    ("Roche/paired.sff", "sff", None),
    ("Roche/paired.sff", "sff-trim", None),
    ]
for filename, format, alphabet in tests:
    assert format in _FormatToIndexedDict
    def funct(fn,fmt,alpha):
        f = lambda x : x.simple_check(fn, fmt, alpha)
        f.__doc__ = "Index %s file %s" % (fmt, fn)
        return f
    setattr(IndexDictTests, "test_%s_%s" \
            % (filename.replace("/","_").replace(".","_"), format),
            funct(filename, format, alphabet))
    del funct

    if format in SeqIO._BinaryFormats:
        continue

    def funct(fn,fmt,alpha):
        f = lambda x : x.get_raw_check(fn, fmt, alpha)
        f.__doc__ = "Index %s file %s get_raw" % (fmt, fn)
        return f
    setattr(IndexDictTests, "test_%s_%s_get_raw" \
            % (filename.replace("/","_").replace(".","_"), format),
            funct(filename, format, alphabet))
    del funct

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
