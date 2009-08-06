# Copyright 2009 by Peter Cock.  All rights reserved.
# Parts copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.Alphabet import generic_protein, generic_nucleotide

def title_to_ids(title):
    """Function to convert a title into the id, name, and description.

    This is just a quick-n-dirty implementation, and is definately not meant
    to handle every FASTA title line case.
    """
    # first split the id information from the description
    # the first item is the id info block, the rest is the description
    all_info = title.split(" ")
    id_info = all_info[0]
    rest = all_info[1:]
    descr = " ".join(rest)

    # now extract the ids from the id block
    # gi|5690369|gb|AF158246.1|AF158246
    id_info_items = id_info.split("|")
    assert id_info_items[2] in ["gb", "emb", "dbj", "pdb"], title
    id = id_info_items[3] # the id with version info
    name = id_info_items[4] # the id without version info

    return id, name, descr

def read_single_with_titles(filename, alphabet) :
    global title_to_ids
    iterator = FastaIterator(open(filename), alphabet, title_to_ids)
    record = iterator.next()
    try :
        second = iterator.next()
    except StopIteration :
        second = None
    assert record is not None and second is None
    return record

def read_title_and_seq(filename) :
    """Crude parser that gets the first record from a FASTA file."""
    handle = open(filename)
    title = handle.readline().rstrip()
    assert title.startswith(">")
    seq = ""
    for line in handle :
        if line.startswith(">") : break
        seq += line.strip()
    handle.close()
    return title[1:], seq


class TitleFunctions(unittest.TestCase) :
    """Cunning unit test where methods are added at run time."""
    def simple_check(self, filename, alphabet) :
        title, seq = read_title_and_seq(filename) #crude parser
        #First check using Bio.SeqIO.FastaIO directly with title function,
        record = read_single_with_titles(filename, alphabet)
        idn, name, descr = title_to_ids(title)
        self.assertEqual(record.id, idn)
        self.assertEqual(record.name, name)
        self.assertEqual(record.description, descr)
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.seq.alphabet, alphabet)
        #Now check using Bio.SeqIO (default settings)
        record = SeqIO.read(open(filename), "fasta", alphabet)
        self.assertEqual(record.id, title.split()[0])
        self.assertEqual(record.name, title.split()[0])
        self.assertEqual(record.description, title)
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.seq.alphabet, alphabet)
        #Uncomment this for tested the methods are calling the right files:
        #print "{%s done}" % filename,

nucleic_files = ['lupine.nu', 'elderberry.nu', 'phlox.nu', 'centaurea.nu',
                 'wisteria.nu', 'sweetpea.nu', 'lavender.nu']

amino_files = ['aster.pro', 'rosemary.pro', 'rose.pro', 'loveliesbleeding.pro']

for filename in nucleic_files :
    name = filename.split(".")[0]
    def funct(fn) :
        f = lambda x : x.simple_check(fn, generic_nucleotide)
        f.__doc__ = "Checking nucleotide file %s" % fn
        return f
    setattr(TitleFunctions, "test_nuc_%s"%name, funct("Nucleic/%s"%filename))
    del funct

for filename in amino_files :
    name = filename.split(".")[0]
    def funct(fn) :
        f = lambda x : x.simple_check(fn, generic_nucleotide)
        f.__doc__ = "Checking protein file %s" % fn
        return f
    setattr(TitleFunctions, "test_pro_%s"%name, funct("Amino/%s"%filename))
    del funct

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
