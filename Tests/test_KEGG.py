"""Tests the basic functionality of the KEGG parsers.
"""

import os

from Bio.KEGG import Enzyme
from Bio.KEGG import Compound
from Bio.KEGG import Map
from Bio.Pathway import Reaction
from Bio.Pathway import System

test_KEGG_Enzyme_files   = ["enzyme.sample", "enzyme.irregular"]
test_KEGG_Compound_files = ["compound.sample", "compound.irregular"] 
test_KEGG_Map_files      = ["map00950.rea"]

def t_KEGG_Enzyme(testfiles):
    """Tests Bio.KEGG.Enzyme functionality."""
    for file in testfiles:
        fh = open(os.path.join("KEGG", file))
        print "Testing Bio.KEGG.Enzyme on " + file + "\n\n"
        records = Enzyme.Iterator(fh, Enzyme.Parser(debug_level=0))
        while 1:
            record = records.next()
            if record is not None:
                print record
            else:
                break
        print "\n"

def t_KEGG_Compound(testfiles):
    """Tests Bio.KEGG.Compound functionality."""
    for file in testfiles:
        fh = open(os.path.join("KEGG", file))
        print "Testing Bio.KEGG.Compound on " + file + "\n\n"
        records = Compound.Iterator(fh, Compound.Parser(debug_level=0))
        while 1:
            record = records.next()
            if record is not None:
                print record
            else:
                break
        print "\n"    

def t_KEGG_Map(testfiles):
    """Tests Bio.KEGG.Map functionality."""
    for file in testfiles:
        fh = open(os.path.join("KEGG", file))
        print "Testing Bio.KEGG.Map on " + file + "\n\n"
        reactions = Map.Iterator(fh, Map.Parser(debug_level=0))
        system = System()
        while 1:
            r = reactions.next()
            if r is not None:
                system.add_reaction(r)
            else:
                break
        # sort the reaction output by the string names, so that the
        # output will be consistent between python versions
        def str_cmp(first, second):
            return cmp(str(first), str(second))
        rxs = system.reactions()
        rxs.sort(str_cmp)
        for x in rxs:
            print str(x)


t_KEGG_Enzyme(test_KEGG_Enzyme_files)
t_KEGG_Compound(test_KEGG_Compound_files)
t_KEGG_Map(test_KEGG_Map_files)


