"""Tests the basic functionality of the KEGG parsers.
"""

import os

from Bio.KEGG import Enzyme
from Bio.KEGG import Compound

test_KEGG_Enzyme_files   = ["enzyme.sample", "enzyme.irregular"]
test_KEGG_Compound_files = ["compound.sample", "compound.irregular"] 

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

t_KEGG_Enzyme(test_KEGG_Enzyme_files)
t_KEGG_Compound(test_KEGG_Compound_files)



