# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# Revisions copyright 2007 by Michiel de Hoon. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests the basic functionality of the KEGG parsers.
"""

from __future__ import print_function

import os

from Bio.KEGG import Enzyme
from Bio.KEGG import Compound
from Bio.KEGG import Map
from Bio.Pathway import System
from Bio import KEGG

test_KEGG_Enzyme_files   = ["enzyme.sample", "enzyme.irregular", "enzyme.new"]
test_KEGG_Compound_files = ["compound.sample", "compound.irregular"]
test_KEGG_Map_files      = ["map00950.rea"]


def t_KEGG_Enzyme(testfiles):
    """Tests Bio.KEGG.Enzyme functionality."""
    for file in testfiles:
        fh = open(os.path.join("KEGG", file))
        print("Testing Bio.KEGG.Enzyme on " + file + "\n\n")
        records = Enzyme.parse(fh)
        for record in records:
            print(record)
        print("\n")
        fh.close()


def t_KEGG_Compound(testfiles):
    """Tests Bio.KEGG.Compound functionality."""
    for file in testfiles:
        fh = open(os.path.join("KEGG", file))
        print("Testing Bio.KEGG.Compound on " + file + "\n\n")
        records = Compound.parse(fh)
        for record in records:
            print(record)
        print("\n")
        fh.close()


def t_KEGG_Map(testfiles):
    """Tests Bio.KEGG.Map functionality."""
    for file in testfiles:
        fh = open(os.path.join("KEGG", file))
        print("Testing Bio.KEGG.Map on " + file + "\n\n")
        reactions = Map.parse(fh)
        system = System()
        for reaction in reactions:
            system.add_reaction(reaction)
        # sort the reaction output by the string names, so that the
        # output will be consistent between python versions
        # def str_cmp(first, second):
        #    return cmp(str(first), str(second))
        rxs = system.reactions()
        # sort: key instead of compare function (for py3 support)
        #  The function str_cmp above can be removed if the
        #  solution below proves resilient
        rxs.sort(key=lambda x:str(x))
        for x in rxs:
            print(str(x))
        print("\n")
        fh.close()

def t_KEGG_Query():
    """Tests Bio.KEGG API Wrapper"""
    print "Testing Bio.KEGG.query\n\n"

    # info tests
    resp = KEGG.info("kegg")
    resp.read()
    print resp.url

    resp = KEGG.info("pathway")
    resp.read()
    print resp.url

    # list tests
    resp = KEGG.list_("pathway")
    resp.read()
    print resp.url

    resp = KEGG.list_("pathway", "hsa")
    resp.read()
    print resp.url

    resp = KEGG.list_("organism")
    resp.read()
    print resp.url

    resp = KEGG.list_("hsa")
    resp.read()
    print resp.url

    resp = KEGG.list_("T01001")
    resp.read()
    print resp.url

    resp = KEGG.list_("hsa:10458+ece:Z5100")
    resp.read()
    print resp.url

    resp = KEGG.list_(["hsa:10458", "ece:Z5100"])
    resp.read()
    print resp.url

    resp = KEGG.list_("cpd:C01290+gl:G00092")
    resp.read()
    print resp.url

    resp = KEGG.list_(["cpd:C01290", "gl:G00092"])
    resp.read()
    print resp.url

    resp = KEGG.list_("C01290+G00092")
    resp.read()
    print resp.url

    resp = KEGG.list_(["C01290", "G00092"])
    resp.read()
    print resp.url

    # find tests
    resp = KEGG.find("genes", "shiga+toxin")
    resp.read()
    print resp.url

    resp = KEGG.find("genes", ["shiga", "toxin"])
    resp.read()
    print resp.url

    resp = KEGG.find("compound", "C7H10O5", "formula")
    resp.read()
    print resp.url

    resp = KEGG.find("compound", "O5C7", "formula")
    resp.read()
    print resp.url

    resp = KEGG.find("compound", "174.05", "exact_mass")
    resp.read()
    print resp.url

    resp = KEGG.find("compound", "300-310", "mol_weight")
    resp.read()
    print resp.url

    # get tests
    resp = KEGG.get("cpd:C01290+gl:G00092")
    resp.read()
    print resp.url

    resp = KEGG.get(["cpd:C01290", "gl:G00092"])
    resp.read()
    print resp.url

    resp = KEGG.get("C01290+G00092")
    resp.read()
    print resp.url

    resp = KEGG.get(["C01290", "G00092"])
    resp.read()
    print resp.url

    resp = KEGG.get("hsa:10458+ece:Z5100")
    resp.read()
    print resp.url

    resp = KEGG.get(["hsa:10458", "ece:Z5100"])
    resp.read()
    print resp.url

    resp = KEGG.get("hsa:10458+ece:Z5100", "aaseq")
    resp.read()
    print resp.url

    resp = KEGG.get(["hsa:10458", "ece:Z5100"], "aaseq")
    resp.read()
    print resp.url

    resp = KEGG.get("hsa05130", "image")
    resp.read()
    print resp.url

    # conv tests
    resp = KEGG.conv("eco", "ncbi-geneid")
    resp.read()
    print resp.url

    resp = KEGG.conv("ncbi-geneid", "eco")
    resp.read()
    print resp.url

    resp = KEGG.conv("ncbi-gi", "hsa:10458+ece:Z5100")
    resp.read()
    print resp.url

    resp = KEGG.conv("ncbi-gi", ["hsa:10458", "ece:Z5100"])
    resp.read()
    print resp.url

    # link tests
    resp = KEGG.link("pathway", "hsa")
    resp.read()
    print resp.url

    resp = KEGG.link("hsa", "pathway")
    resp.read()
    print resp.url

    resp = KEGG.link("pathway", "hsa:10458+ece:Z5100")
    resp.read()
    print resp.url

    resp = KEGG.link("pathway", ["hsa:10458", "ece:Z5100"])
    resp.read()
    print resp.url


t_KEGG_Enzyme(test_KEGG_Enzyme_files)
t_KEGG_Compound(test_KEGG_Compound_files)
t_KEGG_Map(test_KEGG_Map_files)
t_KEGG_Query()
