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
from Bio.KEGG import REST
from Bio.Pathway import System

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
    print("Testing Bio.KEGG.query\n\n")

    # info tests
    resp = REST.kegg_info("kegg")
    resp.read()
    print(resp.url)

    resp = REST.kegg_info("pathway")
    resp.read()
    print(resp.url)

    # list tests
    resp = REST.kegg_list("pathway")
    resp.read()
    print(resp.url)

    resp = REST.kegg_list("pathway", "hsa")
    resp.read()
    print(resp.url)

    resp = REST.kegg_list("organism")
    resp.read()
    print(resp.url)

    resp = REST.kegg_list("hsa")
    resp.read()
    print(resp.url)

    resp = REST.kegg_list("T01001")
    resp.read()
    print(resp.url)

    resp = REST.kegg_list("hsa:10458+ece:Z5100")
    resp.read()
    print(resp.url)

    resp = REST.kegg_list(["hsa:10458", "ece:Z5100"])
    resp.read()
    print(resp.url)

    resp = REST.kegg_list("cpd:C01290+gl:G00092")
    resp.read()
    print(resp.url)

    resp = REST.kegg_list(["cpd:C01290", "gl:G00092"])
    resp.read()
    print(resp.url)

    resp = REST.kegg_list("C01290+G00092")
    resp.read()
    print(resp.url)

    resp = REST.kegg_list(["C01290", "G00092"])
    resp.read()
    print(resp.url)

    # find tests
    resp = REST.kegg_find("genes", "shiga+toxin")
    resp.read()
    print(resp.url)

    resp = REST.kegg_find("genes", ["shiga", "toxin"])
    resp.read()
    print(resp.url)

    resp = REST.kegg_find("compound", "C7H10O5", "formula")
    resp.read()
    print(resp.url)

    resp = REST.kegg_find("compound", "O5C7", "formula")
    resp.read()
    print(resp.url)

    resp = REST.kegg_find("compound", "174.05", "exact_mass")
    resp.read()
    print(resp.url)

    resp = REST.kegg_find("compound", "300-310", "mol_weight")
    resp.read()
    print(resp.url)

    # get tests
    resp = REST.kegg_get("cpd:C01290+gl:G00092")
    resp.read()
    print(resp.url)

    resp = REST.kegg_get(["cpd:C01290", "gl:G00092"])
    resp.read()
    print(resp.url)

    resp = REST.kegg_get("C01290+G00092")
    resp.read()
    print(resp.url)

    resp = REST.kegg_get(["C01290", "G00092"])
    resp.read()
    print(resp.url)

    resp = REST.kegg_get("hsa:10458+ece:Z5100")
    resp.read()
    print(resp.url)

    resp = REST.kegg_get(["hsa:10458", "ece:Z5100"])
    resp.read()
    print(resp.url)

    resp = REST.kegg_get("hsa:10458+ece:Z5100", "aaseq")
    resp.read()
    print(resp.url)

    resp = REST.kegg_get(["hsa:10458", "ece:Z5100"], "aaseq")
    resp.read()
    print(resp.url)

    resp = REST.kegg_get("hsa05130", "image")
    resp.read()
    print(resp.url)

    # conv tests
    resp = REST.kegg_conv("eco", "ncbi-geneid")
    resp.read()
    print(resp.url)

    resp = REST.kegg_conv("ncbi-geneid", "eco")
    resp.read()
    print(resp.url)

    resp = REST.kegg_conv("ncbi-gi", "hsa:10458+ece:Z5100")
    resp.read()
    print(resp.url)

    resp = REST.kegg_conv("ncbi-gi", ["hsa:10458", "ece:Z5100"])
    resp.read()
    print(resp.url)

    # link tests
    resp = REST.kegg_link("pathway", "hsa")
    resp.read()
    print(resp.url)

    resp = REST.kegg_link("hsa", "pathway")
    resp.read()
    print(resp.url)

    resp = REST.kegg_link("pathway", "hsa:10458+ece:Z5100")
    resp.read()
    print(resp.url)

    resp = REST.kegg_link("pathway", ["hsa:10458", "ece:Z5100"])
    resp.read()
    print(resp.url)


t_KEGG_Enzyme(test_KEGG_Enzyme_files)
t_KEGG_Compound(test_KEGG_Compound_files)
t_KEGG_Map(test_KEGG_Map_files)
t_KEGG_Query()
