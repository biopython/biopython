# Copyright 2014 by Kevin Wu.
# Copyright 2014 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for online functionality of the KEGG module."""

# Builtins
import os
import unittest
import tempfile

import requires_internet
requires_internet.check()

from Bio._py3k import _as_string

from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.REST import *


class KEGGTests(unittest.TestCase):
    """Tests for KEGG REST API."""

    def test_info_kegg(self):
        with kegg_info("kegg") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/info/kegg")

    def test_info_pathway(self):
        with kegg_info("pathway") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/info/pathway")

    def test_list_pathway(self):
        with kegg_list("pathway") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/pathway")

    def test_pathway_hsa(self):
        with kegg_list("pathway", "hsa") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/pathway/hsa")

    def test_list_organism(self):
        with kegg_list("organism") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/organism")

    def test_list_hsa(self):
        with kegg_list("hsa") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/hsa")

    def test_list_T01001(self):
        with kegg_list("T01001") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/T01001")

    def test_list_hsa_10458_plus_ece_Z5100(self):
        with kegg_list("hsa:10458+ece:Z5100") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/hsa:10458+ece:Z5100")

    def test_list_hsa_10458_list_ece_Z5100(self):
        with kegg_list(["hsa:10458", "ece:Z5100"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/hsa:10458+ece:Z5100")

    def test_list_cpd_C01290_plus_gl_G0009(self):
        with kegg_list("cpd:C01290+gl:G00092") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/cpd:C01290+gl:G00092")

    def test_list_cpd_C01290_list_gl_G0009(self):
        with kegg_list(["cpd:C01290", "gl:G00092"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/cpd:C01290+gl:G00092")

    def test_list_C01290_plus_G00092(self):
        with kegg_list("C01290+G00092") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/C01290+G00092")

    def test_list_C01290_list_G00092(self):
        with kegg_list(["C01290", "G00092"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/list/C01290+G00092")

    def test_find_genes_shiga_plus_toxin(self):
        with kegg_find("genes", "shiga+toxin") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/find/genes/shiga+toxin")

    def test_find_genes_shiga_list_toxin(self):
        with kegg_find("genes", ["shiga", "toxin"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/find/genes/shiga+toxin")

    def test_find_compound_C7H10O5_formula(self):
        with kegg_find("compound", "C7H10O5", "formula") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/find/compound/C7H10O5/formula")

    def test_find_compound_O5C7_formula(self):
        with kegg_find("compound", "O5C7", "formula") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/find/compound/O5C7/formula")

    def test_find_compound_exact_mass(self):
        with kegg_find("compound", "174.05", "exact_mass") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/find/compound/174.05/exact_mass")

    def test_find_compound_weight(self):
        with kegg_find("compound", "300-310", "mol_weight") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/find/compound/300-310/mol_weight")

    def test_get_cpd_C01290_plus_gl_G00092(self):
        with kegg_get("cpd:C01290+gl:G00092") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/get/cpd:C01290+gl:G00092")

    def test_get_cpd_C01290_list_gl_G00092(self):
        with kegg_get(["cpd:C01290", "gl:G00092"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/get/cpd:C01290+gl:G00092")

    def test_get_C01290_plus_G00092(self):
        with kegg_get(["C01290+G00092"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/get/C01290+G00092")

    def test_get_C01290_list_G00092(self):
        with kegg_get(["C01290", "G00092"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/get/C01290+G00092")

    def test_get_hsa_10458_plus_ece_Z5100(self):
        with kegg_get("hsa:10458+ece:Z5100") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/get/hsa:10458+ece:Z5100")

    def test_get_hsa_10458_list_ece_Z5100(self):
        with kegg_get(["hsa:10458", "ece:Z5100"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/get/hsa:10458+ece:Z5100")

    def test_get_hsa_10458_plus_ece_Z5100_as_aaseq(self):
        with kegg_get("hsa:10458+ece:Z5100", "aaseq") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/get/hsa:10458+ece:Z5100/aaseq")

    def test_get_hsa_10458_list_ece_Z5100_as_aaseq(self):
        with kegg_get(["hsa:10458", "ece:Z5100"], "aaseq") as h:
            data = _as_string(h.read())
            self.assertEqual(data.count(">"), 2)
            self.assertEqual(h.url, "http://rest.kegg.jp/get/hsa:10458+ece:Z5100/aaseq")

    def test_get_hsa05130_image(self):
        with kegg_get("hsa05130", "image") as h:
            data = h.read()
            self.assertEqual(data[:4], b"\x89PNG")
            self.assertEqual(h.url, "http://rest.kegg.jp/get/hsa05130/image")

    def test_conv_eco_ncbi_geneid(self):
        with kegg_conv("eco", "ncbi-geneid") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/conv/eco/ncbi-geneid")

    def test_conv_ncbi_geneid_eco(self):
        with kegg_conv("ncbi-geneid", "eco") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/conv/ncbi-geneid/eco")

    def test_conv_ncbi_gi_hsa_10458_plus_ece_Z5100(self):
        with kegg_conv("ncbi-gi", "hsa:10458+ece:Z5100") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/conv/ncbi-gi/hsa:10458+ece:Z5100")

    def test_conv_ncbi_gi_hsa_10458_list_ece_Z5100(self):
        with kegg_conv("ncbi-gi", ["hsa:10458", "ece:Z5100"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/conv/ncbi-gi/hsa:10458+ece:Z5100")

    def test_link_pathway_hsa(self):
        with kegg_link("pathway", "hsa") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/link/pathway/hsa")

    def test_link_hsa_pathway(self):
        with kegg_link("hsa", "pathway") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/link/hsa/pathway")

    def test_pathway_hsa_10458_plus_ece_Z5100(self):
        with kegg_link("pathway", "hsa:10458+ece:Z5100") as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/link/pathway/hsa:10458+ece:Z5100")

    def test_pathway_hsa_10458_list_ece_Z5100(self):
        with kegg_link("pathway", ["hsa:10458", "ece:Z5100"]) as h:
            h.read()
            self.assertEqual(h.url, "http://rest.kegg.jp/link/pathway/hsa:10458+ece:Z5100")


class KGMLPathwayTests(unittest.TestCase):
    """Tests with metabolic maps."""

    def test_parse_remote_pathway(self):
        """Download a KEGG pathway from the KEGG server and write KGML."""
        # Download the KEGG ko03070 pathway as a filehandle
        with kegg_get("ko03070", "kgml") as h:
            pathway = KGML_parser.read(h)
        self.assertEqual(pathway.name, "path:ko03070")


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
