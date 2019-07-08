# Copyright 2014 by Kevin Wu.
# Copyright 2014 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for online functionality of the KEGG module."""

# Builtins
import unittest

from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.REST import kegg_conv, kegg_find, kegg_get
from Bio.KEGG.REST import kegg_info, kegg_link, kegg_list

from Bio import SeqIO

import requires_internet
requires_internet.check()

# TODO - revert to using with statements once we drop
# Python 2.6 and 2.7, see http://bugs.python.org/issue12487


class KEGGTests(unittest.TestCase):
    """Tests for KEGG REST API."""

    def test_info_kegg(self):
        h = kegg_info("kegg")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/info/kegg")
        h.close()

    def test_info_pathway(self):
        h = kegg_info("pathway")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/info/pathway")
        h.close()

    def test_list_pathway(self):
        h = kegg_list("pathway")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/pathway")
        h.close()

    def test_pathway_hsa(self):
        h = kegg_list("pathway", "hsa")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/pathway/hsa")
        h.close()

    def test_list_organism(self):
        h = kegg_list("organism")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/organism")
        h.close()

    def test_list_hsa(self):
        h = kegg_list("hsa")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/hsa")
        h.close()

    def test_list_T01001(self):
        h = kegg_list("T01001")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/T01001")
        h.close()

    def test_list_hsa_10458_plus_ece_Z5100(self):
        h = kegg_list("hsa:10458+ece:Z5100")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/hsa:10458+ece:Z5100")
        h.close()

    def test_list_hsa_10458_list_ece_Z5100(self):
        h = kegg_list(["hsa:10458", "ece:Z5100"])
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/hsa:10458+ece:Z5100")
        h.close()

    def test_list_cpd_C01290_plus_gl_G0009(self):
        h = kegg_list("cpd:C01290+gl:G00092")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/cpd:C01290+gl:G00092")
        h.close()

    def test_list_cpd_C01290_list_gl_G0009(self):
        h = kegg_list(["cpd:C01290", "gl:G00092"])
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/cpd:C01290+gl:G00092")
        h.close()

    def test_list_C01290_plus_G00092(self):
        h = kegg_list("C01290+G00092")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/C01290+G00092")
        h.close()

    def test_list_C01290_list_G00092(self):
        h = kegg_list(["C01290", "G00092"])
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/list/C01290+G00092")
        h.close()

    def test_find_genes_shiga_plus_toxin(self):
        h = kegg_find("genes", "shiga+toxin")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/find/genes/shiga+toxin")
        h.close()

    def test_find_genes_shiga_list_toxin(self):
        h = kegg_find("genes", ["shiga", "toxin"])
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/find/genes/shiga+toxin")
        h.close()

    def test_find_compound_C7H10O5_formula(self):
        h = kegg_find("compound", "C7H10O5", "formula")
        h.read()
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/find/compound/C7H10O5/formula")
        h.close()

    def test_find_compound_O5C7_formula(self):
        h = kegg_find("compound", "O5C7", "formula")
        h.read()
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/find/compound/O5C7/formula")
        h.close()

    def test_find_compound_exact_mass(self):
        h = kegg_find("compound", "174.05", "exact_mass")
        h.read()
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/find/compound/174.05/exact_mass")
        h.close()

    def test_find_compound_weight(self):
        h = kegg_find("compound", "300-310", "mol_weight")
        h.read()
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/find/compound/300-310/mol_weight")
        h.close()

    def test_get_cpd_C01290_plus_gl_G00092(self):
        h = kegg_get("cpd:C01290+gl:G00092")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/get/cpd:C01290+gl:G00092")
        h.close()

    def test_get_cpd_C01290_list_gl_G00092(self):
        h = kegg_get(["cpd:C01290", "gl:G00092"])
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/get/cpd:C01290+gl:G00092")
        h.close()

    def test_get_C01290_plus_G00092(self):
        h = kegg_get(["C01290+G00092"])
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/get/C01290+G00092")
        h.close()

    def test_get_C01290_list_G00092(self):
        h = kegg_get(["C01290", "G00092"])
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/get/C01290+G00092")
        h.close()

    def test_get_hsa_10458_plus_ece_Z5100(self):
        h = kegg_get("hsa:10458+ece:Z5100")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/get/hsa:10458+ece:Z5100")
        h.close()

    def test_get_hsa_10458_list_ece_Z5100(self):
        h = kegg_get(["hsa:10458", "ece:Z5100"])
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/get/hsa:10458+ece:Z5100")
        h.close()

    def test_get_hsa_10458_plus_ece_Z5100_as_aaseq(self):
        h = kegg_get("hsa:10458+ece:Z5100", "aaseq")
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/get/hsa:10458+ece:Z5100/aaseq")
        data = SeqIO.parse(h, 'fasta')
        self.assertEqual(len(list(data)), 2)
        h.close()

    def test_get_hsa_10458_list_ece_Z5100_as_aaseq(self):
        h = kegg_get(["hsa:10458", "ece:Z5100"], "aaseq")
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/get/hsa:10458+ece:Z5100/aaseq")
        data = SeqIO.parse(h, 'fasta')
        self.assertEqual(len(list(data)), 2)
        h.close()

    def test_get_hsa_10458_plus_ece_Z5100_as_ntseq(self):
        h = kegg_get("hsa:10458+ece:Z5100", "ntseq")
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/get/hsa:10458+ece:Z5100/ntseq")
        data = SeqIO.parse(h, 'fasta')
        self.assertEqual(len(list(data)), 2)
        h.close()

    def test_get_hsa_10458_list_ece_Z5100_as_ntseq(self):
        h = kegg_get(["hsa:10458", "ece:Z5100"], "ntseq")
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/get/hsa:10458+ece:Z5100/ntseq")
        data = SeqIO.parse(h, 'fasta')
        self.assertEqual(len(list(data)), 2)
        h.close()

    def test_get_hsa05130_image(self):
        h = kegg_get("hsa05130", "image")
        data = h.read()
        self.assertEqual(data[:4], b"\x89PNG")
        self.assertEqual(h.url, "http://rest.kegg.jp/get/hsa05130/image")
        h.close()

    def test_conv_eco_ncbi_geneid(self):
        h = kegg_conv("eco", "ncbi-geneid")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/conv/eco/ncbi-geneid")
        h.close()

    def test_conv_ncbi_geneid_eco(self):
        h = kegg_conv("ncbi-geneid", "eco")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/conv/ncbi-geneid/eco")
        h.close()

    def test_conv_ncbi_gi_hsa_10458_plus_ece_Z5100(self):
        h = kegg_conv("ncbi-gi", "hsa:10458+ece:Z5100")
        h.read()
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/conv/ncbi-gi/hsa:10458+ece:Z5100")
        h.close()

    def test_conv_ncbi_gi_hsa_10458_list_ece_Z5100(self):
        h = kegg_conv("ncbi-gi", ["hsa:10458", "ece:Z5100"])
        h.read()
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/conv/ncbi-gi/hsa:10458+ece:Z5100")
        h.close()

    def test_link_pathway_hsa(self):
        h = kegg_link("pathway", "hsa")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/link/pathway/hsa")
        h.close()

    def test_link_hsa_pathway(self):
        h = kegg_link("hsa", "pathway")
        h.read()
        self.assertEqual(h.url, "http://rest.kegg.jp/link/hsa/pathway")
        h.close()

    def test_pathway_hsa_10458_plus_ece_Z5100(self):
        h = kegg_link("pathway", "hsa:10458+ece:Z5100")
        h.read()
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/link/pathway/hsa:10458+ece:Z5100")
        h.close()

    def test_pathway_hsa_10458_list_ece_Z5100(self):
        h = kegg_link("pathway", ["hsa:10458", "ece:Z5100"])
        h.read()
        self.assertEqual(h.url,
                         "http://rest.kegg.jp/link/pathway/hsa:10458+ece:Z5100")
        h.close()


class KGMLPathwayTests(unittest.TestCase):
    """Tests with metabolic maps."""

    def test_parse_remote_pathway(self):
        """Download a KEGG pathway from the KEGG server and parse KGML."""
        h = kegg_get("ko03070", "kgml")
        pathway = KGML_parser.read(h)
        self.assertEqual(pathway.name, "path:ko03070")
        h.close()


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
