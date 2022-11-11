# Copyright 2014 by Kevin Wu.
# Copyright 2014 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for online functionality of the KEGG module."""

# Builtins
import io
import unittest

from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.REST import kegg_conv, kegg_find, kegg_get
from Bio.KEGG.REST import kegg_info, kegg_link, kegg_list

from Bio import SeqIO

import requires_internet

requires_internet.check()


class KEGGTests(unittest.TestCase):
    """Tests for KEGG REST API."""

    def test_info_kegg(self):
        with kegg_info("kegg") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/info/kegg")

    def test_info_pathway(self):
        with kegg_info("pathway") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/info/pathway")

    def test_list_pathway(self):
        with kegg_list("pathway") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/pathway")

    def test_pathway_hsa(self):
        with kegg_list("pathway", "hsa") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/pathway/hsa")

    def test_list_organism(self):
        with kegg_list("organism") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/organism")

    def test_list_hsa(self):
        with kegg_list("hsa") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/hsa")

    def test_list_T01001(self):
        with kegg_list("T01001") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/T01001")

    def test_list_hsa_10458_plus_ece_Z5100(self):
        with kegg_list("hsa:10458+ece:Z5100") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/hsa:10458+ece:Z5100")

    def test_list_hsa_10458_list_ece_Z5100(self):
        with kegg_list(["hsa:10458", "ece:Z5100"]) as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/hsa:10458+ece:Z5100")

    def test_list_cpd_C01290_plus_gl_G0009(self):
        with kegg_list("cpd:C01290+gl:G00092") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/cpd:C01290+gl:G00092")

    def test_list_cpd_C01290_list_gl_G0009(self):
        with kegg_list(["cpd:C01290", "gl:G00092"]) as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/cpd:C01290+gl:G00092")

    def test_list_C01290_plus_G00092(self):
        with kegg_list("C01290+G00092") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/C01290+G00092")

    def test_list_C01290_list_G00092(self):
        with kegg_list(["C01290", "G00092"]) as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/list/C01290+G00092")

    def test_find_genes_shiga_plus_toxin(self):
        with kegg_find("genes", "shiga+toxin") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/find/genes/shiga+toxin")

    def test_find_genes_shiga_list_toxin(self):
        with kegg_find("genes", ["shiga", "toxin"]) as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/find/genes/shiga+toxin")

    def test_find_compound_C7H10O5_formula(self):
        with kegg_find("compound", "C7H10O5", "formula") as handle:
            handle.read()
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/find/compound/C7H10O5/formula"
        )

    def test_find_compound_O5C7_formula(self):
        with kegg_find("compound", "O5C7", "formula") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/find/compound/O5C7/formula")

    def test_find_compound_exact_mass(self):
        with kegg_find("compound", "174.05", "exact_mass") as handle:
            handle.read()
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/find/compound/174.05/exact_mass"
        )

    def test_find_compound_weight(self):
        with kegg_find("compound", "300-310", "mol_weight") as handle:
            handle.read()
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/find/compound/300-310/mol_weight"
        )

    def test_get_br_ko00002(self):
        with kegg_get("br:ko00002", "json") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/get/br:ko00002/json")

    def test_get_cpd_C01290_plus_gl_G00092(self):
        with kegg_get("cpd:C01290+gl:G00092") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/get/cpd:C01290+gl:G00092")

    def test_get_cpd_C01290_list_gl_G00092(self):
        with kegg_get(["cpd:C01290", "gl:G00092"]) as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/get/cpd:C01290+gl:G00092")

    def test_get_C01290_plus_G00092(self):
        with kegg_get(["C01290+G00092"]) as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/get/C01290+G00092")

    def test_get_C01290_list_G00092(self):
        with kegg_get(["C01290", "G00092"]) as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/get/C01290+G00092")

    def test_get_hsa_10458_plus_ece_Z5100(self):
        with kegg_get("hsa:10458+ece:Z5100") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/get/hsa:10458+ece:Z5100")

    def test_get_hsa_10458_list_ece_Z5100(self):
        with kegg_get(["hsa:10458", "ece:Z5100"]) as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/get/hsa:10458+ece:Z5100")

    def test_get_hsa_10458_plus_ece_Z5100_as_aaseq(self):
        with kegg_get("hsa:10458+ece:Z5100", "aaseq") as handle:
            data = SeqIO.parse(handle, "fasta")
            self.assertEqual(len(list(data)), 2)
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/get/hsa:10458+ece:Z5100/aaseq"
        )

    def test_get_hsa_10458_list_ece_Z5100_as_aaseq(self):
        with kegg_get(["hsa:10458", "ece:Z5100"], "aaseq") as handle:
            data = SeqIO.parse(handle, "fasta")
            self.assertEqual(len(list(data)), 2)
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/get/hsa:10458+ece:Z5100/aaseq"
        )

    def test_get_hsa_10458_plus_ece_Z5100_as_ntseq(self):
        with kegg_get("hsa:10458+ece:Z5100", "ntseq") as handle:
            data = SeqIO.parse(handle, "fasta")
            self.assertEqual(len(list(data)), 2)
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/get/hsa:10458+ece:Z5100/ntseq"
        )

    def test_get_hsa_10458_list_ece_Z5100_as_ntseq(self):
        with kegg_get(["hsa:10458", "ece:Z5100"], "ntseq") as handle:
            data = SeqIO.parse(handle, "fasta")
            self.assertEqual(len(list(data)), 2)
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/get/hsa:10458+ece:Z5100/ntseq"
        )

    def test_get_hsa05130_image(self):
        with kegg_get("hsa05130", "image") as handle:
            data = handle.read()
        self.assertEqual(data[:4], b"\x89PNG")
        self.assertEqual(handle.url, "https://rest.kegg.jp/get/hsa05130/image")

    def test_conv_eco_ncbi_geneid(self):
        with kegg_conv("eco", "ncbi-geneid") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/conv/eco/ncbi-geneid")

    def test_conv_ncbi_geneid_eco(self):
        with kegg_conv("ncbi-geneid", "eco") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/conv/ncbi-geneid/eco")

    def test_conv_ncbi_gi_hsa_10458_plus_ece_Z5100(self):
        with kegg_conv("ncbi-gi", "hsa:10458+ece:Z5100") as handle:
            handle.read()
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/conv/ncbi-gi/hsa:10458+ece:Z5100"
        )

    def test_conv_ncbi_gi_hsa_10458_list_ece_Z5100(self):
        with kegg_conv("ncbi-gi", ["hsa:10458", "ece:Z5100"]) as handle:
            handle.read()
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/conv/ncbi-gi/hsa:10458+ece:Z5100"
        )

    def test_link_pathway_hsa(self):
        with kegg_link("pathway", "hsa") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/link/pathway/hsa")

    def test_link_hsa_pathway(self):
        with kegg_link("hsa", "pathway") as handle:
            handle.read()
        self.assertEqual(handle.url, "https://rest.kegg.jp/link/hsa/pathway")

    def test_pathway_hsa_10458_plus_ece_Z5100(self):
        with kegg_link("pathway", "hsa:10458+ece:Z5100") as handle:
            handle.read()
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/link/pathway/hsa:10458+ece:Z5100"
        )

    def test_pathway_hsa_10458_list_ece_Z5100(self):
        with kegg_link("pathway", ["hsa:10458", "ece:Z5100"]) as handle:
            handle.read()
        self.assertEqual(
            handle.url, "https://rest.kegg.jp/link/pathway/hsa:10458+ece:Z5100"
        )


class KGMLPathwayTests(unittest.TestCase):
    """Tests with metabolic maps."""

    def test_parse_remote_pathway(self):
        """Download a KEGG pathway from the KEGG server and parse KGML."""
        with kegg_get("ko03070", "kgml") as handle:
            pathway = KGML_parser.read(handle)
        self.assertEqual(pathway.name, "path:ko03070")

    def test_parser_roundtrip(self):
        """Download a KEGG pathway, write local KGML and check roundtrip."""
        with kegg_get("ko00680", "kgml") as remote_handle:
            pathway = KGML_parser.read(remote_handle)

        with io.StringIO(pathway.get_KGML()) as local_handle:
            roundtrip = KGML_parser.read(local_handle)

        self.assertEqual(pathway.name, roundtrip.name)
        self.assertEqual(len(pathway.relations), len(roundtrip.relations))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
