# Copyright 2008 by Peter Cock.
# Revision copyright 2017 by Markus Piotrowski.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for CodonTable module."""

import unittest

from Bio.Data import IUPACData
from Bio.Data.CodonTable import (
    generic_by_id,
    generic_by_name,
    ambiguous_generic_by_id,
    ambiguous_generic_by_name,
    ambiguous_rna_by_id,
    ambiguous_dna_by_id,
    ambiguous_dna_by_name,
    ambiguous_rna_by_name,
    unambiguous_rna_by_id,
    unambiguous_rna_by_name,
    unambiguous_dna_by_id,
    unambiguous_dna_by_name,
)
from Bio.Data.CodonTable import list_ambiguous_codons, list_possible_proteins
from Bio.Data.CodonTable import TranslationError

exception_list = []
ids = list(unambiguous_dna_by_id)

# Find codon tables with 'ambiguous' stop codons (coding both for stop and
# an amino acid) and put them into exception_list
for id in ids:
    table = unambiguous_dna_by_id[id]
    for codon in table.stop_codons:
        if codon not in table.forward_table:
            continue
        else:
            exception_list.append(table.id)
            break


class BasicSanityTests(unittest.TestCase):
    """Basic tests."""

    def test_number_of_tables(self):
        """Check if we have the same number of tables for each type."""
        self.assertTrue(
            len(unambiguous_dna_by_id)
            == len(unambiguous_rna_by_id)
            == len(generic_by_id)
            == len(ambiguous_dna_by_id)
            == len(ambiguous_rna_by_id)
            == len(ambiguous_generic_by_id)
        )
        self.assertTrue(
            len(unambiguous_dna_by_name)
            == len(unambiguous_rna_by_name)
            == len(generic_by_name)
            == len(ambiguous_dna_by_name)
            == len(ambiguous_rna_by_name)
            == len(ambiguous_generic_by_name)
        )

    def test_complete_tables(self):
        """Check if all unambiguous codon tables have all entries.

        For DNA/RNA tables each of the 64 codons must be either in
        the forward_table or in the stop_codons list.
        For DNA+RNA (generic) tables the number of possible codons is 101.

        There must be at least one start codon and one stop codon.

        The back_tables must have 21 entries (20 amino acids and one
        stop symbol).
        """
        for id in ids:
            nuc_table = generic_by_id[id]
            dna_table = unambiguous_dna_by_id[id]
            rna_table = unambiguous_rna_by_id[id]

            if id not in exception_list:
                self.assertEqual(
                    len(dna_table.forward_table) + len(dna_table.stop_codons), 64
                )
                self.assertEqual(
                    len(rna_table.forward_table) + len(rna_table.stop_codons), 64
                )
                self.assertEqual(
                    len(nuc_table.forward_table) + len(nuc_table.stop_codons), 101
                )
                self.assertTrue(dna_table.stop_codons)
                self.assertTrue(dna_table.start_codons)
                self.assertTrue(rna_table.stop_codons)
                self.assertTrue(rna_table.start_codons)
                self.assertTrue(nuc_table.start_codons)
                self.assertTrue(nuc_table.stop_codons)
                self.assertEqual(len(dna_table.back_table), 21)
                self.assertEqual(len(rna_table.back_table), 21)
                self.assertEqual(len(nuc_table.back_table), 21)

    def test_ambiguous_tables(self):
        """Check if all IDs and all names are present in ambiguous tables."""
        for key, val in generic_by_name.items():
            self.assertIn(key, ambiguous_generic_by_name[key].names)
        for key, val in generic_by_id.items():
            self.assertEqual(ambiguous_generic_by_id[key].id, key)


class AmbiguousCodonsTests(unittest.TestCase):
    """Tests for ambiguous codons."""

    def test_list_ambiguous_codons(self):
        """Check if stop codons are properly extended."""
        self.assertEqual(
            list_ambiguous_codons(["TGA", "TAA"], IUPACData.ambiguous_dna_values),
            ["TGA", "TAA", "TRA"],
        )
        self.assertEqual(
            list_ambiguous_codons(["TAG", "TGA"], IUPACData.ambiguous_dna_values),
            ["TAG", "TGA"],
        )
        self.assertEqual(
            list_ambiguous_codons(["TAG", "TAA"], IUPACData.ambiguous_dna_values),
            ["TAG", "TAA", "TAR"],
        )
        self.assertEqual(
            list_ambiguous_codons(["UAG", "UAA"], IUPACData.ambiguous_rna_values),
            ["UAG", "UAA", "UAR"],
        )
        self.assertEqual(
            list_ambiguous_codons(
                ["TGA", "TAA", "TAG"], IUPACData.ambiguous_dna_values
            ),
            ["TGA", "TAA", "TAG", "TAR", "TRA"],
        )

    def test_coding(self):
        """Check a few ambiguous codons for correct coding."""
        for id in ids:
            amb_dna = ambiguous_dna_by_id[id]
            amb_rna = ambiguous_rna_by_id[id]
            amb_nuc = ambiguous_generic_by_id[id]

            self.assertEqual(amb_rna.forward_table["GUU"], "V")
            self.assertEqual(amb_rna.forward_table["GUN"], "V")
            self.assertEqual(amb_dna.forward_table["GTT"], "V")
            self.assertEqual(amb_dna.forward_table["GTN"], "V")
            self.assertEqual(amb_rna.forward_table["ACN"], "T")
            self.assertEqual(amb_nuc.forward_table["GUU"], "V")
            self.assertEqual(amb_nuc.forward_table["GUN"], "V")
            self.assertEqual(amb_nuc.forward_table["GTT"], "V")
            self.assertEqual(amb_nuc.forward_table["GTN"], "V")
            # And finally something evil, an RNA-DNA mixture:
            self.assertEqual(amb_nuc.forward_table["UTU"], "F")
            if id != 23:
                self.assertEqual(amb_rna.forward_table["UUN"], "X")  # F or L
                self.assertEqual(amb_dna.forward_table["TTN"], "X")  # F or L
                self.assertEqual(amb_nuc.forward_table.get("TTN"), "X")
                self.assertEqual(amb_nuc.forward_table["UUN"], "X")  # F or L
                self.assertEqual(amb_nuc.forward_table["TTN"], "X")  # F or L
                self.assertEqual(amb_nuc.forward_table["UTN"], "X")  # F or L

    def test_stop_codons(self):
        """Test various ambiguous codons as stop codon.

        Stop codons should not appear in forward tables. This should give a
        KeyError. If an ambiguous codon may code for both (stop codon and
        amino acid), this should raise a TranslationError.
        """
        for id in ids:
            rna = unambiguous_rna_by_id[id]
            amb_dna = ambiguous_dna_by_id[id]
            amb_rna = ambiguous_rna_by_id[id]
            amb_nuc = ambiguous_generic_by_id[id]

            # R = A or G, so URR = UAA or UGA / TRA = TAA or TGA = stop codons
            if (
                "UAA" in amb_rna.stop_codons
                and "UGA" in amb_rna.stop_codons
                and id not in (28, 32)
            ):
                self.assertEqual(amb_dna.forward_table.get("TRA", "X"), "X")
                with self.assertRaises(KeyError):
                    amb_dna.forward_table["TRA"]
                    amb_rna.forward_table["URA"]
                    amb_nuc.forward_table["URA"]
                self.assertIn("URA", amb_nuc.stop_codons)
                self.assertIn("URA", amb_rna.stop_codons)
                self.assertIn("TRA", amb_nuc.stop_codons)
                self.assertIn("TRA", amb_dna.stop_codons)

            if (
                "UAG" in rna.stop_codons
                and "UAA" in rna.stop_codons
                and "UGA" in rna.stop_codons
                and id not in (28, 32)
            ):
                with self.assertRaises(KeyError):
                    amb_dna.forward_table["TAR"]
                    amb_rna.forward_table["UAR"]
                    amb_nuc.forward_table["UAR"]
                with self.assertRaises(TranslationError):
                    amb_nuc.forward_table["URR"]
                self.assertIn("UAR", amb_nuc.stop_codons)
                self.assertIn("UAR", amb_rna.stop_codons)
                self.assertIn("TAR", amb_nuc.stop_codons)
                self.assertIn("TAR", amb_dna.stop_codons)
                self.assertIn("URA", amb_nuc.stop_codons)
                self.assertIn("URA", amb_rna.stop_codons)
                self.assertIn("TRA", amb_nuc.stop_codons)
                self.assertIn("TRA", amb_dna.stop_codons)

    def test_start_codons(self):
        """Test various ambiguous codons as start codon."""
        for id in ids:
            rna = unambiguous_rna_by_id[id]
            amb_dna = ambiguous_dna_by_id[id]
            amb_rna = ambiguous_rna_by_id[id]

            if (
                "UUG" in rna.start_codons
                and "CUG" in rna.start_codons
                and "AUG" in rna.start_codons
                and "UUG" not in rna.start_codons
            ):
                self.assertNotIn("NUG", amb_rna.start_codons)
                self.assertNotIn("RUG", amb_rna.start_codons)
                self.assertNotIn("WUG", amb_rna.start_codons)
                self.assertNotIn("KUG", amb_rna.start_codons)
                self.assertNotIn("SUG", amb_rna.start_codons)
                self.assertNotIn("DUG", amb_rna.start_codons)

                self.assertNotIn("NTG", amb_dna.start_codons)
                self.assertNotIn("RTG", amb_dna.start_codons)
                self.assertNotIn("WTG", amb_dna.start_codons)
                self.assertNotIn("KTG", amb_dna.start_codons)
                self.assertNotIn("STG", amb_dna.start_codons)
                self.assertNotIn("DTG", amb_dna.start_codons)


class SingleTableTests(unittest.TestCase):
    """Several test for individual codon tables.

    In general we want to check:
        - that we can call the codon table by its name(s) or id
        - the number of start and stop codons
        - deviations from the standard code
        - that stop codons do not appear in the forward_table (except in
          tables 27, 28, 31).
    """

    def test_table01(self):
        """Check table 1: Standard."""
        self.assertEqual(ambiguous_dna_by_id[1].names, ["Standard", "SGC0"])
        self.assertEqual(
            ambiguous_dna_by_name["Standard"].stop_codons,
            ambiguous_dna_by_id[1].stop_codons,
        )
        self.assertEqual(
            generic_by_id[1].start_codons, ["TTG", "UUG", "CTG", "CUG", "ATG", "AUG"]
        )
        self.assertEqual(len(unambiguous_dna_by_id[1].start_codons), 3)
        self.assertEqual(len(unambiguous_dna_by_id[1].stop_codons), 3)

    def test_table02(self):
        """Check table 2: Vertebrate Mitochondrial.

        Table 2 Vertebrate Mitochondrial has TAA and TAG -> TAR,
        plus AGA and AGG -> AGR as stop codons.
        """
        self.assertEqual(generic_by_name["Vertebrate Mitochondrial"].id, 2)
        self.assertEqual(generic_by_name["SGC1"].id, 2)
        self.assertIn("SGC1", generic_by_id[2].names)
        self.assertIn("AGR", ambiguous_dna_by_id[2].stop_codons)
        self.assertIn("TAR", ambiguous_dna_by_id[2].stop_codons)
        self.assertIn("AGR", ambiguous_rna_by_id[2].stop_codons)
        self.assertIn("UAR", ambiguous_rna_by_id[2].stop_codons)
        self.assertIn("AGR", ambiguous_generic_by_id[2].stop_codons)
        self.assertIn("UAR", ambiguous_generic_by_id[2].stop_codons)
        self.assertIn("TAR", ambiguous_generic_by_id[2].stop_codons)
        self.assertEqual(len(unambiguous_dna_by_id[2].start_codons), 5)
        self.assertEqual(len(unambiguous_dna_by_id[2].stop_codons), 4)

    def test_table03(self):
        """Check table 3: Yeast Mitochondrial.

        Start codons ATG and ATA -> ATR and start codon GTG (since ver. 4.4)
        Stop codons TAA and TAG -> TAR
        TGA codes for W (instead of stop) and CTN codes for T (instead of L).
        """
        self.assertEqual(generic_by_name["Yeast Mitochondrial"].id, 3)
        self.assertEqual(generic_by_name["SGC2"].id, 3)
        self.assertIn("SGC2", generic_by_id[3].names)
        self.assertEqual(len(unambiguous_dna_by_id[3].start_codons), 3)
        self.assertEqual(len(unambiguous_dna_by_id[3].stop_codons), 2)
        self.assertIn("ATR", ambiguous_dna_by_id[3].start_codons)
        self.assertIn("GTG", unambiguous_dna_by_id[3].start_codons)
        self.assertIn("TAR", ambiguous_dna_by_id[3].stop_codons)
        self.assertNotIn("TGA", ambiguous_dna_by_id[3].stop_codons)
        self.assertEqual(generic_by_id[3].forward_table["UGA"], "W")
        self.assertEqual(ambiguous_rna_by_id[3].forward_table["CUN"], "T")

    def test_table04(self):
        """Check table 4: Mold Mitochondrial and others.

        Stop codons TAA and TAG -> TAR
        TGA codes for W (instead of stop).
        """
        self.assertEqual(generic_by_name["Mold Mitochondrial"].id, 4)
        self.assertEqual(generic_by_name["Mycoplasma"].id, 4)
        self.assertIn("SGC3", generic_by_id[4].names)
        self.assertEqual(len(unambiguous_dna_by_id[4].start_codons), 8)
        self.assertEqual(len(unambiguous_dna_by_id[4].stop_codons), 2)
        self.assertIn("ATN", ambiguous_dna_by_id[4].start_codons)
        self.assertIn("TAR", ambiguous_dna_by_id[4].stop_codons)
        self.assertNotIn("TGA", ambiguous_dna_by_id[4].stop_codons)
        self.assertEqual(ambiguous_rna_by_id[4].forward_table["UGA"], "W")

    def test_table05(self):
        """Check table 5: Invertebrate Mitochondrial.

        Stop codons TAA and TAG -> TAR
        TGA codes for W (instead of stop), AGR codes for S (instead of R),
        ATA for M (instead of I).
        """
        self.assertEqual(generic_by_name["Invertebrate Mitochondrial"].id, 5)
        self.assertEqual(generic_by_name["SGC4"].id, 5)
        self.assertIn("SGC4", generic_by_id[5].names)
        self.assertEqual(len(unambiguous_dna_by_id[5].start_codons), 6)
        self.assertEqual(len(unambiguous_dna_by_id[5].stop_codons), 2)
        self.assertIn("ATN", ambiguous_dna_by_id[5].start_codons)
        self.assertIn("KTG", ambiguous_dna_by_id[5].start_codons)
        self.assertIn("TAR", ambiguous_dna_by_id[5].stop_codons)
        self.assertNotIn("TGA", ambiguous_dna_by_id[5].stop_codons)
        self.assertEqual(ambiguous_rna_by_id[5].forward_table["UGA"], "W")
        self.assertEqual(ambiguous_dna_by_id[5].forward_table["AGR"], "S")
        self.assertEqual(generic_by_id[5].forward_table["AUA"], "M")

    def test_table06(self):
        """Check table 6: Ciliate and Other Nuclear.

        Only one stop codon TGA. TAA and TAG code for Q.
        """
        dna_table = unambiguous_dna_by_id[6]
        nuc_table = generic_by_id[6]
        amb_rna_table = ambiguous_rna_by_id[6]
        amb_nuc_table = ambiguous_generic_by_id[6]

        self.assertEqual(generic_by_name["Ciliate Nuclear"].id, 6)
        self.assertEqual(generic_by_name["Hexamita Nuclear"].id, 6)
        self.assertIn("SGC5", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 1)
        self.assertNotIn("UAR", amb_rna_table.stop_codons)
        self.assertEqual(amb_nuc_table.forward_table["UAR"], "Q")

    def test_table09(self):
        """Check table 9: Echinoderm and Flatworm Mitochondrial.

        Stop codons TAA and TAG -> TAR
        TGA codes for W (instead of stop), AGR codes for S (instead of R),
        AAA for N (instead of K).
        """
        dna_table = unambiguous_dna_by_id[9]
        rna_table = unambiguous_rna_by_id[9]
        nuc_table = generic_by_id[9]
        amb_rna_table = ambiguous_rna_by_id[9]
        amb_nuc_table = ambiguous_generic_by_id[9]

        self.assertEqual(generic_by_name["Echinoderm Mitochondrial"].id, 9)
        self.assertEqual(generic_by_name["Flatworm Mitochondrial"].id, 9)
        self.assertIn("SGC8", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 2)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("UAR", amb_rna_table.stop_codons)
        self.assertNotIn("TGA", dna_table.stop_codons)
        self.assertEqual(nuc_table.forward_table["UGA"], "W")
        self.assertEqual(amb_nuc_table.forward_table["AGR"], "S")
        self.assertEqual(rna_table.forward_table["AAA"], "N")

    def test_table10(self):
        """Check table 10: Euplotid Nuclear.

        Stop codons TAA and TAG -> TAR
        TGA codes for C (instead of stop).
        """
        dna_table = unambiguous_dna_by_id[10]
        nuc_table = generic_by_id[10]
        amb_rna_table = ambiguous_rna_by_id[10]

        self.assertEqual(generic_by_name["Euplotid Nuclear"].id, 10)
        self.assertIn("SGC9", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("UAR", amb_rna_table.stop_codons)
        self.assertNotIn("TGA", dna_table.stop_codons)
        self.assertEqual(nuc_table.forward_table["UGA"], "C")

    def test_table11(self):
        """Check table 11: Bacterial, Archaeal and Plant Plastid."""
        dna_table = unambiguous_dna_by_id[11]
        nuc_table = generic_by_id[11]

        self.assertEqual(generic_by_name["Bacterial"].id, 11)
        self.assertEqual(generic_by_name["Archaeal"].id, 11)
        self.assertIn("Plant Plastid", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 7)
        self.assertEqual(len(dna_table.stop_codons), 3)

    def test_table12(self):
        """Check table 12: Alternative Yeast Nuclear.

        CTG codes for S (instead of L).
        """
        dna_table = unambiguous_dna_by_id[12]
        nuc_table = generic_by_id[12]

        self.assertEqual(generic_by_name["Alternative Yeast Nuclear"].id, 12)
        self.assertIn("Alternative Yeast Nuclear", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 2)
        self.assertEqual(len(dna_table.stop_codons), 3)
        self.assertEqual(nuc_table.forward_table["CUG"], "S")

    def test_table13(self):
        """Check table 13: Ascidian Mitochondrial.

        Stop codons TAA and TAG -> TAR
        TGA codes for W (instead of stop), AGR codes for G (instead of R),
        ATA for M (instead of I).
        """
        dna_table = unambiguous_dna_by_id[13]
        nuc_table = generic_by_id[13]
        amb_rna_table = ambiguous_rna_by_id[13]
        amb_nuc_table = ambiguous_generic_by_id[13]

        self.assertEqual(generic_by_name["Ascidian Mitochondrial"].id, 13)
        self.assertIn("Ascidian Mitochondrial", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 4)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("UAR", amb_rna_table.stop_codons)
        self.assertNotIn("TGA", dna_table.stop_codons)
        self.assertEqual(nuc_table.forward_table["UGA"], "W")
        self.assertEqual(amb_nuc_table.forward_table["AGR"], "G")
        self.assertEqual(amb_rna_table.forward_table["AUR"], "M")

    def test_table14(self):
        """Check table 14: Alternative Flatworm Mitochondrial.

        Only one stop codon TAG. TAA and TGA code for Y and W, respecitively
        (instead of stop). AGR codes for S (instead of R), AAA for N (instead
        of K).
        """
        dna_table = unambiguous_dna_by_id[14]
        rna_table = unambiguous_rna_by_id[14]
        nuc_table = generic_by_id[14]
        amb_rna_table = ambiguous_rna_by_id[14]
        amb_nuc_table = ambiguous_generic_by_id[14]

        self.assertEqual(generic_by_name["Alternative Flatworm Mitochondrial"].id, 14)
        self.assertIn("Alternative Flatworm Mitochondrial", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 1)
        self.assertNotIn("URA", amb_rna_table.stop_codons)
        self.assertEqual(nuc_table.forward_table["UAA"], "Y")
        self.assertEqual(rna_table.forward_table["UGA"], "W")
        self.assertEqual(amb_nuc_table.forward_table["AGR"], "S")
        self.assertEqual(rna_table.forward_table["AAA"], "N")

    def test_table16(self):
        """Check table 16: Chlorophycean Mitochondrial.

        Stop codons TAA and TGA -> TRA
        TAG codes for L (instead of stop).
        """
        dna_table = unambiguous_dna_by_id[16]
        nuc_table = generic_by_id[16]
        amb_rna_table = ambiguous_rna_by_id[16]

        self.assertEqual(generic_by_name["Chlorophycean Mitochondrial"].id, 16)
        self.assertIn("Chlorophycean Mitochondrial", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("URA", amb_rna_table.stop_codons)
        self.assertNotIn("TAG", dna_table.stop_codons)
        self.assertEqual(nuc_table.forward_table["UAG"], "L")

    def test_table21(self):
        """Check table 21: Trematode Mitochondrial.

        Stop codons TAA and TAG -> TAR
        TGA codes for W (instead of stop), ATA codes for M (instead of I),
        AGR codes for S (instead of R), AAA for N (instead of K).
        """
        dna_table = unambiguous_dna_by_id[21]
        rna_table = unambiguous_rna_by_id[21]
        nuc_table = generic_by_id[21]
        amb_rna_table = ambiguous_rna_by_id[21]
        amb_nuc_table = ambiguous_generic_by_id[21]

        self.assertEqual(generic_by_name["Trematode Mitochondrial"].id, 21)
        self.assertIn("Trematode Mitochondrial", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 2)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("UAR", amb_rna_table.stop_codons)
        self.assertNotIn("TGA", dna_table.stop_codons)
        self.assertEqual(rna_table.forward_table["AUA"], "M")
        self.assertEqual(nuc_table.forward_table["UGA"], "W")
        self.assertEqual(amb_nuc_table.forward_table["AGR"], "S")
        self.assertEqual(rna_table.forward_table["AAA"], "N")

    def test_table22(self):
        """Check table 22: Scenedesmus obliquus Mitochondrial.

        Stop codons TAA, TCA and TGA -> TVA
        TAG codes for L (instead of stop).
        """
        dna_table = unambiguous_dna_by_id[22]
        nuc_table = generic_by_id[22]
        amb_rna_table = ambiguous_rna_by_id[22]

        self.assertEqual(generic_by_name["Scenedesmus obliquus Mitochondrial"].id, 22)
        self.assertIn("Scenedesmus obliquus Mitochondrial", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 3)
        self.assertIn("UVA", amb_rna_table.stop_codons)
        self.assertNotIn("TAG", dna_table.stop_codons)
        self.assertEqual(nuc_table.forward_table["UAG"], "L")

    def test_table23(self):
        """Check table 9: Thraustochytrium Mitochondrial.

        TTA codes for stop (instead of L).
        """
        dna_table = unambiguous_dna_by_id[23]
        nuc_table = generic_by_id[23]
        amb_rna_table = ambiguous_rna_by_id[23]

        self.assertEqual(generic_by_name["Thraustochytrium Mitochondrial"].id, 23)
        self.assertIn("Thraustochytrium Mitochondrial", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 3)
        self.assertEqual(len(dna_table.stop_codons), 4)
        self.assertIn("UUA", amb_rna_table.stop_codons)
        self.assertNotIn("TTA", dna_table.forward_table.keys())

    def test_table24(self):
        """Check table 24: Pterobranchia Mitochondrial.

        Stop codons TAA and TAG -> TAR
        TGA codes for W (instead of stop), AGA codes for S (instead of R),
        AGG for K (instead of R).
        """
        dna_table = unambiguous_dna_by_id[24]
        nuc_table = generic_by_id[24]
        amb_rna_table = ambiguous_rna_by_id[24]
        amb_nuc_table = ambiguous_generic_by_id[24]

        self.assertEqual(generic_by_name["Pterobranchia Mitochondrial"].id, 24)
        self.assertIn("Pterobranchia Mitochondrial", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 4)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("UAR", amb_rna_table.stop_codons)
        self.assertNotIn("TGA", dna_table.stop_codons)
        self.assertEqual(nuc_table.forward_table["UGA"], "W")
        self.assertEqual(amb_nuc_table.forward_table["AGA"], "S")
        self.assertEqual(dna_table.forward_table["AGG"], "K")

    def test_table25(self):
        """Check table 25: Candidate Division SR1 and Gracilibacteria.

        Stop codons TAA and TAG -> TAR
        TGA codes for G (instead of stop).
        """
        dna_table = unambiguous_dna_by_id[25]
        nuc_table = generic_by_id[25]
        amb_rna_table = ambiguous_rna_by_id[25]

        self.assertEqual(generic_by_name["Candidate Division SR1"].id, 25)
        self.assertEqual(generic_by_name["Gracilibacteria"].id, 25)
        self.assertIn("Candidate Division SR1", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 3)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("UAR", amb_rna_table.stop_codons)
        self.assertNotIn("TGA", dna_table.stop_codons)
        self.assertEqual(nuc_table.forward_table["UGA"], "G")

    def test_table26(self):
        """Check table 26: Pachysolen tannophilus Nuclear.

        CTG codes for A (instead of L).
        """
        dna_table = unambiguous_dna_by_id[26]
        nuc_table = generic_by_id[26]

        self.assertEqual(generic_by_name["Pachysolen tannophilus Nuclear"].id, 26)
        self.assertIn("Pachysolen tannophilus Nuclear", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 2)
        self.assertEqual(len(dna_table.stop_codons), 3)
        self.assertEqual(nuc_table.forward_table["CTG"], "A")

    def test_table27(self):
        """Check table 27: Karyorelict Nuclear.

        NOTE: This code has no unambiguous stop codon! TGA codes for either
        stop or W, dependent on the context.
        TAR codes for Q (instead of stop)
        """
        dna_table = unambiguous_dna_by_id[27]
        nuc_table = generic_by_id[27]
        amb_rna_table = ambiguous_rna_by_id[27]
        amb_nuc_table = ambiguous_generic_by_id[27]

        self.assertEqual(generic_by_name["Karyorelict Nuclear"].id, 27)
        self.assertIn("Karyorelict Nuclear", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 1)
        self.assertNotIn("UAR", amb_rna_table.stop_codons)
        self.assertEqual(amb_nuc_table.forward_table["UAR"], "Q")
        self.assertEqual(nuc_table.forward_table["UGA"], "W")
        self.assertIn("TGA", dna_table.stop_codons)

    def test_table28(self):
        """Check table 28: Condylostoma Nuclear.

        NOTE: This code has no unambiguous stop codon! TAR codes for either
        stop or Q, TGA for stop or W, dependent on the context.
        """
        dna_table = unambiguous_dna_by_id[28]
        nuc_table = generic_by_id[28]
        amb_dna_table = ambiguous_dna_by_id[28]
        amb_rna_table = ambiguous_rna_by_id[28]

        self.assertEqual(generic_by_name["Condylostoma Nuclear"].id, 28)
        self.assertIn("Condylostoma Nuclear", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 3)
        self.assertIn("UAR", amb_rna_table.stop_codons)
        self.assertIn("TGA", amb_dna_table.stop_codons)
        self.assertEqual(amb_dna_table.forward_table["TAR"], "Q")
        self.assertEqual(nuc_table.forward_table["UGA"], "W")

    def test_table29(self):
        """Check table 29: Mesodinium Nuclear.

        TAR codes for Y (instead of stop).
        """
        dna_table = unambiguous_dna_by_id[29]
        nuc_table = generic_by_id[29]
        amb_rna_table = ambiguous_rna_by_id[29]
        amb_nuc_table = ambiguous_generic_by_id[29]

        self.assertEqual(generic_by_name["Mesodinium Nuclear"].id, 29)
        self.assertIn("Mesodinium Nuclear", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 1)
        self.assertNotIn("UAR", amb_rna_table.stop_codons)
        self.assertEqual(amb_nuc_table.forward_table["UAR"], "Y")

    def test_table30(self):
        """Check table 30: Peritrich Nuclear.

        TAR codes for E (instead of stop).
        """
        dna_table = unambiguous_dna_by_id[30]
        nuc_table = generic_by_id[30]
        amb_rna_table = ambiguous_rna_by_id[30]
        amb_nuc_table = ambiguous_generic_by_id[30]

        self.assertEqual(generic_by_name["Peritrich Nuclear"].id, 30)
        self.assertIn("Peritrich Nuclear", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 1)
        self.assertNotIn("UAR", amb_rna_table.stop_codons)
        self.assertEqual(amb_nuc_table.forward_table["UAR"], "E")

    def test_table31(self):
        """Check table 31: Blastocrithidia Nuclear.

        NOTE: This code has no unambiguous stop codon! TAR codes for either
        stop or E, dependent on the context.
        TGA codes for W (instead of stop)
        """
        dna_table = unambiguous_dna_by_id[31]
        nuc_table = generic_by_id[31]
        amb_dna_table = ambiguous_dna_by_id[31]
        amb_rna_table = ambiguous_rna_by_id[31]

        self.assertEqual(generic_by_name["Blastocrithidia Nuclear"].id, 31)
        self.assertIn("Blastocrithidia Nuclear", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 1)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("UAR", amb_rna_table.stop_codons)
        self.assertNotIn("TGA", amb_dna_table.stop_codons)
        self.assertEqual(amb_dna_table.forward_table["TAR"], "E")
        self.assertEqual(nuc_table.forward_table["UGA"], "W")

    def test_table32(self):
        """Check table 32: Balanophoraceae Plastid.

        In v4.4, TAG and TGA were simple stop codons, while TAA
        coded for either stop or W, dependent on the context.

        In v4.5, TAA and TGA were just simple stop codons. TAG
        was not a stop codon any more. TAA was just a stop codon.
        """
        dna_table = unambiguous_dna_by_id[32]
        nuc_table = generic_by_id[32]
        amb_dna_table = ambiguous_dna_by_id[32]
        amb_rna_table = ambiguous_rna_by_id[32]

        self.assertEqual(generic_by_name["Balanophoraceae Plastid"].id, 32)
        self.assertIn("Balanophoraceae Plastid", nuc_table.names)
        self.assertEqual(len(dna_table.start_codons), 7)
        for codon in ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"):
            self.assertIn(codon, dna_table.start_codons)
        self.assertEqual(len(dna_table.stop_codons), 2)
        self.assertIn("URA", amb_rna_table.stop_codons)
        self.assertNotIn("UAA", nuc_table.forward_table)


class ErrorConditions(unittest.TestCase):
    """Tests for module specific errors."""

    def test_list_possible_proteins(self):
        """Raise errors in list_possible proteins."""
        table = unambiguous_dna_by_id[1]
        amb_values = {"T": "T", "G": "G", "A": "A", "R": ("A", "G")}
        with self.assertRaises(TranslationError):
            # Can be stop or amino acid:
            codon = ["T", "R", "R"]
            list_possible_proteins(codon, table.forward_table, amb_values)
        with self.assertRaises(KeyError):
            # Is a stop codon:
            codon = ["T", "G", "A"]
            list_possible_proteins(codon, table.forward_table, amb_values)

    def test_ambiguous_forward_table(self):
        """Raise errors in AmbiguousForwardTable."""
        table = ambiguous_dna_by_id[1]
        self.assertIsNone(table.forward_table.get("ZZZ"))
        with self.assertRaises(KeyError):
            table.forward_table["ZZZ"]  # KeyError it's a stop codon
            table.forward_table["TGA"]  # KeyError stop codon
        with self.assertRaises(TranslationError):
            table.forward_table["WWW"]  # Translation error does not code


class PrintTable(unittest.TestCase):
    """Test for __str__ in CodonTable."""

    def test_print_table(self):
        """Test output of __str__ function."""
        table = generic_by_id[1]
        output = table.__str__()

        expected_output = """Table 1 Standard, SGC0

  |  U      |  C      |  A      |  G      |
--+---------+---------+---------+---------+--
U | UUU F   | UCU S   | UAU Y   | UGU C   | U
U | UUC F   | UCC S   | UAC Y   | UGC C   | C
U | UUA L   | UCA S   | UAA Stop| UGA Stop| A
U | UUG L(s)| UCG S   | UAG Stop| UGG W   | G
--+---------+---------+---------+---------+--
C | CUU L   | CCU P   | CAU H   | CGU R   | U
C | CUC L   | CCC P   | CAC H   | CGC R   | C
C | CUA L   | CCA P   | CAA Q   | CGA R   | A
C | CUG L(s)| CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | AUU I   | ACU T   | AAU N   | AGU S   | U
A | AUC I   | ACC T   | AAC N   | AGC S   | C
A | AUA I   | ACA T   | AAA K   | AGA R   | A
A | AUG M(s)| ACG T   | AAG K   | AGG R   | G
--+---------+---------+---------+---------+--
G | GUU V   | GCU A   | GAU D   | GGU G   | U
G | GUC V   | GCC A   | GAC D   | GGC G   | C
G | GUA V   | GCA A   | GAA E   | GGA G   | A
G | GUG V   | GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--"""

        self.assertEqual(output, expected_output)

        table = unambiguous_dna_by_id[1]
        table.id = ""
        output = table.__str__()

        expected_output = """Table ID unknown Standard, SGC0

  |  T      |  C      |  A      |  G      |
--+---------+---------+---------+---------+--
T | TTT F   | TCT S   | TAT Y   | TGT C   | T
T | TTC F   | TCC S   | TAC Y   | TGC C   | C
T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
--+---------+---------+---------+---------+--
C | CTT L   | CCT P   | CAT H   | CGT R   | T
C | CTC L   | CCC P   | CAC H   | CGC R   | C
C | CTA L   | CCA P   | CAA Q   | CGA R   | A
C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | ATT I   | ACT T   | AAT N   | AGT S   | T
A | ATC I   | ACC T   | AAC N   | AGC S   | C
A | ATA I   | ACA T   | AAA K   | AGA R   | A
A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
--+---------+---------+---------+---------+--
G | GTT V   | GCT A   | GAT D   | GGT G   | T
G | GTC V   | GCC A   | GAC D   | GGC G   | C
G | GTA V   | GCA A   | GAA E   | GGA G   | A
G | GTG V   | GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--"""

        self.assertEqual(output, expected_output)
        # We need to set table.id to the correct value, otherwise
        # following tests may fail!
        table.id = 1


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
