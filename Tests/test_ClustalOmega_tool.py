# Copyright 2008-2011 by Peter Cock.  All rights reserved.
# Revisions copyright 2012 by Christian Brueffer.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for ClustalOmega tool."""

import os
import unittest
import warnings

from subprocess import getoutput

from Bio import MissingExternalDependencyError
from Bio import SeqIO
from Bio import Align

from Bio import BiopythonDeprecationWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=BiopythonDeprecationWarning)
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio.Application import ApplicationError

#################################################################

# Try to avoid problems when the OS is in another language
os.environ["LANG"] = "C"

clustalo_exe = None
try:
    output = getoutput("clustalo --help")
    if output.startswith("Clustal Omega"):
        clustalo_exe = "clustalo"
except FileNotFoundError:
    pass

if not clustalo_exe:
    raise MissingExternalDependencyError(
        "Install clustalo if you want to use Clustal Omega from Biopython."
    )


class ClustalOmegaTestCase(unittest.TestCase):
    def setUp(self):
        self.files_to_clean = set()

    def tearDown(self):
        for filename in self.files_to_clean:
            if os.path.isfile(filename):
                os.remove(filename)

    def standard_test_procedure(self, cline):
        """Shared test procedure used by all tests."""
        # Overwrite existing files.
        cline.force = True

        # Mark output files for later cleanup.
        self.add_file_to_clean(cline.outfile)
        if cline.guidetree_out:
            self.add_file_to_clean(cline.guidetree_out)

        input_records = SeqIO.to_dict(SeqIO.parse(cline.infile, "fasta"))
        self.assertEqual(str(eval(repr(cline))), str(cline))
        output, error = cline()
        self.assertTrue(not output or output.strip().startswith("CLUSTAL"))

        # Test if ClustalOmega executed successfully.
        self.assertTrue(
            error.strip() == ""
            or error.startswith(
                (
                    "WARNING: Sequence type is DNA.",
                    "WARNING: DNA alignment is still experimental.",
                )
            )
        )

        # TODO - Try and parse this with Bio.Nexus?
        if cline.guidetree_out:
            self.assertTrue(os.path.isfile(cline.guidetree_out))

    def add_file_to_clean(self, filename):
        """Add a file for deferred removal by the tearDown routine."""
        self.files_to_clean.add(filename)


#################################################################


class ClustalOmegaTestErrorConditions(ClustalOmegaTestCase):
    def test_empty_file(self):
        """Test an empty file."""
        input_file = "does_not_exist.fasta"
        self.assertFalse(os.path.isfile(input_file))
        cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            message = str(err)
            self.assertTrue(
                "Cannot open sequence file" in message
                or "Cannot open input file" in message
                or "Non-zero return code" in message,
                message,
            )
        else:
            self.fail(f"Should have failed, returned:\n{stdout}\n{stderr}")

    def test_single_sequence(self):
        """Test an input file containing a single sequence."""
        input_file = "Fasta/f001"
        self.assertTrue(os.path.isfile(input_file))
        self.assertEqual(len(list(SeqIO.parse(input_file, "fasta"))), 1)
        cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            self.assertIn("contains 1 sequence, nothing to align", str(err))
        else:
            self.fail(f"Should have failed, returned:\n{stdout}\n{stderr}")

    def test_invalid_format(self):
        """Test an input file in an invalid format."""
        input_file = "Medline/pubmed_result1.txt"
        self.assertTrue(os.path.isfile(input_file))
        cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
        with self.assertRaises(ApplicationError) as cm:
            stdout, stderr = cline()
            self.fail(f"Should have failed, returned:\n{stdout}\n{stderr}")
        err = str(cm.exception)
        # Ideally we'd catch the return code and raise the specific
        # error for "invalid format".
        self.assertIn("Can't determine format of sequence file", err)


#################################################################


class ClustalOmegaTestNormalConditions(ClustalOmegaTestCase):
    def test_simple_fasta(self):
        """Test a simple fasta file."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp_test.aln"

        cline = ClustalOmegaCommandline(
            clustalo_exe, infile=input_file, outfile=output_file, outfmt="clustal"
        )

        self.standard_test_procedure(cline)
        alignment = Align.read(cline.outfile, "clustal")
        self.assertEqual(
            str(alignment),
            """\
gi|134891         0 GATCCCTACCCTTNCCGTTGGTCTCTNTCGCTGACTCGAGGCACCTAACATCCATTCACA
                  0 ---------..-........|......|....|......|..............|.----
gi|129628         0 ---------MP-VVVVASSKGGAGKSTTAVVLGTELAHKGVPVTMLDCDPNRSLTI----

gi|134891        60 CCCAACACAGGCCAGCGACTTCTGGGGCTCAGCCACAGACATGGTTTGTNACTNTTGAGC
                 60 -----.|.||.......|....|-------------------......|.......||..
gi|129628        46 -----WANAGEVPENITALSDVT-------------------ESSIVKTIKQHDVDGAVV

gi|134891       120 TTCTGTTCCTAGAGAATCCTAGAGGCTTGATTGGCCCAGGCTGCTGTNTGTNCTGGAGG-
                120 ...--------..|.|......|..............|...|..............|..-
gi|129628        82 IVD--------LEGVASRMVSRAISQADLVLIPMRPKALDATIGAQSLQLIAEEEEAIDR

gi|134891       179 -CAAAGAATCCCTACCTCCTAGGGGTGAAAGGAAATNAAAATGGAAAGTTCTTGTAGCGC
                180 -.|.|...|....|.......|........|------------...........||....
gi|129628       134 KIAHAVVFTMVSPAIRSHEYTGIKASLIENG------------VEIIEPPLVERTAYSAL

gi|134891       238 AAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTTC 285
                240 ...|..........|..........|.|-----.|.....|.|..-- 287
gi|129628       182 FQFGGNLHSMKSKQGNMAAAIENAEAFA-----MAIFKKLTEALR-- 222
""",
        )
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"],
            "                    *      *    *      *              *           * **       *    *                         *       **               * *      *              *   *              *     * *   *    *       *        *                       **       *          *          * *      *     * *    ",
        )

    def test_properties(self):
        """Test setting options via properties."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp_test.aln"

        cline = ClustalOmegaCommandline(clustalo_exe)
        cline.infile = input_file
        cline.outfile = output_file
        cline.outfmt = "clustal"

        self.standard_test_procedure(cline)
        alignment = Align.read(cline.outfile, "clustal")
        self.assertEqual(
            str(alignment),
            """\
gi|134891         0 GATCCCTACCCTTNCCGTTGGTCTCTNTCGCTGACTCGAGGCACCTAACATCCATTCACA
                  0 ---------..-........|......|....|......|..............|.----
gi|129628         0 ---------MP-VVVVASSKGGAGKSTTAVVLGTELAHKGVPVTMLDCDPNRSLTI----

gi|134891        60 CCCAACACAGGCCAGCGACTTCTGGGGCTCAGCCACAGACATGGTTTGTNACTNTTGAGC
                 60 -----.|.||.......|....|-------------------......|.......||..
gi|129628        46 -----WANAGEVPENITALSDVT-------------------ESSIVKTIKQHDVDGAVV

gi|134891       120 TTCTGTTCCTAGAGAATCCTAGAGGCTTGATTGGCCCAGGCTGCTGTNTGTNCTGGAGG-
                120 ...--------..|.|......|..............|...|..............|..-
gi|129628        82 IVD--------LEGVASRMVSRAISQADLVLIPMRPKALDATIGAQSLQLIAEEEEAIDR

gi|134891       179 -CAAAGAATCCCTACCTCCTAGGGGTGAAAGGAAATNAAAATGGAAAGTTCTTGTAGCGC
                180 -.|.|...|....|.......|........|------------...........||....
gi|129628       134 KIAHAVVFTMVSPAIRSHEYTGIKASLIENG------------VEIIEPPLVERTAYSAL

gi|134891       238 AAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTTC 285
                240 ...|..........|..........|.|-----.|.....|.|..-- 287
gi|129628       182 FQFGGNLHSMKSKQGNMAAAIENAEAFA-----MAIFKKLTEALR-- 222
""",
        )
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"],
            "                    *      *    *      *              *           * **       *    *                         *       **               * *      *              *   *              *     * *   *    *       *        *                       **       *          *          * *      *     * *    ",
        )

    def test_input_filename_with_space(self):
        """Test an input filename containing a space."""
        input_file = "Clustalw/temp horses.fasta"
        with open(input_file, "w") as handle:
            SeqIO.write(SeqIO.parse("Phylip/hennigian.phy", "phylip"), handle, "fasta")
        output_file = "temp_test.aln"

        cline = ClustalOmegaCommandline(
            clustalo_exe, infile=input_file, outfile=output_file, outfmt="clustal"
        )

        self.add_file_to_clean(input_file)
        self.standard_test_procedure(cline)
        alignment = Align.read(cline.outfile, "clustal")
        self.assertEqual(
            str(alignment),
            """\
A                 0 -CACACACAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAA 40
B                 0 -CACACAACAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAA 40
C                 0 -CACAACAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAA 40
D                 0 -CAACAAAACAAAAAAAAACAAAAAAAAAAAAAAAAAAAAA 40
E                 0 -CAACAAAAACAAAAAAAACAAAAAAAAAAAAAAAAAAAAA 40
F                 0 ACAAAAAAAACACACAAAACAAAAAAAAAAAAAAAAAAAA- 40
G                 0 ACAAAAAAAACACAACAAACAAAAAAAAAAAAAAAAAAAA- 40
H                 0 ACAAAAAAAACAACAAAAACAAAAAAAAAAAAAAAAAAAA- 40
I                 0 ACAAAAAAAAACAAAACAACAAAAAAAAAAAAAAAAAAAA- 40
J                 0 ACAAAAAAAAACAAAAACACAAAAAAAAAAAAAAAAAAAA- 40
""",
        )
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"],
            " **               ********************** ",
        )

    def test_output_filename_with_spaces(self):
        """Test an output filename containing spaces."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp with spaces.aln"

        cline = ClustalOmegaCommandline(
            clustalo_exe, infile=input_file, outfile=output_file, outfmt="clustal"
        )
        self.standard_test_procedure(cline)
        alignment = Align.read(cline.outfile, "clustal")
        self.assertEqual(
            str(alignment),
            """\
gi|134891         0 GATCCCTACCCTTNCCGTTGGTCTCTNTCGCTGACTCGAGGCACCTAACATCCATTCACA
                  0 ---------..-........|......|....|......|..............|.----
gi|129628         0 ---------MP-VVVVASSKGGAGKSTTAVVLGTELAHKGVPVTMLDCDPNRSLTI----

gi|134891        60 CCCAACACAGGCCAGCGACTTCTGGGGCTCAGCCACAGACATGGTTTGTNACTNTTGAGC
                 60 -----.|.||.......|....|-------------------......|.......||..
gi|129628        46 -----WANAGEVPENITALSDVT-------------------ESSIVKTIKQHDVDGAVV

gi|134891       120 TTCTGTTCCTAGAGAATCCTAGAGGCTTGATTGGCCCAGGCTGCTGTNTGTNCTGGAGG-
                120 ...--------..|.|......|..............|...|..............|..-
gi|129628        82 IVD--------LEGVASRMVSRAISQADLVLIPMRPKALDATIGAQSLQLIAEEEEAIDR

gi|134891       179 -CAAAGAATCCCTACCTCCTAGGGGTGAAAGGAAATNAAAATGGAAAGTTCTTGTAGCGC
                180 -.|.|...|....|.......|........|------------...........||....
gi|129628       134 KIAHAVVFTMVSPAIRSHEYTGIKASLIENG------------VEIIEPPLVERTAYSAL

gi|134891       238 AAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTTC 285
                240 ...|..........|..........|.|-----.|.....|.|..-- 287
gi|129628       182 FQFGGNLHSMKSKQGNMAAAIENAEAFA-----MAIFKKLTEALR-- 222
""",
        )
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"],
            "                    *      *    *      *              *           * **       *    *                         *       **               * *      *              *   *              *     * *   *    *       *        *                       **       *          *          * *      *     * *    ",
        )

    def test_large_fasta_file(self):
        """Test a large fasta input file."""
        # Create a large input file by converting another example file
        # (See Bug 2804, this will produce so much output on stdout that
        # subprocess could suffer a deadlock and hang).  Using all the
        # records should show the deadlock but is very slow - just thirty
        # seems to lockup on Mac OS X, even 20 on Linux (without the fix).
        input_file = "temp_cw_prot.fasta"
        records = list(SeqIO.parse("NBRF/Cw_prot.pir", "pir"))[:40]
        with open(input_file, "w") as handle:
            SeqIO.write(records, handle, "fasta")
        del handle, records
        output_file = "temp_cw_prot.aln"

        cline = ClustalOmegaCommandline(
            clustalo_exe, infile=input_file, outfile=output_file, outfmt="clustal"
        )

        self.add_file_to_clean(input_file)
        self.standard_test_procedure(cline)
        alignment = Align.read(cline.outfile, "clustal")

    def test_newtree_files(self):
        """Test requesting a guide tree."""
        input_file = "Fasta/f002"
        output_file = "temp_test.aln"
        newtree_file = "temp_test.dnd"
        alignment_text = """\
gi|134891         0 CGGACCAGACGGACACAGGGAGAAGCTAGTTTCTTTCATGTGATTGANATNATGACTCTA
gi|134891         0 ---------CGGAGCCAGCGAGCATAT---------------------------------
gi|159293         0 ------------------------------------------------------------

gi|134891        60 CTCCTAAAAGGGAAAAANCAATATCCTTGTTTACAGAAGAGAAACAAACAAGCCCCACTC
gi|134891        18 ----------------------------------------------------GCTGCATG
gi|159293         0 --------------------------------------------GATCAAATCTGCACTG

gi|134891       120 AGCTCAGTCACAGGAGAGANCACAGAAAGTCTTAGGATCATGANCTCTGAA-AAAAAGAG
gi|134891        26 -------------------------AGGACCTTTCTATCTTACATTATGGC-TGGGAATC
gi|159293        16 TGTCTACATATAGGAAAGGTCCTGGTGTGTGCTAATGTTCCCAATGCAGGACTTGAGGAA

gi|134891       179 AAACCTTATCTTTNCTTTGTGGTTCCTTTAAACACACTCACACACACTTGGTCAGAGATG
gi|134891        60 TTACTCTTTCATCTG-------ATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATC
gi|159293        76 GAGCTCTGTTATATGTTTCCATTTCTCTTTATCAAAGATAACCAAACCTTATGGCCCTT-

gi|134891       239 CTGTGCTTCTTGGAAGCAAGGNCTCAAAGGCAAGGTGCACGC----------AGAGGGAC
gi|134891       113 CTGGTTTTACAGATGTGAAACTT----TCAAGAGATTTACTGACTTTCCTAGAATA----
gi|159293       135 ---ATAACAATGGAGGCACTGGCTGCCTCTTAATTTTCAATCATGGACCTAAAGAAGTAC

gi|134891       289 GTTTGA--GTCTGGGATGAAGCATGTNCGTATTATTTATATGATGGAATTTCACGTTTTT
gi|134891       165 --------GT--------------TTCTCTACTGGAAACCTGATGCTTTTATAAGCCATT
gi|159293       192 TCTGAAGGGTCTCAACAATGCCAGGTGGGGACAGATATACTCAGAGATTATCCAGGTCTG

gi|134891       347 ATGTNAAGCNTGACAACACCAGGCAGGTATGAGAGGA-AAGCAAGGCCCGTCCATNGCTG
gi|134891       203 GTGATTAGGATGACTGTTACAGGCTTAGCTTTGTGTGAAANCCAGTCACCTTT------C
gi|159293       252 CCTCCCAGCGAGCC-----------TGGA------GT-ACACCAGACCCTCCTAGAGAAA

gi|134891       406 TCCGTACNCTTACGGNTTGCTTGTNGGAGNCATTTNGGTATTGTTTGTTGTAANANCCAA
gi|134891       257 TCCTAGGTAATGAGTAGTGCTGTTCATATTACTNT-------AAGTTCTATAGCATACTT
gi|159293       294 TCTGTT------------------------------------ATAATTTACCACCCACTT

gi|134891       466 AANGGGCTTTGGNNTGGNAAAA----GGGCAGANNGGGGGGGTTGGTGTNGTTTTTTGG-
gi|134891       310 GCNATCCTTTANCCATGCTTATCATANGTACCATTTGAGGAATTGNTT-----TGCCCTT
gi|159293       318 ATCCACCTTTAAACTTGGGGAA----GGNNGCN------TTTCAAATTAAATTTAATCNT

gi|134891       521 GGGGANNNTTTNGATTTGG-------TNCCGGGNTTTNGTTTNCCNCGGNACCGGNTTTT
gi|134891       365 TTG-GGTTTNTTNTTGGTAA--ANNNTTCCCGGGTGGGGGNGGTNNNGAAA---------
gi|159293       368 NGGGGGNTTTTAAACTTTAACCCTTTTNCCNTTNTNGGGGTNGGNANTTGNCCCCNTTAA

gi|134891       574 GGTTGGGGNCCATTTNTGNGGGGCNTTGGNGTTNCNTTNCCCNNNTNNGANTGGTTTNA
gi|134891       413 -----------------------------------------------------------
gi|159293       428 AGGGGGNNCCCCT-NCNNGGGGGAATAA-AACAA----------NTTNNTTT--TTT--

gi|134891       633
gi|134891       413
gi|159293       471
"""
        clustal_consensus = "                                                                                                                      *                                 *    *          *              *  * *  *           *   **   ** *       * *  *         *            *     *              *  *  *             *               **               *    *         * *     *     *   *       **   * *                        *  * ** * *           **                                              *        *        ****      *   *      *                  *      *        *     * *               * **    *   *                                                                                "

        cline = ClustalOmegaCommandline(
            clustalo_exe,
            infile=input_file,
            outfile=output_file,
            guidetree_out=newtree_file,
            outfmt="clustal",
        )

        self.standard_test_procedure(cline)
        alignment = Align.read(cline.outfile, "clustal")
        self.assertEqual(str(alignment), alignment_text)
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"], clustal_consensus
        )

        cline.guidetree_out = "temp with space.dnd"
        self.standard_test_procedure(cline)
        alignment = Align.read(cline.outfile, "clustal")
        self.assertEqual(str(alignment), alignment_text)
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"], clustal_consensus
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
