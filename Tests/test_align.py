# Copyright 2000-2001 by Brad Chapman.  All rights reserved.
# Revisions copyright 2007-2003 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Test alignment stuff.

Right now we've got tests for:

- Reading and Writing clustal format
- Reading and Writing fasta format
- Converting between formats

"""

# standard library
import os
import unittest
import warnings
from io import StringIO

# biopython
from Bio import BiopythonDeprecationWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio import AlignIO, Align
from Bio.Align import MultipleSeqAlignment, Alignment
from Bio import motifs


class TestBasics(unittest.TestCase):
    def test_empty_alignment(self):
        """Very simple tests on an empty alignment."""
        alignment = MultipleSeqAlignment([])
        self.assertEqual(alignment.get_alignment_length(), 0)
        self.assertEqual(len(alignment), 0)
        alignment = alignment.alignment  # new-style Alignment object
        self.assertEqual(alignment.length, 0)
        self.assertEqual(len(alignment), 0)

    def test_basic_alignment(self):
        """Basic tests on a simple alignment of three sequences."""
        msa = MultipleSeqAlignment([])
        letters = "AbcDefGhiJklMnoPqrStuVwxYz"
        msa.append(SeqRecord(Seq(letters), id="mixed"))
        msa.append(SeqRecord(Seq(letters.lower()), id="lower"))
        msa.append(SeqRecord(Seq(letters.upper()), id="upper"))
        msa.append(SeqRecord(Seq(letters), id="duplicate"))
        del msa[3]
        self.assertEqual(msa.get_alignment_length(), 26)
        self.assertEqual(len(msa), 3)
        self.assertEqual(msa[0].seq, letters)
        self.assertEqual(msa[1].seq, letters.lower())
        self.assertEqual(msa[2].seq, letters.upper())
        self.assertEqual(msa[0].id, "mixed")
        self.assertEqual(msa[1].id, "lower")
        self.assertEqual(msa[2].id, "upper")
        for col, letter in enumerate(letters):
            self.assertEqual(msa[:, col], letter + letter.lower() + letter.upper())
        # Check row extractions:
        self.assertEqual(msa[0].id, "mixed")
        self.assertEqual(msa[-1].id, "upper")
        # Check sub-alignment extraction by row slicing:
        self.assertIsInstance(msa[::-1], MultipleSeqAlignment)
        self.assertEqual(msa[::-1][0].id, "upper")
        self.assertEqual(msa[::-1][2].id, "mixed")
        # create a new-style Alignment object
        alignment = msa.alignment
        self.assertEqual(alignment.shape, (3, 26))
        self.assertEqual(len(alignment), 3)
        self.assertEqual(alignment.sequences[0].seq, letters)
        self.assertEqual(alignment.sequences[1].seq, letters.lower())
        self.assertEqual(alignment.sequences[2].seq, letters.upper())
        self.assertEqual(alignment.sequences[0].id, "mixed")
        self.assertEqual(alignment.sequences[1].id, "lower")
        self.assertEqual(alignment.sequences[2].id, "upper")
        for col, letter in enumerate(letters):
            self.assertEqual(
                alignment[:, col], letter + letter.lower() + letter.upper()
            )
        # Check row extractions:
        self.assertEqual(alignment[0], letters)
        self.assertEqual(alignment[-1], letters.upper())
        # Check sub-alignment extraction by row slicing:
        self.assertIsInstance(alignment[::-1], Alignment)
        self.assertEqual(alignment[::-1].sequences[0].id, "upper")
        self.assertEqual(alignment[::-1].sequences[2].id, "mixed")


class TestReading(unittest.TestCase):
    def test_read_clustal1(self):
        """Parse an alignment file and get an alignment object."""
        opuntia_clustal_header = """\
CLUSTAL X (1.81) multiple sequence alignment


"""
        opuntia_clustal_body = """\
gi|6273285|gb|AF191659.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273284|gb|AF191658.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273287|gb|AF191661.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273286|gb|AF191660.1|AF191      TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273290|gb|AF191664.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273289|gb|AF191663.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273291|gb|AF191665.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
                                    ******* **** *************************************

gi|6273285|gb|AF191659.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273284|gb|AF191658.1|AF191      TATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273287|gb|AF191661.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
gi|6273286|gb|AF191660.1|AF191      TATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
gi|6273290|gb|AF191664.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273289|gb|AF191663.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273291|gb|AF191665.1|AF191      TATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
                                    ******          ********  **** ********* *********

gi|6273285|gb|AF191659.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGT
gi|6273284|gb|AF191658.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273287|gb|AF191661.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273286|gb|AF191660.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273290|gb|AF191664.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273289|gb|AF191663.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTAT
gi|6273291|gb|AF191665.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
                                    ************************************ *********** *

gi|6273285|gb|AF191659.1|AF191      ACCAGA
gi|6273284|gb|AF191658.1|AF191      ACCAGA
gi|6273287|gb|AF191661.1|AF191      ACCAGA
gi|6273286|gb|AF191660.1|AF191      ACCAGA
gi|6273290|gb|AF191664.1|AF191      ACCAGA
gi|6273289|gb|AF191663.1|AF191      ACCAGA
gi|6273291|gb|AF191665.1|AF191      ACCAGA
                                    ******


"""  # noqa : W291
        opuntia_fasta = """\
>gi|6273285|gb|AF191659.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA
>gi|6273284|gb|AF191658.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273287|gb|AF191661.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273286|gb|AF191660.1|AF191
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273290|gb|AF191664.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273289|gb|AF191663.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTATACCAGA
>gi|6273291|gb|AF191665.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
"""
        opuntia_fasta_oneline = """\
>gi|6273285|gb|AF191659.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA
>gi|6273284|gb|AF191658.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273287|gb|AF191661.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273286|gb|AF191660.1|AF191
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273290|gb|AF191664.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273289|gb|AF191663.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA
>gi|6273291|gb|AF191665.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
"""
        opuntia_fasta_oneline_with_description = """\
>gi|6273285|gb|AF191659.1|AF191 gi|6273285|gb|AF191659.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA
>gi|6273284|gb|AF191658.1|AF191 gi|6273284|gb|AF191658.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273287|gb|AF191661.1|AF191 gi|6273287|gb|AF191661.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273286|gb|AF191660.1|AF191 gi|6273286|gb|AF191660.1|AF191
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273290|gb|AF191664.1|AF191 gi|6273290|gb|AF191664.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273289|gb|AF191663.1|AF191 gi|6273289|gb|AF191663.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA
>gi|6273291|gb|AF191665.1|AF191 gi|6273291|gb|AF191665.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
"""
        path = os.path.join(os.getcwd(), "Clustalw", "opuntia.aln")
        msa = AlignIO.read(path, "clustal")
        opuntia_clustal = opuntia_clustal_header + opuntia_clustal_body
        self.assertEqual(format(msa, "clustal"), opuntia_clustal)
        self.assertEqual(format(msa, "fasta"), opuntia_fasta)
        # create a new-style Alignment object
        alignment = msa.alignment
        self.assertEqual(format(alignment, "clustal"), opuntia_clustal_body)
        # New-style Alignment objects generate FASTA format with the sequence
        # on one line. Also, the clustal parser in Bio.AlignIO generates
        # SeqRecords with an (identical) ID and a description; the clustal
        # parser in Bio.Align generates SeqRecords with an ID only.
        self.assertEqual(
            format(alignment, "fasta"), opuntia_fasta_oneline_with_description
        )
        alignment = Align.read(path, "clustal")
        self.assertEqual(format(alignment, "fasta"), opuntia_fasta_oneline)

    def test_read_clustal2(self):
        """Parse an alignment file and get an alignment object."""
        clustalw_clustal_header = """\
CLUSTAL X (1.81) multiple sequence alignment


"""
        clustalw_clustal_body = """\
gi|4959044|gb|AAD34209.1|AF069      MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNN
gi|671626|emb|CAA85685.1|           ---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFR
                                              * *: ::    :.   :*  :  :. : . :*  ::   .

gi|4959044|gb|AAD34209.1|AF069      LLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDW
gi|671626|emb|CAA85685.1|           VTPQPG-----------------VPPEEAGAAVAAESSTGT---------
                                    :   **                  **:...   *.*** ..         

gi|4959044|gb|AAD34209.1|AF069      LNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQT
gi|671626|emb|CAA85685.1|           WTTVWTDGLTSLDRYKG-----RCYHIEPVPG------------------
                                     .:*   * *: .* :*        : :* .*                  

gi|4959044|gb|AAD34209.1|AF069      SENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTE
gi|671626|emb|CAA85685.1|           -EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIP
                                     *::.  .    .:: :*..*  :* .*   .. .  :    .  :    

gi|4959044|gb|AAD34209.1|AF069      VPTTRAQRRARSRSPEHRRTRARAERSMSPLQPTSEIPRRAPTLEQSSEN
gi|671626|emb|CAA85685.1|           VAYVKTFQGPPHGIQVERDKLNKYGRPLLGCTIKPKLGLSAKNYGRAVYE
                                    *. .:: : .      .* .  :  *.:     ..::   * .  ::  :

gi|4959044|gb|AAD34209.1|AF069      EPEGSSRTRHHVTLRQQISGPELLGRGLFAASGSRNPSQGTSSSDTGSNS
gi|671626|emb|CAA85685.1|           CLRGGLDFTKDDENVNSQPFMRWRDRFLFCAEAIYKAQAETGEIKGHYLN
                                      .*.    :.    :. .  .  .* **.*..  :..  *.. .    .

gi|4959044|gb|AAD34209.1|AF069      ESSGSGQRPPTIVLDLQVRRVRPGEYRQRDSIASRTRSRSQAPNNTVTYE
gi|671626|emb|CAA85685.1|           ATAG-----------------------TCEEMIKRAIFARELGVPIVMHD
                                     ::*                         :.: .*:    :     * ::

gi|4959044|gb|AAD34209.1|AF069      SERGGFRRTFSRSERAGVRTYVSTIRIPIRRILNTGLSETTSVAIQTMLR
gi|671626|emb|CAA85685.1|           YLTGGFTANTSLAHYCRDNGLLLHIHRAMHAVIDRQKNHGMHFRVLAKAL
                                       ***  . * :. .  .  :  *: .:: :::   ..   . : :   

gi|4959044|gb|AAD34209.1|AF069      QIMTGFGELSYFMYSDSDSEPSASVSSRNVERVESRNGRGSSGGGNSSGS
gi|671626|emb|CAA85685.1|           RLSGGDHIHSGTVVGKLEGERDITLGFVDLLRDDFIEKDRSRGIYFTQDW
                                    ::  *    *  : .. :.* . ::.  :: * :  :   * *   :.. 

gi|4959044|gb|AAD34209.1|AF069      SSSSSPSPSSSGESSESSSKMFEGSSEGGSSGPSRKDGRHRAPVTFDESG
gi|671626|emb|CAA85685.1|           VSLPGVIPVASG-----------------------------GIHVWHMPA
                                     * ..  * :**                             .  .:. ..

gi|4959044|gb|AAD34209.1|AF069      SLPFFSLAQFFLLNEDDEDQPRGLTKEQIDNLAMRSFGENDALKTCSVCI
gi|671626|emb|CAA85685.1|           LTEIFGDDSVLQFGGGTLGHPWGNAPGAVANRVA-----------VEACV
                                       :*.  ..: :. .  .:* * :   : * .             ..*:

gi|4959044|gb|AAD34209.1|AF069      TEYTEGDKLRKLPCSHEFHVHCIDRWLSE-NSTCPICRRAVLSSGNRESV
gi|671626|emb|CAA85685.1|           KARNEG---RDLAAEGNAIIREACKWSPELAAACEVWKEIKFEFPAMD--
                                    .  .**   *.*... :  ::   :* .*  ::* : :.  :.    :  

gi|4959044|gb|AAD34209.1|AF069      V
gi|671626|emb|CAA85685.1|           -
                                     


"""  # noqa : W291

        path = os.path.join(os.curdir, "Clustalw", "clustalw.aln")
        msa = AlignIO.read(path, "clustal")
        clustalw_clustal = clustalw_clustal_header + clustalw_clustal_body
        self.assertEqual(format(msa, "clustal"), clustalw_clustal)
        # create a new-style Alignment object
        alignment = msa.alignment
        self.assertEqual(format(alignment, "clustal"), clustalw_clustal_body)

    def test_read_write_clustal(self):
        """Test the base alignment stuff."""
        path = os.path.join(os.getcwd(), "Clustalw", "opuntia.aln")
        msa = AlignIO.read(path, "clustal")
        self.assertEqual(len(msa), 7)
        seq_record = msa[0]
        self.assertEqual(seq_record.description, "gi|6273285|gb|AF191659.1|AF191")
        self.assertEqual(
            seq_record.seq,
            Seq(
                "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA"
            ),
        )
        seq_record = msa[1]
        self.assertEqual(seq_record.description, "gi|6273284|gb|AF191658.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = msa[2]
        self.assertEqual(seq_record.description, "gi|6273287|gb|AF191661.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = msa[3]
        self.assertEqual(seq_record.description, "gi|6273286|gb|AF191660.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = msa[4]
        self.assertEqual(seq_record.description, "gi|6273290|gb|AF191664.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = msa[5]
        self.assertEqual(seq_record.description, "gi|6273289|gb|AF191663.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA",
        )
        seq_record = msa[6]
        self.assertEqual(seq_record.description, "gi|6273291|gb|AF191665.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(msa.get_alignment_length(), 156)
        align_info = AlignInfo.SummaryInfo(msa)
        with self.assertWarns(BiopythonDeprecationWarning):
            consensus = align_info.dumb_consensus(ambiguous="N")
        self.assertIsInstance(consensus, Seq)
        self.assertEqual(
            consensus,
            "TATACATTAAAGNAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTNCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        with self.assertWarns(BiopythonDeprecationWarning):
            dictionary = align_info.replacement_dictionary(
                skip_chars=None, letters="ACGT"
            )
        self.assertEqual(len(dictionary), 16)
        self.assertAlmostEqual(dictionary[("A", "A")], 1395.0, places=1)
        self.assertAlmostEqual(dictionary[("A", "C")], 3.0, places=1)
        self.assertAlmostEqual(dictionary[("A", "G")], 13.0, places=1)
        self.assertAlmostEqual(dictionary[("A", "T")], 6.0, places=1)
        self.assertAlmostEqual(dictionary[("C", "A")], 3.0, places=1)
        self.assertAlmostEqual(dictionary[("C", "C")], 271.0, places=1)
        self.assertAlmostEqual(dictionary[("C", "G")], 0, places=1)
        self.assertAlmostEqual(dictionary[("C", "T")], 16.0, places=1)
        self.assertAlmostEqual(dictionary[("G", "A")], 5.0, places=1)
        self.assertAlmostEqual(dictionary[("G", "C")], 0, places=1)
        self.assertAlmostEqual(dictionary[("G", "G")], 480.0, places=1)
        self.assertAlmostEqual(dictionary[("G", "T")], 0, places=1)
        self.assertAlmostEqual(dictionary[("T", "A")], 6.0, places=1)
        self.assertAlmostEqual(dictionary[("T", "C")], 12.0, places=1)
        self.assertAlmostEqual(dictionary[("T", "G")], 0, places=1)
        self.assertAlmostEqual(dictionary[("T", "T")], 874.0, places=1)
        alignment = msa.alignment
        dictionary = alignment.substitutions
        self.assertEqual(len(dictionary), 4)
        self.assertEqual(dictionary.shape, (4, 4))
        self.assertEqual(len(dictionary.keys()), 16)
        self.assertAlmostEqual(dictionary[("A", "A")], 1395)
        self.assertAlmostEqual(dictionary[("A", "C")], 3)
        self.assertAlmostEqual(dictionary[("A", "G")], 13)
        self.assertAlmostEqual(dictionary[("A", "T")], 6)
        self.assertAlmostEqual(dictionary[("C", "A")], 3)
        self.assertAlmostEqual(dictionary[("C", "C")], 271)
        self.assertAlmostEqual(dictionary[("C", "G")], 0)
        self.assertAlmostEqual(dictionary[("C", "T")], 16)
        self.assertAlmostEqual(dictionary[("G", "A")], 5)
        self.assertAlmostEqual(dictionary[("G", "C")], 0)
        self.assertAlmostEqual(dictionary[("G", "G")], 480)
        self.assertAlmostEqual(dictionary[("G", "T")], 0)
        self.assertAlmostEqual(dictionary[("T", "A")], 6)
        self.assertAlmostEqual(dictionary[("T", "C")], 12)
        self.assertAlmostEqual(dictionary[("T", "G")], 0)
        self.assertAlmostEqual(dictionary[("T", "T")], 874)
        with self.assertWarns(BiopythonDeprecationWarning):
            matrix = align_info.pos_specific_score_matrix(consensus, ["N", "-"])

        motif = motifs.Motif("ACGT", alignment)
        counts = motif.counts
        for i in range(alignment.length):
            for letter in "ACGT":
                self.assertAlmostEqual(counts[letter][i], matrix[i][letter])
        self.assertEqual(counts.calculate_consensus(identity=0.7), consensus)
        self.assertEqual(
            str(matrix),
            """\
    A   C   G   T
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  1.0 0.0 0.0 6.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
N  4.0 0.0 3.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
C  0.0 7.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
C  0.0 7.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 4.0
A  4.0 0.0 0.0 0.0
T  0.0 0.0 0.0 3.0
A  3.0 0.0 0.0 0.0
T  0.0 0.0 0.0 1.0
A  1.0 0.0 0.0 0.0
T  0.0 0.0 0.0 1.0
A  1.0 0.0 0.0 0.0
T  0.0 0.0 0.0 1.0
A  1.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
C  1.0 6.0 0.0 0.0
A  6.0 0.0 0.0 1.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
N  0.0 3.0 0.0 4.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 2.0 0.0 5.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
T  0.0 1.0 0.0 6.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
T  0.0 0.0 0.0 7.0
G  1.0 0.0 6.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
""",
        )

        with self.assertWarns(BiopythonDeprecationWarning):
            matrix = align_info.pos_specific_score_matrix(chars_to_ignore=["N", "-"])

        alignment = msa.alignment
        motif = motifs.Motif("ACGT", alignment)
        counts = motif.counts
        for i in range(alignment.length):
            for letter in "ACGT":
                self.assertAlmostEqual(counts[letter][i], matrix[i][letter])

        self.assertEqual(
            str(matrix),
            """\
    A   C   G   T
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  1.0 0.0 0.0 6.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
X  4.0 0.0 3.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
C  0.0 7.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
C  0.0 7.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 4.0
A  4.0 0.0 0.0 0.0
T  0.0 0.0 0.0 3.0
A  3.0 0.0 0.0 0.0
T  0.0 0.0 0.0 1.0
A  1.0 0.0 0.0 0.0
T  0.0 0.0 0.0 1.0
A  1.0 0.0 0.0 0.0
T  0.0 0.0 0.0 1.0
A  1.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
C  1.0 6.0 0.0 0.0
A  6.0 0.0 0.0 1.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
X  0.0 3.0 0.0 4.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 2.0 0.0 5.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
T  0.0 1.0 0.0 6.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
T  0.0 0.0 0.0 7.0
G  1.0 0.0 6.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
""",
        )

        second_seq = msa[1].seq
        with self.assertWarns(BiopythonDeprecationWarning):
            matrix = align_info.pos_specific_score_matrix(second_seq, ["N", "-"])

        alignment = msa.alignment
        motif = motifs.Motif("ACGT", alignment)
        counts = motif.counts
        for i in range(alignment.length):
            for letter in "ACGT":
                self.assertAlmostEqual(counts[letter][i], matrix[i][letter])

        self.assertEqual(
            str(matrix),
            """\
    A   C   G   T
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  1.0 0.0 0.0 6.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  4.0 0.0 3.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
C  0.0 7.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
G  0.0 0.0 7.0 0.0
C  0.0 7.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 4.0
A  4.0 0.0 0.0 0.0
-  0.0 0.0 0.0 3.0
-  3.0 0.0 0.0 0.0
-  0.0 0.0 0.0 1.0
-  1.0 0.0 0.0 0.0
-  0.0 0.0 0.0 1.0
-  1.0 0.0 0.0 0.0
-  0.0 0.0 0.0 1.0
-  1.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
C  1.0 6.0 0.0 0.0
A  6.0 0.0 0.0 1.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
T  0.0 3.0 0.0 4.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
C  0.0 2.0 0.0 5.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
C  0.0 7.0 0.0 0.0
T  0.0 1.0 0.0 6.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
T  0.0 0.0 0.0 7.0
G  1.0 0.0 6.0 0.0
T  0.0 0.0 0.0 7.0
A  7.0 0.0 0.0 0.0
C  0.0 7.0 0.0 0.0
C  0.0 7.0 0.0 0.0
A  7.0 0.0 0.0 0.0
G  0.0 0.0 7.0 0.0
A  7.0 0.0 0.0 0.0
""",
        )
        e_freq_table = {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}
        with self.assertWarns(BiopythonDeprecationWarning):
            value = align_info.information_content(
                5, 50, chars_to_ignore=["N"], e_freq_table=e_freq_table
            )
        self.assertAlmostEqual(value, 88.42309908538343)  # MultipleSeqAlignment
        value = sum(motif[5:50].relative_entropy)
        self.assertAlmostEqual(value, 88.42309908538343)  # Alignment
        with self.assertWarns(BiopythonDeprecationWarning):
            value = align_info.information_content(
                e_freq_table=e_freq_table, chars_to_ignore=["N", "-"]
            )
        self.assertAlmostEqual(value, 306.2080592664532)  # MultipleSeqAlignment
        relative_entropy = motif.relative_entropy
        value = sum(relative_entropy)
        self.assertAlmostEqual(value, 306.2080592664532)  # Alignment
        self.assertEqual(align_info.get_column(1), "AAAAAAA")
        self.assertAlmostEqual(align_info.ic_vector[1], 2.00)
        self.assertEqual(align_info.get_column(7), "TTTATTT")
        self.assertAlmostEqual(align_info.ic_vector[7], 1.4083272214176725)
        self.assertAlmostEqual(relative_entropy[0], 2.0)
        self.assertAlmostEqual(relative_entropy[1], 2.0)
        self.assertAlmostEqual(relative_entropy[2], 2.0)
        self.assertAlmostEqual(relative_entropy[3], 2.0)
        self.assertAlmostEqual(relative_entropy[4], 2.0)
        self.assertAlmostEqual(relative_entropy[5], 2.0)
        self.assertAlmostEqual(relative_entropy[6], 2.0)
        self.assertAlmostEqual(relative_entropy[7], 1.4083272214176723)
        self.assertAlmostEqual(relative_entropy[8], 2.0)
        self.assertAlmostEqual(relative_entropy[9], 2.0)
        self.assertAlmostEqual(relative_entropy[10], 2.0)
        self.assertAlmostEqual(relative_entropy[11], 2.0)
        self.assertAlmostEqual(relative_entropy[12], 1.0147718639657484)
        self.assertAlmostEqual(relative_entropy[13], 2.0)
        self.assertAlmostEqual(relative_entropy[14], 2.0)
        self.assertAlmostEqual(relative_entropy[15], 2.0)
        self.assertAlmostEqual(relative_entropy[16], 2.0)
        self.assertAlmostEqual(relative_entropy[17], 2.0)
        self.assertAlmostEqual(relative_entropy[18], 2.0)
        self.assertAlmostEqual(relative_entropy[19], 2.0)
        self.assertAlmostEqual(relative_entropy[20], 2.0)
        self.assertAlmostEqual(relative_entropy[21], 2.0)
        self.assertAlmostEqual(relative_entropy[22], 2.0)
        self.assertAlmostEqual(relative_entropy[23], 2.0)
        self.assertAlmostEqual(relative_entropy[24], 2.0)
        self.assertAlmostEqual(relative_entropy[25], 2.0)
        self.assertAlmostEqual(relative_entropy[26], 2.0)
        self.assertAlmostEqual(relative_entropy[27], 2.0)
        self.assertAlmostEqual(relative_entropy[28], 2.0)
        self.assertAlmostEqual(relative_entropy[29], 2.0)
        self.assertAlmostEqual(relative_entropy[30], 2.0)
        self.assertAlmostEqual(relative_entropy[31], 2.0)
        self.assertAlmostEqual(relative_entropy[32], 2.0)
        self.assertAlmostEqual(relative_entropy[33], 2.0)
        self.assertAlmostEqual(relative_entropy[34], 2.0)
        self.assertAlmostEqual(relative_entropy[35], 2.0)
        self.assertAlmostEqual(relative_entropy[36], 2.0)
        self.assertAlmostEqual(relative_entropy[37], 2.0)
        self.assertAlmostEqual(relative_entropy[38], 2.0)
        self.assertAlmostEqual(relative_entropy[39], 2.0)
        self.assertAlmostEqual(relative_entropy[40], 2.0)
        self.assertAlmostEqual(relative_entropy[41], 2.0)
        self.assertAlmostEqual(relative_entropy[42], 2.0)
        self.assertAlmostEqual(relative_entropy[43], 2.0)
        self.assertAlmostEqual(relative_entropy[44], 2.0)
        self.assertAlmostEqual(relative_entropy[45], 2.0)
        self.assertAlmostEqual(relative_entropy[46], 2.0)
        self.assertAlmostEqual(relative_entropy[47], 2.0)
        self.assertAlmostEqual(relative_entropy[48], 2.0)
        self.assertAlmostEqual(relative_entropy[49], 2.0)
        self.assertAlmostEqual(relative_entropy[50], 2.0)
        self.assertAlmostEqual(relative_entropy[51], 2.0)
        self.assertAlmostEqual(relative_entropy[52], 2.0)
        self.assertAlmostEqual(relative_entropy[53], 2.0)
        self.assertAlmostEqual(relative_entropy[54], 2.0)
        self.assertAlmostEqual(relative_entropy[55], 2.0)
        self.assertAlmostEqual(relative_entropy[56], 2.0)
        self.assertAlmostEqual(relative_entropy[57], 2.0)
        self.assertAlmostEqual(relative_entropy[58], 2.0)
        self.assertAlmostEqual(relative_entropy[59], 2.0)
        self.assertAlmostEqual(relative_entropy[60], 2.0)
        self.assertAlmostEqual(relative_entropy[61], 2.0)
        self.assertAlmostEqual(relative_entropy[62], 2.0)
        self.assertAlmostEqual(relative_entropy[63], 2.0)
        self.assertAlmostEqual(relative_entropy[64], 2.0)
        self.assertAlmostEqual(relative_entropy[65], 2.0)
        self.assertAlmostEqual(relative_entropy[66], 2.0)
        self.assertAlmostEqual(relative_entropy[67], 2.0)
        self.assertAlmostEqual(relative_entropy[68], 2.0)
        self.assertAlmostEqual(relative_entropy[69], 2.0)
        self.assertAlmostEqual(relative_entropy[70], 2.0)
        self.assertAlmostEqual(relative_entropy[71], 2.0)
        self.assertAlmostEqual(relative_entropy[72], 2.0)
        self.assertAlmostEqual(relative_entropy[73], 2.0)
        self.assertAlmostEqual(relative_entropy[74], 1.4083272214176723)
        self.assertAlmostEqual(relative_entropy[75], 1.4083272214176723)
        self.assertAlmostEqual(relative_entropy[76], 2.0)
        self.assertAlmostEqual(relative_entropy[77], 2.0)
        self.assertAlmostEqual(relative_entropy[78], 2.0)
        self.assertAlmostEqual(relative_entropy[79], 2.0)
        self.assertAlmostEqual(relative_entropy[80], 1.0147718639657484)
        self.assertAlmostEqual(relative_entropy[81], 2.0)
        self.assertAlmostEqual(relative_entropy[82], 2.0)
        self.assertAlmostEqual(relative_entropy[83], 2.0)
        self.assertAlmostEqual(relative_entropy[84], 2.0)
        self.assertAlmostEqual(relative_entropy[85], 2.0)
        self.assertAlmostEqual(relative_entropy[86], 2.0)
        self.assertAlmostEqual(relative_entropy[87], 2.0)
        self.assertAlmostEqual(relative_entropy[88], 2.0)
        self.assertAlmostEqual(relative_entropy[89], 2.0)
        self.assertAlmostEqual(relative_entropy[90], 1.136879431433369)
        self.assertAlmostEqual(relative_entropy[91], 2.0)
        self.assertAlmostEqual(relative_entropy[92], 2.0)
        self.assertAlmostEqual(relative_entropy[93], 2.0)
        self.assertAlmostEqual(relative_entropy[94], 2.0)
        self.assertAlmostEqual(relative_entropy[95], 2.0)
        self.assertAlmostEqual(relative_entropy[96], 2.0)
        self.assertAlmostEqual(relative_entropy[97], 2.0)
        self.assertAlmostEqual(relative_entropy[98], 2.0)
        self.assertAlmostEqual(relative_entropy[99], 2.0)
        self.assertAlmostEqual(relative_entropy[100], 2.0)
        self.assertAlmostEqual(relative_entropy[101], 2.0)
        self.assertAlmostEqual(relative_entropy[102], 2.0)
        self.assertAlmostEqual(relative_entropy[103], 2.0)
        self.assertAlmostEqual(relative_entropy[104], 2.0)
        self.assertAlmostEqual(relative_entropy[105], 2.0)
        self.assertAlmostEqual(relative_entropy[106], 2.0)
        self.assertAlmostEqual(relative_entropy[107], 2.0)
        self.assertAlmostEqual(relative_entropy[108], 2.0)
        self.assertAlmostEqual(relative_entropy[109], 2.0)
        self.assertAlmostEqual(relative_entropy[110], 2.0)
        self.assertAlmostEqual(relative_entropy[111], 2.0)
        self.assertAlmostEqual(relative_entropy[112], 2.0)
        self.assertAlmostEqual(relative_entropy[113], 2.0)
        self.assertAlmostEqual(relative_entropy[114], 2.0)
        self.assertAlmostEqual(relative_entropy[115], 2.0)
        self.assertAlmostEqual(relative_entropy[116], 2.0)
        self.assertAlmostEqual(relative_entropy[117], 2.0)
        self.assertAlmostEqual(relative_entropy[118], 2.0)
        self.assertAlmostEqual(relative_entropy[119], 2.0)
        self.assertAlmostEqual(relative_entropy[120], 2.0)
        self.assertAlmostEqual(relative_entropy[121], 2.0)
        self.assertAlmostEqual(relative_entropy[122], 2.0)
        self.assertAlmostEqual(relative_entropy[123], 2.0)
        self.assertAlmostEqual(relative_entropy[124], 2.0)
        self.assertAlmostEqual(relative_entropy[125], 2.0)
        self.assertAlmostEqual(relative_entropy[126], 2.0)
        self.assertAlmostEqual(relative_entropy[127], 2.0)
        self.assertAlmostEqual(relative_entropy[128], 2.0)
        self.assertAlmostEqual(relative_entropy[129], 2.0)
        self.assertAlmostEqual(relative_entropy[130], 2.0)
        self.assertAlmostEqual(relative_entropy[131], 2.0)
        self.assertAlmostEqual(relative_entropy[132], 2.0)
        self.assertAlmostEqual(relative_entropy[133], 2.0)
        self.assertAlmostEqual(relative_entropy[134], 2.0)
        self.assertAlmostEqual(relative_entropy[135], 2.0)
        self.assertAlmostEqual(relative_entropy[136], 1.4083272214176723)
        self.assertAlmostEqual(relative_entropy[137], 2.0)
        self.assertAlmostEqual(relative_entropy[138], 2.0)
        self.assertAlmostEqual(relative_entropy[139], 2.0)
        self.assertAlmostEqual(relative_entropy[140], 2.0)
        self.assertAlmostEqual(relative_entropy[141], 2.0)
        self.assertAlmostEqual(relative_entropy[142], 2.0)
        self.assertAlmostEqual(relative_entropy[143], 2.0)
        self.assertAlmostEqual(relative_entropy[144], 2.0)
        self.assertAlmostEqual(relative_entropy[145], 2.0)
        self.assertAlmostEqual(relative_entropy[146], 2.0)
        self.assertAlmostEqual(relative_entropy[147], 2.0)
        self.assertAlmostEqual(relative_entropy[148], 1.4083272214176723)
        self.assertAlmostEqual(relative_entropy[149], 2.0)
        self.assertAlmostEqual(relative_entropy[150], 2.0)
        self.assertAlmostEqual(relative_entropy[151], 2.0)
        self.assertAlmostEqual(relative_entropy[152], 2.0)
        self.assertAlmostEqual(relative_entropy[153], 2.0)
        self.assertAlmostEqual(relative_entropy[154], 2.0)
        self.assertAlmostEqual(relative_entropy[155], 2.0)
        handle = StringIO()
        with self.assertWarns(BiopythonDeprecationWarning):
            AlignInfo.print_info_content(align_info, fout=handle)
        self.assertEqual(
            handle.getvalue(),
            """\
0 T 2.000
1 A 2.000
2 T 2.000
3 A 2.000
4 C 2.000
5 A 2.000
6 T 2.000
7 T 1.408
8 A 2.000
9 A 2.000
10 A 2.000
11 G 2.000
12 A 1.015
13 A 2.000
14 G 2.000
15 G 2.000
16 G 2.000
17 G 2.000
18 G 2.000
19 A 2.000
20 T 2.000
21 G 2.000
22 C 2.000
23 G 2.000
24 G 2.000
25 A 2.000
26 T 2.000
27 A 2.000
28 A 2.000
29 A 2.000
30 T 2.000
31 G 2.000
32 G 2.000
33 A 2.000
34 A 2.000
35 A 2.000
36 G 2.000
37 G 2.000
38 C 2.000
39 G 2.000
40 A 2.000
41 A 2.000
42 A 2.000
43 G 2.000
44 A 2.000
45 A 2.000
46 A 2.000
47 G 2.000
48 A 2.000
49 A 2.000
50 T 2.000
51 A 2.000
52 T 2.000
53 A 2.000
54 T 2.000
55 A 2.000
56 - 2.000
57 - 2.000
58 - 2.000
59 - 2.000
60 - 2.000
61 - 2.000
62 - 2.000
63 - 2.000
64 - 2.000
65 - 2.000
66 A 2.000
67 T 2.000
68 A 2.000
69 T 2.000
70 A 2.000
71 T 2.000
72 T 2.000
73 T 2.000
74 C 1.408
75 A 1.408
76 A 2.000
77 A 2.000
78 T 2.000
79 T 2.000
80 T 1.015
81 C 2.000
82 C 2.000
83 T 2.000
84 T 2.000
85 A 2.000
86 T 2.000
87 A 2.000
88 T 2.000
89 A 2.000
90 C 1.137
91 C 2.000
92 C 2.000
93 A 2.000
94 A 2.000
95 A 2.000
96 T 2.000
97 A 2.000
98 T 2.000
99 A 2.000
100 A 2.000
101 A 2.000
102 A 2.000
103 A 2.000
104 T 2.000
105 A 2.000
106 T 2.000
107 C 2.000
108 T 2.000
109 A 2.000
110 A 2.000
111 T 2.000
112 A 2.000
113 A 2.000
114 A 2.000
115 T 2.000
116 T 2.000
117 A 2.000
118 G 2.000
119 A 2.000
120 T 2.000
121 G 2.000
122 A 2.000
123 A 2.000
124 T 2.000
125 A 2.000
126 T 2.000
127 C 2.000
128 A 2.000
129 A 2.000
130 A 2.000
131 G 2.000
132 A 2.000
133 A 2.000
134 T 2.000
135 C 2.000
136 C 1.408
137 A 2.000
138 T 2.000
139 T 2.000
140 G 2.000
141 A 2.000
142 T 2.000
143 T 2.000
144 T 2.000
145 A 2.000
146 G 2.000
147 T 2.000
148 G 1.408
149 T 2.000
150 A 2.000
151 C 2.000
152 C 2.000
153 A 2.000
154 G 2.000
155 A 2.000
""",
        )
        # create a new-style Alignment object
        del seq_record
        del align_info
        del consensus
        del dictionary
        del matrix
        del second_seq
        del e_freq_table
        del value
        del handle
        alignment = msa.alignment
        self.assertEqual(len(alignment), 7)
        seq_record = alignment.sequences[0]
        self.assertEqual(seq_record.description, "gi|6273285|gb|AF191659.1|AF191")
        self.assertEqual(
            alignment[0],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(
            seq_record.seq,
            Seq(
                "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA"
            ),
        )
        seq_record = alignment.sequences[1]
        self.assertEqual(seq_record.description, "gi|6273284|gb|AF191658.1|AF191")
        self.assertEqual(
            alignment[1],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATAATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = alignment.sequences[2]
        self.assertEqual(seq_record.description, "gi|6273287|gb|AF191661.1|AF191")
        self.assertEqual(
            alignment[2],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = alignment.sequences[3]
        self.assertEqual(seq_record.description, "gi|6273286|gb|AF191660.1|AF191")
        self.assertEqual(
            alignment[3],
            "TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(
            seq_record.seq,
            "TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = alignment.sequences[4]
        self.assertEqual(seq_record.description, "gi|6273290|gb|AF191664.1|AF191")
        self.assertEqual(
            alignment[4],
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = alignment.sequences[5]
        self.assertEqual(seq_record.description, "gi|6273289|gb|AF191663.1|AF191")
        self.assertEqual(
            alignment[5],
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA",
        )
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA",
        )
        seq_record = alignment.sequences[6]
        self.assertEqual(seq_record.description, "gi|6273291|gb|AF191665.1|AF191")
        self.assertEqual(
            alignment[6],
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(alignment.shape, (7, 156))
        substitutions = alignment.substitutions
        self.assertEqual(len(substitutions), 4)
        self.assertEqual(substitutions.shape, (4, 4))
        self.assertAlmostEqual(substitutions[("A", "A")], 1395)
        self.assertAlmostEqual(substitutions[("A", "C")], 3)
        self.assertAlmostEqual(substitutions[("A", "G")], 13)
        self.assertAlmostEqual(substitutions[("A", "T")], 6)
        self.assertAlmostEqual(substitutions[("C", "A")], 3)
        self.assertAlmostEqual(substitutions[("C", "C")], 271)
        self.assertAlmostEqual(substitutions[("C", "G")], 0)
        self.assertAlmostEqual(substitutions[("C", "T")], 16)
        self.assertAlmostEqual(substitutions[("G", "A")], 5)
        self.assertAlmostEqual(substitutions[("G", "C")], 0)
        self.assertAlmostEqual(substitutions[("G", "G")], 480)
        self.assertAlmostEqual(substitutions[("G", "T")], 0)
        self.assertAlmostEqual(substitutions[("T", "A")], 6)
        self.assertAlmostEqual(substitutions[("T", "C")], 12)
        self.assertAlmostEqual(substitutions[("T", "G")], 0)
        self.assertAlmostEqual(substitutions[("T", "T")], 874)
        motif = motifs.Motif(alphabet="ACGT", alignment=alignment)
        self.assertEqual(
            motif.consensus,
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(
            motif.degenerate_consensus,
            "TATACATTAAAGRAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTYCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        matrix = motif.counts
        self.assertEqual(
            str(matrix),
            """\
        0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26     27     28     29     30     31     32     33     34     35     36     37     38     39     40     41     42     43     44     45     46     47     48     49     50     51     52     53     54     55     56     57     58     59     60     61     62     63     64     65     66     67     68     69     70     71     72     73     74     75     76     77     78     79     80     81     82     83     84     85     86     87     88     89     90     91     92     93     94     95     96     97     98     99    100    101    102    103    104    105    106    107    108    109    110    111    112    113    114    115    116    117    118    119    120    121    122    123    124    125    126    127    128    129    130    131    132    133    134    135    136    137    138    139    140    141    142    143    144    145    146    147    148    149    150    151    152    153    154    155
A:   0.00   7.00   0.00   7.00   0.00   7.00   0.00   1.00   7.00   7.00   7.00   0.00   4.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   1.00   6.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   1.00   0.00   7.00   0.00   0.00   7.00   0.00   7.00
C:   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   3.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   2.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   3.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00
T:   7.00   0.00   7.00   0.00   0.00   0.00   7.00   6.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   1.00   0.00   0.00   7.00   7.00   4.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   5.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   7.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00
""",
        )
        self.assertEqual(
            format(motif, "transfac"),
            """\
P0      A      C      G      T
01      0      0      0      7      T
02      7      0      0      0      A
03      0      0      0      7      T
04      7      0      0      0      A
05      0      7      0      0      C
06      7      0      0      0      A
07      0      0      0      7      T
08      1      0      0      6      T
09      7      0      0      0      A
10      7      0      0      0      A
11      7      0      0      0      A
12      0      0      7      0      G
13      4      0      3      0      R
14      7      0      0      0      A
15      0      0      7      0      G
16      0      0      7      0      G
17      0      0      7      0      G
18      0      0      7      0      G
19      0      0      7      0      G
20      7      0      0      0      A
21      0      0      0      7      T
22      0      0      7      0      G
23      0      7      0      0      C
24      0      0      7      0      G
25      0      0      7      0      G
26      7      0      0      0      A
27      0      0      0      7      T
28      7      0      0      0      A
29      7      0      0      0      A
30      7      0      0      0      A
31      0      0      0      7      T
32      0      0      7      0      G
33      0      0      7      0      G
34      7      0      0      0      A
35      7      0      0      0      A
36      7      0      0      0      A
37      0      0      7      0      G
38      0      0      7      0      G
39      0      7      0      0      C
40      0      0      7      0      G
41      7      0      0      0      A
42      7      0      0      0      A
43      7      0      0      0      A
44      0      0      7      0      G
45      7      0      0      0      A
46      7      0      0      0      A
47      7      0      0      0      A
48      0      0      7      0      G
49      7      0      0      0      A
50      7      0      0      0      A
51      0      0      0      7      T
52      7      0      0      0      A
53      0      0      0      7      T
54      7      0      0      0      A
55      0      0      0      7      T
56      7      0      0      0      A
57      0      0      0      4      T
58      4      0      0      0      A
59      0      0      0      3      T
60      3      0      0      0      A
61      0      0      0      1      T
62      1      0      0      0      A
63      0      0      0      1      T
64      1      0      0      0      A
65      0      0      0      1      T
66      1      0      0      0      A
67      7      0      0      0      A
68      0      0      0      7      T
69      7      0      0      0      A
70      0      0      0      7      T
71      7      0      0      0      A
72      0      0      0      7      T
73      0      0      0      7      T
74      0      0      0      7      T
75      1      6      0      0      C
76      6      0      0      1      A
77      7      0      0      0      A
78      7      0      0      0      A
79      0      0      0      7      T
80      0      0      0      7      T
81      0      3      0      4      Y
82      0      7      0      0      C
83      0      7      0      0      C
84      0      0      0      7      T
85      0      0      0      7      T
86      7      0      0      0      A
87      0      0      0      7      T
88      7      0      0      0      A
89      0      0      0      7      T
90      7      0      0      0      A
91      0      2      0      5      T
92      0      7      0      0      C
93      0      7      0      0      C
94      7      0      0      0      A
95      7      0      0      0      A
96      7      0      0      0      A
97      0      0      0      7      T
98      7      0      0      0      A
99      0      0      0      7      T
100      7      0      0      0      A
101      7      0      0      0      A
102      7      0      0      0      A
103      7      0      0      0      A
104      7      0      0      0      A
105      0      0      0      7      T
106      7      0      0      0      A
107      0      0      0      7      T
108      0      7      0      0      C
109      0      0      0      7      T
110      7      0      0      0      A
111      7      0      0      0      A
112      0      0      0      7      T
113      7      0      0      0      A
114      7      0      0      0      A
115      7      0      0      0      A
116      0      0      0      7      T
117      0      0      0      7      T
118      7      0      0      0      A
119      0      0      7      0      G
120      7      0      0      0      A
121      0      0      0      7      T
122      0      0      7      0      G
123      7      0      0      0      A
124      7      0      0      0      A
125      0      0      0      7      T
126      7      0      0      0      A
127      0      0      0      7      T
128      0      7      0      0      C
129      7      0      0      0      A
130      7      0      0      0      A
131      7      0      0      0      A
132      0      0      7      0      G
133      7      0      0      0      A
134      7      0      0      0      A
135      0      0      0      7      T
136      0      7      0      0      C
137      0      1      0      6      T
138      7      0      0      0      A
139      0      0      0      7      T
140      0      0      0      7      T
141      0      0      7      0      G
142      7      0      0      0      A
143      0      0      0      7      T
144      0      0      0      7      T
145      0      0      0      7      T
146      7      0      0      0      A
147      0      0      7      0      G
148      0      0      0      7      T
149      1      0      6      0      G
150      0      0      0      7      T
151      7      0      0      0      A
152      0      7      0      0      C
153      0      7      0      0      C
154      7      0      0      0      A
155      0      0      7      0      G
156      7      0      0      0      A
XX
//
""",
        )
        self.assertAlmostEqual(sum(motif[5:50].relative_entropy), 88.42309908538343)
        relative_entropy = motif.relative_entropy
        self.assertAlmostEqual(sum(relative_entropy[5:50]), 88.42309908538343)
        self.assertAlmostEqual(sum(relative_entropy), 306.20805926645323)
        self.assertEqual(alignment[:, 1], "AAAAAAA")
        self.assertAlmostEqual(motif.relative_entropy[1], 2.0)
        self.assertEqual(alignment[:, 7], "TTTATTT")
        self.assertAlmostEqual(relative_entropy[7], 1.4083272214176723)

    def test_read_fasta(self):
        path = os.path.join(os.curdir, "Quality", "example.fasta")
        msa = AlignIO.read(path, "fasta")
        self.assertEqual(len(msa), 3)
        seq_record = msa[0]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_413_324")
        self.assertEqual(seq_record.seq, "CCCTTCTTGTCTTCAGCGTTTCTCC")
        seq_record = msa[1]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_540_792")
        self.assertEqual(seq_record.seq, "TTGGCAGGCCAAGGCCGATGGATCA")
        seq_record = msa[2]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_443_348")
        self.assertEqual(seq_record.seq, "GTTGCTTCTGGCGTGGGTGGGGGGG")
        self.assertEqual(msa.get_alignment_length(), 25)
        align_info = AlignInfo.SummaryInfo(msa)
        with self.assertWarns(BiopythonDeprecationWarning):
            consensus = align_info.dumb_consensus(ambiguous="N", threshold=0.6)
        self.assertIsInstance(consensus, Seq)
        self.assertEqual(consensus, "NTNGCNTNNNNNGNNGGNTGGNTCN")
        self.assertEqual(
            str(msa),
            """\
Alignment with 3 rows and 25 columns
CCCTTCTTGTCTTCAGCGTTTCTCC EAS54_6_R1_2_1_413_324
TTGGCAGGCCAAGGCCGATGGATCA EAS54_6_R1_2_1_540_792
GTTGCTTCTGGCGTGGGTGGGGGGG EAS54_6_R1_2_1_443_348""",
        )
        alignment = msa.alignment
        self.assertEqual(len(alignment), 3)
        seq_record = alignment.sequences[0]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_413_324")
        self.assertEqual(seq_record.seq, "CCCTTCTTGTCTTCAGCGTTTCTCC")
        seq_record = alignment.sequences[1]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_540_792")
        self.assertEqual(seq_record.seq, "TTGGCAGGCCAAGGCCGATGGATCA")
        seq_record = alignment.sequences[2]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_443_348")
        self.assertEqual(seq_record.seq, "GTTGCTTCTGGCGTGGGTGGGGGGG")
        self.assertEqual(alignment.length, 25)
        motif = motifs.Motif(alphabet="ACGT", alignment=alignment)
        self.assertEqual(motif.consensus, "CTCGCATCCCAAGCAGGATGGATCA")
        self.assertEqual(motif.degenerate_consensus, "BYBKYHKBBBVHKBVSSDKKKVKSV")
        self.assertEqual(
            str(alignment),
            """\
EAS54_6_R         0 CCCTTCTTGTCTTCAGCGTTTCTCC 25
EAS54_6_R         0 TTGGCAGGCCAAGGCCGATGGATCA 25
EAS54_6_R         0 GTTGCTTCTGGCGTGGGTGGGGGGG 25
""",
        )
        self.assertEqual(motif.counts.calculate_consensus(identity=0.6), consensus)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
