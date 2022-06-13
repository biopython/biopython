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
from io import StringIO

# biopython
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


opuntia_clustal = """\
CLUSTAL X (1.81) multiple sequence alignment


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

clustalw_clustal = """\
CLUSTAL X (1.81) multiple sequence alignment


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


class TestBasics(unittest.TestCase):
    def test_empty_alignment(self):
        """Very simple tests on an empty alignment."""
        alignment = MultipleSeqAlignment([])
        self.assertEqual(alignment.get_alignment_length(), 0)
        self.assertEqual(len(alignment), 0)

    def test_basic_alignment(self):
        """Basic tests on a simple alignment of three sequences."""
        alignment = MultipleSeqAlignment([])
        letters = "AbcDefGhiJklMnoPqrStuVwxYz"
        alignment.append(SeqRecord(Seq(letters), id="mixed"))
        alignment.append(SeqRecord(Seq(letters.lower()), id="lower"))
        alignment.append(SeqRecord(Seq(letters.upper()), id="upper"))
        self.assertEqual(alignment.get_alignment_length(), 26)
        self.assertEqual(len(alignment), 3)
        self.assertEqual(alignment[0].seq, letters)
        self.assertEqual(alignment[1].seq, letters.lower())
        self.assertEqual(alignment[2].seq, letters.upper())
        self.assertEqual(alignment[0].id, "mixed")
        self.assertEqual(alignment[1].id, "lower")
        self.assertEqual(alignment[2].id, "upper")
        for (col, letter) in enumerate(letters):
            self.assertEqual(
                alignment[:, col], letter + letter.lower() + letter.upper()
            )
        # Check row extractions:
        self.assertEqual(alignment[0].id, "mixed")
        self.assertEqual(alignment[-1].id, "upper")
        # Check sub-alignment extraction by row slicing:
        self.assertIsInstance(alignment[::-1], MultipleSeqAlignment)
        self.assertEqual(alignment[::-1][0].id, "upper")
        self.assertEqual(alignment[::-1][2].id, "mixed")


class TestReading(unittest.TestCase):
    def test_read_clustal1(self):
        """Parse an alignment file and get an alignment object."""
        path = os.path.join(os.getcwd(), "Clustalw", "opuntia.aln")
        alignment = AlignIO.read(path, "clustal")
        self.assertEqual(format(alignment, "clustal"), opuntia_clustal)

    def test_read_clustal2(self):
        """Parse an alignment file and get an alignment object."""
        path = os.path.join(os.curdir, "Clustalw", "clustalw.aln")
        alignment = AlignIO.read(path, "clustal")
        self.assertEqual(format(alignment, "clustal"), clustalw_clustal)

    def test_read_write_clustal(self):
        """Test the base alignment stuff."""
        path = os.path.join(os.getcwd(), "Clustalw", "opuntia.aln")
        alignment = AlignIO.read(path, "clustal")
        self.assertEqual(len(alignment), 7)
        seq_record = alignment[0]
        self.assertEqual(seq_record.description, "gi|6273285|gb|AF191659.1|AF191")
        self.assertEqual(
            seq_record.seq,
            Seq(
                "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA"
            ),
        )
        seq_record = alignment[1]
        self.assertEqual(seq_record.description, "gi|6273284|gb|AF191658.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = alignment[2]
        self.assertEqual(seq_record.description, "gi|6273287|gb|AF191661.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = alignment[3]
        self.assertEqual(seq_record.description, "gi|6273286|gb|AF191660.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = alignment[4]
        self.assertEqual(seq_record.description, "gi|6273290|gb|AF191664.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        seq_record = alignment[5]
        self.assertEqual(seq_record.description, "gi|6273289|gb|AF191663.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA",
        )
        seq_record = alignment[6]
        self.assertEqual(seq_record.description, "gi|6273291|gb|AF191665.1|AF191")
        self.assertEqual(
            seq_record.seq,
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        self.assertEqual(alignment.get_alignment_length(), 156)
        align_info = AlignInfo.SummaryInfo(alignment)
        consensus = align_info.dumb_consensus()
        self.assertIsInstance(consensus, Seq)
        self.assertEqual(
            consensus,
            "TATACATTAAAGXAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTXCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
        )
        dictionary = align_info.replacement_dictionary(skip_chars=None, letters="ACGT")
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
        matrix = align_info.pos_specific_score_matrix(consensus, ["N", "-"])
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

        matrix = align_info.pos_specific_score_matrix(chars_to_ignore=["N", "-"])
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

        second_seq = alignment[1].seq
        matrix = align_info.pos_specific_score_matrix(second_seq, ["N", "-"])
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
        value = align_info.information_content(
            5, 50, chars_to_ignore=["N"], e_freq_table=e_freq_table
        )
        self.assertAlmostEqual(value, 88.42, places=2)
        value = align_info.information_content(
            e_freq_table=e_freq_table, chars_to_ignore=["N"]
        )
        self.assertAlmostEqual(value, 287.55, places=2)
        self.assertEqual(align_info.get_column(1), "AAAAAAA")
        self.assertAlmostEqual(align_info.ic_vector[1], 2.00, places=2)
        self.assertEqual(align_info.get_column(7), "TTTATTT")
        self.assertAlmostEqual(align_info.ic_vector[7], 1.41, places=2)
        handle = StringIO()
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
56 - 0.682
57 - 0.682
58 - 0.333
59 - 0.333
60 - -0.115
61 - -0.115
62 - -0.115
63 - -0.115
64 - -0.115
65 - -0.115
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

    def test_read_fasta(self):
        path = os.path.join(os.curdir, "Quality", "example.fasta")
        alignment = AlignIO.read(path, "fasta")
        self.assertEqual(len(alignment), 3)
        seq_record = alignment[0]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_413_324")
        self.assertEqual(seq_record.seq, "CCCTTCTTGTCTTCAGCGTTTCTCC")
        seq_record = alignment[1]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_540_792")
        self.assertEqual(seq_record.seq, "TTGGCAGGCCAAGGCCGATGGATCA")
        seq_record = alignment[2]
        self.assertEqual(seq_record.description, "EAS54_6_R1_2_1_443_348")
        self.assertEqual(seq_record.seq, "GTTGCTTCTGGCGTGGGTGGGGGGG")
        self.assertEqual(alignment.get_alignment_length(), 25)
        align_info = AlignInfo.SummaryInfo(alignment)
        consensus = align_info.dumb_consensus(ambiguous="N", threshold=0.6)
        self.assertIsInstance(consensus, Seq)
        self.assertEqual(consensus, "NTNGCNTNNNNNGNNGGNTGGNTCN")
        self.assertEqual(
            str(alignment),
            """\
Alignment with 3 rows and 25 columns
CCCTTCTTGTCTTCAGCGTTTCTCC EAS54_6_R1_2_1_413_324
TTGGCAGGCCAAGGCCGATGGATCA EAS54_6_R1_2_1_540_792
GTTGCTTCTGGCGTGGGTGGGGGGG EAS54_6_R1_2_1_443_348""",
        )

    def test_format_conversion(self):
        """Parse the alignment file and get an alignment object."""
        path = os.path.join(os.curdir, "Clustalw", "opuntia.aln")
        alignment = AlignIO.read(path, "clustal")
        self.assertEqual(format(alignment, "fasta"), opuntia_fasta)
        self.assertEqual(format(alignment, "clustal"), opuntia_clustal)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
