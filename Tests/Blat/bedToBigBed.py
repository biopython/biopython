"""Hello."""


import sys
import unittest
import tempfile

from Bio import Align
from Bio.Align import bigbed


for i, argument in enumerate(sys.argv):
    if argument == "--large":
        large = True
        sys.argv.pop(i)
        break
else:
    large = False


class BinaryTestBaseClass(unittest.TestCase):
    def assertBinaryEqual(self, file1, file2):
        blocksize = 1024
        n = 0
        while True:
            data1 = file1.read(blocksize)
            data2 = file2.read(blocksize)
            if data1 == b"" and data2 == b"":
                return
            if data1 == data2:
                n += len(data1)
                continue
            n1 = n + len(data1)
            n2 = n + len(data2)
            if n1 < n2:
                return self.fail(f"unequal file sizes: {n1} bytes vs >= {n2} bytes")
            if n1 > n2:
                return self.fail(f"unequal file sizes: >= {n1} bytes vs {n2} bytes")
            for i, (c1, c2) in enumerate(zip(data1, data2)):
                if c1 != c2:
                    return self.fail(f"bytes at position {n1+i} differ: {c1} vs {c2}")


@unittest.skipUnless(large is True, "large file; use --large to run")
class TestAlign_big(BinaryTestBaseClass):
    def test_a_compressed(self):
        # bedToBigBed -as=bed12.as ucsc.bed hg38.chrom.sizes ucsc.bb
        with open("bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_b_uncompressed(self):
        # bedToBigBed -as=bed12.as -unc ucsc.bed hg38.chrom.sizes ucsc.unc.bb
        with open("bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.unc.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=False,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_c_bed3(self):
        # cut -f 1-3 ucsc.bed > ucsc.bed3.bed
        # bedToBigBed -as=bed3.as -type=bed3 ucsc.bed3.bed hg38.chrom.sizes ucsc.bed3.bb
        with open("bed3.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed3.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=3,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_d_bed4(self):
        # cut -f 1-4 ucsc.bed > ucsc.bed4.bed
        # bedToBigBed -as=bed4.as -type=bed4 ucsc.bed4.bed hg38.chrom.sizes ucsc.bed4.bb
        with open("bed4.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed4.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=4,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_e_bed5(self):
        # cut -f 1-5 ucsc.bed > ucsc.bed5.bed
        # bedToBigBed -as=bed5.as -type=bed5 ucsc.bed5.bed hg38.chrom.sizes ucsc.bed5.bb
        with open("bed5.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed5.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=5,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_f_bed6(self):
        # cut -f 1-6 ucsc.bed > ucsc.bed6.bed
        # bedToBigBed -as=bed6.as -type=bed6 ucsc.bed6.bed hg38.chrom.sizes ucsc.bed6.bb
        with open("bed6.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed6.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=6,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_g_bed7(self):
        # cut -f 1-7 ucsc.bed > ucsc.bed7.bed
        # bedToBigBed -as=bed7.as -type=bed7 ucsc.bed7.bed hg38.chrom.sizes ucsc.bed7.bb
        with open("bed7.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed7.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=7,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_h_bed8(self):
        # cut -f 1-8 ucsc.bed > ucsc.bed8.bed
        # bedToBigBed -as=bed8.as -type=bed8 ucsc.bed8.bed hg38.chrom.sizes ucsc.bed8.bb
        with open("bed8.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed8.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=8,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_i_bed9(self):
        # cut -f 1-9 ucsc.bed > ucsc.bed9.bed
        # bedToBigBed -as=bed9.as -type=bed9 ucsc.bed9.bed hg38.chrom.sizes ucsc.bed9.bb
        with open("bed9.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed9.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=9,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_j_extraindex(self):
        # bedToBigBed -as=bed12.as -extraIndex=name ucsc.bed hg38.chrom.sizes ucsc.indexed.bb
        with open("bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.indexed.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                extraIndex=["name"],
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_k_anogam(self):
        # bedToBigBed -as=bed12.as anoGam3.bed anoGam3.chrom.sizes anoGam3.bb
        with open("bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "anoGam3.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_l_ailmel(self):
        # bedToBigBed -as=bed12.as ailMel1.bed ailMel1.chrom.sizes ailMel1.bb
        with open("bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ailMel1.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_m_bisbis(self):
        # bedToBigBed -as=bed12.as bisBis1.bed bisBis1.chrom.sizes bisBis1.bb
        with open("bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "bisBis1.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
