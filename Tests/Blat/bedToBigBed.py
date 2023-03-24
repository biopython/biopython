"""Hello."""

import sys


from Bio import Align
from Bio.Align import bigbed


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.bigbed."
    ) from None


if False:
    # bedToBigBed -as=bed12.as anoGam3.bed anoGam3.chrom.sizes anoGam3.bb
    bigBedFileName = "anoGam3.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed12.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open(bigBedFileName, "rb")
    correct = stream.read()
    assert data == correct
    print("test 11 ok")
    # return len(data)
    sys.exit(1)


if False:
    # bedToBigBed -as=bed12.as ailMel1.bed ailMel1.chrom.sizes ailMel1.bb
    bigBedFileName = "ailMel1.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed12.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open(bigBedFileName, "rb")
    correct = stream.read()
    assert data == correct
    print("test 12 ok")
    # return len(data)


if False:
    # bedToBigBed -as=bed12.as bisBis1.bed bisBis1.chrom.sizes bisBis1.bb
    bigBedFileName = "bisBis1.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed12.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open(bigBedFileName, "rb")
    correct = stream.read()
    assert data == correct
    print("test 13 ok")
    # return len(data)


def test1():
    # bedToBigBed -as=bed12.as ucsc.bed hg38.chrom.sizes ucsc.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed12.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open(bigBedFileName, "rb")
    correct = stream.read()
    assert data == correct
    print("test 1 ok")
    return len(data)


def test2():
    # bedToBigBed -as=bed12.as -unc ucsc.bed hg38.chrom.sizes ucsc.unc.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed12.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        declaration=declaration,
        compress=False,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open("ucsc.unc.bb", "rb")
    correct = stream.read()
    assert data == correct
    print("test 2 ok")
    return len(data)


def test3():
    # cut -f 1-3 ucsc.bed > ucsc.bed3.bed
    # bedToBigBed -as=bed3.as -type=bed3 ucsc.bed3.bed hg38.chrom.sizes ucsc.bed3.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed3.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        bedN=3,
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open("ucsc.bed3.bb", "rb")
    correct = stream.read()
    assert data == correct
    print("test 3 ok")
    return len(data)


def test4():
    # cut -f 1-4 ucsc.bed > ucsc.bed4.bed
    # bedToBigBed -as=bed4.as -type=bed4 ucsc.bed4.bed hg38.chrom.sizes ucsc.bed4.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed4.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        bedN=4,
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open("ucsc.bed4.bb", "rb")
    correct = stream.read()
    assert data == correct
    print("test 4 ok")
    return len(data)


def test5():
    # cut -f 1-5 ucsc.bed > ucsc.bed5.bed
    # bedToBigBed -as=bed5.as -type=bed5 ucsc.bed5.bed hg38.chrom.sizes ucsc.bed5.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed5.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        bedN=5,
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open("ucsc.bed5.bb", "rb")
    correct = stream.read()
    assert data == correct
    print("test 5 ok")
    return len(data)


def test6():
    # cut -f 1-6 ucsc.bed > ucsc.bed6.bed
    # bedToBigBed -as=bed6.as -type=bed6 ucsc.bed6.bed hg38.chrom.sizes ucsc.bed6.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed6.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        bedN=6,
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open("ucsc.bed6.bb", "rb")
    correct = stream.read()
    assert data == correct
    print("test 6 ok")
    return len(data)


def test7():
    # cut -f 1-7 ucsc.bed > ucsc.bed7.bed
    # bedToBigBed -as=bed7.as -type=bed7 ucsc.bed7.bed hg38.chrom.sizes ucsc.bed7.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed7.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        bedN=7,
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open("ucsc.bed7.bb", "rb")
    correct = stream.read()
    assert data == correct
    print("test 7 ok")
    return len(data)


def test8():
    # cut -f 1-8 ucsc.bed > ucsc.bed8.bed
    # bedToBigBed -as=bed8.as -type=bed8 ucsc.bed8.bed hg38.chrom.sizes ucsc.bed8.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed8.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        bedN=8,
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open("ucsc.bed8.bb", "rb")
    correct = stream.read()
    assert data == correct
    print("test 8 ok")
    return len(data)


def test9():
    # cut -f 1-9 ucsc.bed > ucsc.bed9.bed
    # bedToBigBed -as=bed9.as -type=bed9 ucsc.bed9.bed hg38.chrom.sizes ucsc.bed9.bb
    bigBedFileName = "ucsc.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed9.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        bedN=9,
        declaration=declaration,
        compress=True,
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open("ucsc.bed9.bb", "rb")
    correct = stream.read()
    assert data == correct
    print("test 9 ok")
    return len(data)


def test10():
    # bedToBigBed -as=bed12.as -extraIndex=name ucsc.bed hg38.chrom.sizes ucsc.indexed.bb
    bigBedFileName = "ucsc.indexed.bb"
    alignments = Align.parse(bigBedFileName, "bigbed")
    data = open("bed12.as").read()
    declaration = bigbed.AutoSQLTable.from_string(data)
    output = open("test.bb", "wb")
    Align.write(
        alignments,
        output,
        "bigbed",
        declaration=declaration,
        extraIndex=["name"],
    )
    output.close()
    stream = open("test.bb", "rb")
    data = stream.read()
    stream.close()
    stream = open(bigBedFileName, "rb")
    correct = stream.read()
    assert data == correct
    print("test 10 ok")
    return len(data)


size1 = test1()
size2 = test2()
size3 = test3()
size4 = test4()
size5 = test5()
size6 = test6()
size7 = test7()
size8 = test8()
size9 = test9()
size10 = test10()
