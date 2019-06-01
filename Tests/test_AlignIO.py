# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for AlignIO module."""

from __future__ import print_function


import unittest


import os
from Bio._py3k import StringIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

test_write_read_alignment_formats = sorted(AlignIO._FormatToWriter)
test_write_read_align_with_seq_count = test_write_read_alignment_formats \
                                     + ["fasta", "tab"]

# test_files is a list of tuples containing:
# - string:  file format
# - integer: number of sequences per alignment
# - integer: number of alignments
# - string:  relative filename
#
# Most of the input files are also used by test_SeqIO.py,
# and by other additional tests as noted below.
test_files = [
# Following examples are also used in test_Clustalw.py
    ("clustal", 2, 1, 'Clustalw/cw02.aln'),
    ("clustal", 7, 1, 'Clustalw/opuntia.aln'),
    ("clustal", 5, 1, 'Clustalw/hedgehog.aln'),
    ("clustal", 2, 1, 'Clustalw/odd_consensus.aln'),
    ("clustal", 20, 1, 'Clustalw/protein.aln'),  # Used in the tutorial
    ("clustal", 20, 1, 'Clustalw/promals3d.aln'),  # Nonstandard header
# Following examples are also used in test_GFF.py
    ("fasta", 3, 1, 'GFF/multi.fna'),  # Trivial nucleotide alignment
# Following example is also used in test_Nexus.py
    ("nexus", 9, 1, 'Nexus/test_Nexus_input.nex'),
    ("nexus", 2, 1, 'Nexus/codonposset.nex'),
    ("stockholm", 2, 1, 'Stockholm/simple.sth'),
    ("stockholm", 6, 1, 'Stockholm/funny.sth'),
    ("phylip", 6, 1, 'Phylip/reference_dna.phy'),
    ("phylip", 6, 1, 'Phylip/reference_dna2.phy'),
    ("phylip", 10, 1, 'Phylip/hennigian.phy'),
    ("phylip", 10, 1, 'Phylip/horses.phy'),
    ("phylip", 10, 1, 'Phylip/random.phy'),
    ("phylip", 3, 1, 'Phylip/interlaced.phy'),
    ("phylip", 4, 1, 'Phylip/interlaced2.phy'),
    ("phylip-relaxed", 12, 1, 'ExtendedPhylip/primates.phyx'),
    ("phylip-sequential", 3, 1, 'Phylip/sequential.phy'),
    ("phylip-sequential", 4, 1, 'Phylip/sequential2.phy'),
    ("emboss", 4, 1, 'Emboss/alignret.txt'),
    ("emboss", 2, 5, 'Emboss/needle.txt'),
    ("emboss", 2, 1, 'Emboss/needle_asis.txt'),
    ("emboss", 2, 1, 'Emboss/water.txt'),
    ("emboss", 2, 1, 'Emboss/water2.txt'),
    ("emboss", 2, 1, 'Emboss/matcher_simple.txt'),
    ("emboss", 2, 5, 'Emboss/matcher_pair.txt'),
    ("emboss", 2, 1, 'Emboss/emboss_pair_aln_full_blank_line.txt'),
    ("fasta-m10", 2, 4, 'Fasta/output001.m10'),
    ("fasta-m10", 2, 6, 'Fasta/output002.m10'),
    ("fasta-m10", 2, 3, 'Fasta/output003.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output004.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output005.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output006.m10'),
    ("fasta-m10", 2, 9, 'Fasta/output007.m10'),
    ("fasta-m10", 2, 12, 'Fasta/output008.m10'),
    ("ig", 16, 1, 'IntelliGenetics/VIF_mase-pro.txt'),
    ("pir", 2, 1, 'NBRF/clustalw.pir'),
    ("maf", 3, 2, 'MAF/humor.maf'),
    ("maf", None, 3, "MAF/bug2453.maf"),  # Have 5, 5, 4 sequences
    ("maf", None, 3, "MAF/ucsc_test.maf"),  # Have 5, 5, 4 sequences
    ("maf", None, 48, "MAF/ucsc_mm9_chr10.maf"),
    ("mauve", None, 5, 'Mauve/simple.xmfa'),
    ]


def str_summary(text, max_len=40):
    """Cuts *text* if too long for summaries."""
    if len(text) <= max_len:
        return text
    else:
        return text[:max_len - 4] + "..." + text[-3:]


def alignment_summary(alignment, vertical_threshold=5):
    """Return a concise summary of an Alignment object as a string."""
    lines = []
    alignment_len = alignment.get_alignment_length()
    rec_count = len(alignment)
    if rec_count < vertical_threshold:
        # Show each sequence row horizontally
        for record in alignment:
            line = "  %s %s" % (str_summary(str(record.seq)), record.id)
            lines.append(line)
        for key, value in alignment.column_annotations.items():
            if isinstance(value, str):
                line = "  %s %s" % (str_summary(str(value)), key)
                lines.append(line)
    else:
        # Show each sequence row vertically
        for i in range(min(5, alignment_len)):
            line = "  " + str_summary(alignment[:, i]) + " alignment column %i" % i
            lines.append(line)
        if alignment_len > 5:
            i = alignment_len - 1
            line = "  " + str_summary("|" * rec_count) + " ..."
            lines.append(line)
            line = "  " + str_summary(alignment[:, i]) + " alignment column %i" % i
            lines.append(line)
    return "\n".join(lines)


def check_simple_write_read(alignments, indent=" "):
    # print(indent+"Checking we can write and then read back these alignments")
    for format in test_write_read_align_with_seq_count:
        records_per_alignment = len(alignments[0])
        for a in alignments:
            if records_per_alignment != len(a):
                records_per_alignment = None
        # Can we expect this format to work?
        if not records_per_alignment \
                and format not in test_write_read_alignment_formats:
            continue

        print(indent + "Checking can write/read as '%s' format" % format)

        # Going to write to a handle...
        handle = StringIO()

        try:
            c = AlignIO.write(alignments, handle=handle, format=format)
            assert c == len(alignments)
        except ValueError as e:
            # This is often expected to happen, for example when we try and
            # write sequences of different lengths to an alignment file.
            print(indent + "Failed: %s" % str(e))
            # Carry on to the next format:
            continue

        # First, try with the seq_count
        if records_per_alignment:
            handle.flush()
            handle.seek(0)
            try:
                alignments2 = list(AlignIO.parse(
                    handle=handle,
                    format=format,
                    seq_count=records_per_alignment))
            except ValueError as e:
                # This is BAD.  We can't read our own output.
                # I want to see the output when called from the test harness,
                # run_tests.py (which can be funny about new lines on Windows)
                handle.seek(0)
                raise ValueError("%s\n\n%s\n\n%s" % (
                    str(e), repr(handle.read()), repr(alignments2)))
            simple_alignment_comparison(alignments, alignments2, format)

        if format in test_write_read_alignment_formats:
            # Don't need the seq_count
            handle.flush()
            handle.seek(0)
            try:
                alignments2 = list(AlignIO.parse(handle=handle, format=format))
            except ValueError as e:
                # This is BAD.  We can't read our own output.
                # I want to see the output when called from the test harness,
                # run_tests.py (which can be funny about new lines on Windows)
                handle.seek(0)
                raise ValueError("%s\n\n%s\n\n%s" % (
                    str(e), repr(handle.read()), repr(alignments2)))
            simple_alignment_comparison(alignments, alignments2, format)

        if len(alignments) > 1:
            # Try writing just one Alignment (not a list)
            handle = StringIO()
            AlignIO.write(alignments[0:1], handle, format)
            assert handle.getvalue() == alignments[0].format(format)


def simple_alignment_comparison(alignments, alignments2, format):
    assert len(alignments) == len(alignments2)
    for a1, a2 in zip(alignments, alignments2):
        assert a1.get_alignment_length() == a2.get_alignment_length()
        assert len(a1) == len(a2)
        for r1, r2 in zip(a1, a2):
            # Check the bare minimum (ID and sequence) as
            # many formats can't store more than that.

            # Check the sequence
            assert str(r1.seq) == str(r2.seq)

            # Beware of different quirks and limitations in the
            # valid character sets and the identifier lengths!
            if format in ["phylip", "phylip-sequential"]:
                assert r1.id.replace("[", "").replace("]", "")[:10] == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format == "phylip-relaxed":
                assert r1.id.replace(" ", "").replace(':', '|') == r2.id, \
                        "'%s' vs '%s'" % (r1.id, r2.id)
            elif format == "clustal":
                assert r1.id.replace(" ", "_")[:30] == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format in ["stockholm", "maf"]:
                assert r1.id.replace(" ", "_") == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format == "fasta":
                assert r1.id.split()[0] == r2.id
            else:
                assert r1.id == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)


# Main tests...
for (t_format, t_per, t_count, t_filename) in test_files[:0]:
    print("Testing reading %s format file %s with %i alignments"
          % (t_format, t_filename, t_count))
    assert os.path.isfile(t_filename), t_filename

    # Try as an iterator using handle
    with open(t_filename, "r") as handle:
        alignments = list(AlignIO.parse(handle, format=t_format))
    msg = "Found %i alignments but expected %i" % (len(alignments), t_count)
    assert len(alignments) == t_count, msg
    if t_per is not None:
        for alignment in alignments:
            assert len(alignment) == t_per, \
                "Expected %i records per alignment, got %i" \
                % (t_per, len(alignment))

    # Try using the iterator with a for loop and a filename not handle
    alignments2 = []
    for record in AlignIO.parse(t_filename, format=t_format):
        alignments2.append(record)
    assert len(alignments2) == t_count

    # Try using the iterator with the next() method
    alignments3 = []
    seq_iterator = AlignIO.parse(t_filename, format=t_format)
    while True:
        try:
            record = next(seq_iterator)
        except StopIteration:
            break
        assert record is not None, "Should raise StopIteration not return None"
        alignments3.append(record)

    # Try a mixture of next() and list (a torture test!)
    seq_iterator = AlignIO.parse(t_filename, format=t_format)
    try:
        record = next(seq_iterator)
    except StopIteration:
        record = None
    if record is not None:
        alignments4 = [record]
        alignments4.extend(list(seq_iterator))
    else:
        alignments4 = []
    assert len(alignments4) == t_count

    # Try a mixture of next() and for loop (a torture test!)
    seq_iterator = AlignIO.parse(t_filename, format=t_format)
    try:
        record = next(seq_iterator)
    except StopIteration:
        record = None
    if record is not None:
        alignments5 = [record]
        for record in seq_iterator:
            alignments5.append(record)
    else:
        alignments5 = []
    assert len(alignments5) == t_count

    # Check Bio.AlignIO.read(...)
    if t_count == 1:
        with open(t_filename) as handle:
            alignment = AlignIO.read(handle, format=t_format)
        assert isinstance(alignment, MultipleSeqAlignment)
    else:
        try:
            with open(t_filename) as handle:
                alignment = AlignIO.read(handle, t_format)
            raise ValueError("Bio.AlignIO.read(...) should have failed")
        except ValueError:
            # Expected to fail
            pass

    # Show the alignment
    for i, alignment in enumerate(alignments):
        if i < 3 or i + 1 == t_count:
            print(" Alignment %i, with %i sequences of length %i"
                  % (i,
                     len(alignment),
                     alignment.get_alignment_length()))
            print(alignment_summary(alignment))
        elif i == 3:
            print(" ...")

    # Check AlignInfo.SummaryInfo likes the alignment
    summary = AlignInfo.SummaryInfo(alignment)
    dumb_consensus = summary.dumb_consensus()
    # gap_consensus = summary.gap_consensus()
    if t_format != "nexus":
        # Hack for bug 2535
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        try:
            info_content = summary.information_content()
        except ValueError as err:
            if str(err) != "Error in alphabet: not Nucleotide or Protein, supply expected frequencies":
                raise err

    if t_count == 1 and t_format not in ["nexus", "emboss", "fasta-m10"]:
        # print(" Trying to read a triple concatenation of the input file")
        with open(t_filename, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        assert len(list(AlignIO.parse(handle=handle, format=t_format, seq_count=t_per))) == 3
        handle.close()

    # Some alignment file formats have magic characters which mean
    # use the letter in this position in the first sequence.
    # They should all have been converted by the parser, but if
    # not reversing the record order might expose an error.  Maybe.
    alignments.reverse()
    check_simple_write_read(alignments)

print("Finished tested reading files")


class TestAlignIO_exceptions(unittest.TestCase):

    t_formats = list(AlignIO._FormatToWriter) + list(SeqIO._FormatToWriter)

    def test_phylip_reject_duplicate(self):
        """Check that writing duplicated IDs after truncation fails for PHYLIP."""
        handle = StringIO()
        sequences = [SeqRecord(Seq('AAAA'), id='longsequencename1'),
                     SeqRecord(Seq('AAAA'), id='longsequencename2'),
                     SeqRecord(Seq('AAAA'), id='other_sequence'), ]
        alignment = MultipleSeqAlignment(sequences)
        with self.assertRaises(ValueError) as cm:
            AlignIO.write(alignment, handle, 'phylip')
        self.assertEqual("Repeated name 'longsequen' (originally 'longsequencename2'), possibly due to truncation", str(cm.exception))

    def test_parsing_empty_files(self):
        """Check that parsing an empty file returns an empty list."""
        for t_format in AlignIO._FormatToIterator:
            handle = StringIO()
            alignments = list(AlignIO.parse(handle, t_format))
            self.assertEqual(alignments, [])

    def test_writing_empty_files(self):
        """Check that writers can cope with no alignments."""
        for t_format in self.t_formats:
            handle = StringIO()
            number = AlignIO.write([], handle, t_format)
            self.assertEqual(number, 0)

    def test_writing_not_alignments(self):
        """Check that writers reject records that are not alignments."""
        records = list(AlignIO.read("Clustalw/opuntia.aln", "clustal"))
        for t_format in self.t_formats:
            handle = StringIO()
            self.assertRaises(Exception, AlignIO.write, [records], handle, t_format)


class TestAlignIO_reading(unittest.TestCase):

    def check_simple_write_read(self, alignments, indent=" "):
        for format in test_write_read_align_with_seq_count:
            records_per_alignment = len(alignments[0])
            for a in alignments:
                if records_per_alignment != len(a):
                    records_per_alignment = None
            # Can we expect this format to work?
            if not records_per_alignment \
                    and format not in test_write_read_alignment_formats:
                continue
    
            # Going to write to a handle...
            handle = StringIO()

            if format == 'nexus':
                with self.assertRaises(ValueError) as cm:
                    c = AlignIO.write(alignments, handle=handle, format=format)
                self.assertEqual("We can only write one Alignment to a Nexus file.", str(cm.exception))
                continue
            c = AlignIO.write(alignments, handle=handle, format=format)
            self.assertEqual(c, len(alignments))

            # First, try with the seq_count
            if records_per_alignment:
                handle.flush()
                handle.seek(0)
                alignments2 = list(AlignIO.parse(
                    handle=handle,
                    format=format,
                    seq_count=records_per_alignment))
                simple_alignment_comparison(alignments, alignments2, format)
    
            if format in test_write_read_alignment_formats:
                # Don't need the seq_count
                handle.flush()
                handle.seek(0)
                alignments2 = list(AlignIO.parse(handle=handle, format=format))
                simple_alignment_comparison(alignments, alignments2, format)
    
            if len(alignments) > 1:
                # Try writing just one Alignment (not a list)
                handle = StringIO()
                AlignIO.write(alignments[0:1], handle, format)
                self.assertEqual(handle.getvalue(), alignments[0].format(format))
    
    def test_reading_alignments_clustal(self):
        path = 'Clustalw/cw02.aln'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="clustal"))
            self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments2 = []
        for record in AlignIO.parse(path, format="clustal"):
            alignments2.append(record)
        self.assertEqual(len(alignments2), 1)
    
        # Try using the iterator with the next() method
        alignments3 = []
        seq_iterator = AlignIO.parse(path, format="clustal")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments3.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="clustal")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.get_alignment_length(), 601)
        self.assertEqual(alignment_summary(alignment), """\
  MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVN...SVV gi|4959044|gb|AAD34209.1|AF069
  ---------MSPQTETKASVGFKAGVKEYKLTYYTP...--- gi|671626|emb|CAA85685.1|
            * *: ::    :.   :*  :  :. ...    clustal_consensus""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=2))), 3)
        handle.close()
        path = 'Clustalw/opuntia.aln'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="clustal"))
            self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 7)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments2 = []
        for record in AlignIO.parse(path, format="clustal"):
            alignments2.append(record)
        self.assertEqual(len(alignments2), 1)
    
        # Try using the iterator with the next() method
        alignments3 = []
        seq_iterator = AlignIO.parse(path, format="clustal")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments3.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="clustal")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignment), 7)
        self.assertEqual(alignment.get_alignment_length(), 156)
        self.assertEqual(alignment_summary(alignment), """\
  TTTTTTT alignment column 0
  AAAAAAA alignment column 1
  TTTTTTT alignment column 2
  AAAAAAA alignment column 3
  CCCCCCC alignment column 4
  ||||||| ...
  AAAAAAA alignment column 155""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=7))), 3)
        handle.close()
        path = 'Clustalw/hedgehog.aln'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="clustal"))
            self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 5)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments2 = []
        for record in AlignIO.parse(path, format="clustal"):
            alignments2.append(record)
        self.assertEqual(len(alignments2), 1)
    
        # Try using the iterator with the next() method
        alignments3 = []
        seq_iterator = AlignIO.parse(path, format="clustal")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments3.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="clustal")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignment), 5)
        self.assertEqual(alignment.get_alignment_length(), 447)
        self.assertEqual(alignment_summary(alignment), """\
  M---- alignment column 0
  F---- alignment column 1
  N---- alignment column 2
  L---- alignment column 3
  V---- alignment column 4
  ||||| ...
  ---SS alignment column 446""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=5))), 3)
        handle.close()
        path = 'Clustalw/odd_consensus.aln'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="clustal"))
            self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments2 = []
        for record in AlignIO.parse(path, format="clustal"):
            alignments2.append(record)
        self.assertEqual(len(alignments2), 1)
    
        # Try using the iterator with the next() method
        alignments3 = []
        seq_iterator = AlignIO.parse(path, format="clustal")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments3.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="clustal")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.get_alignment_length(), 687)
        self.assertEqual(alignment_summary(alignment), """\
  ------------------------------------...TAG AT3G20900.1-CDS
  ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGT...TAG AT3G20900.1-SEQ
                                      ...*** clustal_consensus""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=2))), 3)
        handle.close()
        path = 'Clustalw/protein.aln'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="clustal"))
            self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 20)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments2 = []
        for record in AlignIO.parse(path, format="clustal"):
            alignments2.append(record)
        self.assertEqual(len(alignments2), 1)
    
        # Try using the iterator with the next() method
        alignments3 = []
        seq_iterator = AlignIO.parse(path, format="clustal")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments3.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="clustal")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignment), 20)
        self.assertEqual(alignment.get_alignment_length(), 411)
        self.assertEqual(alignment_summary(alignment), """\
  -M------------------ alignment column 0
  -T------------------ alignment column 1
  -V------------------ alignment column 2
  -L-----------------M alignment column 3
  -E---------------MMS alignment column 4
  |||||||||||||||||||| ...
  -------------------T alignment column 410""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=20))), 3)
        handle.close()
        path = 'Clustalw/promals3d.aln'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="clustal"))
            self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 20)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments2 = []
        for record in AlignIO.parse(path, format="clustal"):
            alignments2.append(record)
        self.assertEqual(len(alignments2), 1)
    
        # Try using the iterator with the next() method
        alignments3 = []
        seq_iterator = AlignIO.parse(path, format="clustal")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments3.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="clustal")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="clustal")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignment), 20)
        self.assertEqual(alignment.get_alignment_length(), 414)
        self.assertEqual(alignment_summary(alignment), """\
  MMMMMMMMMMMMMMMM-M-- alignment column 0
  -----------------T-- alignment column 1
  -----------------V-- alignment column 2
  -----------------L-- alignment column 3
  -S---------------E-- alignment column 4
  |||||||||||||||||||| ...
  -T------------------ alignment column 413""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=20))), 3)
        handle.close()
    
    def test_reading_alignments_fasta(self):
        path = 'GFF/multi.fna'  # Trivial nucleotide alignment
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="fasta"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 3)
       
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="fasta"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
      
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="fasta")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
     
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="fasta")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="fasta")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="fasta")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignment), 3)
        self.assertEqual(alignment.get_alignment_length(), 8)
        self.assertEqual(alignment_summary(alignment), """\
  ACGTCGCG test1
  GGGGCCCC test2
  AAACACAC test3""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="fasta", seq_count=3))), 3)
        handle.close()
    
    def test_reading_alignments_nexus(self):
        path = 'Nexus/test_Nexus_input.nex'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="nexus"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 9)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="nexus"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="nexus")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="nexus")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="nexus")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="nexus")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignment), 9)
        self.assertEqual(alignment.get_alignment_length(), 48)
        self.assertEqual(alignment_summary(alignment), """\
  AAAAAAAAc alignment column 0
  -----c?tc alignment column 1
  CCCCCCCCc alignment column 2
  --c-?a-tc alignment column 3
  GGGGGGGGc alignment column 4
  ||||||||| ...
  tt--?ag?c alignment column 47""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        path = 'Nexus/codonposset.nex'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="nexus"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="nexus"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="nexus")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="nexus")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="nexus")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="nexus")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.get_alignment_length(), 22)
        self.assertEqual(alignment_summary(alignment), """\
  AAAAAGGCATTGTGGTGGGAAT Aegotheles
  ?????????TTGTGGTGGGAAT Aerodramus""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()

    def test_reading_alignments_stockholm(self):
        path = 'Stockholm/simple.sth'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="stockholm"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="stockholm"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments3 = []
        seq_iterator = AlignIO.parse(path, format="stockholm")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments3.append(record)

        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="stockholm")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="stockholm")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="stockholm")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.get_alignment_length(), 104)
        self.assertEqual(alignment_summary(alignment), """\
  UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCA...UGU AP001509.1
  AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGA...GAU AE007476.1
  .................<<<<<<<<...<<<<<<<....... secondary_structure""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="stockholm", seq_count=2))), 3)
        handle.close()
        path = 'Stockholm/funny.sth'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="stockholm"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 6)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="stockholm"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments3 = []
        seq_iterator = AlignIO.parse(path, format="stockholm")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments3.append(record)

        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="stockholm")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="stockholm")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="stockholm")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignment), 6)
        self.assertEqual(alignment.get_alignment_length(), 43)
        self.assertEqual(alignment_summary(alignment), """\
  MMMEEE alignment column 0
  TQIVVV alignment column 1
  CHEMMM alignment column 2
  RVALLL alignment column 3
  ASDTTT alignment column 4
  |||||| ...
  SYSEEE alignment column 42""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="stockholm", seq_count=6))), 3)
        handle.close()
        
    def test_reading_alignments_phylip(self):
        path = 'Phylip/reference_dna.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 6)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 6)
        self.assertEqual(alignments[0].get_alignment_length(), 13)
        self.assertEqual(alignment_summary(alignments[0]), """\
  CCTTCG alignment column 0
  GGAAAG alignment column 1
  ATAAAC alignment column 2
  TTTTAA alignment column 3
  GAGGAG alignment column 4
  |||||| ...
  CTTTTC alignment column 12""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=6))), 3)
        handle.close()
        path = 'Phylip/reference_dna2.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 6)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 6)
        self.assertEqual(alignments[0].get_alignment_length(), 39)
        self.assertEqual(alignment_summary(alignments[0]), """\
  CCTTCG alignment column 0
  GGAAAG alignment column 1
  ATAAAC alignment column 2
  TTTTAA alignment column 3
  GAGGAG alignment column 4
  |||||| ...
  CTTTTC alignment column 38""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=6))), 3)
        handle.close()
        path = 'Phylip/hennigian.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 10)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 10)
        self.assertEqual(alignments[0].get_alignment_length(), 40)
        self.assertEqual(alignment_summary(alignments[0]), """\
  CCCCCAAAAA alignment column 0
  AAAAACCCCC alignment column 1
  CCCAAAAAAA alignment column 2
  AAACCAAAAA alignment column 3
  CCAAAAAAAA alignment column 4
  |||||||||| ...
  AAAAAAAAAA alignment column 39""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=10))), 3)
        handle.close()
        path = 'Phylip/horses.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 10)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 10)
        self.assertEqual(alignments[0].get_alignment_length(), 40)
        self.assertEqual(alignment_summary(alignments[0]), """\
  AACCCCCCCC alignment column 0
  AAAACCCCCC alignment column 1
  AAAAAAAAAC alignment column 2
  ACAAAAAAAA alignment column 3
  ACACCCCCCC alignment column 4
  |||||||||| ...
  AAAAAAAAAA alignment column 39""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=10))), 3)
        handle.close()
        path = 'Phylip/random.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 10)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 10)
        self.assertEqual(alignments[0].get_alignment_length(), 40)
        self.assertEqual(alignment_summary(alignments[0]), """\
  CAAAACAAAC alignment column 0
  AACAACCACC alignment column 1
  CAAAACAAAA alignment column 2
  ACAACACACA alignment column 3
  CCAAAACCAA alignment column 4
  |||||||||| ...
  AAAAAAAAAA alignment column 39""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=10))), 3)
        handle.close()
        path = 'Phylip/interlaced.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 3)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 3)
        self.assertEqual(alignments[0].get_alignment_length(), 384)
        self.assertEqual(alignment_summary(alignments[0]), """\
  -----MKVILLFVLAVFTVFVSS-------------...I-- CYS1_DICDI
  MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPV...VAA ALEU_HORVU
  ------MWATLPLLCAGAWLLGV--------PVCGA...PLV CATH_HUMAN""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=3))), 3)
        handle.close()
        path = 'Phylip/interlaced2.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 4)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 4)
        self.assertEqual(alignments[0].get_alignment_length(), 131)
        self.assertEqual(alignment_summary(alignments[0]), """\
  TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTG...SHE IXI_234
  TSPASIRPPAGPSSR---------RPSPPGPRRPTG...SHE IXI_235
  TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPG...SHE IXI_236
  TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT-...SHE IXI_237""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=4))), 3)
        handle.close()
        path = 'ExtendedPhylip/primates.phyx'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip-relaxed"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 12)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip-relaxed"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip-relaxed")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip-relaxed")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip-relaxed")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip-relaxed")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 12)
        self.assertEqual(alignments[0].get_alignment_length(), 898)
        self.assertEqual(alignment_summary(alignments[0]), """\
  AAAAAAAAAAAA alignment column 0
  AAAAAAAAAAAA alignment column 1
  GGGGGGGGGGGG alignment column 2
  TCCCCCCCCCCC alignment column 3
  TTTTTTTTTTTT alignment column 4
  |||||||||||| ...
  TTTTTTTTTTTT alignment column 897""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip-relaxed", seq_count=12))), 3)
        handle.close()
        path = 'Phylip/sequential.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip-sequential"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 3)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip-sequential"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip-sequential")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip-sequential")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip-sequential")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip-sequential")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 3)
        self.assertEqual(alignments[0].get_alignment_length(), 384)
        self.assertEqual(alignment_summary(alignments[0]), """\
  -----MKVILLFVLAVFTVFVSS-------------...I-- CYS1_DICDI
  MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPV...VAA ALEU_HORVU
  ------MWATLPLLCAGAWLLGV--------PVCGA...PLV CATH_HUMAN""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip-sequential", seq_count=3))), 3)
        handle.close()
        path = 'Phylip/sequential2.phy'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="phylip-sequential"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 4)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="phylip-sequential"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="phylip-sequential")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip-sequential")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="phylip-sequential")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="phylip-sequential")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        # Show the alignment
        self.assertEqual(len(alignments[0]), 4)
        self.assertEqual(alignments[0].get_alignment_length(), 131)
        self.assertEqual(alignment_summary(alignments[0]), """\
  TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTG...SHE IXI_234
  TSPASIRPPAGPSSR---------RPSPPGPRRPTG...SHE IXI_235
  TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPG...SHE IXI_236
  TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT-...SHE IXI_237""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip-sequential", seq_count=4))), 3)
        handle.close()
        
    def test_reading_alignments_emboss(self):
        path = 'Emboss/alignret.txt'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="emboss"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 4)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="emboss"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="emboss")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="emboss")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        alignment = alignments[0]
        self.assertEqual(len(alignment), 4)
        self.assertEqual(alignment.get_alignment_length(), 131)
        self.assertEqual(alignment_summary(alignment), """\
  TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTG...SHE IXI_234
  TSPASIRPPAGPSSR---------RPSPPGPRRPTG...SHE IXI_235
  TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPG...SHE IXI_236
  TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT-...SHE IXI_237""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        path = 'Emboss/needle.txt'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="emboss"))
        self.assertEqual(len(alignments), 5)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="emboss"):
            alignments.append(record)
        self.assertEqual(len(alignments), 5)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="emboss")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 5)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 5)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 5)
    
        # Check Bio.AlignIO.read(...)
        self.assertRaises(ValueError, AlignIO.read, handle, "emboss")
        # Show the alignment
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 124)
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 119)
        self.assertEqual(len(alignments[2]), 2)
        self.assertEqual(alignments[2].get_alignment_length(), 120)
        self.assertEqual(len(alignments[3]), 2)
        self.assertEqual(alignments[3].get_alignment_length(), 118)
        self.assertEqual(len(alignments[4]), 2)
        self.assertEqual(alignments[4].get_alignment_length(), 125)
        self.assertEqual(alignment_summary(alignments[0]), """\
  KILIVDD----QYGIRILLNEVFNKEGYQTFQAANG...--- ref_rec
  -VLLADDHALVRRGFRLMLED--DPEIEIVAEAGDG...GET gi|94968718|receiver""")
        self.assertEqual(alignment_summary(alignments[1]), """\
  KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQAL...--- ref_rec
  -ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRAL...LKR gi|94968761|receiver""")
        self.assertEqual(alignment_summary(alignments[2]), """\
  -KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQA...--- ref_rec
  LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAA...HKE gi|94967506|receiver""")
        self.assertEqual(alignment_summary(alignments[3]), """\
  KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQAL...DAV ref_rec
  -VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEAL...EVL gi|94970045|receiver""")
        self.assertEqual(alignment_summary(alignments[4]), """\
  KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQAL...--- ref_rec
  TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEAL...KRQ gi|94970041|receiver""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)
        path = 'Emboss/needle_asis.txt'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="emboss"))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="emboss"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="emboss")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="emboss")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 3653)
        self.assertEqual(alignment_summary(alignment), """\
  TATTTTTTGGATTTTTTTCTAGATTTTCTAGGTTAT...GAA asis
  ------------------------------------...GAA asis""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        path = 'Emboss/water.txt'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="emboss"))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="emboss"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="emboss")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="emboss")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        # Show the alignment
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 131)
        self.assertEqual(alignment_summary(alignments[0]), """\
  TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTG...SHE IXI_234
  TSPASIRPPAGPSSR---------RPSPPGPRRPTG...SHE IXI_235""")
    
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        path = 'Emboss/water2.txt'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="emboss"))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="emboss"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="emboss")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="emboss")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 18)
        self.assertEqual(alignment_summary(alignments[0]), """\
  CGTTTGAGT-CTGGGATG asis
  CGTTTGAGTACTGGGATG asis""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        path = 'Emboss/matcher_simple.txt'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="emboss"))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="emboss"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="emboss")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="emboss")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.get_alignment_length(), 16)
        self.assertEqual(alignment_summary(alignment), """\
  GPPPQSPDENRAGESS AF069992_1
  GVPPEEAGAAVAAESS CAA85685.1""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        path ='Emboss/matcher_pair.txt'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="emboss"))
        self.assertEqual(len(alignments), 5)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="emboss"):
            alignments.append(record)
        self.assertEqual(len(alignments), 5)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="emboss")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 5)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 5)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, "emboss")

        self.assertEquals(len(alignments[0]), 2)
        self.assertEquals(alignments[0].get_alignment_length(), 145)
        self.assertEquals(len(alignments[1]), 2)
        self.assertEquals(alignments[1].get_alignment_length(), 13)
        self.assertEquals(len(alignments[2]), 2)
        self.assertEquals(alignments[2].get_alignment_length(), 18)
        self.assertEquals(len(alignments[3]), 2)
        self.assertEquals(alignments[3].get_alignment_length(), 10)
        self.assertEquals(len(alignments[4]), 2)
        self.assertEquals(alignments[4].get_alignment_length(), 10)
        self.assertEquals(alignment_summary(alignments[0]), """\
  LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFP...SKY HBA_HUMAN
  LTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYP...HKY HBB_HUMAN""")
        self.assertEquals(alignment_summary(alignments[1]), """\
  KKVADALTNAVAH HBA_HUMAN
  QKVVAGVANALAH HBB_HUMAN""")
        self.assertEquals(alignment_summary(alignments[2]), """\
  KLRVDPVNFKLLSHCLLV HBA_HUMAN
  KVNVDEVGGEALGRLLVV HBB_HUMAN""")
        self.assertEquals(alignment_summary(alignments[3]), """\
  LSALSDLHAH HBA_HUMAN
  LGAFSDGLAH HBB_HUMAN""")
        self.assertEquals(alignment_summary(alignments[4]), """\
  VKAAWGKVGA HBA_HUMAN
  VQAAYQKVVA HBB_HUMAN""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)
        path = 'Emboss/emboss_pair_aln_full_blank_line.txt'
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format="emboss"))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse(path, format="emboss"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse(path, format="emboss")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and for loop (a torture test!)
        seq_iterator = AlignIO.parse(path, format="emboss")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format="emboss")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignment.get_alignment_length(), 1450)
        self.assertEqual(alignment_summary(alignment), """\
  GGCAGGTGCATAGCTTGAGCCTAGGAGTTCAAGTCC...AAA hg38_chrX_131691529_131830643_47210_48660
  G--------------------------TTCAAGGCC...AAA mm10_chrX_50555743_50635321_27140_27743""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_fasta_m10(self):
        with open('Fasta/output001.m10', "r") as handle:
            alignments = list(AlignIO.parse(handle, format='fasta-m10'))
        self.assertEqual(len(alignments), 4)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
        alignments = []
        for record in AlignIO.parse('Fasta/output001.m10', format='fasta-m10'):
            alignments.append(record)
        self.assertEqual(len(alignments), 4)
        alignments = []
        seq_iterator = AlignIO.parse('Fasta/output001.m10', format='fasta-m10')
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('Fasta/output001.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 4)
        seq_iterator = AlignIO.parse('Fasta/output001.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 4)
        with open('Fasta/output001.m10') as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, 'fasta-m10')
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 108)
        self.assertEqual(alignment_summary(alignments[0]), """\
  SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECG...LSR gi|10955263|ref|NP_052604.1|
  AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFG...LSR gi|152973457|ref|YP_001338508.1|""")
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 64)
        self.assertEqual(alignment_summary(alignments[1]), """\
  AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MR...AIT gi|10955263|ref|NP_052604.1|
  ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIM...TLT gi|152973588|ref|YP_001338639.1|""")
        self.assertEqual(len(alignments[2]), 2)
        self.assertEqual(alignments[2].get_alignment_length(), 38)
        self.assertEqual(alignment_summary(alignments[2]), """\
  MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE gi|10955264|ref|NP_052605.1|
  IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ gi|152973462|ref|YP_001338513.1|""")
        self.assertEqual(len(alignments[3]), 2)
        self.assertEqual(alignments[3].get_alignment_length(), 43)
        self.assertEqual(alignment_summary(alignments[3]), """\
  SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEM...TYC gi|10955265|ref|NP_052606.1|
  SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLI...EYC gi|152973545|ref|YP_001338596.1|""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)

        with open('Fasta/output002.m10', "r") as handle:
            alignments = list(AlignIO.parse(handle, format='fasta-m10'))
        self.assertEqual(len(alignments), 6)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
        alignments = []
        for record in AlignIO.parse('Fasta/output002.m10', format='fasta-m10'):
            alignments.append(record)
        self.assertEqual(len(alignments), 6)
        alignments = []
        seq_iterator = AlignIO.parse('Fasta/output002.m10', format='fasta-m10')
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('Fasta/output002.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 6)
        seq_iterator = AlignIO.parse('Fasta/output002.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 6)
    
        with open('Fasta/output002.m10') as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, 'fasta-m10')
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 88)
        self.assertEqual(alignment_summary(alignments[0]), """\
  SGSNTRRRAISRPVR--LTAEEDQEIRKRAAECG-K...AEV gi|10955263|ref|NP_052604.1|
  SQRSTRRKPENQPTRVILFNKPYDVLPQFTDEAGRK...VQV gi|162139799|ref|NP_309634.2|""")
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 53)
        self.assertEqual(alignment_summary(alignments[1]), """\
  EIRKRAAECGKTVSGFLRAAA-LGKKV----NSLTD...KKL gi|10955263|ref|NP_052604.1|
  EIKPRGTSKGEAIAAFMQEAPFIGRTPVFLGDDLTD...VKI gi|15831859|ref|NP_310632.1|""")
        self.assertEqual(len(alignments[2]), 2)
        self.assertEqual(alignments[2].get_alignment_length(), 92)
        self.assertEqual(alignment_summary(alignments[2]), """\
  SEFFSKIESDLKKKKSKGDVFFDLIIPNG-----GK...ATS gi|10955264|ref|NP_052605.1|
  TELNSELAKAMKVDAQRG-AFVSQVLPNSSAAKAGI...QSS gi|15829419|ref|NP_308192.1|""")
        self.assertEqual(len(alignments[5]), 2)
        self.assertEqual(alignments[5].get_alignment_length(), 157)
        self.assertEqual(alignment_summary(alignments[5]), """\
  QYIMTTSNGDRVRAKIYKRGSIQFQGKYLQIASLIN...REI gi|10955265|ref|NP_052606.1|
  EFIRLLSDHDQFEKDQISELTVAANALKLEVAK--N...KKV gi|15833861|ref|NP_312634.1|""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)
        with open('Fasta/output003.m10', "r") as handle:
            alignments = list(AlignIO.parse(handle, format='fasta-m10'))
        self.assertEqual(len(alignments), 3)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
        alignments = []
        for record in AlignIO.parse('Fasta/output003.m10', format='fasta-m10'):
            alignments.append(record)
        self.assertEqual(len(alignments), 3)
        alignments = []
        seq_iterator = AlignIO.parse('Fasta/output003.m10', format='fasta-m10')
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('Fasta/output003.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 3)
        seq_iterator = AlignIO.parse('Fasta/output003.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 3)
        with open('Fasta/output003.m10') as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, 'fasta-m10')
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 55)
        self.assertEqual(alignment_summary(alignments[0]), """\
  ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGG...IAA gi|152973837|ref|YP_001338874.1|
  VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGK...LGA gi|10955263|ref|NP_052604.1|""")
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 22)
        self.assertEqual(alignment_summary(alignments[1]), """\
  DDAEHLFRTLSSR-LDALQDGN gi|152973840|ref|YP_001338877.1|
  DDRANLFEFLSEEGITITEDNN gi|10955265|ref|NP_052606.1|""")
        self.assertEqual(len(alignments[2]), 2)
        self.assertEqual(alignments[2].get_alignment_length(), 63)
        self.assertEqual(alignment_summary(alignments[2]), """\
  VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIR...QAM gi|152973841|ref|YP_001338878.1|
  VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSAS...KGI gi|10955264|ref|NP_052605.1|""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)
        with open('Fasta/output004.m10', "r") as handle:
            alignments = list(AlignIO.parse(handle, format='fasta-m10'))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
        alignments = []
        for record in AlignIO.parse('Fasta/output004.m10', format='fasta-m10'):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        alignments = []
        seq_iterator = AlignIO.parse('Fasta/output004.m10', format='fasta-m10')
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('Fasta/output004.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
        seq_iterator = AlignIO.parse('Fasta/output004.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        with open('Fasta/output004.m10') as handle:
            alignment = AlignIO.read(handle, format='fasta-m10')
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 102)
        self.assertEqual(alignment_summary(alignments[0]), """\
  AAAAAAGATAAAAAATATCAAATAGAAGCAATAAAA...TCA ref|NC_002127.1|:c1351-971
  AGAGAAAATAAAACAAGTAATAAAATATTAATGGAA...ACA ref|NC_002695.1|:1970775-1971404""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)
        with open('Fasta/output005.m10', "r") as handle:
            alignments = list(AlignIO.parse(handle, format='fasta-m10'))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
        alignments = []
        for record in AlignIO.parse('Fasta/output005.m10', format='fasta-m10'):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        alignments = []
        seq_iterator = AlignIO.parse('Fasta/output005.m10', format='fasta-m10')
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('Fasta/output005.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
        seq_iterator = AlignIO.parse('Fasta/output005.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        with open('Fasta/output005.m10') as handle:
            alignment = AlignIO.read(handle, format='fasta-m10')
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 110)
        self.assertEqual(alignment_summary(alignments[0]), """\
  IKNKDKTLFIVYAT-DIYSPSEFFSKIESDLKKKKS...LSK gi|10955264|ref|NP_052605.1|
  IKDELPVAFCSWASLDLECEVKYINDVTSLYAKDWM...MSE gi|10955282|ref|NP_052623.1|""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)
        with open('Fasta/output006.m10', "r") as handle:
            alignments = list(AlignIO.parse(handle, format='fasta-m10'))
        self.assertEqual(len(alignments), 1)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
        alignments = []
        for record in AlignIO.parse('Fasta/output006.m10', format='fasta-m10'):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        alignments = []
        seq_iterator = AlignIO.parse('Fasta/output006.m10', format='fasta-m10')
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('Fasta/output006.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
        seq_iterator = AlignIO.parse('Fasta/output006.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        with open('Fasta/output006.m10') as handle:
            alignment = AlignIO.read(handle, format='fasta-m10')
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 131)
        self.assertEqual(alignment_summary(alignments[0]), """\
  GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACA...CAT query
  GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACA...CAT gi|116660610|gb|EG558221.1|EG558221""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)

        with open('Fasta/output007.m10', "r") as handle:
            alignments = list(AlignIO.parse(handle, format='fasta-m10'))
        self.assertEqual(len(alignments), 9)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
        alignments = []
        for record in AlignIO.parse('Fasta/output007.m10', format='fasta-m10'):
            alignments.append(record)
        self.assertEqual(len(alignments), 9)
        alignments = []
        seq_iterator = AlignIO.parse('Fasta/output007.m10', format='fasta-m10')
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('Fasta/output007.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 9)
        seq_iterator = AlignIO.parse('Fasta/output007.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 9)
        with open('Fasta/output007.m10') as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, 'fasta-m10')
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 108)
        self.assertEqual(alignment_summary(alignments[0]), """\
  SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECG...LSR gi|10955263|ref|NP_052604.1|
  AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFG...LSR gi|152973457|ref|YP_001338508.1|""")
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 64)
        self.assertEqual(alignment_summary(alignments[1]), """\
  AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MR...AIT gi|10955263|ref|NP_052604.1|
  ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIM...TLT gi|152973588|ref|YP_001338639.1|""")
        self.assertEqual(len(alignments[2]), 2)
        self.assertEqual(alignments[2].get_alignment_length(), 45)
        self.assertEqual(alignment_summary(alignments[2]), """\
  EIRKRAAECGKTVSGFLRAAA-----LGKKVNSLTD...VMR gi|10955263|ref|NP_052604.1|
  ELVKLIADMGISVRALLRKNVEPYEELGLEEDKFTD...MLQ gi|152973480|ref|YP_001338531.1|""")
        self.assertEqual(len(alignments[8]), 2)
        self.assertEqual(alignments[8].get_alignment_length(), 64)
        self.assertEqual(alignment_summary(alignments[8]), """\
  ISGTYKGIDFLIKLMPSGGNTTIGRASGQNNTYFDE...FSD gi|10955265|ref|NP_052606.1|
  IDGVITAFD-LRTGMNISKDKVVAQIQGMDPVW---...YPD gi|152973505|ref|YP_001338556.1|""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)
        with open('Fasta/output008.m10', "r") as handle:
            alignments = list(AlignIO.parse(handle, format='fasta-m10'))
        self.assertEqual(len(alignments), 12)
        for alignment in alignments:
            self.assertEqual(len(alignment), 2)
        alignments = []
        for record in AlignIO.parse('Fasta/output008.m10', format='fasta-m10'):
            alignments.append(record)
        self.assertEqual(len(alignments), 12)
        alignments = []
        seq_iterator = AlignIO.parse('Fasta/output008.m10', format='fasta-m10')
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('Fasta/output008.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 12)
        seq_iterator = AlignIO.parse('Fasta/output008.m10', format='fasta-m10')
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 12)
        with open('Fasta/output008.m10') as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, 'fasta-m10')
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 65)
        self.assertEqual(alignment_summary(alignments[0]), """\
  LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH-...QML sp|Q9NSY1|BMP2K_HUMAN
  IPHQLPHALRHRPAQEAAHASQLHPAQPGCGQPLHG...GLL gi|283855822|gb|GQ290312.1|""")
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 201)
        self.assertEqual(alignment_summary(alignments[1]), """\
  GPEIL---LGQ-GPPQQPPQQHRVLQQLQQGDWRLQ...NRS sp|Q9NSY1|BMP2K_HUMAN
  GPELLRALLQQNGCGTQPLRVPTVLPG*AMAVLHAG...QKS gi|57163782|ref|NM_001009242.1|""")
        self.assertEqual(len(alignments[2]), 2)
        self.assertEqual(alignments[2].get_alignment_length(), 348)
        self.assertEqual(alignment_summary(alignments[2]), """\
  MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQ...APA sp|P08100|OPSD_HUMAN
  MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQ...APA gi|57163782|ref|NM_001009242.1|""")
        self.assertEqual(len(alignments[11]), 2)
        self.assertEqual(alignments[11].get_alignment_length(), 31)
        self.assertEqual(alignment_summary(alignments[11]), """\
  AQQQESATTQKAEKEVTRMVIIMVIAFLICW sp|P08100|OPSD_HUMAN
  SQQIRNATTMMMTMRVTSFSAFWVVADSCCW gi|283855822|gb|GQ290312.1|""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        check_simple_write_read(alignments)

    def test_reading_alignments_ig(self):
        with open('IntelliGenetics/VIF_mase-pro.txt', "r") as handle:
            alignments = list(AlignIO.parse(handle, format="ig"))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 16)
        alignments = []
        for record in AlignIO.parse('IntelliGenetics/VIF_mase-pro.txt', format="ig"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        alignments = []
        seq_iterator = AlignIO.parse('IntelliGenetics/VIF_mase-pro.txt', format="ig")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('IntelliGenetics/VIF_mase-pro.txt', format="ig")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
    
        seq_iterator = AlignIO.parse('IntelliGenetics/VIF_mase-pro.txt', format="ig")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        with open('IntelliGenetics/VIF_mase-pro.txt') as handle:
            alignment = AlignIO.read(handle, format="ig")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
    

        self.assertEqual(len(alignments[0]), 16)
        self.assertEqual(alignments[0].get_alignment_length(), 298)
        self.assertEqual(alignment_summary(alignments[0]), """\
  MMMMMMMMMMMMMMMM alignment column 0
  EEEEEEETEEEENEEE alignment column 1
  NNNNNNNAEEEEQRKK alignment column 2
  --------DEEEEE-- alignment column 3
  --------KKKKKK-- alignment column 4
  |||||||||||||||| ...
  HHHHHHH-AAAAL-R- alignment column 297""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open('IntelliGenetics/VIF_mase-pro.txt', "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="ig", seq_count=16))), 3)
        handle.close()
    
    def test_reading_alignments_pir(self):
        with open('NBRF/clustalw.pir', "r") as handle:
            alignments = list(AlignIO.parse(handle, format="pir"))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
    
        # Try using the iterator with a for loop and a filename not handle
        alignments = []
        for record in AlignIO.parse('NBRF/clustalw.pir', format="pir"):
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try using the iterator with the next() method
        alignments = []
        seq_iterator = AlignIO.parse('NBRF/clustalw.pir', format="pir")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
    
        # Try a mixture of next() and list (a torture test!)
        seq_iterator = AlignIO.parse('NBRF/clustalw.pir', format="pir")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 1)
        seq_iterator = AlignIO.parse('NBRF/clustalw.pir', format="pir")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 1)
        with open('NBRF/clustalw.pir') as handle:
            alignment = AlignIO.read(handle, format="pir")
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 2527)
        self.assertEqual(alignment_summary(alignments[0]), """\
  ------------------------------------...--- 804Angiostrongylus_cantonensis
  ------------------------------------...--- 815Parelaphostrongylus_odocoil""")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        try:
            info_content = summary.information_content()
        except ValueError as err:
            if str(err) != "Error in alphabet: not Nucleotide or Protein, supply expected frequencies":
                raise err

        with open('NBRF/clustalw.pir', "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="pir", seq_count=2))), 3)
        handle.close()
    
    def test_reading_alignments_maf(self):
        with open('MAF/humor.maf', "r") as handle:
            alignments = list(AlignIO.parse(handle, format="maf"))
        self.assertEqual(len(alignments), 2)
        alignments = []
        for record in AlignIO.parse('MAF/humor.maf', format="maf"):
            alignments.append(record)
        self.assertEqual(len(alignments), 2)
        alignments = []
        seq_iterator = AlignIO.parse('MAF/humor.maf', format="maf")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse('MAF/humor.maf', format="maf")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 2)
        seq_iterator = AlignIO.parse('MAF/humor.maf', format="maf")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 2)
        with open('MAF/humor.maf') as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, "maf")
        self.assertEqual(len(alignments[0]), 3)
        self.assertEqual(alignments[0].get_alignment_length(), 5486)
        self.assertEqual(alignment_summary(alignments[0]), """\
  gcacagcctttactccctgactgcgtttatattctg...CCG NM_006987
  gcacagcctttactccctgactgcgtttatattctg...TTG mm3
  gcacagcctttactccctgactgcgtttatattctg...CCG rn3""")
        self.assertEqual(len(alignments[1]), 3)
        self.assertEqual(alignments[1].get_alignment_length(), 5753)
        self.assertEqual(alignment_summary(alignments[1]), """\
  tttgtccatgttggtcaggctggtctcgaactcccc...GGT NM_018289
  tttgtccatgttggtcaggctggtctcgaactcccc...GGT mm3
  tttgtccatgttggtcaggctggtctcgaactcccc...GGT rn3""")
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignments[1])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        alignments.reverse()
        check_simple_write_read(alignments)
        with open("MAF/bug2453.maf", "r") as handle:
            alignments = list(AlignIO.parse(handle, format="maf"))
        self.assertEqual(len(alignments), 3)
        alignments = []
        for record in AlignIO.parse("MAF/bug2453.maf", format="maf"):
            alignments.append(record)
        self.assertEqual(len(alignments), 3)
        alignments = []
        seq_iterator = AlignIO.parse("MAF/bug2453.maf", format="maf")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse("MAF/bug2453.maf", format="maf")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 3)
        seq_iterator = AlignIO.parse("MAF/bug2453.maf", format="maf")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 3)
        with open("MAF/bug2453.maf") as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, "maf")
        self.assertEqual(len(alignments[0]), 5)
        self.assertEqual(alignments[0].get_alignment_length(), 42)
        self.assertEqual(alignment_summary(alignments[0]), """\
  AAA-- alignment column 0
  AAAAA alignment column 1
  AAAAA alignment column 2
  ---T- alignment column 3
  GGGGG alignment column 4
  ||||| ...
  GGGGG alignment column 41""")
        self.assertEqual(len(alignments[1]), 5)
        self.assertEqual(alignments[1].get_alignment_length(), 6)
        self.assertEqual(alignment_summary(alignments[1]), """\
  TTTTt alignment column 0
  AAAAa alignment column 1
  AAAAa alignment column 2
  AAAAg alignment column 3
  GGGGg alignment column 4
  ||||| ...
  AAAAa alignment column 5""")
        self.assertEqual(len(alignments[2]), 4)
        self.assertEqual(alignments[2].get_alignment_length(), 13)
        self.assertEqual(alignment_summary(alignments[2]), """\
  gcagctgaaaaca hg16.chr7
  gcagctgaaaaca panTro1.chr6
  gcagctgaaaaca baboon
  ACAGCTGAAAATA mm4.chr6""")
        summary = AlignInfo.SummaryInfo(alignments[1])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        alignments.reverse()
        check_simple_write_read(alignments)
        with open("MAF/ucsc_test.maf", "r") as handle:
            alignments = list(AlignIO.parse(handle, format="maf"))
        self.assertEqual(len(alignments), 3)
        alignments = []
        for record in AlignIO.parse("MAF/ucsc_test.maf", format="maf"):
            alignments.append(record)
        self.assertEqual(len(alignments), 3)
        alignments = []
        seq_iterator = AlignIO.parse("MAF/ucsc_test.maf", format="maf")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse("MAF/ucsc_test.maf", format="maf")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 3)
        seq_iterator = AlignIO.parse("MAF/ucsc_test.maf", format="maf")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 3)
        with open("MAF/ucsc_test.maf") as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, "maf")
        self.assertEqual(len(alignments[0]), 5)
        self.assertEqual(alignments[0].get_alignment_length(), 42)
        self.assertEqual(alignment_summary(alignments[0]), """\
  AAA-- alignment column 0
  AAAAA alignment column 1
  AAAAA alignment column 2
  ---T- alignment column 3
  GGGGG alignment column 4
  ||||| ...
  GGGGG alignment column 41""")
        self.assertEqual(len(alignments[1]), 5)
        self.assertEqual(alignments[1].get_alignment_length(), 6)
        self.assertEqual(alignment_summary(alignments[1]), """\
  TTTTt alignment column 0
  AAAAa alignment column 1
  AAAAa alignment column 2
  AAAAg alignment column 3
  GGGGg alignment column 4
  ||||| ...
  AAAAa alignment column 5""")
        self.assertEqual(len(alignments[2]), 4)
        self.assertEqual(alignments[2].get_alignment_length(), 13)
        self.assertEqual(alignment_summary(alignments[2]), """\
  gcagctgaaaaca hg16.chr7
  gcagctgaaaaca panTro1.chr6
  gcagctgaaaaca baboon
  ACAGCTGAAAATA mm4.chr6""")
        summary = AlignInfo.SummaryInfo(alignments[2])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        try:
            info_content = summary.information_content()
        except ValueError as err:
            if str(err) != "Error in alphabet: not Nucleotide or Protein, supply expected frequencies":
                raise err
        alignments.reverse()
        check_simple_write_read(alignments)
        with open("MAF/ucsc_mm9_chr10.maf", "r") as handle:
            alignments = list(AlignIO.parse(handle, format="maf"))
        self.assertEqual(len(alignments), 48)
        alignments = []
        for record in AlignIO.parse("MAF/ucsc_mm9_chr10.maf", format="maf"):
            alignments.append(record)
        self.assertEqual(len(alignments), 48)
        alignments = []
        seq_iterator = AlignIO.parse("MAF/ucsc_mm9_chr10.maf", format="maf")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        seq_iterator = AlignIO.parse("MAF/ucsc_mm9_chr10.maf", format="maf")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 48)
        seq_iterator = AlignIO.parse("MAF/ucsc_mm9_chr10.maf", format="maf")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 48)
        with open("MAF/ucsc_mm9_chr10.maf") as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, "maf")
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 164)
        self.assertEqual(alignment_summary(alignments[0]), """\
  TCATAGGTATTTATTTTTAAATATGGTTTGCTTTAT...GTT mm9.chr10
  TCACAGATATTTACTATTAAATATGGTTTGTTATAT...GTT oryCun1.scaffold_133159""")
        self.assertEqual(len(alignments[1]), 4)
        self.assertEqual(alignments[1].get_alignment_length(), 466)
        self.assertEqual(alignment_summary(alignments[1]), """\
  AGTCTTTCCAATGGGACCTGTGAGTCCTAACTATGC...CTG mm9.chr10
  AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTC...TTC ponAbe2.chr6
  AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTC...TTC panTro2.chr6
  AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTC...TTC hg18.chr6""")
        self.assertEqual(len(alignments[2]), 5)
        self.assertEqual(alignments[2].get_alignment_length(), 127)
        self.assertEqual(alignment_summary(alignments[2]), """\
  TTTTT alignment column 0
  GGGGG alignment column 1
  GGGGG alignment column 2
  GGGGG alignment column 3
  TTTTC alignment column 4
  ||||| ...
  CCCCC alignment column 126""")
        self.assertEqual(len(alignments[47]), 6)
        self.assertEqual(alignments[47].get_alignment_length(), 46)
        self.assertEqual(alignment_summary(alignments[47]), """\
  TTTTTT alignment column 0
  GGGGGG alignment column 1
  TTTTTT alignment column 2
  TTTTTT alignment column 3
  TGGGAT alignment column 4
  |||||| ...
  tTTTT- alignment column 45""")
        summary = AlignInfo.SummaryInfo(alignments[47])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        alignments.reverse()
        check_simple_write_read(alignments)
 
    def test_reading_alignments_mauve(self):
        with open('Mauve/simple.xmfa', "r") as handle:
            alignments = list(AlignIO.parse(handle, format="mauve"))
        self.assertEqual(len(alignments), 5)
        alignments = []
        for record in AlignIO.parse('Mauve/simple.xmfa', format="mauve"):
            alignments.append(record)
        self.assertEqual(len(alignments), 5)
        alignments = []
        seq_iterator = AlignIO.parse('Mauve/simple.xmfa', format="mauve")
        while True:
            try:
                record = next(seq_iterator)
            except StopIteration:
                break
            self.assertIsNotNone(record)
            alignments.append(record)
        self.assertEqual(len(alignments), 5)
        seq_iterator = AlignIO.parse('Mauve/simple.xmfa', format="mauve")
        record = next(seq_iterator)
        alignments = [record]
        alignments.extend(list(seq_iterator))
        self.assertEqual(len(alignments), 5)
        seq_iterator = AlignIO.parse('Mauve/simple.xmfa', format="mauve")
        record = next(seq_iterator)
        alignments = [record]
        for record in seq_iterator:
            alignments.append(record)
        self.assertEqual(len(alignments), 5)
        with open('Mauve/simple.xmfa') as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, "mauve")
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 5670)
        self.assertEqual(alignment_summary(alignments[0]), """\
  ATATTAGGTTTTTACCTACCCAGGAAAAGCCAACCA...AAT 1/0-5670
  ATATTAGGTTTTTACCTACCCAGGAAAAGCCAACCA...AAT 2/0-5670""")
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 4420)
        self.assertEqual(alignment_summary(alignments[1]), """\
  GAACATCAGCACCTGAGTTGCTAAAGTCATTTAGAG...CTC 1/5670-9940
  GAACATCAGCACCTGAGTTGCTAAAGTCATTTAGAG...CTC 2/7140-11410""")
        self.assertEqual(len(alignments[2]), 1)
        self.assertEqual(alignments[2].get_alignment_length(), 4970)
        self.assertEqual(alignment_summary(alignments[2]), """\
  TCTACCAACCACCACAGACATCAATCACTTCTGCTG...GAC 1/9940-14910""")
        self.assertEqual(len(alignments[3]), 1)
        self.assertEqual(alignments[3].get_alignment_length(), 1470)
        self.assertEqual(len(alignments[4]), 1)
        self.assertEqual(alignments[4].get_alignment_length(), 1470)
        self.assertEqual(alignment_summary(alignments[4]), """\
  ATTCGCACATAAGAATGTACCTTGCTGTAATTTATA...ATA 2/11410-12880""")

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignments[4])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)
    
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
