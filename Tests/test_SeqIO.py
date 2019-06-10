# Copyright 2007-2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SeqIO module."""

from Bio._py3k import basestring

import unittest
import warnings


try:
    from StringIO import StringIO  # Python 2
    # Can't use cStringIO, quoting the documentation,
    #   "Unlike the StringIO module, this module is not able to accept
    #    Unicode strings that cannot be encoded as plain ASCII strings."
    # Therefore can't use from Bio._py3k import StringIO
except ImportError:
    from io import StringIO  # Python 3
from io import BytesIO

from Bio import BiopythonWarning, BiopythonParserWarning
from Bio import SeqIO
from Bio import AlignIO
from Bio.AlignIO import PhylipIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UnknownSeq
from Bio import Alphabet
from Bio.Align import MultipleSeqAlignment


# TODO - Check that desired warnings are issued. Used to do that by capturing
# warnings to stdout and verifying via the print-and-compare check. However,
# there was some frustrating cross-platform inconsistency I couldn't resolve.

protein_alphas = [Alphabet.generic_protein]
dna_alphas = [Alphabet.generic_dna]
rna_alphas = [Alphabet.generic_rna]
nucleotide_alphas = [Alphabet.generic_nucleotide,
                     Alphabet.Gapped(Alphabet.generic_nucleotide)]
no_alpha_formats = {
    'clustal',
    'emboss',
    'fasta', 'fasta-2line',
    'fastq', 'fastq-illumina', 'fastq-solexa',
    'ig',
    'phylip', 'phylip-relaxed', 'phylip-sequential',
    'qual',
    'stockholm',
    'tab',
}
possible_unknown_seq_formats = {
    "embl",
    "genbank", "gb",
    "imgt",
    "qual",
}

# List of formats including alignment only file formats we can read AND write.
# The list is initially hard coded to preserve the original order of the unit
# test output, with any new formats added since appended to the end.
test_write_read_alignment_formats = ["fasta", "clustal", "phylip", "stockholm",
                                     "phylip-relaxed"]
for format in sorted(SeqIO._FormatToWriter):
    if format not in test_write_read_alignment_formats:
        test_write_read_alignment_formats.append(format)
for format in sorted(AlignIO._FormatToWriter):
    if format not in test_write_read_alignment_formats:
        test_write_read_alignment_formats.append(format)
test_write_read_alignment_formats.remove("gb")  # an alias for genbank
test_write_read_alignment_formats.remove("fastq-sanger")  # an alias for fastq


class ForwardOnlyHandle(object):
    """Mimic a network handle without seek and tell methods etc."""

    def __init__(self, handle):
        """Initialize."""
        self._handle = handle

    def __iter__(self):
        """Iterate."""
        return iter(self._handle)

    def read(self, length=None):
        if length is None:
            return self._handle.read()
        else:
            return self._handle.read(length)

    def readline(self):
        return self._handle.readline()

    def readlines(self):
        return self._handle.readlines()

    def close(self):
        return self._handle.close()


def col_summary(col_text):
    if len(col_text) < 65:
        return col_text
    else:
        return col_text[:60] + "..." + col_text[-5:]


def alignment_summary(alignment, index=" "):
    """Return a concise summary of an Alignment object as a string."""
    answer = []
    alignment_len = alignment.get_alignment_length()
    rec_count = len(alignment)
    for i in range(min(5, alignment_len)):
        answer.append(index + col_summary(alignment[:, i]) +
                      " alignment column %i" % i)
    if alignment_len > 5:
        i = alignment_len - 1
        answer.append(index + col_summary("|" * rec_count) +
                      " ...")
        answer.append(index + col_summary(alignment[:, i]) +
                      " alignment column %i" % i)
    return "\n".join(answer)


class TestSeqIO(unittest.TestCase):

    def setUp(self):
        self.addTypeEqualityFunc(SeqRecord, self.compare_record)

    def compare_record(self, record_one, record_two, msg=None):
        """Attempt strict SeqRecord comparison."""
        if not isinstance(record_one, SeqRecord):
            self.failureException(msg)
        if not isinstance(record_two, SeqRecord):
            self.failureException(msg)
        if record_one.seq is None:
            self.failureException(msg)
        if record_two.seq is None:
            self.failureException(msg)
        if record_one.id != record_two.id:
            self.failureException(msg)
        if record_one.name != record_two.name:
            self.failureException(msg)
        if record_one.description != record_two.description:
            self.failureException(msg)
        if len(record_one) != len(record_two):
            self.failureException(msg)
        if isinstance(record_one.seq, UnknownSeq) \
                and isinstance(record_two.seq, UnknownSeq):
            # Jython didn't like us comparing the string of very long UnknownSeq
            # object (out of heap memory error)
            if record_one.seq._character != record_two.seq._character:
                self.failureException(msg)
        elif str(record_one.seq) != str(record_two.seq):
            self.failureException(msg)
        # TODO - check features and annotation (see code for BioSQL tests)
        for key in set(record_one.letter_annotations).intersection(
                record_two.letter_annotations):
            if record_one.letter_annotations[key] != \
               record_two.letter_annotations[key]:
                self.failureException(msg)

    def check_simple_write_read(self, records, t_format, t_count, messages):
        messages = iter(messages)
        for format in test_write_read_alignment_formats:
            if format not in possible_unknown_seq_formats \
                    and isinstance(records[0].seq, UnknownSeq) \
                    and len(records[0].seq) > 100:
                # Skipping for speed.  Some of the unknown sequences are
                # rather long, and it seems a bit pointless to record them.
                continue
            # Going to write to a handle...
            if format in SeqIO._BinaryFormats:
                handle = BytesIO()
            else:
                handle = StringIO()

            try:
                with warnings.catch_warnings():
                    # e.g. data loss
                    warnings.simplefilter("ignore", BiopythonWarning)
                    c = SeqIO.write(
                        sequences=records, handle=handle, format=format)
                self.assertEqual(c, len(records))
            except (ValueError, TypeError) as e:
                message = next(messages)
                self.assertEqual(str(e), message)
                if records[0].seq.alphabet.letters is not None:
                    self.assertNotEqual(format, t_format,
                                        "Should be able to re-write "
                                        "in the original format!")
                # Carry on to the next format:
                continue

            handle.flush()
            handle.seek(0)
            # Now ready to read back from the handle...
            try:
                records2 = list(SeqIO.parse(handle=handle, format=format))
            except ValueError as e:
                # This is BAD.  We can't read our own output.
                # I want to see the output when called from the test harness,
                # run_tests.py (which can be funny about new lines on Windows)
                handle.seek(0)
                raise ValueError("%s\n\n%s\n\n%s"
                                 % (str(e), repr(handle.read()), repr(records)))

            self.assertEqual(len(records2), t_count)
            for r1, r2 in zip(records, records2):
                # Check the bare minimum (ID and sequence) as
                # many formats can't store more than that.
                self.assertEqual(len(r1), len(r2))

                # Check the sequence
                if format in ["gb", "genbank", "embl", "imgt"]:
                    # The GenBank/EMBL parsers will convert to upper case.
                    if isinstance(r1.seq, UnknownSeq) \
                            and isinstance(r2.seq, UnknownSeq):
                        # Jython didn't like us comparing the string of very long
                        # UnknownSeq object (out of heap memory error)
                        self.assertEqual(r1.seq._character.upper(),
                                         r2.seq._character)
                    else:
                        self.assertEqual(str(r1.seq).upper(), str(r2.seq))
                elif format == "qual":
                    self.assertIsInstance(r2.seq, UnknownSeq)
                    self.assertEqual(len(r2), len(r1))
                elif format == "nib":
                    self.assertEqual(str(r1.seq).upper(), str(r2.seq))
                else:
                    self.assertEqual(str(r1.seq), str(r2.seq))
                # Beware of different quirks and limitations in the
                # valid character sets and the identifier lengths!
                if format in ["phylip", "phylip-sequential"]:
                    self.assertEqual(PhylipIO.sanitize_name(r1.id, 10), r2.id,
                                     "'%s' vs '%s'" % (r1.id, r2.id))
                elif format == "phylip-relaxed":
                    self.assertEqual(PhylipIO.sanitize_name(r1.id), r2.id,
                                     "'%s' vs '%s'" % (r1.id, r2.id))
                elif format == "clustal":
                    self.assertEqual(r1.id.replace(" ", "_")[:30], r2.id,
                                     "'%s' vs '%s'" % (r1.id, r2.id))
                elif format == "stockholm":
                    r1_id = r1.id.replace(" ", "_")
                    if "start" in r1.annotations and "end" in r1.annotations:
                        suffix = "/%d-%d" % (r1.annotations["start"],
                                             r1.annotations["end"])
                        if not r1_id.endswith(suffix):
                            r1_id += suffix

                    self.assertEqual(r1_id, r2.id,
                                     "'%s' vs '%s'" % (r1.id, r2.id))
                elif format == "maf":
                    self.assertEqual(r1.id.replace(" ", "_"), r2.id,
                                     "'%s' vs '%s'" % (r1.id, r2.id))
                elif format in ["fasta", "fasta-2line"]:
                    self.assertEqual(r1.id.split()[0], r2.id)
                elif format == 'nib':
                    self.assertEqual(r2.id, '<unknown id>')
                else:
                    self.assertEqual(r1.id, r2.id,
                                     "'%s' vs '%s'" % (r1.id, r2.id))

            if len(records) > 1:
                # Try writing just one record (passing a SeqRecord, not a list)
                if format in SeqIO._BinaryFormats:
                    handle = BytesIO()
                else:
                    handle = StringIO()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", BiopythonWarning)
                    SeqIO.write(records[0], handle, format)
                    self.assertEqual(handle.getvalue(), records[0].format(format))

    def perform_test(self, t_format, t_alignment, t_filename, t_count,
                     expected_records, expected_lengths, expected_alignment,
                     expected_messages):
        if t_format in SeqIO._BinaryFormats:
            mode = "rb"
        else:
            mode = "r"

        with warnings.catch_warnings():
            # e.g. BiopythonParserWarning: Dropping bond qualifier in feature
            # location
            # e.g. First line is not a 'HEADER'; can't determine PDB ID
            warnings.simplefilter("ignore", BiopythonParserWarning)

            # e.g. WARNING: Chain C is discontinuous at line 2645
            warnings.simplefilter("ignore", PDBConstructionWarning)

            # Try as an iterator using handle
            h = open(t_filename, mode)
            records = list(SeqIO.parse(handle=h, format=t_format))
            h.close()
            self.assertEqual(len(records), t_count,
                             "Found %i records but expected %i"
                             % (len(records), t_count))

            # Try using the iterator with a for loop, and a filename not handle
            records2 = []
            for record in SeqIO.parse(t_filename, format=t_format):
                records2.append(record)
            self.assertEqual(len(records2), t_count)

            # Try using the iterator with the next() method
            records3 = []
            h = open(t_filename, mode)
            seq_iterator = SeqIO.parse(handle=h, format=t_format)
            while True:
                try:
                    record = next(seq_iterator)
                except StopIteration:
                    break
                self.assertIsNotNone(record,
                                     "Should raise StopIteration, "
                                     "not return None")
                records3.append(record)
            h.close()
            self.assertEqual(len(records3), t_count)

            # Try a mixture of next() and list (a torture test!)
            h = open(t_filename, mode)
            seq_iterator = SeqIO.parse(handle=h, format=t_format)
            try:
                record = next(seq_iterator)
            except StopIteration:
                record = None
            if record is not None:
                records4 = [record]
                records4.extend(list(seq_iterator))
            else:
                records4 = []
            h.close()
            self.assertEqual(len(records4), t_count)

            # Try a mixture of next() and for loop (a torture test!)
            # with a forward-only-handle
            if t_format == "abi":
                # Temp hack
                h = open(t_filename, mode)
            else:
                h = ForwardOnlyHandle(open(t_filename, mode))
            seq_iterator = SeqIO.parse(h, format=t_format)
            try:
                record = next(seq_iterator)
            except StopIteration:
                record = None
            if record is not None:
                records5 = [record]
                for record in seq_iterator:
                    records5.append(record)
            else:
                records5 = []
            h.close()
            self.assertEqual(len(records5), t_count)

            for i in range(t_count):
                record = records[i]

                # Check returned expected object type
                self.assertIsInstance(record, SeqRecord)
                if t_format in possible_unknown_seq_formats:
                    if (not isinstance(record.seq, Seq)) and \
                       (not isinstance(record.seq, UnknownSeq)):
                        self.failureException("Expected a Seq or UnknownSeq object")
                else:
                    self.assertIsInstance(record.seq, Seq)
                self.assertIsInstance(record.id, basestring)
                self.assertIsInstance(record.name, basestring)
                self.assertIsInstance(record.description, basestring)
                self.assertTrue(record.id)

                if "accessions" in record.annotations:
                    accs = record.annotations["accessions"]
                    # Check for blanks, or entries with leading/trailing spaces
                    for acc in accs:
                        self.assertTrue(acc,
                                        "Bad accession in annotations: %s"
                                        % repr(acc))
                        self.assertEqual(acc, acc.strip(),
                                         "Bad accession in annotations: %s"
                                         % repr(acc))
                    self.assertEqual(len(set(accs)), len(accs),
                                     "Repeated accession in annotations: %s"
                                     % repr(accs))
                for ref in record.dbxrefs:
                    self.assertTrue(ref,
                                    "Bad cross reference in dbxrefs: %s"
                                    % repr(ref))
                    self.assertEqual(ref, ref.strip(),
                                     "Bad cross reference in dbxrefs: %s"
                                     % repr(ref))
                self.assertEqual(len(set(record.dbxrefs)), len(record.dbxrefs),
                                 "Repeated cross reference in dbxrefs: %s"
                                 % repr(record.dbxrefs))

                # Check the lists obtained by the different methods agree
                self.assertEqual(record, records2[i])
                self.assertEqual(record, records3[i])
                self.assertEqual(record, records4[i])
                self.assertEqual(record, records5[i])

                if i == t_count - 1:
                    i = -1
                if i < 3:
                    self.assertEqual(record.id, expected_records[i].id)
                    self.assertEqual(record.name, expected_records[i].name)
                    if record.seq is None:
                        length = None
                        seq = record.seq
                    else:
                        seq = str(record.seq)
                        length = len(seq)
                        if length > 50:
                            seq = str(seq[:40]) + "..." + str(seq[-7:])
                    self.assertEqual(seq, expected_records[i].seq)
                    self.assertEqual(length, expected_lengths[i])

            # Check Bio.SeqIO.read(...)
            if t_count == 1:
                record = SeqIO.read(t_filename, format=t_format)
                self.assertIsInstance(record, SeqRecord)
            else:
                self.assertRaises(ValueError, SeqIO.read, t_filename, t_format)

            # Check alphabets
            for record in records:
                base_alpha = Alphabet._get_base_alphabet(record.seq.alphabet)
                if isinstance(base_alpha, Alphabet.SingleLetterAlphabet):
                    if t_format in no_alpha_formats:
                        # Too harsh?
                        self.assertIs(base_alpha, Alphabet.single_letter_alphabet)
                else:
                    base_alpha = None
            if base_alpha is None:
                good = []
                bad = []
                given_alpha = None
            elif isinstance(base_alpha, Alphabet.ProteinAlphabet):
                good = protein_alphas
                bad = dna_alphas + rna_alphas + nucleotide_alphas
            elif isinstance(base_alpha, Alphabet.RNAAlphabet):
                good = nucleotide_alphas + rna_alphas
                bad = protein_alphas + dna_alphas
            elif isinstance(base_alpha, Alphabet.DNAAlphabet):
                good = nucleotide_alphas + dna_alphas
                bad = protein_alphas + rna_alphas
            elif isinstance(base_alpha, Alphabet.NucleotideAlphabet):
                good = nucleotide_alphas
                bad = protein_alphas
            else:
                self.assertIn(t_format, no_alpha_formats,
                              "Got %s from %s file"
                              % (repr(base_alpha), t_format))
                good = protein_alphas + dna_alphas + rna_alphas + nucleotide_alphas
                bad = []
            for given_alpha in good:
                # These should all work...
                given_base = Alphabet._get_base_alphabet(given_alpha)
                for record in SeqIO.parse(t_filename, t_format, given_alpha):
                    base_alpha = Alphabet._get_base_alphabet(record.seq.alphabet)
                    self.assertIsInstance(base_alpha, given_base.__class__)
                    self.assertEqual(base_alpha, given_base)
                if t_count == 1:
                    h = open(t_filename, mode)
                    record = SeqIO.read(h, t_format, given_alpha)
                    h.close()
                    self.assertIsInstance(base_alpha, given_base.__class__)
                    self.assertEqual(base_alpha, given_base)
            for given_alpha in bad:
                # These should all fail...
                with open(t_filename, mode) as h:
                    records6 = SeqIO.parse(h, t_format, given_alpha)
                    self.assertRaises(ValueError, next, records6,
                                      "Forcing wrong alphabet, %s, "
                                      "should fail (%s)"
                                      % (repr(given_alpha), t_filename))
            del good, bad, given_alpha, base_alpha

            if t_alignment:
                alignment = MultipleSeqAlignment(SeqIO.parse(
                        handle=t_filename, format=t_format))
                self.assertEqual(len(alignment), t_count)
                alignment_len = alignment.get_alignment_length()

                # Check the record order agrees, and double check the
                # sequence lengths all agree too.
                for i in range(t_count):
                    self.assertEqual(records[i], alignment[i])
                    self.assertEqual(len(records[i].seq), alignment_len)

                self.assertEqual(str(alignment_summary(alignment)),
                                 expected_alignment)

        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        records.reverse()
        self.check_simple_write_read(records, t_format, t_count,
                                     expected_messages)

    def test_sff1(self):
        records = [SeqRecord(Seq("tcagGGTCTACATGTTGGTTAACCCGTACTGATTTGAATT...GGGCTTa"), "E3MFGYR02JWQ7T", "E3MFGYR02JWQ7T"),
                   SeqRecord(Seq("tcagTTTTTTTTGGAAAGGAAAACGGACGTACTCATAGAT...AAATGcc"), "E3MFGYR02JA6IL", "E3MFGYR02JA6IL"),
                   SeqRecord(Seq("tcagAAAGACAAGTGGTATCAACGCAGAGTGGCCATTACG...gagaacg"), "E3MFGYR02JHD4H", "E3MFGYR02JHD4H"),
                   SeqRecord(Seq("tcagAATCATCCACTTTTTAACGTTTTGTTTTGTTCATCT...nnnnnnn"), "E3MFGYR02F7Z7G", "E3MFGYR02F7Z7G"),
                   ]
        lengths = [265, 271, 310, 219]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "More than one sequence found",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("sff", False, "Roche/E3MFGYR02_random_10_reads.sff", 10, records, lengths, alignment, messages)

    def test_clustal1(self):
        records = [SeqRecord(Seq("MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSE...GNRESVV"), "gi|4959044|gb|AAD34209.1|AF069"),
                   SeqRecord(Seq("---------MSPQTETKASVGFKAGVKEYKLTYYTPEYET...PAMD---"), "gi|671626|emb|CAA85685.1|"),
                   ]
        lengths = [601, 601]
        alignment = """\
 M- alignment column 0
 E- alignment column 1
 N- alignment column 2
 S- alignment column 3
 D- alignment column 4
 || ...
 V- alignment column 600"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("clustal", True, "Clustalw/cw02.aln", 2, records, lengths, alignment, messages)

    def test_clustal2(self):
        records = [SeqRecord(Seq("TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCG...TACCAGA"), "gi|6273285|gb|AF191659.1|AF191"),
                   SeqRecord(Seq("TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCG...TACCAGA"), "gi|6273284|gb|AF191658.1|AF191"),
                   SeqRecord(Seq("TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCG...TACCAGA"), "gi|6273287|gb|AF191661.1|AF191"),
                   SeqRecord(Seq("TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCG...TACCAGA"), "gi|6273291|gb|AF191665.1|AF191"),
                   ]
        lengths = [156, 156, 156, 156]
        alignment = """\
 TTTTTTT alignment column 0
 AAAAAAA alignment column 1
 TTTTTTT alignment column 2
 AAAAAAA alignment column 3
 CCCCCCC alignment column 4
 ||||||| ...
 AAAAAAA alignment column 155"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|6273291|gb|AF191665.1|AF191).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|6273291|gb|AF191665.1|AF191).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|6273291|gb|AF191665.1|AF191).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|6273291|gb|AF191665.1|AF191).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|6273291|gb|AF191665.1|AF191).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("clustal", True, "Clustalw/opuntia.aln", 7, records, lengths, alignment, messages)

    def test_clustal3(self):
        records = [SeqRecord(Seq("MFNLVSGTGGSSCCHRRNCFANRKKFFTMLLIFLLYMVSQ...-------"), "gi|167877390|gb|EDS40773.1|"),
                   SeqRecord(Seq("---------------------MRSASAAALLLAALLVVQA...-------"), "gi|167234445|ref|NP_001107837."),
                   SeqRecord(Seq("-------------MPQR----SLRHQLGMILVFFLLVTSH...D------"), "gi|74100009|gb|AAZ99217.1|"),
                   SeqRecord(Seq("----------------------LAADDQGRLLYSDFLTFL...GMAVKSS"), "gi|56122354|gb|AAV74328.1|"),
                   ]
        lengths = [447, 447, 447, 447]
        alignment = """\
 M---- alignment column 0
 F---- alignment column 1
 N---- alignment column 2
 L---- alignment column 3
 V---- alignment column 4
 ||||| ...
 ---SS alignment column 446"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|56122354|gb|AAV74328.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|56122354|gb|AAV74328.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|56122354|gb|AAV74328.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|56122354|gb|AAV74328.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|56122354|gb|AAV74328.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("clustal", True, "Clustalw/hedgehog.aln", 5, records, lengths, alignment, messages)

    def test_clustal4(self):
        records = [SeqRecord(Seq("----------------------------------------...AGAGTAG"), "AT3G20900.1-CDS"),
                   SeqRecord(Seq("ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAAC...AGAGTAG"), "AT3G20900.1-SEQ"),
                   ]
        lengths = [687, 687]
        alignment = """\
 -A alignment column 0
 -T alignment column 1
 -G alignment column 2
 -A alignment column 3
 -A alignment column 4
 || ...
 GG alignment column 686"""
        messages = ["Repeated name 'AT3G20900.' (originally 'AT3G20900.1-CDS'), possibly due to truncation",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AT3G20900.1-SEQ).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AT3G20900.1-SEQ).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AT3G20900.1-SEQ).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AT3G20900.1-SEQ).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AT3G20900.1-SEQ).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    "Repeated name 'AT3G20900.' (originally 'AT3G20900.1-CDS'), possibly due to truncation",
                    ]
        self.perform_test("clustal", True, "Clustalw/odd_consensus.aln", 2, records, lengths, alignment, messages)

    def test_fasta1(self):
        records = [SeqRecord(Seq("GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTC...GGTTTTT"), "gi|5049839|gb|AI730987.1|AI730987", "gi|5049839|gb|AI730987.1|AI730987"),
                   ]
        lengths = [655]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5049839|gb|AI730987.1|AI730987).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5049839|gb|AI730987.1|AI730987).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5049839|gb|AI730987.1|AI730987).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5049839|gb|AI730987.1|AI730987).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5049839|gb|AI730987.1|AI730987).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/lupine.nu", 1, records, lengths, alignment, messages)

    def test_fasta2(self):
        records = [SeqRecord(Seq("ATGAAGTTAAGCACTCTTCTCATCTTATCTTTTCCTTTCC...GTCGTTT"), "gi|4218935|gb|AF074388.1|AF074388", "gi|4218935|gb|AF074388.1|AF074388"),
                   ]
        lengths = [2050]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4218935|gb|AF074388.1|AF074388).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4218935|gb|AF074388.1|AF074388).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4218935|gb|AF074388.1|AF074388).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4218935|gb|AF074388.1|AF074388).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4218935|gb|AF074388.1|AF074388).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/elderberry.nu", 1, records, lengths, alignment, messages)

    def test_fasta3(self):
        records = [SeqRecord(Seq("TCGAAACCTGCCTAGCAGAACGACCCGCGAACTTGTATTC...CACGACC"), "gi|5052071|gb|AF067555.1|AF067555", "gi|5052071|gb|AF067555.1|AF067555"),
                   ]
        lengths = [623]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5052071|gb|AF067555.1|AF067555).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5052071|gb|AF067555.1|AF067555).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5052071|gb|AF067555.1|AF067555).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5052071|gb|AF067555.1|AF067555).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5052071|gb|AF067555.1|AF067555).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/phlox.nu", 1, records, lengths, alignment, messages)

    def test_fasta4(self):
        records = [SeqRecord(Seq("CCTGTCACTTAACTTTTTGTTCATAAGGTATATATGGGGG...GTTAGAG"), "gi|4104054|gb|AH007193.1|SEG_CVIGS", "gi|4104054|gb|AH007193.1|SEG_CVIGS"),
                   ]
        lengths = [1002]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4104054|gb|AH007193.1|SEG_CVIGS).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4104054|gb|AH007193.1|SEG_CVIGS).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4104054|gb|AH007193.1|SEG_CVIGS).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4104054|gb|AH007193.1|SEG_CVIGS).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4104054|gb|AH007193.1|SEG_CVIGS).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/centaurea.nu", 1, records, lengths, alignment, messages)

    def test_fasta5(self):
        records = [SeqRecord(Seq("GCTCCATTTTTTACACATTTCTATGAACTAATTGGTTCAT...ATGATGA"), "gi|5817701|gb|AF142731.1|AF142731", "gi|5817701|gb|AF142731.1|AF142731"),
                   ]
        lengths = [2551]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5817701|gb|AF142731.1|AF142731).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5817701|gb|AF142731.1|AF142731).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5817701|gb|AF142731.1|AF142731).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5817701|gb|AF142731.1|AF142731).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5817701|gb|AF142731.1|AF142731).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/wisteria.nu", 1, records, lengths, alignment, messages)

    def test_fasta6(self):
        records = [SeqRecord(Seq("CAGGCTGCGCGGTTTCTATTTATGAAGAACAAGGTCCGTA...GTTTGTT"), "gi|3176602|gb|U78617.1|LOU78617", "gi|3176602|gb|U78617.1|LOU78617"),
                   ]
        lengths = [309]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3176602|gb|U78617.1|LOU78617).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3176602|gb|U78617.1|LOU78617).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3176602|gb|U78617.1|LOU78617).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3176602|gb|U78617.1|LOU78617).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3176602|gb|U78617.1|LOU78617).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/sweetpea.nu", 1, records, lengths, alignment, messages)

    def test_fasta7(self):
        records = [SeqRecord(Seq("GGCTCTTAAGTCATGTCTAGGCAGGTGTGCACAAGTTTAG...GTAGGTG"), "gi|5690369|gb|AF158246.1|AF158246", "gi|5690369|gb|AF158246.1|AF158246"),
                   ]
        lengths = [550]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5690369|gb|AF158246.1|AF158246).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5690369|gb|AF158246.1|AF158246).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5690369|gb|AF158246.1|AF158246).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5690369|gb|AF158246.1|AF158246).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|5690369|gb|AF158246.1|AF158246).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/lavender.nu", 1, records, lengths, alignment, messages)

    def test_fasta8(self):
        records = [SeqRecord(Seq("GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLL...FIVGANI"), "gi|3298468|dbj|BAA31520.1|", "gi|3298468|dbj|BAA31520.1|"),
                   ]
        lengths = [107]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/aster.pro", 1, records, lengths, alignment, messages)

    def test_fasta9(self):
        records = [SeqRecord(Seq("GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLL...FIVGANI"), "gi|3298468|dbj|BAA31520.1|", "gi|3298468|dbj|BAA31520.1|"),
                   ]
        lengths = [107]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/aster_no_wrap.pro", 1, records, lengths, alignment, messages)

    def test_fasta_2line1(self):
        records = [SeqRecord(Seq("GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLL...FIVGANI"), "gi|3298468|dbj|BAA31520.1|", "gi|3298468|dbj|BAA31520.1|"),
                   ]
        lengths = [107]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3298468|dbj|BAA31520.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta-2line", False, "Fasta/aster_no_wrap.pro", 1, records, lengths, alignment, messages)

    def test_fasta10(self):
        records = [SeqRecord(Seq("XAGLPVIMCLKSNNHQKYLRYQSDNIQQYGLLQFSADKIL...IELGQNN"), "gi|2781234|pdb|1JLY|B", "gi|2781234|pdb|1JLY|B"),
                   ]
        lengths = [304]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|2781234|pdb|1JLY|B).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|2781234|pdb|1JLY|B).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|2781234|pdb|1JLY|B).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|2781234|pdb|1JLY|B).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|2781234|pdb|1JLY|B).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/loveliesbleeding.pro", 1, records, lengths, alignment, messages)

    def test_fasta11(self):
        records = [SeqRecord(Seq("MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSE...GNRESVV"), "gi|4959044|gb|AAD34209.1|AF069992_1", "gi|4959044|gb|AAD34209.1|AF069992_1"),
                   ]
        lengths = [600]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4959044|gb|AAD34209.1|AF069992_1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4959044|gb|AAD34209.1|AF069992_1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4959044|gb|AAD34209.1|AF069992_1).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4959044|gb|AAD34209.1|AF069992_1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|4959044|gb|AAD34209.1|AF069992_1).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/rose.pro", 1, records, lengths, alignment, messages)

    def test_fasta12(self):
        records = [SeqRecord(Seq("MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAF...FEFPAMD"), "gi|671626|emb|CAA85685.1|", "gi|671626|emb|CAA85685.1|"),
                   ]
        lengths = [473]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|671626|emb|CAA85685.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/rosemary.pro", 1, records, lengths, alignment, messages)

    def test_fasta13(self):
        # Protein
        records = [SeqRecord(Seq("MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAA...YVMFAVA"), "gi|3318709|pdb|1A91|", "gi|3318709|pdb|1A91|"),
                   ]
        lengths = [79]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3318709|pdb|1A91|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3318709|pdb|1A91|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3318709|pdb|1A91|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3318709|pdb|1A91|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|3318709|pdb|1A91|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "Fasta/f001", 1, records, lengths, alignment, messages)

    def test_fasta14(self):
        # DNA
        records = [SeqRecord(Seq("CGGACCAGACGGACACAGGGAGAAGCTAGTTTCTTTCATG...GGTTTNA"), "gi|1348912|gb|G26680|G26680", "gi|1348912|gb|G26680|G26680"),
                   SeqRecord(Seq("CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATC...NNNGAAA"), "gi|1348917|gb|G26685|G26685", "gi|1348917|gb|G26685|G26685"),
                   SeqRecord(Seq("GATCAAATCTGCACTGTGTCTACATATAGGAAAGGTCCTG...NTTTTTT"), "gi|1592936|gb|G29385|G29385", "gi|1592936|gb|G29385|G29385"),
                   ]
        lengths = [633, 413, 471]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|1592936|gb|G29385|G29385).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|1592936|gb|G29385|G29385).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|1592936|gb|G29385|G29385).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|1592936|gb|G29385|G29385).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|1592936|gb|G29385|G29385).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("fasta", False, "Fasta/f002", 3, records, lengths, alignment, messages)

    def test_fasta15(self):
        # Protein with gaps
        records = [SeqRecord(Seq("CPDSINAALICRGEKMSIAIMAGVLEARGH-N--VTVIDP...INIVAIA"), "AK1H_ECOLI/1-378", "AK1H_ECOLI/1-378"),
                   SeqRecord(Seq("-----------------VEDAVKATIDCRGEKLSIAMMKA...SALAQAN"), "AKH_HAEIN/1-382", "AKH_HAEIN/1-382"),
                   ]
        lengths = [378, 382]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AKH_HAEIN/1-382).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AKH_HAEIN/1-382).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AKH_HAEIN/1-382).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AKH_HAEIN/1-382).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AKH_HAEIN/1-382).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("fasta", False, "Fasta/fa01", 2, records, lengths, alignment, messages)

    def test_fasta16(self):
        # FASTA -> Tabbed
        records = [SeqRecord(Seq("TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTG...ACCCCTG"), "gi|45478711|ref|NC_005816.1|", "gi|45478711|ref|NC_005816.1|"),
                   ]
        lengths = [9609]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478711|ref|NC_005816.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478711|ref|NC_005816.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478711|ref|NC_005816.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478711|ref|NC_005816.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478711|ref|NC_005816.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "GenBank/NC_005816.fna", 1, records, lengths, alignment, messages)

    def test_fasta17(self):
        records = [SeqRecord(Seq("ATGGTCACTTTTGAGACAGTTATGGAAATTAAAATCCTGC...GGCGTGA"), "ref|NC_005816.1|:87-1109", "ref|NC_005816.1|:87-1109"),
                   SeqRecord(Seq("GTGATGATGGAACTGCAACATCAACGACTGATGGCGCTCG...TGAGTAA"), "ref|NC_005816.1|:1106-1888", "ref|NC_005816.1|:1106-1888"),
                   SeqRecord(Seq("GTGAACAAACAACAACAAACTGCGCTGAATATGGCGCGAT...AACATAA"), "ref|NC_005816.1|:2925-3119", "ref|NC_005816.1|:2925-3119"),
                   SeqRecord(Seq("TTGGCTGATTTGAAAAAGCTACAGGTTTACGGACCTGAGT...CAAGTAA"), "ref|NC_005816.1|:c8360-8088", "ref|NC_005816.1|:c8360-8088"),
                   ]
        lengths = [1023, 783, 195, 273]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=ref|NC_005816.1|:c8360-8088).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=ref|NC_005816.1|:c8360-8088).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=ref|NC_005816.1|:c8360-8088).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=ref|NC_005816.1|:c8360-8088).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=ref|NC_005816.1|:c8360-8088).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("fasta", False, "GenBank/NC_005816.ffn", 10, records, lengths, alignment, messages)

    def test_fasta18(self):
        records = [SeqRecord(Seq("MVTFETVMEIKILHKQGMSSRAIARELGISRNTVKRYLQA...SFCRGVA"), "gi|45478712|ref|NP_995567.1|", "gi|45478712|ref|NP_995567.1|"),
                   SeqRecord(Seq("MMMELQHQRLMALAGQLQLESLISAAPALSQQAVDQEWSY...IAEANPE"), "gi|45478713|ref|NP_995568.1|", "gi|45478713|ref|NP_995568.1|"),
                   SeqRecord(Seq("MNKQQQTALNMARFIRSQSLILLEKLDALDADEQAAMCER...AESETGT"), "gi|45478714|ref|NP_995569.1|", "gi|45478714|ref|NP_995569.1|"),
                   SeqRecord(Seq("MADLKKLQVYGPELPRPYADTVKGSRYKNMKELRVQFSGR...LNTLESK"), "gi|45478721|ref|NP_995576.1|", "gi|45478721|ref|NP_995576.1|"),
                   ]
        lengths = [340, 260, 64, 90]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("fasta", False, "GenBank/NC_005816.faa", 10, records, lengths, alignment, messages)

    def test_fasta19(self):
        records = [SeqRecord(Seq("MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYTI...YGVKKPK"), "gi|7525080|ref|NP_051037.1|", "gi|7525080|ref|NP_051037.1|"),
                   SeqRecord(Seq("MTAILERRESESLWGRFCNWITSTENRLYIGWFGVLMIPT...EAPSTNG"), "gi|7525013|ref|NP_051039.1|", "gi|7525013|ref|NP_051039.1|"),
                   SeqRecord(Seq("MDKFQGYLEFDGARQQSFLYPLFFREYIYVLAYDHGLNRL...NDLVNHE"), "gi|126022795|ref|NP_051040.2|", "gi|126022795|ref|NP_051040.2|"),
                   SeqRecord(Seq("MAIHLYKTSTPSTRNGAVDSQVKSNPRNNLICGQHHCGKG...ILRRRSK"), "gi|7525099|ref|NP_051123.1|", "gi|7525099|ref|NP_051123.1|"),
                   ]
        lengths = [123, 353, 504, 274]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|7525099|ref|NP_051123.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|7525099|ref|NP_051123.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|7525099|ref|NP_051123.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|7525099|ref|NP_051123.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|7525099|ref|NP_051123.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("fasta", False, "GenBank/NC_000932.faa", 85, records, lengths, alignment, messages)

    def test_tab1(self):
        records = [SeqRecord(Seq("MVTFETVMEIKILHKQGMSSRAIARELGISRNTVKRYLQA...SFCRGVA"), "gi|45478712|ref|NP_995567.1|", "gi|45478712|ref|NP_995567.1|"),
                   SeqRecord(Seq("MMMELQHQRLMALAGQLQLESLISAAPALSQQAVDQEWSY...IAEANPE"), "gi|45478713|ref|NP_995568.1|", "gi|45478713|ref|NP_995568.1|"),
                   SeqRecord(Seq("MNKQQQTALNMARFIRSQSLILLEKLDALDADEQAAMCER...AESETGT"), "gi|45478714|ref|NP_995569.1|", "gi|45478714|ref|NP_995569.1|"),
                   SeqRecord(Seq("MADLKKLQVYGPELPRPYADTVKGSRYKNMKELRVQFSGR...LNTLESK"), "gi|45478721|ref|NP_995576.1|", "gi|45478721|ref|NP_995576.1|"),
                   ]
        lengths = [340, 260, 64, 90]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|45478721|ref|NP_995576.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("tab", False, "GenBank/NC_005816.tsv", 10, records, lengths, alignment, messages)

    def test_fasta20(self):
        # upper case
        records = [SeqRecord(Seq("GGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTC...GTGCTTC"), "gi|9629357|ref|NC_001802.1|", "gi|9629357|ref|NC_001802.1|"),
                   ]
        lengths = [9181]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|NC_001802.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|NC_001802.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|NC_001802.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|NC_001802.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|NC_001802.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "GFF/NC_001802.fna", 1, records, lengths, alignment, messages)

    def test_fasta21(self):
        # lower case
        records = [SeqRecord(Seq("ggtctctctggttagaccagatctgagcctgggagctctc...gtgcttc"), "gi|9629357|ref|nc_001802.1|", "gi|9629357|ref|nc_001802.1|"),
                   ]
        lengths = [9181]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|nc_001802.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|nc_001802.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|nc_001802.1|).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|nc_001802.1|).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|9629357|ref|nc_001802.1|).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", False, "GFF/NC_001802lc.fna", 1, records, lengths, alignment, messages)

    def test_fasta22(self):
        # Trivial nucleotide alignment
        records = [SeqRecord(Seq("ACGTCGCG"), "test1", "test1"),
                   SeqRecord(Seq("GGGGCCCC"), "test2", "test2"),
                   SeqRecord(Seq("AAACACAC"), "test3", "test3"),
                   ]
        lengths = [8, 8, 8]
        alignment = """\
 AGA alignment column 0
 CGA alignment column 1
 GGA alignment column 2
 TGC alignment column 3
 CCA alignment column 4
 ||| ...
 GCC alignment column 7"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=test3).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=test3).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=test3).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=test3).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=test3).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fasta", True, "GFF/multi.fna", 3, records, lengths, alignment, messages)

    def test_fasta23(self):
        # contains blank line
        records = [SeqRecord(Seq("GATCCCTACCCTTNCCGTTGGTCTCTNTCGCTGACTCGAG...TTATTTC"), "gi|1348916|gb|G26684|G26684", "gi|1348916|gb|G26684|G26684"),
                   SeqRecord(Seq("MPVVVVASSKGGAGKSTTAVVLGTELAHKGVPVTMLDCDP...KLTEALR"), "gi|129628|sp|P07175|PARA_AGRTU", "gi|129628|sp|P07175|PARA_AGRTU"),
                   ]
        lengths = [285, 222]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|129628|sp|P07175|PARA_AGRTU).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|129628|sp|P07175|PARA_AGRTU).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|129628|sp|P07175|PARA_AGRTU).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|129628|sp|P07175|PARA_AGRTU).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|129628|sp|P07175|PARA_AGRTU).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("fasta", False, "Registry/seqs.fasta", 2, records, lengths, alignment, messages)

    def test_nexus1(self):
        records = [SeqRecord(Seq("A-C-G-Tc-gtgtgtgctct-t-t------ac-gtgtgtgctct-t-t"), "t1", "t1"),
                   SeqRecord(Seq("A-C-GcTc-gtg-----tct-t-t----acac-gtg-----tct-t-t"), "t2 the name", "t2 the name"),
                   SeqRecord(Seq("A-CcGcTc-gtgtgtgct--------acacac-gtgtgtgct------"), "isn'that [a] strange name?", "isn'that [a] strange name?"),
                   SeqRecord(Seq("cccccccc-cccccccccccNc-ccccccccc-cccccccccccNc-c"), "t9", "t9"),
                   ]
        lengths = [48, 48, 48, 48]
        alignment = """\
 AAAAAAAAc alignment column 0
 -----c?tc alignment column 1
 CCCCCCCCc alignment column 2
 --c-?a-tc alignment column 3
 GGGGGGGGc alignment column 4
 ||||||||| ...
 tt--?ag?c alignment column 47"""
        messages = ["Whitespace not allowed in identifier: one should be punished, for (that)!",
                    "Cannot have spaces in EMBL accession, 'one should be punished, for (that)!'",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=t9).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=t9).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=t9).",
                    "Invalid whitespace in 'one should be punished, for (that)!' for LOCUS line",
                    "Cannot have spaces in EMBL accession, 'one should be punished, for (that)!'",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=t9).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=t9).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("nexus", True, "Nexus/test_Nexus_input.nex", 9, records, lengths, alignment, messages)

    def test_swiss1(self):
        records = [SeqRecord(Seq("MGARGAPSRRRQAGRRLRYLPTGSFPFLLLLLLLCIQLGG...YSDLDFE"), "Q13454", "N33_HUMAN"),
                   ]
        lengths = [348]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=Q13454).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13454).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13454).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13454).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13454).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp001", 1, records, lengths, alignment, messages)

    def test_swiss2(self):
        records = [SeqRecord(Seq("MADQRQRSLSTSGESLYHVLGLDKNATSDDIKKSYRKLAL...YHTDGFN"), "P54101", "CSP_MOUSE"),
                   ]
        lengths = [198]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P54101).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P54101).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P54101).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P54101).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P54101).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp002", 1, records, lengths, alignment, messages)

    def test_swiss3(self):
        records = [SeqRecord(Seq("MDDREDLVYQAKLAEQAERYDEMVESMKKVAGMDVELTVE...DVEDENQ"), "P42655", "143E_HUMAN"),
                   ]
        lengths = [255]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P42655).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P42655).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P42655).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P42655).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P42655).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp003", 1, records, lengths, alignment, messages)

    def test_swiss4(self):
        records = [SeqRecord(Seq("TVKWIEAVALSDILEGDVLGVTVEGKELALYEVEGEIYAT...RVMIDLS"), "P23082", "NDOA_PSEPU"),
                   ]
        lengths = [103]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P23082).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P23082).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P23082).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P23082).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P23082).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp004", 1, records, lengths, alignment, messages)

    def test_swiss5(self):
        records = [SeqRecord(Seq("MNLLLTLLTNTTLALLLVFIAFWLPQLNVYAEKTSPYECG...EGLEWAE"), "P24973", "NU3M_BALPH"),
                   ]
        lengths = [115]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P24973).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P24973).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P24973).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P24973).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P24973).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp005", 1, records, lengths, alignment, messages)

    def test_swiss6(self):
        records = [SeqRecord(Seq("MTPHTHVRGPGDILQLTMAFYGSRALISAVELDLFTLLAG...AIGRKPR"), "P39896", "TCMO_STRGA"),
                   ]
        lengths = [339]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P39896).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P39896).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P39896).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P39896).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P39896).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp006", 1, records, lengths, alignment, messages)

    def test_swiss7(self):
        records = [SeqRecord(Seq("MANAGLQLLGFILAFLGWIGAIVSTALPQWRIYSYAGDNI...SSGKDYV"), "O95832", "CLD1_HUMAN"),
                   ]
        lengths = [211]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=O95832).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=O95832).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=O95832).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=O95832).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=O95832).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp007", 1, records, lengths, alignment, messages)

    def test_swiss8(self):
        records = [SeqRecord(Seq("MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPG...SLTACKV"), "P01892", "1A02_HUMAN"),
                   ]
        lengths = [365]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P01892).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01892).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01892).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01892).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01892).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp008", 1, records, lengths, alignment, messages)

    def test_swiss9(self):
        records = [SeqRecord(Seq("MAPAMEEIRQAQRAEGPAAVLAIGTSTPPNALYQADYPDY...VPIAGAE"), "O23729", "CHS3_BROFI"),
                   ]
        lengths = [394]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=O23729).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=O23729).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=O23729).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=O23729).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=O23729).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp009", 1, records, lengths, alignment, messages)

    def test_swiss10(self):
        records = [SeqRecord(Seq("MDKLDANVSSEEGFGSVEKVVLLTFLSTVILMAILGNLLV...AAQPSDT"), "Q13639", "5H4_HUMAN"),
                   ]
        lengths = [388]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp010", 1, records, lengths, alignment, messages)

    def test_swiss11(self):
        records = [SeqRecord(Seq("MGRRVPALRQLLVLAVLLLKPSQLQSRELSGSRCPEPCDC...PPRALTH"), "P16235", "LSHR_RAT"),
                   ]
        lengths = [700]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P16235).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P16235).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P16235).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P16235).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P16235).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp011", 1, records, lengths, alignment, messages)

    def test_swiss12(self):
        records = [SeqRecord(Seq("MQIFVKTLTGKTITLEVESSDTIDNVKTKIQDKEGIPPDQ...LRLRGGN"), "Q9Y736", "Q9Y736"),
                   ]
        lengths = [153]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=Q9Y736).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q9Y736).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q9Y736).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q9Y736).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q9Y736).",
                    "ncbiTaxID should be of type string or int",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp012", 1, records, lengths, alignment, messages)

    def test_swiss13(self):
        records = [SeqRecord(Seq("MGSKMASASRVVQVVKPHTPLIRFPDRRDNPKPNVSEALR...IQRGGPE"), "P82909", "P82909"),
                   ]
        lengths = [102]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P82909).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P82909).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P82909).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P82909).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P82909).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp013", 1, records, lengths, alignment, messages)

    def test_swiss14(self):
        records = [SeqRecord(Seq("TQSNPNEQNVELNRTSLYWGLLLIFVLAVLFSNYFFN"), "P12166", "PSBL_ORYSA"),
                   ]
        lengths = [37]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P12166).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P12166).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P12166).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P12166).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P12166).",
                    "ncbiTaxID should be of type string or int",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp014", 1, records, lengths, alignment, messages)

    def test_swiss15(self):
        records = [SeqRecord(Seq("MSFQAPRRLLELAGQSLLRDQALAISVLDELPRELFPRLF...FIGPTPC"), "IPI00383150", "IPI00383150.2"),
                   ]
        lengths = [457]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=IPI00383150).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IPI00383150).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IPI00383150).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IPI00383150).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IPI00383150).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp015", 1, records, lengths, alignment, messages)

    def test_swiss16(self):
        records = [SeqRecord(Seq("MMFSGFNADYEASSSRCSSASPAGDSLSYYHSPADSFSSM...SPTLLAL"), "P01100", "FOS_HUMAN"),
                   ]
        lengths = [380]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P01100).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01100).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01100).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01100).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01100).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/sp016", 1, records, lengths, alignment, messages)

    def test_swiss17(self):
        records = [SeqRecord(Seq("ARRERMTAREEASLRTLEGRRRATLLSARQGMMSARGDFL...TKNFGFV"), "Q62671", "EDD_RAT"),
                   ]
        lengths = [920]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=Q62671).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q62671).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q62671).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q62671).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q62671).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "Registry/EDD_RAT.dat", 1, records, lengths, alignment, messages)

    def test_uniprot_xml1(self):
        records = [SeqRecord(Seq("MDLINNKLNIEIQKFCLDLEKKYNINYNNLIDLWFNKEST...CLNDIPI"), "Q91G55", "043L_IIV6"),
                   ]
        lengths = [116]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=Q91G55).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q91G55).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q91G55).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q91G55).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q91G55).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("uniprot-xml", False, "SwissProt/uni001", 1, records, lengths, alignment, messages)

    def test_uniprot_xml2(self):
        records = [SeqRecord(Seq("MDLINNKLNIEIQKFCLDLEKKYNINYNNLIDLWFNKEST...CLNDIPI"), "Q91G55", "043L_IIV6"),
                   SeqRecord(Seq("MSGGSLINSIAINTRIKKIKKSLLQNYTKEKTDMIKILYL...SSEHMTV"), "O55717", "094L_IIV6"),
                   SeqRecord(Seq("MYFYKKYLHFFFVVSKFKFFLKMQVPFGCNMKGLGVLLGL...SLPTYYG"), "P0C9J6", "11011_ASFP4"),
                   ]
        lengths = [116, 118, 302]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P0C9J6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P0C9J6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P0C9J6).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P0C9J6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P0C9J6).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("uniprot-xml", False, "SwissProt/uni002", 3, records, lengths, alignment, messages)

    def test_uniprot_xml3(self):
        records = [SeqRecord(Seq("MDKLDANVSSEEGFGSVEKVVLLTFLSTVILMAILGNLLV...AAQPSDT"), "Q13639", "5HT4R_HUMAN"),
                   ]
        lengths = [388]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("uniprot-xml", False, "SwissProt/Q13639.xml", 1, records, lengths, alignment, messages)

    def test_swiss18(self):
        records = [SeqRecord(Seq("MDKLDANVSSEEGFGSVEKVVLLTFLSTVILMAILGNLLV...AAQPSDT"), "Q13639", "5HT4R_HUMAN"),
                   ]
        lengths = [388]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Q13639).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/Q13639.txt", 1, records, lengths, alignment, messages)

    def test_uniprot_xml4(self):
        records = [SeqRecord(Seq("FIVVVAVNSTLLTINAGDYIFYTDWAWTSFVVFSISQSTM...LNTWTYR"), "H2CNN8", "H2CNN8_9ARCH"),
                   ]
        lengths = [196]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("uniprot-xml", False, "SwissProt/H2CNN8.xml", 1, records, lengths, alignment, messages)

    def test_swiss19(self):
        records = [SeqRecord(Seq("FIVVVAVNSTLLTINAGDYIFYTDWAWTSFVVFSISQSTM...LNTWTYR"), "H2CNN8", "H2CNN8_9ARCH"),
                   ]
        lengths = [196]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=H2CNN8).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/H2CNN8.txt", 1, records, lengths, alignment, messages)

    def test_uniprot_xml5(self):
        records = [SeqRecord(Seq("MTMAAAQGKLSPDAIDNEVISNGSAKDYLDPPPAPLVDAG...SNYDAAV"), "F2CXE6", "F2CXE6_HORVD"),
                   ]
        lengths = [291]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("uniprot-xml", False, "SwissProt/F2CXE6.xml", 1, records, lengths, alignment, messages)

    def test_swiss20(self):
        records = [SeqRecord(Seq("MTMAAAQGKLSPDAIDNEVISNGSAKDYLDPPPAPLVDAG...SNYDAAV"), "F2CXE6", "F2CXE6_HORVD"),
                   ]
        lengths = [291]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=F2CXE6).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("swiss", False, "SwissProt/F2CXE6.txt", 1, records, lengths, alignment, messages)

    def test_genbank1(self):
        records = [SeqRecord(Seq("GGCAAGATGGCGCCGGTGGGGGTGGAGAAGAAGCTGCTGC...AAAAAAA"), "NM_006141.1", "NM_006141"),
                   ]
        lengths = [1622]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NM_006141.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NM_006141.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NM_006141.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NM_006141.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NM_006141.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/noref.gb", 1, records, lengths, alignment, messages)

    def test_genbank2(self):
        records = [SeqRecord(Seq("AACAAAACACACATCAAAAACGATTTTACAAGAAAAAAAT...AAAAAAA"), "X55053.1", "ATCOR66M"),
                   SeqRecord(Seq("ATTTGGCCTATAAATATAAACCCTTAAGCCCACATATCTT...AATTATA"), "X62281.1", "ATKIN2"),
                   SeqRecord(Seq("AAAAAAACACAACAAAACTCAATAAATAAACAAATGGCAG...AAGCTTC"), "M81224.1", "BNAKINI"),
                   SeqRecord(Seq("ATGGCAGACAACAAGCAGAGCTTCCAAGCCGGTCAAGCCG...CAAGTAG"), "AF297471.1", "AF297471"),
                   ]
        lengths = [513, 880, 441, 497]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AF297471.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AF297471.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AF297471.1).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AF297471.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AF297471.1).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("genbank", False, "GenBank/cor6_6.gb", 6, records, lengths, alignment, messages)

    def test_genbank3(self):
        records = [SeqRecord(Seq("CACAGGCCCAGAGCCACTCCTGCCTACAGGTTCTGAGGGC...AAAAAAA"), "AL109817.1", "IRO125195"),
                   ]
        lengths = [1326]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=AL109817.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL109817.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL109817.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL109817.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL109817.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/iro.gb", 1, records, lengths, alignment, messages)

    def test_genbank4(self):
        records = [SeqRecord(Seq("GATCATGCATGCACTCCAGCCTGGGACAAGAGCGAAACTC...GTTTGCA"), "U05344.1", "HUGLUT1"),
                   ]
        lengths = [741]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=U05344.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U05344.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U05344.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U05344.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U05344.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/pri1.gb", 1, records, lengths, alignment, messages)

    def test_genbank5(self):
        records = [SeqRecord(Seq("AAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTT...GAAGCTT"), "AC007323.5", "AC007323"),
                   ]
        lengths = [86436]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=AC007323.5).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AC007323.5).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AC007323.5).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AC007323.5).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AC007323.5).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/arab1.gb", 1, records, lengths, alignment, messages)

    def test_genbank6(self):
        records = [SeqRecord(Seq("MNNRWILHAAFLLCFSTTALSINYKQLQLQERTNIRKCQE...LTRNFQN"), "NP_034640.1", "NP_034640"),
                   ]
        lengths = [182]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/protein_refseq.gb", 1, records, lengths, alignment, messages)

    def test_genbank7(self):
        records = [SeqRecord(Seq("MNNRWILHAAFLLCFSTTALSINYKQLQLQERTNIRKCQE...LTRNFQN"), "NP_034640.1", "NP_034640"),
                   ]
        lengths = [182]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_034640.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/protein_refseq2.gb", 1, records, lengths, alignment, messages)

    def test_genbank8(self):
        records = [SeqRecord(Seq("TCCAGGGGATTCACGCGCAATATGTTTCCCTCGCTCGTCT...TCGATTG"), "AL138972.1", "DMBR25B3"),
                   ]
        lengths = [154329]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=AL138972.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL138972.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL138972.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL138972.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL138972.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/extra_keywords.gb", 1, records, lengths, alignment, messages)

    def test_genbank9(self):
        records = [SeqRecord(Seq("GAATTCAGATAGAATGTAGACAAGAGGGATGGTGAGGAAA...CAAAGGC"), "U18266.1", "HSTMPO1"),
                   ]
        lengths = [2509]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=U18266.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U18266.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U18266.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U18266.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U18266.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/one_of.gb", 1, records, lengths, alignment, messages)

    def test_genbank10(self):
        # contig, no sequence
        records = [SeqRecord(Seq("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNN"), "NT_019265.6", "NT_019265"),
                   ]
        lengths = [1250660]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NT_019265.6).",
                    ]
        self.perform_test("genbank", False, "GenBank/NT_019265.gb", 1, records, lengths, alignment, messages)

    def test_genbank11(self):
        records = [SeqRecord(Seq("TTAATTAACTGTCTTCGATTGCGTTTAATTGACGGTTTTC...TCAGCGC"), "NC_002678.1", "NC_002678"),
                   ]
        lengths = [180]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NC_002678.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_002678.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_002678.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_002678.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_002678.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/origin_line.gb", 1, records, lengths, alignment, messages)

    def test_genbank12(self):
        records = [SeqRecord(Seq("MEECWVTEIANGSKDGLDSNPMKDYMILSGPQKTAVAVLC...DLDLSDC"), "NP_001832.1", "NP_001832"),
                   ]
        lengths = [360]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NP_001832.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_001832.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_001832.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_001832.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_001832.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/blank_seq.gb", 1, records, lengths, alignment, messages)

    def test_genbank13(self):
        records = [SeqRecord(Seq("VKDGYIVDDRNCTYFCGRNAYCNEECTKLKGESGYCQWAS...KGPGRCN"), "P01485", "SCX3_BUTOC"),
                   ]
        lengths = [64]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=P01485).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01485).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01485).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01485).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=P01485).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/dbsource_wrap.gb", 1, records, lengths, alignment, messages)

    def test_genbank14(self):
        # See also AE017046.embl
        records = [SeqRecord(Seq("TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTG...ACCCCTG"), "NC_005816.1", "NC_005816"),
                   ]
        lengths = [9609]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NC_005816.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_005816.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_005816.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_005816.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_005816.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/NC_005816.gb", 1, records, lengths, alignment, messages)

    def test_genbank15(self):
        records = [SeqRecord(Seq("ATGGGCGAACGACGGGAATTGAACCCGCGATGGTGAATTC...GGGCATC"), "NC_000932.1", "NC_000932"),
                   ]
        lengths = [154478]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NC_000932.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_000932.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_000932.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_000932.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_000932.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/NC_000932.gb", 1, records, lengths, alignment, messages)

    def test_genbank16(self):
        """Test parsing Genbank file from Vector NTI with an odd LOCUS line."""
        records = [SeqRecord(Seq("GCTAGCGGAGTGTATACTGGCTTACTATGTTGGCACTGAT...GCCCATG"), "pBAD30", "pBAD30"),
                   ]
        lengths = [4923]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=pBAD30).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=pBAD30).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=pBAD30).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=pBAD30).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=pBAD30).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/pBAD30.gb", 1, records, lengths, alignment, messages)

    def test_genbank17(self):
        records = [SeqRecord(Seq("ATGTCTGGCAACCAGTATACTGAGGAAGTTATGGAGGGAG...GGATTAA"), "AB000048.1", "AB000048"),
                   SeqRecord(Seq("ATGTCTGGCAACCAGTATACTGAGGAAGTTATGGAGGGAG...GGATTAA"), "AB000049.1", "AB000049"),
                   SeqRecord(Seq("ATGAGTGATGGAGCAGTTCAACCAGACGGTGGTCAACCTG...ATATTAA"), "AB000050.1", "AB000050"),
                   ]
        lengths = [2007, 2007, 1755]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AB000050.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AB000050.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AB000050.1).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AB000050.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AB000050.1).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        # This example is a truncated copy of gbvrl1.seq from
        # ftp://ftp.ncbi.nih.gov/genbank/gbvrl1.seq.gz
        # including an NCBI header, and the first three records.
        self.perform_test("genbank", False, "GenBank/gbvrl1_start.seq", 3, records, lengths, alignment, messages)

    def test_genbank18(self):
        records = [SeqRecord(Seq("GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCG...ACCTGCA"), "NC_001422.1", "NC_001422"),
                   ]
        lengths = [5386]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NC_001422.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_001422.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_001422.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_001422.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NC_001422.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GFF/NC_001422.gbk", 1, records, lengths, alignment, messages)

    def test_genbank19(self):
        records = [SeqRecord(Seq("MKVKVLSLLVPALLVAGAANAAEVYNKDGNKLDLYGKVDG...LGLVYQF"), "NP_416719.1", "NP_416719"),
                   ]
        lengths = [367]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=NP_416719.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_416719.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_416719.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_416719.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NP_416719.1).",
                    "Missing SFF flow information",
                    ]
        # Generated with Entrez.efetch("protein", id="16130152",
        # rettype="gbwithparts")
        self.perform_test("genbank", False, "GenBank/NP_416719.gbwithparts", 1, records, lengths, alignment, messages)

    def test_genbank20(self):
        """Test parsing GenPept file with nasty bond locations."""
        records = [SeqRecord(Seq("AYTTFSATKNDQLKEPMFFGQPVQVARYDQQKYDIFEKLI...DLSNFQL"), "1MRR_A", "1MRR_A"),
                   ]
        lengths = [375]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=1MRR_A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1MRR_A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1MRR_A).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1MRR_A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1MRR_A).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("genbank", False, "GenBank/1MRR_A.gp", 1, records, lengths, alignment, messages)

    def test_genbank21(self):
        # This and the next one are a pair, and should be roughly equivalent.
        records = [SeqRecord(Seq("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNN"), "DS830848.1", "DS830848"),
                   ]
        lengths = [1311]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=DS830848.1).",
                    ]
        self.perform_test("genbank", False, "GenBank/DS830848.gb", 1, records, lengths, alignment, messages)

    def test_embl1(self):
        records = [SeqRecord(Seq("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNN"), "DS830848.1", "DS830848"),
                   ]
        lengths = [1311]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=DS830848.1).",
                    ]
        self.perform_test("embl", False, "EMBL/DS830848.embl", 1, records, lengths, alignment, messages)

    def test_embl2(self):
        records = [SeqRecord(Seq("CLARIIRYFYNAKA"), "A00022.1", "A00022"),
                   SeqRecord(Seq("CIARIIRYFYNAKA"), "A00028.1", "A00028"),
                   SeqRecord(Seq("RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRA...ERTCGGA"), "A00031.1", "A00031"),
                   SeqRecord(Seq("MAIGTLEATTLIRGRAMTTVLPSPELIASFVDIVGPGNAL...LREDRGE"), "CQ797900.1", "CQ797900"),
                   ]
        lengths = [14, 14, 58, 496]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CQ797900.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CQ797900.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CQ797900.1).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CQ797900.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CQ797900.1).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("embl", False, "EMBL/epo_prt_selection.embl", 9, records, lengths, alignment, messages)

    def test_embl3(self):
        records = [SeqRecord(Seq("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX...XXXXXXX"), "NRP00000001", "NRP00000001"),
                   SeqRecord(Seq("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX...XXXXXXX"), "NRP00000002", "NRP00000002"),
                   SeqRecord(Seq("XXXXXXXXXXXXXXXXXXXXXXXXX"), "NRP00210944", "NRP00210944"),
                   SeqRecord(Seq("XXXXXXXXXXXXXXXXXXXXXXXXX"), "NRP00210945", "NRP00210945"),
                   ]
        lengths = [358, 65, 25, 25]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NRP00210945).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NRP00210945).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NRP00210945).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NRP00210945).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=NRP00210945).",
                    "Sequence type is UnknownSeq but SeqXML requires sequence",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("embl", False, "EMBL/patents.embl", 4, records, lengths, alignment, messages)

    def test_embl4(self):
        records = [SeqRecord(Seq("AAACAAACCAAATATGGATTTTATTGTAGCCATATTTGCT...AAAAAAA"), "X56734.1", "X56734"),
                   ]
        lengths = [1859]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=X56734.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=X56734.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=X56734.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=X56734.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=X56734.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("embl", False, "EMBL/TRBG361.embl", 1, records, lengths, alignment, messages)

    def test_embl5(self):
        records = [SeqRecord(Seq("GCCGAGCTGACCCAGTCTCCATCCTCCCTGTCTGCATCTG...TATGAGA"), "DD231055.1", "DD231055"),
                   ]
        lengths = [315]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("embl", False, "EMBL/DD231055_edited.embl", 1, records, lengths, alignment, messages)

    def test_embl6(self):
        records = [SeqRecord(Seq("GCCGAGCTGACCCAGTCTCCATCCTCCCTGTCTGCATCTG...TATGAGA"), "DD231055.1", "DD231055"),
                   ]
        lengths = [315]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DD231055.1).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("embl", False, "EMBL/DD231055_edited2.embl", 1, records, lengths, alignment, messages)

    def test_embl7(self):
        records = [SeqRecord(Seq("GATCAGTAGACCCAGCGACAGCAGGGCGGGGCCCAGCAGG...CGAGCAT"), "AL031232", "SC10H5"),
                   ]
        lengths = [4870]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=AL031232).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL031232).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL031232).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL031232).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AL031232).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("embl", False, "EMBL/SC10H5.embl", 1, records, lengths, alignment, messages)

    def test_embl8(self):
        records = [SeqRecord(Seq("CAATTACTGCAATGCCCTCGTAATTAAGTGAATTTACAAT...CATCACC"), "U87107.1", "U87107"),
                   ]
        lengths = [8840]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=U87107.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U87107.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U87107.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U87107.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=U87107.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("embl", False, "EMBL/U87107.embl", 1, records, lengths, alignment, messages)

    def test_embl9(self):
        records = [SeqRecord(Seq("ATGAAGCTTTTAGTGCGCCTAGCACCGATCCTGCTGCTAG...GGTGTGA"), "AAA03323.1", "AAA03323"),
                   ]
        lengths = [1545]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=AAA03323.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AAA03323.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AAA03323.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AAA03323.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AAA03323.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("embl", False, "EMBL/AAA03323.embl", 1, records, lengths, alignment, messages)

    def test_embl10(self):
        records = [SeqRecord(Seq("TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTG...ACCCCTG"), "AE017046.1", "AE017046"),
                   ]
        lengths = [9609]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=AE017046.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE017046.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE017046.1).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE017046.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE017046.1).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("embl", False, "EMBL/AE017046.embl", 1, records, lengths, alignment, messages)

    def test_embl11(self):
        records = [SeqRecord(Seq("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNN"), "AJ229040.1", "AJ229040"),
                   SeqRecord(Seq("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNN"), "AL954800.2", "AL954800"),
                   ]
        lengths = [958952, 87191216]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=AL954800.2).",
                    ]
        self.perform_test("embl", False, "EMBL/Human_contigs.embl", 2, records, lengths, alignment, messages)

    def test_embl12(self):
        records = [SeqRecord(Seq("DVVMTQTPLSLPVSLGDQASISCRSSQSLVHRNGNTYLHW...AGTKLEI"), "DI500001", "DI500001"),
                   SeqRecord(Seq("GSSSSSSSSSSSSSSSSXYSCFWKTCT"), "DI500002", "DI500002"),
                   SeqRecord(Seq("GSSNRAT"), "DI500003", "DI500003"),
                   SeqRecord(Seq("MKRSKRFAVLAQRPVNQDGLIGEWPEEGLIAMDSPFDPVS...HIDLVRE"), "DI500020", "DI500020"),
                   ]
        lengths = [111, 27, 7, 754]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DI500020).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DI500020).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DI500020).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DI500020).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=DI500020).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("embl", False, "EMBL/kipo_prt_sample.embl", 20, records, lengths, alignment, messages)

    def test_embl13(self):
        """Test parsing embl file with wrapped locations and unspecified type."""
        records = [SeqRecord(Seq("CGACTTTCCACTGCCCTCTACGCCCGCGCAATGGGTCGTA...CTACGTT"), "Test", "Tester"),
                   ]
        lengths = [120]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Test).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Test).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Test).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Test).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Test).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("embl", False, "EMBL/location_wrap.embl", 1, records, lengths, alignment, messages)

    def test_embl14(self):
        """Test parsing file with features over-indented for EMBL."""
        records = [SeqRecord(Seq("GATTGATCAATGCAGGCTGTTATGACTCAGGAATCTGCAC...CACATCA"), "A04195", "A04195"),
                   ]
        lengths = [51]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("embl", False, "EMBL/A04195.imgt", 1, records, lengths, alignment, messages)

    def test_imgt1(self):
        """Test parsing file with features over-indented for EMBL."""
        records = [SeqRecord(Seq("GATTGATCAATGCAGGCTGTTATGACTCAGGAATCTGCAC...CACATCA"), "A04195", "A04195"),
                   ]
        lengths = [51]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=A04195).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("imgt", False, "EMBL/A04195.imgt", 1, records, lengths, alignment, messages)

    def test_imgt2(self):
        records = [SeqRecord(Seq("CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGG...ATTAAAA"), "HLA00001.1", "HLA00001"),
                   SeqRecord(Seq("GATTGGGGAGTCCCAGCCTTGGGGATTCCCCAACTCCGCA...TGACCCT"), "HLA02169.1", "HLA02169"),
                   SeqRecord(Seq("ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCT...AGTGTGA"), "HLA14798.1", "HLA14798"),
                   SeqRecord(Seq("GCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCG...AGATGGG"), "HLA03131.1", "HLA03131"),
                   ]
        lengths = [3503, 3291, 2903, 822]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA03131.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA03131.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA03131.1).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA03131.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA03131.1).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("imgt", False, "EMBL/hla_3260_sample.imgt", 8, records, lengths, alignment, messages)

    def test_stockholm1(self):
        records = [SeqRecord(Seq("UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAU...UUAAUGU"), "AP001509.1", "AP001509.1"),
                   SeqRecord(Seq("AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-U...GUGAGAU"), "AE007476.1", "AE007476.1"),
                   ]
        lengths = [104, 104]
        alignment = """\
 UA alignment column 0
 UA alignment column 1
 AA alignment column 2
 AA alignment column 3
 UU alignment column 4
 || ...
 UU alignment column 103"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE007476.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE007476.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE007476.1).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE007476.1).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=AE007476.1).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("stockholm", True, "Stockholm/simple.sth", 2, records, lengths, alignment, messages)

    def test_stockholm2(self):
        records = [SeqRecord(Seq("MTCRAQLIAVPRASSLAE--AIACAQKM----RVSRVPVYERS"), "O83071/192-246", "O83071"),
                   SeqRecord(Seq("MQHVSAPVFVFECTRLAY--VQHKLRAH----SRAVAIVLDEY"), "O83071/259-312", "O83071"),
                   SeqRecord(Seq("MIEADKVAHVQVGNNLEH--ALLVLTKT----GYTAIPVLDPS"), "O31698/18-71", "O31698"),
                   SeqRecord(Seq("EVMLTDIPRLHINDPIMK--GFGMVINN------GFVCVENDE"), "363253|refseq_protein.50.proto_past_mitoc_micro_vira|gi|94986659|ref|YP_594592.1|awsonia_intraceuaris_PHE/MN1-00", "363253|refseq_protein.50.proto_past_mitoc_micro_vira|gi|94986659|ref|YP_594592.1|awsonia_intraceuaris_PHE/MN1-00"),
                   ]
        lengths = [43, 43, 43, 43]
        alignment = """\
 MMMEEE alignment column 0
 TQIVVV alignment column 1
 CHEMMM alignment column 2
 RVALLL alignment column 3
 ASDTTT alignment column 4
 |||||| ...
 SYSEEE alignment column 42"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=363253|refseq_protein.50.proto_past_mitoc_micro_vira|gi|94986659|ref|YP_594592.1|awsonia_intraceuaris_PHE/MN1-00).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=363253|refseq_protein.50.proto_past_mitoc_micro_vira|gi|94986659|ref|YP_594592.1|awsonia_intraceuaris_PHE/MN1-00).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=363253|refseq_protein.50.proto_past_mitoc_micro_vira|gi|94986659|ref|YP_594592.1|awsonia_intraceuaris_PHE/MN1-00).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=363253|refseq_protein.50.proto_past_mitoc_micro_vira|gi|94986659|ref|YP_594592.1|awsonia_intraceuaris_PHE/MN1-00).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=363253|refseq_protein.50.proto_past_mitoc_micro_vira|gi|94986659|ref|YP_594592.1|awsonia_intraceuaris_PHE/MN1-00).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("stockholm", True, "Stockholm/funny.sth", 6, records, lengths, alignment, messages)

    def test_phylip1(self):
        records = [SeqRecord(Seq("CGATGCTTACCGC"), "Archaeopt", "Archaeopt"),
                   SeqRecord(Seq("CGTTACTCGTTGT"), "Hesperorni", "Hesperorni"),
                   SeqRecord(Seq("TAATGTTAATTGT"), "Baluchithe", "Baluchithe"),
                   SeqRecord(Seq("GGCAGCCAATCAC"), "B.subtilis", "B.subtilis"),
                   ]
        lengths = [13, 13, 13, 13]
        alignment = """\
 CCTTCG alignment column 0
 GGAAAG alignment column 1
 ATAAAC alignment column 2
 TTTTAA alignment column 3
 GAGGAG alignment column 4
 |||||| ...
 CTTTTC alignment column 12"""
        messages = ["Whitespace not allowed in identifier: B. virgini",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("phylip", True, "Phylip/reference_dna.phy", 6, records, lengths, alignment, messages)

    def test_phylip2(self):
        records = [SeqRecord(Seq("CGATGCTTACCGCCGATGCTTACCGCCGATGCTTACCGC"), "Archaeopt", "Archaeopt"),
                   SeqRecord(Seq("CGTTACTCGTTGTCGTTACTCGTTGTCGTTACTCGTTGT"), "Hesperorni", "Hesperorni"),
                   SeqRecord(Seq("TAATGTTAATTGTTAATGTTAATTGTTAATGTTAATTGT"), "Baluchithe", "Baluchithe"),
                   SeqRecord(Seq("GGCAGCCAATCACGGCAGCCAATCACGGCAGCCAATCAC"), "B.subtilis", "B.subtilis"),
                   ]
        lengths = [39, 39, 39, 39]
        alignment = """\
 CCTTCG alignment column 0
 GGAAAG alignment column 1
 ATAAAC alignment column 2
 TTTTAA alignment column 3
 GAGGAG alignment column 4
 |||||| ...
 CTTTTC alignment column 38"""
        messages = ["Whitespace not allowed in identifier: B. virgini",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=B.subtilis).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("phylip", True, "Phylip/reference_dna2.phy", 6, records, lengths, alignment, messages)

    def test_phylip3(self):
        records = [SeqRecord(Seq("CACACACAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAA"), "A", "A"),
                   SeqRecord(Seq("CACACAACAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAA"), "B", "B"),
                   SeqRecord(Seq("CACAACAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAA"), "C", "C"),
                   SeqRecord(Seq("ACAAAAAAAAACAAAAACACAAAAAAAAAAAAAAAAAAAA"), "J", "J"),
                   ]
        lengths = [40, 40, 40, 40]
        alignment = """\
 CCCCCAAAAA alignment column 0
 AAAAACCCCC alignment column 1
 CCCAAAAAAA alignment column 2
 AAACCAAAAA alignment column 3
 CCAAAAAAAA alignment column 4
 |||||||||| ...
 AAAAAAAAAA alignment column 39"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("phylip", True, "Phylip/hennigian.phy", 10, records, lengths, alignment, messages)

    def test_phylip4(self):
        records = [SeqRecord(Seq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), "Mesohippus", "Mesohippus"),
                   SeqRecord(Seq("AAACCCCCCCAAAAAAAAACAAAAAAAAAAAAAAAAAAAA"), "Hypohippus", "Hypohippus"),
                   SeqRecord(Seq("CAAAAAAAAAAAAAAAACACAAAAAAAAAAAAAAAAAAAA"), "Archaeohip", "Archaeohip"),
                   SeqRecord(Seq("CCCACCCCCCCCCACACCCCAAAAAAAAAAAAAAAAAAAA"), "Pliohippus", "Pliohippus"),
                   ]
        lengths = [40, 40, 40, 40]
        alignment = """\
 AACCCCCCCC alignment column 0
 AAAACCCCCC alignment column 1
 AAAAAAAAAC alignment column 2
 ACAAAAAAAA alignment column 3
 ACACCCCCCC alignment column 4
 |||||||||| ...
 AAAAAAAAAA alignment column 39"""
        messages = ["Whitespace not allowed in identifier: M. secundu",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Pliohippus).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Pliohippus).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Pliohippus).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Pliohippus).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=Pliohippus).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("phylip", True, "Phylip/horses.phy", 10, records, lengths, alignment, messages)

    def test_phylip5(self):
        records = [SeqRecord(Seq("CACACAACCAAACAAACCACAAAAAAAAAAAAAAAAAAAA"), "A", "A"),
                   SeqRecord(Seq("AAACCACACACACAAACCCAAAAAAAAAAAAAAAAAAAAA"), "B", "B"),
                   SeqRecord(Seq("ACAAAACCAAACCACCCACAAAAAAAAAAAAAAAAAAAAA"), "C", "C"),
                   SeqRecord(Seq("CCAAAAACACCCAACCCAACAAAAAAAAAAAAAAAAAAAA"), "J", "J"),
                   ]
        lengths = [40, 40, 40, 40]
        alignment = """\
 CAAAACAAAC alignment column 0
 AACAACCACC alignment column 1
 CAAAACAAAA alignment column 2
 ACAACACACA alignment column 3
 CCAAAACCAA alignment column 4
 |||||||||| ...
 AAAAAAAAAA alignment column 39"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=J).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("phylip", True, "Phylip/random.phy", 10, records, lengths, alignment, messages)

    def test_phylip6(self):
        records = [SeqRecord(Seq("-----MKVILLFVLAVFTVFVSS---------------RG...STSII--"), "CYS1_DICDI", "CYS1_DICDI"),
                   SeqRecord(Seq("MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRA...SYPVVAA"), "ALEU_HORVU", "ALEU_HORVU"),
                   SeqRecord(Seq("------MWATLPLLCAGAWLLGV--------PVCGAAELS...SYPIPLV"), "CATH_HUMAN", "CATH_HUMAN"),
                   ]
        lengths = [384, 384, 384]
        alignment = """\
 -M- alignment column 0
 -A- alignment column 1
 -H- alignment column 2
 -A- alignment column 3
 -R- alignment column 4
 ||| ...
 -AV alignment column 383"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CATH_HUMAN).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CATH_HUMAN).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CATH_HUMAN).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CATH_HUMAN).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CATH_HUMAN).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("phylip", True, "Phylip/interlaced.phy", 3, records, lengths, alignment, messages)

    def test_phylip7(self):
        records = [SeqRecord(Seq("TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCC...AGDRSHE"), "IXI_234", "IXI_234"),
                   SeqRecord(Seq("TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCC...AGDRSHE"), "IXI_235", "IXI_235"),
                   SeqRecord(Seq("TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPGRPCC...AGDRSHE"), "IXI_236", "IXI_236"),
                   SeqRecord(Seq("TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----C...AGDRSHE"), "IXI_237", "IXI_237"),
                   ]
        lengths = [131, 131, 131, 131]
        alignment = """\
 TTTT alignment column 0
 SSSS alignment column 1
 PPPP alignment column 2
 AAAA alignment column 3
 SSSS alignment column 4
 |||| ...
 EEEE alignment column 130"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("phylip", True, "Phylip/interlaced2.phy", 4, records, lengths, alignment, messages)

    def test_emboss1(self):
        records = [SeqRecord(Seq("TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCC...AGDRSHE"), "IXI_234"),
                   SeqRecord(Seq("TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCC...AGDRSHE"), "IXI_235"),
                   SeqRecord(Seq("TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPGRPCC...AGDRSHE"), "IXI_236"),
                   SeqRecord(Seq("TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----C...AGDRSHE"), "IXI_237"),
                   ]
        lengths = [131, 131, 131, 131]
        alignment = """\
 TTTT alignment column 0
 SSSS alignment column 1
 PPPP alignment column 2
 AAAA alignment column 3
 SSSS alignment column 4
 |||| ...
 EEEE alignment column 130"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_237).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("emboss", True, "Emboss/alignret.txt", 4, records, lengths, alignment, messages)

    def test_emboss2(self):
        records = [SeqRecord(Seq("KILIVDD----QYGIRILLNEVFNKEGYQTFQAANGLQAL...-------"), "ref_rec"),
                   SeqRecord(Seq("-VLLADDHALVRRGFRLMLED--DPEIEIVAEAGDGAQAV...RVANGET"), "gi|94968718|receiver"),
                   SeqRecord(Seq("KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVT...-------"), "ref_rec"),
                   SeqRecord(Seq("TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVR...AVLQKRQ"), "gi|94970041|receiver"),
                   ]
        lengths = [124, 124, 119, 125]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|94970041|receiver).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|94970041|receiver).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|94970041|receiver).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|94970041|receiver).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=gi|94970041|receiver).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("emboss", False, "Emboss/needle.txt", 10, records, lengths, alignment, messages)

    def test_emboss3(self):
        records = [SeqRecord(Seq("TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCC...AGDRSHE"), "IXI_234"),
                   SeqRecord(Seq("TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCC...AGDRSHE"), "IXI_235"),
                   ]
        lengths = [131, 131]
        alignment = """\
 TT alignment column 0
 SS alignment column 1
 PP alignment column 2
 AA alignment column 3
 SS alignment column 4
 || ...
 EE alignment column 130"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_235).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_235).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_235).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_235).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=IXI_235).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("emboss", True, "Emboss/water.txt", 2, records, lengths, alignment, messages)

    def test_phd1(self):
        records = [SeqRecord(Seq("ctccgtcggaacatcatcggatcctatcacagagtttttg...aagcgtg"), "34_222_(80-A03-19).b.ab1", "34_222_(80-A03-19).b.ab1"),
                   SeqRecord(Seq("cgggatcccacctgatccgaggtcaacctgaaaaaatatg...agccaag"), "425_103_(81-A03-19).g.ab1", "425_103_(81-A03-19).g.ab1"),
                   SeqRecord(Seq("acataaatcaaattactnaccaacacacaaaccngtctcg...tgctttn"), "425_7_(71-A03-19).b.ab1", "425_7_(71-A03-19).b.ab1"),
                   ]
        lengths = [876, 862, 1280]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "More than one sequence found",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("phd", False, "Phd/phd1", 3, records, lengths, alignment, messages)

    def test_phd2(self):
        records = [SeqRecord(Seq("actttggtcgcctgcaggtaccggtccgngattcccgggt...ggtgaga"), "ML4924R", "ML4924R"),
                   ]
        lengths = [180]
        alignment = None
        messages = ["Missing SFF flow information",
                    ]
        self.perform_test("phd", False, "Phd/phd2", 1, records, lengths, alignment, messages)

    def test_phd3(self):
        records = [SeqRecord(Seq("gccaatcaggtttctctgcaagcccctttagcagctgagc"), "HWI-EAS94_4_1_1_537_446", "HWI-EAS94_4_1_1_537_446"),
                   SeqRecord(Seq("gccatggcacatatatgaaggtcagaggacaacttgctgt"), "HWI-EAS94_4_1_1_602_99", "HWI-EAS94_4_1_1_602_99"),
                   ]
        lengths = [40, 40]
        alignment = None
        messages = ["Repeated name 'HWI-EAS94_' (originally 'HWI-EAS94_4_1_1_537_446'), possibly due to truncation",
                    "More than one sequence found",
                    "Missing SFF flow information",
                    "Repeated name 'HWI-EAS94_' (originally 'HWI-EAS94_4_1_1_537_446'), possibly due to truncation",
                    ]
        self.perform_test("phd", False, "Phd/phd_solexa", 2, records, lengths, alignment, messages)

    def test_phd4(self):
        records = [SeqRecord(Seq("ggggatgaaagggatctcggtggtaggtga"), "EBE03TV04IHLTF.77-243", "EBE03TV04IHLTF.77-243"),
                   ]
        lengths = [30]
        alignment = None
        messages = ["Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "Missing SFF flow information",
                    ]
        self.perform_test("phd", False, "Phd/phd_454", 1, records, lengths, alignment, messages)

    def test_ace1(self):
        records = [SeqRecord(Seq("aatacgGGATTGCCCTAGTAACGGCGAGTGAAGCGGCAAC...CTAGtac"), "Contig1", "Contig1"),
                   SeqRecord(Seq("cacggatgatagcttcgcgacactagcttttcagctaacc...cttgtag"), "Contig2", "Contig2"),
                   ]
        lengths = [856, 3296]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "More than one sequence found",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("ace", False, "Ace/contig1.ace", 2, records, lengths, alignment, messages)

    def test_ace2(self):
        records = [SeqRecord(Seq("agccccgggccgtggggttccttgagcactcccaaagttc...gggtttg"), "Contig1", "Contig1"),
                   ]
        lengths = [1475]
        alignment = None
        messages = ["Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "Missing SFF flow information",
                    ]
        self.perform_test("ace", False, "Ace/consed_sample.ace", 1, records, lengths, alignment, messages)

    def test_ace3(self):
        records = [SeqRecord(Seq("AGTTTTAGTTTTCCTCTGAAGCAAGCACACCTTCCCTTTC...TCACATT"), "Contig1", "Contig1"),
                   ]
        lengths = [1222]
        alignment = None
        messages = ["Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "Missing SFF flow information",
                    ]
        self.perform_test("ace", False, "Ace/seq.cap.ace", 1, records, lengths, alignment, messages)

    def test_ig1(self):
        records = [SeqRecord(Seq("ATGGAGCCAGTAGATCCTAACCTAGAGCCCTGGAAACACC...ATTCGCT"), "A_U455", "A_U455"),
                   SeqRecord(Seq("ATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATC...CGATTAG"), "B_HXB2R", "B_HXB2R"),
                   SeqRecord(Seq("%CAGGAAGTCAGCCTAAAACTCCTTGTACTAAGTGTTTTG...AGATTAA"), "C_UG268A", "C_UG268A"),
                   SeqRecord(Seq("ATGTCCTCAACGGACCAGATATGCCAGACACAGAGGGTAC...GAATCTT"), "SYK_SYK", "SYK_SYK"),
                   ]
        lengths = [303, 306, 267, 330]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK_SYK).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK_SYK).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK_SYK).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK_SYK).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK_SYK).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("ig", False, "IntelliGenetics/TAT_mase_nuc.txt", 17, records, lengths, alignment, messages)

    def test_ig2(self):
        records = [SeqRecord(Seq("MEN--RW-QVMIVWQVDRMRIRTWKSLVKHHMYRSKKA-K...-----GH"), "most-likely", "most-likely"),
                   SeqRecord(Seq("MEN--RW-QVMIVWQVDRMKIRTWNSLVKHHMYVSKKA-Q...-----RH"), "U455", "U455"),
                   SeqRecord(Seq("MEN--RW-QVMIVWQVDRMRIRTWKSLVKHHMYVSGKA-R...-----GH"), "HXB2R", "HXB2R"),
                   SeqRecord(Seq("MEK--EW-IVVPTWRMTPRQIDRLQHIIKTHKYKSKELEK...-------"), "SYK", "SYK"),
                   ]
        lengths = [298, 298, 298, 298]
        alignment = """\
 MMMMMMMMMMMMMMMM alignment column 0
 EEEEEEETEEEENEEE alignment column 1
 NNNNNNNAEEEEQRKK alignment column 2
 --------DEEEEE-- alignment column 3
 --------KKKKKK-- alignment column 4
 |||||||||||||||| ...
 HHHHHHH-AAAAL-R- alignment column 297"""
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=SYK).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("ig", True, "IntelliGenetics/VIF_mase-pro.txt", 16, records, lengths, alignment, messages)

    def test_ig3(self):
        """Test parsing a MASE alignment with sequence O_ANT70 being shorter."""
        records = [SeqRecord(Seq("ATGc?tcattt?ga??t?ttagcaaTaa?agcattaatag...?gA?ctg"), "VPU_CONSENSUS", "VPU_CONSENSUS"),
                   SeqRecord(Seq("ATGACACCTTTGGAAATCTGGGCAATAACAGGGCTGATAG...-AATTTG"), "A_U455", "A_U455"),
                   SeqRecord(Seq("ATGCAATCTTTACAAATATTAGCAATAGTATCATTAGTAG...-GATCTG"), "B_SF2", "B_SF2"),
                   SeqRecord(Seq("ATGACTAATATATTTGAGTATGCTTTT-------------...-GACGAA"), "CPZANT", "CPZANT"),
                   ]
        lengths = [294, 294, 294, 294]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CPZANT).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CPZANT).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CPZANT).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CPZANT).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=CPZANT).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("ig", False, "IntelliGenetics/vpu_nucaligned.txt", 9, records, lengths, alignment, messages)

    def test_pir1(self):
        records = [SeqRecord(Seq("ATGCTGGTCATGGCGCCCCGAACCGTCCTCCTGCTGCTCT...AGCTTGA"), "HLA:HLA00132", "HLA:HLA00132"),
                   SeqRecord(Seq("ATGCTGGTCATGGCGCCCCGAACCGTCCTCCTGCTGCTCT...AAGAGTT"), "HLA:HLA00133", "HLA:HLA00133"),
                   SeqRecord(Seq("GCTCCCACTCCATGAGGTATTTCTACACCTCCGTGTCCCG...CGCGCTG"), "HLA:HLA00134", "HLA:HLA00134"),
                   SeqRecord(Seq("ATGCGGGTCACGGCGCCCCGAACCCTCCTCCTGCTGCTCT...CGCGCGG"), "HLA:HLA01135", "HLA:HLA01135"),
                   ]
        lengths = [1089, 1009, 546, 619]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01135).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01135).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01135).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01135).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01135).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("pir", False, "NBRF/B_nuc.pir", 444, records, lengths, alignment, messages)

    def test_pir2(self):
        records = [SeqRecord(Seq("MRVMAPRTLILLLSGALALTETWACSHSMKYFFTSVSRPG...SLIACKA"), "HLA:HLA00401", "HLA:HLA00401"),
                   SeqRecord(Seq("MRVMAPRTLILLLSGALALTETWACSHSMKYFFTSVSRPG...SLIACKA"), "HLA:HLA00402", "HLA:HLA00402"),
                   SeqRecord(Seq("MRVMAPRTLILLLSGALALTETWACSHSMKYFFTSVSRPG...SLIACKA"), "HLA:HLA01075", "HLA:HLA01075"),
                   SeqRecord(Seq("MRVMAPRALLLLLSGGLALTETWACSHSMRYFDTAVSRPG...SLIACKA"), "HLA:HLA00484", "HLA:HLA00484"),
                   ]
        lengths = [366, 366, 366, 366]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00484).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00484).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00484).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00484).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00484).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("pir", False, "NBRF/Cw_prot.pir", 111, records, lengths, alignment, messages)

    def test_pir3(self):
        records = [SeqRecord(Seq("ATGGGTCATGAACAGAACCAAGGAGCTGCGCTGCTACAGA...TGACTGA"), "HLA:HLA00485", "HLA:HLA00485"),
                   SeqRecord(Seq("CTCCTACTCCAATGTGGCCAGATGACCTGCAAAACCACAC...TATTGGG"), "HLA:HLA00486", "HLA:HLA00486"),
                   SeqRecord(Seq("GGGTTTCCTATCGCTGAAGTGTTCACGCTGAAGCCCCTGG...CTATTGG"), "HLA:HLA00487", "HLA:HLA00487"),
                   SeqRecord(Seq("GGGTTTCCTATCGCTGAAGTGTTCACGCTGAAGCCCCTGG...CTATTGG"), "HLA:HLA00488", "HLA:HLA00488"),
                   ]
        lengths = [786, 564, 279, 279]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00488).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00488).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00488).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00488).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA00488).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("pir", False, "NBRF/DMA_nuc.pir", 4, records, lengths, alignment, messages)

    def test_pir4(self):
        records = [SeqRecord(Seq("MITFLPLLLGLSLGCTGAGGFVAHVESTCLLDDAGTPKDF...SEGWHIS"), "HLA:HLA00489", "HLA:HLA00489"),
                   SeqRecord(Seq("PPSVQVAKTTPFNTREPVMLACYVWGFYPAEVTITWRKNG...EPILRDW"), "HLA:HLA00490", "HLA:HLA00490"),
                   SeqRecord(Seq("PPSVQVAKTTPFNTREPVMLACYVWGFYPAEVTITWRKNG...EPILRDW"), "HLA:HLA00491", "HLA:HLA00491"),
                   SeqRecord(Seq("GFVAHVESTCLLDDAGTPKDFTYCIFFNKDLLTCWDPEEN...EPILRDW"), "HLA:HLA01083", "HLA:HLA01083"),
                   ]
        lengths = [263, 94, 94, 188]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01083).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01083).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01083).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01083).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=HLA:HLA01083).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("pir", False, "NBRF/DMB_prot.pir", 6, records, lengths, alignment, messages)

    def test_pir5(self):
        records = [SeqRecord(Seq("----------------------------------------...-------"), "804Angiostrongylus_cantonensis", "804Angiostrongylus_cantonensis"),
                   SeqRecord(Seq("----------------------------------------...-------"), "815Parelaphostrongylus_odocoil", "815Parelaphostrongylus_odocoil"),
                   ]
        lengths = [2527, 2527]
        alignment = """\
 -- alignment column 0
 -- alignment column 1
 -- alignment column 2
 -- alignment column 3
 -- alignment column 4
 || ...
 -- alignment column 2526"""
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=815Parelaphostrongylus_odocoil).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=815Parelaphostrongylus_odocoil).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=815Parelaphostrongylus_odocoil).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=815Parelaphostrongylus_odocoil).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=815Parelaphostrongylus_odocoil).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("pir", True, "NBRF/clustalw.pir", 2, records, lengths, alignment, messages)

    def test_fasta24(self):
        records = [SeqRecord(Seq("CCCTTCTTGTCTTCAGCGTTTCTCC"), "EAS54_6_R1_2_1_413_324", "EAS54_6_R1_2_1_413_324"),
                   SeqRecord(Seq("TTGGCAGGCCAAGGCCGATGGATCA"), "EAS54_6_R1_2_1_540_792", "EAS54_6_R1_2_1_540_792"),
                   SeqRecord(Seq("GTTGCTTCTGGCGTGGGTGGGGGGG"), "EAS54_6_R1_2_1_443_348", "EAS54_6_R1_2_1_443_348"),
                   ]
        lengths = [25, 25, 25]
        alignment = """\
 CTG alignment column 0
 CTT alignment column 1
 CGT alignment column 2
 TGG alignment column 3
 TCC alignment column 4
 ||| ...
 CAG alignment column 24"""
        messages = ["Repeated name 'EAS54_6_R1' (originally 'EAS54_6_R1_2_1_540_792'), possibly due to truncation",
                    "Need a DNA, RNA or Protein alphabet",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=EAS54_6_R1_2_1_443_348).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=EAS54_6_R1_2_1_443_348).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=EAS54_6_R1_2_1_443_348).",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=EAS54_6_R1_2_1_443_348).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=EAS54_6_R1_2_1_443_348).",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    "Repeated name 'EAS54_6_R1' (originally 'EAS54_6_R1_2_1_540_792'), possibly due to truncation",
                    ]
        self.perform_test("fasta", True, "Quality/example.fasta", 3, records, lengths, alignment, messages)

    def test_qual1(self):
        records = [SeqRecord(Seq("?????????????????????????"), "EAS54_6_R1_2_1_413_324", "EAS54_6_R1_2_1_413_324"),
                   SeqRecord(Seq("?????????????????????????"), "EAS54_6_R1_2_1_540_792", "EAS54_6_R1_2_1_540_792"),
                   SeqRecord(Seq("?????????????????????????"), "EAS54_6_R1_2_1_443_348", "EAS54_6_R1_2_1_443_348"),
                   ]
        lengths = [25, 25, 25]
        alignment = None
        messages = ["Repeated name 'EAS54_6_R1' (originally 'EAS54_6_R1_2_1_540_792'), possibly due to truncation",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "Sequence type is UnknownSeq but SeqXML requires sequence",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    "Repeated name 'EAS54_6_R1' (originally 'EAS54_6_R1_2_1_540_792'), possibly due to truncation",
                    ]
        self.perform_test("qual", False, "Quality/example.qual", 3, records, lengths, alignment, messages)

    def test_fastq1(self):
        records = [SeqRecord(Seq("CCCTTCTTGTCTTCAGCGTTTCTCC"), "EAS54_6_R1_2_1_413_324", "EAS54_6_R1_2_1_413_324"),
                   SeqRecord(Seq("TTGGCAGGCCAAGGCCGATGGATCA"), "EAS54_6_R1_2_1_540_792", "EAS54_6_R1_2_1_540_792"),
                   SeqRecord(Seq("GTTGCTTCTGGCGTGGGTGGGGGGG"), "EAS54_6_R1_2_1_443_348", "EAS54_6_R1_2_1_443_348"),
                   ]
        lengths = [25, 25, 25]
        alignment = """\
 CTG alignment column 0
 CTT alignment column 1
 CGT alignment column 2
 TGG alignment column 3
 TCC alignment column 4
 ||| ...
 CAG alignment column 24"""
        messages = ["Repeated name 'EAS54_6_R1' (originally 'EAS54_6_R1_2_1_540_792'), possibly due to truncation",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    "Repeated name 'EAS54_6_R1' (originally 'EAS54_6_R1_2_1_540_792'), possibly due to truncation",
                    ]
        self.perform_test("fastq", True, "Quality/example.fastq", 3, records, lengths, alignment, messages)

    def test_fastq2(self):
        records = [SeqRecord(Seq("CCCTTCTTGTCTTCAGCGTTTCTCC"), "EAS54_6_R1_2_1_413_324", "EAS54_6_R1_2_1_413_324"),
                   SeqRecord(Seq("TTGGCAGGCCAAGGCCGATGGATCA"), "EAS54_6_R1_2_1_540_792", "EAS54_6_R1_2_1_540_792"),
                   SeqRecord(Seq("GTTGCTTCTGGCGTGGGTGGGGGGG"), "EAS54_6_R1_2_1_443_348", "EAS54_6_R1_2_1_443_348"),
                   ]
        lengths = [25, 25, 25]
        alignment = """\
 CTG alignment column 0
 CTT alignment column 1
 CGT alignment column 2
 TGG alignment column 3
 TCC alignment column 4
 ||| ...
 CAG alignment column 24"""
        messages = ["Repeated name 'EAS54_6_R1' (originally 'EAS54_6_R1_2_1_540_792'), possibly due to truncation",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    "Repeated name 'EAS54_6_R1' (originally 'EAS54_6_R1_2_1_540_792'), possibly due to truncation",
                    ]
        self.perform_test("fastq", True, "Quality/example_dos.fastq", 3, records, lengths, alignment, messages)

    def test_fastq3(self):
        records = [SeqRecord(Seq("TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA"), "071113_EAS56_0053:1:1:998:236", "071113_EAS56_0053:1:1:998:236"),
                   SeqRecord(Seq("ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG"), "071113_EAS56_0053:1:1:182:712", "071113_EAS56_0053:1:1:182:712"),
                   SeqRecord(Seq("TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT"), "071113_EAS56_0053:1:1:153:10", "071113_EAS56_0053:1:1:153:10"),
                   SeqRecord(Seq("TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA"), "071113_EAS56_0053:1:3:990:501", "071113_EAS56_0053:1:3:990:501"),
                   ]
        lengths = [36, 36, 36, 36]
        alignment = """\
 TATT alignment column 0
 TCGG alignment column 1
 TCTG alignment column 2
 CCTG alignment column 3
 TACA alignment column 4
 |||| ...
 AGTA alignment column 35"""
        messages = ["Repeated name '071113_EAS' (originally '071113_EAS56_0053:1:1:153:10'), possibly due to truncation",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    "Repeated name '071113_EAS' (originally '071113_EAS56_0053:1:1:153:10'), possibly due to truncation",
                    ]
        self.perform_test("fastq", True, "Quality/tricky.fastq", 4, records, lengths, alignment, messages)

    def test_fastq4(self):
        records = [SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN"), "Test", "Test"),
                   ]
        lengths = [41]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fastq", False, "Quality/sanger_faked.fastq", 1, records, lengths, alignment, messages)

    def test_fastq5(self):
        records = [SeqRecord(Seq("ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG...GACTGAN"), "Test", "Test"),
                   ]
        lengths = [94]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fastq", False, "Quality/sanger_93.fastq", 1, records, lengths, alignment, messages)

    def test_fastq_illumina1(self):
        records = [SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN"), "Test", "Test"),
                   ]
        lengths = [41]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fastq-illumina", False, "Quality/illumina_faked.fastq", 1, records, lengths, alignment, messages)

    def test_fastq_solexa1(self):
        records = [SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN"), "slxa_0001_1_0001_01", "slxa_0001_1_0001_01"),
                   ]
        lengths = [46]
        alignment = None
        messages = ["Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    ]
        self.perform_test("fastq-solexa", False, "Quality/solexa_faked.fastq", 1, records, lengths, alignment, messages)

    def test_fastq_solexa2(self):
        records = [SeqRecord(Seq("GATGTGCAATACCTTTGTAGAGGAA"), "SLXA-B3_649_FC8437_R1_1_1_610_79", "SLXA-B3_649_FC8437_R1_1_1_610_79"),
                   SeqRecord(Seq("GGTTTGAGAAAGAGAAATGAGATAA"), "SLXA-B3_649_FC8437_R1_1_1_397_389", "SLXA-B3_649_FC8437_R1_1_1_397_389"),
                   SeqRecord(Seq("GAGGGTGTTGATCATGATGATGGCG"), "SLXA-B3_649_FC8437_R1_1_1_850_123", "SLXA-B3_649_FC8437_R1_1_1_850_123"),
                   SeqRecord(Seq("GTATTATTTAATGGCATACACTCAA"), "SLXA-B3_649_FC8437_R1_1_1_183_714", "SLXA-B3_649_FC8437_R1_1_1_183_714"),
                   ]
        lengths = [25, 25, 25, 25]
        alignment = """\
 GGGGG alignment column 0
 AGAGT alignment column 1
 TTGAA alignment column 2
 GTGAT alignment column 3
 TTGAT alignment column 4
 ||||| ...
 AAGGA alignment column 24"""
        messages = ["Repeated name 'SLXA-B3_64' (originally 'SLXA-B3_649_FC8437_R1_1_1_362_549'), possibly due to truncation",
                    "Need a DNA, RNA or Protein alphabet",
                    "Need a Nucleotide or Protein alphabet",
                    "Need a DNA, RNA or Protein alphabet",
                    "More than one sequence found",
                    "Need a DNA, RNA or Protein alphabet",
                    "Missing SFF flow information",
                    "Need a DNA, RNA or Protein alphabet",
                    "Repeated name 'SLXA-B3_64' (originally 'SLXA-B3_649_FC8437_R1_1_1_362_549'), possibly due to truncation",
                    ]
        self.perform_test("fastq-solexa", True, "Quality/solexa_example.fastq", 5, records, lengths, alignment, messages)

    def test_seqxml1(self):
        records = [SeqRecord(Seq("CTTTGATTCACCATTTACTGGGAGCCCACGGCGATCTGGG...TCCCTGA"), "ENSMUSG00000076441"),
                   SeqRecord(Seq("ACGTMRWSYKVHDBXN.-"), "fake1"),
                   SeqRecord(Seq("AAGGCGTTAAACCC"), "fake2"),
                   SeqRecord(Seq("G"), "minimal"),
                   ]
        lengths = [2460, 18, 14, 1]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=minimal).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=minimal).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=minimal).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=minimal).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=minimal).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("seqxml", False, "SeqXML/dna_example.xml", 4, records, lengths, alignment, messages)

    def test_seqxml2(self):
        records = [SeqRecord(Seq("UGCACUGUGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGG...CCUAAAG"), "gga-let-7a-1"),
                   SeqRecord(Seq("ACGUMRWSYKVHDBXN.-"), "fake1"),
                   SeqRecord(Seq("GGAGAU"), "fake2"),
                   SeqRecord(Seq("UGUGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCAUACCCGCAAC"), "empty description"),
                   ]
        lengths = [90, 18, 6, 48]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Cannot have spaces in EMBL accession, 'empty description'",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=empty description).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=empty description).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=empty description).",
                    "Invalid whitespace in 'empty description' for LOCUS line",
                    "Cannot have spaces in EMBL accession, 'empty description'",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=empty description).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=empty description).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("seqxml", False, "SeqXML/rna_example.xml", 5, records, lengths, alignment, messages)

    def test_seqxml3(self):
        records = [SeqRecord(Seq("MSSKGSVVLAYSGGLDTSCILVWLKEQGYDVIAYLANIGQ...QSKVTAK"), "ENSMUSP00000099904"),
                   SeqRecord(Seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ.-*"), "fake1"),
                   SeqRecord(Seq("GAKKVFIEDVSKEFVEEFIWPAVQSSALYE"), "fake2"),
                   SeqRecord(Seq("PPPWAKKVFIEDVIIAGSKEFVEEFIWPAVQSE"), "UniprotProtein"),
                   ]
        lengths = [412, 29, 30, 33]
        alignment = None
        messages = ["Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=UniprotProtein).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=UniprotProtein).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=UniprotProtein).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=UniprotProtein).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=UniprotProtein).",
                    "Missing SFF flow information",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    "Sequences must all be the same length",
                    ]
        self.perform_test("seqxml", False, "SeqXML/protein_example.xml", 5, records, lengths, alignment, messages)

    def test_abi1(self):
        records = [SeqRecord(Seq("TGATNTTNACNNTTTTGAANCANTGAGTTAATAGCAATNC...NNNNNNG"), "D11F", "310"),
                   ]
        lengths = [868]
        alignment = None
        messages = ["Missing SFF flow information",
                    ]
        self.perform_test("abi", False, "Abi/310.ab1", 1, records, lengths, alignment, messages)

    def test_abi2(self):
        records = [SeqRecord(Seq("CAAGATTGCATTCATGATCTACGATTACTAGCGATTCCAG...CCTTTTA"), "16S_S2_1387R", "3100"),
                   ]
        lengths = [795]
        alignment = None
        messages = ["Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "Missing SFF flow information",
                    ]
        self.perform_test("abi", False, "Abi/3100.ab1", 1, records, lengths, alignment, messages)

    def test_abi3(self):
        records = [SeqRecord(Seq("GGGCGAGCKYYAYATTTTGGCAAGAATTGAGCTCTATGGC...ACCTTTC"), "226032_C-ME-18_pCAGseqF", "3730"),
                   ]
        lengths = [1165]
        alignment = None
        messages = ["Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "Missing SFF flow information",
                    ]
        self.perform_test("abi", False, "Abi/3730.ab1", 1, records, lengths, alignment, messages)

    def test_pdb_atom1(self):
        records = [SeqRecord(Seq("MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLL...MMTACQG"), "1A8O:A"),
                   ]
        lengths = [70]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "source should be of type string",
                    "Missing SFF flow information",
                    ]
        self.perform_test("pdb-atom", False, "PDB/1A8O.pdb", 1, records, lengths, alignment, messages)

    def test_pdb_atom2(self):
        records = [SeqRecord(Seq("LVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:A"),
                   SeqRecord(Seq("LVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:B"),
                   SeqRecord(Seq("LVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:C"),
                   SeqRecord(Seq("LVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:E"),
                   ]
        lengths = [26, 26, 26, 26]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "source should be of type string",
                    "Missing SFF flow information",
                    ]
        self.perform_test("pdb-atom", False, "PDB/2BEG.pdb", 5, records, lengths, alignment, messages)

    def test_pdb_atom3(self):
        records = [SeqRecord(Seq("MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEA...LNYIPNR"), "????:A"),
                   ]
        lengths = [51]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=????:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=????:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=????:A).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=????:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=????:A).",
                    "source should be of type string",
                    "Missing SFF flow information",
                    ]
        self.perform_test("pdb-atom", False, "PDB/1LCD.pdb", 1, records, lengths, alignment, messages)

    def test_pdb_seqres1(self):
        records = [SeqRecord(Seq("MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLL...MMTACQG"), "1A8O:A", "1A8O:A"),
                   ]
        lengths = [70]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("pdb-seqres", False, "PDB/1A8O.pdb", 1, records, lengths, alignment, messages)

    def test_pdb_seqres2(self):
        records = [SeqRecord(Seq("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:A", "2BEG:A"),
                   SeqRecord(Seq("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:B", "2BEG:B"),
                   SeqRecord(Seq("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:C", "2BEG:C"),
                   SeqRecord(Seq("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:E", "2BEG:E"),
                   ]
        lengths = [42, 42, 42, 42]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("pdb-seqres", False, "PDB/2BEG.pdb", 5, records, lengths, alignment, messages)

    def test_cif_atom1(self):
        records = [SeqRecord(Seq("MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLL...MMTACQG"), "1A8O:A"),
                   ]
        lengths = [70]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("cif-atom", False, "PDB/1A8O.cif", 1, records, lengths, alignment, messages)

    def test_cif_atom2(self):
        records = [SeqRecord(Seq("LVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:A"),
                   SeqRecord(Seq("LVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:B"),
                   SeqRecord(Seq("LVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:C"),
                   SeqRecord(Seq("LVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:E"),
                   ]
        lengths = [26, 26, 26, 26]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("cif-atom", False, "PDB/2BEG.cif", 5, records, lengths, alignment, messages)

    def test_cif_seqres1(self):
        records = [SeqRecord(Seq("MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLL...MMTACQG"), "1A8O:A", "1A8O:A"),
                   ]
        lengths = [70]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=1A8O:A).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("cif-seqres", False, "PDB/1A8O.cif", 1, records, lengths, alignment, messages)

    def test_cif_seqres2(self):
        records = [SeqRecord(Seq("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:A", "2BEG:A"),
                   SeqRecord(Seq("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:B", "2BEG:B"),
                   SeqRecord(Seq("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:C", "2BEG:C"),
                   SeqRecord(Seq("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"), "2BEG:E", "2BEG:E"),
                   ]
        lengths = [42, 42, 42, 42]
        alignment = None
        messages = ["No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "More than one sequence found",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "No suitable quality scores found in letter_annotations of SeqRecord (id=2BEG:E).",
                    "Missing SFF flow information",
                    ]
        self.perform_test("cif-seqres", False, "PDB/2BEG.cif", 5, records, lengths, alignment, messages)

    def test_empty_file(self):
        """Check parsers can cope with an empty file."""
        for t_format in SeqIO._FormatToIterator:
            if t_format in SeqIO._BinaryFormats or \
               t_format in ("uniprot-xml", "pdb-seqres", "pdb-atom", "cif-atom"):
                # Not allowed empty SFF files.
                continue
            handle = StringIO()
            records = list(SeqIO.parse(handle, t_format))
            self.assertEqual(len(records), 0)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
