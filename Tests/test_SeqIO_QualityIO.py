# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Additional unit tests for Bio.SeqIO.QualityIO (covering FASTQ and QUAL)."""
import os
import unittest
import warnings
from Bio.Alphabet import generic_dna
from Bio.SeqIO import QualityIO
from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, MutableSeq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO
from Bio.Data.IUPACData import ambiguous_dna_letters, ambiguous_rna_letters

def truncation_expected(format) :
    if format in ["fastq-solexa", "fastq-illumina"] :
        return 62
    elif format in ["fastq", "fastq-sanger"] :
        return 93
    else :
        assert format in ["fasta", "qual"]
        return None

#Top level function as this makes it easier to use for debugging:
def write_read(filename, in_format, out_format) :
    records = list(SeqIO.parse(open(filename),in_format))
    #Write it out...
    handle = StringIO()
    SeqIO.write(records, handle, out_format)
    handle.seek(0)
    #Now load it back and check it agrees,
    records2 = list(SeqIO.parse(handle,out_format))
    compare_records(records, records2, truncation_expected(out_format))

def compare_record(old, new, truncate=None) :
    """Quality aware SeqRecord comparision.

    This will check the mapping between Solexa and PHRED scores.
    It knows to ignore UnknownSeq objects for string matching (i.e. QUAL files).
    """
    if old.id != new.id :
        raise ValueError("'%s' vs '%s' " % (old.id, new.id))
    if old.description != new.description \
    and (old.id+" "+old.description).strip() != new.description :
        raise ValueError("'%s' vs '%s' " % (old.description, new.description))
    if len(old.seq) != len(new.seq) :
        raise ValueError("%i vs %i" % (len(old.seq), len(new.seq)))
    if isinstance(old.seq, UnknownSeq) or isinstance(new.seq, UnknownSeq) :
        pass
    elif str(old.seq) != str(new.seq) :
        if len(old.seq) < 200 :
            raise ValueError("'%s' vs '%s'" % (old.seq, new.seq))
        else :
            raise ValueError("'%s...' vs '%s...'" % (old.seq[:100], new.seq[:100]))
    if "phred_quality" in old.letter_annotations \
    and "phred_quality" in new.letter_annotations \
    and old.letter_annotations["phred_quality"] != new.letter_annotations["phred_quality"] :
        if truncate and [min(q,truncate) for q in old.letter_annotations["phred_quality"]] == \
                        [min(q,truncate) for q in new.letter_annotations["phred_quality"]] :
            pass
        else :
            raise ValuerError("Mismatch in phred_quality")
    if "solexa_quality" in old.letter_annotations \
    and "solexa_quality" in new.letter_annotations \
    and old.letter_annotations["solexa_quality"] != new.letter_annotations["solexa_quality"] :
        if truncate and [min(q,truncate) for q in old.letter_annotations["solexa_quality"]] == \
                        [min(q,truncate) for q in new.letter_annotations["solexa_quality"]] :
            pass
        else :
            raise ValueError("Mismatch in phred_quality")
    if "phred_quality" in old.letter_annotations \
    and "solexa_quality" in new.letter_annotations :
        #Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
        #Assume "old" is the original, and "new" has been converted.
        converted = [round(QualityIO.solexa_quality_from_phred(q)) \
                     for q in old.letter_annotations["phred_quality"]]
        if truncate :
            converted = [min(q,truncate) for q in converted]
        if converted != new.letter_annotations["solexa_quality"] :
            print
            print old.letter_annotations["phred_quality"]
            print converted
            print new.letter_annotations["solexa_quality"]
            raise ValueError("Mismatch in phred_quality vs solexa_quality")
    if "solexa_quality" in old.letter_annotations \
    and "phred_quality" in new.letter_annotations :
        #Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
        #Assume "old" is the original, and "new" has been converted.
        converted = [round(QualityIO.phred_quality_from_solexa(q)) \
                     for q in old.letter_annotations["solexa_quality"]]
        if truncate :
            converted = [min(q,truncate) for q in converted]
        if converted != new.letter_annotations["phred_quality"] :
            print old.letter_annotations["solexa_quality"]
            print converted
            print new.letter_annotations["phred_quality"]
            raise ValueError("Mismatch in solexa_quality vs phred_quality")
    return True

def compare_records(old_list, new_list, truncate_qual=None) :
    """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list) :
        raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list) :
        if not compare_record(old,new,truncate_qual) :
            return False
    return True


class TestFastqErrors(unittest.TestCase) :
    """Test reject invalid FASTQ files."""
    def setUp(self):
        warnings.resetwarnings()

    def check_fails(self, filename, good_count, formats=None, raw=True):
        if not formats :
            formats = ["fastq-sanger", "fastq-solexa", "fastq-illumina"]
        for format in formats :
            handle = open(filename, "rU")
            records = SeqIO.parse(handle, format)
            for i in range(good_count) :
                record = records.next() #Make sure no errors!
                self.assert_(isinstance(record, SeqRecord))
            self.assertRaises(ValueError, records.next)
            handle.close()

    def check_general_fails(self, filename, good_count) :
        handle = open(filename, "rU")
        tuples = QualityIO.FastqGeneralIterator(handle)
        for i in range(good_count) :
            title, seq, qual = tuples.next() #Make sure no errors!
        self.assertRaises(ValueError, tuples.next)
        handle.close()

    def check_general_passes(self, filename, record_count) :
        handle = open(filename, "rU")
        tuples = QualityIO.FastqGeneralIterator(handle)
        #This "raw" parser doesn't check the ASCII characters which means
        #certain invalid FASTQ files will get parsed without errors.
        count = 0
        for title, seq, qual in tuples :
            self.assertEqual(len(seq), len(qual))
            count += 1
        self.assertEqual(count, record_count)
        handle.close()

    def test_space(self):
        """Reject FASTQ with spaces in seq/qual"""
        self.check_fails("Quality/error_spaces.fastq", 0)
        self.check_general_fails("Quality/error_spaces.fastq", 0)

    def test_tabs(self):
        """Reject FASTQ with tabs in seq/qual"""
        self.check_fails("Quality/error_tabs.fastq", 0)
        self.check_general_fails("Quality/error_tabs.fastq", 0)

    def test_no_qual(self):
        """Reject FASTQ with missing qualities"""
        self.check_fails("Quality/error_no_qual.fastq", 0)
        self.check_general_fails("Quality/error_no_qual.fastq", 0)

    def test_long_qual(self):
        """Reject FASTQ with longer qual than seq"""
        self.check_fails("Quality/error_long_qual.fastq", 3)
        self.check_general_fails("Quality/error_long_qual.fastq", 3)

    def test_short_qual(self):
        """Reject FASTQ with shorted qual than seq"""
        self.check_fails("Quality/error_short_qual.fastq", 2)
        self.check_general_fails("Quality/error_short_qual.fastq", 2)

    def test_diff_ids(self):
        """Reject FASTQ where + and @ identifers disagree"""
        self.check_fails("Quality/error_diff_ids.fastq", 2)
        self.check_general_fails("Quality/error_diff_ids.fastq", 2)

    def test_trunc_at_seq(self):
        """Reject FASTQ truncated at the sequence"""
        self.check_fails("Quality/error_trunc_at_seq.fastq", 4)
        self.check_general_fails("Quality/error_trunc_at_seq.fastq", 4)

    def test_trunc_at_seq(self):
        """Reject FASTQ truncated at the plus line"""
        self.check_fails("Quality/error_trunc_at_plus.fastq", 4)
        self.check_general_fails("Quality/error_trunc_at_plus.fastq", 4)

    def test_trunc_at_seq(self):
        """Reject FASTQ truncated at the quality"""
        self.check_fails("Quality/error_trunc_at_qual.fastq", 4)
        self.check_general_fails("Quality/error_trunc_at_qual.fastq", 4)

    def test_qual_null(self):
        """Reject FASTQ with null (ASCII 0) in the quality"""
        self.check_fails("Quality/error_qual_null.fastq", 0)
        self.check_general_passes("Quality/error_qual_null.fastq", 5)

    def test_qual_tab(self):
        """Reject FASTQ with tab (ASCII 9) in the quality"""
        self.check_fails("Quality/error_qual_tab.fastq", 4)
        self.check_general_passes("Quality/error_qual_tab.fastq", 5)

    def test_qual_vtab(self):
        """Reject FASTQ with vertical tab (ASCII 11) in quality"""
        self.check_fails("Quality/error_qual_vtab.fastq", 0)
        self.check_general_passes("Quality/error_qual_vtab.fastq", 5)

    def test_qual_escape(self):
        """Reject FASTQ with escape (ASCII 27) in quality"""
        self.check_fails("Quality/error_qual_escape.fastq", 4)
        self.check_general_passes("Quality/error_qual_escape.fastq", 5)

    def test_qual_unit_sep(self):
        """Reject FASTQ with unit sep (ASCII 31) in quality"""
        self.check_fails("Quality/error_qual_unit_sep.fastq", 2)
        self.check_general_passes("Quality/error_qual_unit_sep.fastq", 5)

    def test_qual_space(self):
        """Reject FASTQ with space (ASCII 32) in the quality"""
        self.check_fails("Quality/error_qual_space.fastq", 3)
        self.check_general_passes("Quality/error_qual_space.fastq", 5)

    def test_qual_del(self):
        """Reject FASTQ with delete (ASCI 127) in quality"""
        self.check_fails("Quality/error_qual_del.fastq", 3)
        self.check_general_passes("Quality/error_qual_del.fastq", 5)

class TestQual(unittest.TestCase):
    """Tests with QUAL files."""
    def setUp(self):
        warnings.resetwarnings()

    def test_paired(self):
        """Check FASTQ parsing matches FASTA+QUAL parsing"""
        records1 = list(\
            QualityIO.PairedFastaQualIterator(open("Quality/example.fasta"),
                                              open("Quality/example.qual")))
        records2 = list(SeqIO.parse(open("Quality/example.fastq"),"fastq"))
        self.assert_(compare_records(records1, records2))

    def test_qual(self):
        """Check FASTQ parsing matches QUAL parsing"""
        records1 = list(SeqIO.parse(open("Quality/example.qual"),"qual"))
        records2 = list(SeqIO.parse(open("Quality/example.fastq"),"fastq"))
        #Will ignore the unknown sequences :)
        self.assert_(compare_records(records1, records2))

    def test_qual_out(self):
        """Check FASTQ to QUAL output"""
        records = SeqIO.parse(open("Quality/example.fastq"),"fastq")
        h = StringIO("")
        SeqIO.write(records, h, "qual")
        self.assertEqual(h.getvalue(),open("Quality/example.qual").read())

    def test_fasta(self):
        """Check FASTQ parsing matches FASTA parsing"""
        records1 = list(SeqIO.parse(open("Quality/example.fasta"),"fasta"))
        records2 = list(SeqIO.parse(open("Quality/example.fastq"),"fastq"))
        self.assert_(compare_records(records1, records2))

    def test_fasta_out(self):
        """Check FASTQ to FASTA output"""
        records = SeqIO.parse(open("Quality/example.fastq"),"fastq")
        h = StringIO("")
        SeqIO.write(records, h, "fasta")
        self.assertEqual(h.getvalue(),open("Quality/example.fasta").read())


class TestReadWrite(unittest.TestCase) :
    """Test can read and write back files."""
    def setUp(self):
        warnings.resetwarnings()

    def test_fastq_2000(self) :
        """Read and write back simple example with upper case 2000bp read"""
        data = "@%s\n%s\n+\n%s\n" \
               % ("id descr goes here", "ACGT"*500, "!@a~"*500)
        handle = StringIO("")
        self.assertEqual(1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq"))
        self.assertEqual(data, handle.getvalue())

    def test_fastq_1000(self) :
        """Read and write back simple example with mixed case 1000bp read"""
        data = "@%s\n%s\n+\n%s\n" \
               % ("id descr goes here", "ACGTNncgta"*100, "abcd!!efgh"*100)
        handle = StringIO("")
        self.assertEqual(1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq"))
        self.assertEqual(data, handle.getvalue())

    def test_fastq_dna(self) :
        """Read and write back simple example with ambiguous DNA"""
        #First in upper case...        
        data = "@%s\n%s\n+\n%s\n" \
               % ("id descr goes here",
                  ambiguous_dna_letters.upper(),
                  "".join(chr(33+q) for q in range(len(ambiguous_dna_letters))))
        handle = StringIO("")
        self.assertEqual(1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq"))
        self.assertEqual(data, handle.getvalue())
        #Now in lower case...
        data = "@%s\n%s\n+\n%s\n" \
               % ("id descr goes here",
                  ambiguous_dna_letters.lower(),
                  "".join(chr(33+q) for q in range(len(ambiguous_dna_letters))))
        handle = StringIO("")
        self.assertEqual(1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq"))
        self.assertEqual(data, handle.getvalue())

    def test_fastq_rna(self) :
        """Read and write back simple example with ambiguous RNA"""
        #First in upper case...        
        data = "@%s\n%s\n+\n%s\n" \
               % ("id descr goes here",
                  ambiguous_rna_letters.upper(),
                  "".join(chr(33+q) for q in range(len(ambiguous_rna_letters))))
        handle = StringIO("")
        self.assertEqual(1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq"))
        self.assertEqual(data, handle.getvalue())
        #Now in lower case...
        data = "@%s\n%s\n+\n%s\n" \
               % ("id descr goes here",
                  ambiguous_rna_letters.lower(),
                  "".join(chr(33+q) for q in range(len(ambiguous_rna_letters))))
        handle = StringIO("")
        self.assertEqual(1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq"))
        self.assertEqual(data, handle.getvalue())


class TestWriteRead(unittest.TestCase) :
    """Test can write and read back files."""
    def setUp(self):
        warnings.resetwarnings()

    def test_generated(self) :
        """Write and read back odd SeqRecord objects"""
        record1 = SeqRecord(Seq("ACGT"*500, generic_dna),  id="Test", description="Long "*500,
                           letter_annotations={"phred_quality":[40,30,20,10]*500})
        record2 = SeqRecord(MutableSeq("NGGC"*1000),  id="Mut", description="very "*1000+"long",
                           letter_annotations={"phred_quality":[0,5,5,10]*1000})
        record3 = SeqRecord(UnknownSeq(2000,character="N"),  id="Unk", description="l"+("o"*1000)+"ng",
                           letter_annotations={"phred_quality":[0,1]*1000})
        record4 = SeqRecord(Seq("ACGT"*500),  id="no_descr", description="", name="",
                           letter_annotations={"phred_quality":[40,50,60,62]*500})
        record5 = SeqRecord(Seq("",generic_dna),  id="empty_p", description="(could have been trimmed lots)",
                           letter_annotations={"phred_quality":[]})
        record6 = SeqRecord(Seq(""),  id="empty_s", description="(could have been trimmed lots)",
                           letter_annotations={"solexa_quality":[]})
        record7 = SeqRecord(Seq("ACNN"*500),  id="Test_Sol", description="Long "*500,
                           letter_annotations={"solexa_quality":[40,30,0,-5]*500})
        record8 = SeqRecord(Seq("ACGT"),  id="HighQual", description="With very large qualities that even Sanger FASTQ can't hold!",
                           letter_annotations={"solexa_quality":[0,10,100,1000]})
        #TODO - Record with no identifier?
        records = [record1, record2, record3, record4, record5, record6, record7, record8]
        #TODO - Have a Biopython defined "DataLossWarning?"
        warnings.simplefilter('ignore', UserWarning)
        for format in ["fasta", "fastq", "fastq-solexa", "fastq-illumina", "qual"] :
            handle = StringIO()
            SeqIO.write(records, handle, format)
            handle.seek(0)
            compare_records(records,
                            list(SeqIO.parse(handle, format)),
                            truncation_expected(format))
            
    def test_tricky(self) :
        """Write and read back tricky.fastq"""
        write_read(os.path.join("Quality", "tricky.fastq"), "fastq", "fastq")
        write_read(os.path.join("Quality", "tricky.fastq"), "fastq", "fastq-sanger")
        write_read(os.path.join("Quality", "tricky.fastq"), "fastq-sanger", "fastq")
        write_read(os.path.join("Quality", "tricky.fastq"), "fastq-sanger", "fastq-sanger")
        write_read(os.path.join("Quality", "tricky.fastq"), "fastq", "qual")

    def test_sanger_93(self) :
        """Write and read back sanger_93.fastq"""
        write_read(os.path.join("Quality", "sanger_93.fastq"), "fastq-sanger", "fasta")
        write_read(os.path.join("Quality", "sanger_93.fastq"), "fastq-sanger", "fastq-sanger")
        write_read(os.path.join("Quality", "sanger_93.fastq"), "fastq-sanger", "qual")
        #TODO - Have a Biopython defined "DataLossWarning?"
        #TODO - On Python 2.6+ we can check this warning is really triggered
        warnings.simplefilter('ignore', UserWarning)
        write_read(os.path.join("Quality", "sanger_93.fastq"), "fastq-sanger", "fastq-solexa")
        write_read(os.path.join("Quality", "sanger_93.fastq"), "fastq-sanger", "fastq-illumina")

    def test_sanger_faked(self) :
        """Write and read back sanger_faked.fastq"""
        write_read(os.path.join("Quality", "sanger_faked.fastq"), "fastq-sanger", "fasta")
        write_read(os.path.join("Quality", "sanger_faked.fastq"), "fastq-sanger", "fastq-sanger")
        write_read(os.path.join("Quality", "sanger_faked.fastq"), "fastq-sanger", "fastq-solexa")
        write_read(os.path.join("Quality", "sanger_faked.fastq"), "fastq-sanger", "fastq-illumina")
        write_read(os.path.join("Quality", "sanger_faked.fastq"), "fastq-sanger", "qual")

    def test_example_fasta(self) :
        """Write and read back example.fasta"""
        write_read(os.path.join("Quality", "example.fasta"), "fasta", "fasta")
        #TODO - tests to check can't write FASTQ or QUAL...

    def test_example_fastq(self) :
        """Write and read back example.fastq"""
        write_read(os.path.join("Quality", "example.fastq"), "fastq-sanger", "fasta")
        write_read(os.path.join("Quality", "example.fastq"), "fastq-sanger", "fastq-sanger")
        write_read(os.path.join("Quality", "example.fastq"), "fastq-sanger", "fastq-solexa")
        write_read(os.path.join("Quality", "example.fastq"), "fastq-sanger", "fastq-illumina")
        write_read(os.path.join("Quality", "example.fastq"), "fastq-sanger", "qual")

    def test_example_qual(self) :
        """Write and read back example.qual"""
        write_read(os.path.join("Quality", "example.qual"), "qual", "fasta")
        write_read(os.path.join("Quality", "example.qual"), "qual", "qual")
        write_read(os.path.join("Quality", "example.qual"), "qual", "fastq")
        write_read(os.path.join("Quality", "example.qual"), "qual", "fastq-sanger")
        write_read(os.path.join("Quality", "example.qual"), "qual", "fastq-solexa")
        write_read(os.path.join("Quality", "example.qual"), "qual", "fastq-illumina")

    def test_solexa_faked(self) :
        """Write and read back solexa_faked.fastq"""
        write_read(os.path.join("Quality", "solexa_faked.fastq"), "fastq-solexa", "fasta")
        write_read(os.path.join("Quality", "solexa_faked.fastq"), "fastq-solexa", "fastq-sanger")
        write_read(os.path.join("Quality", "solexa_faked.fastq"), "fastq-solexa", "fastq-solexa")
        write_read(os.path.join("Quality", "solexa_faked.fastq"), "fastq-solexa", "fastq-illumina")
        write_read(os.path.join("Quality", "solexa_faked.fastq"), "fastq-solexa", "qual")

    def test_solexa_example(self) :
        """Write and read back solexa_example.fastq"""
        write_read(os.path.join("Quality", "solexa_example.fastq"), "fastq-solexa", "fasta")
        write_read(os.path.join("Quality", "solexa_example.fastq"), "fastq-solexa", "fastq-sanger")
        write_read(os.path.join("Quality", "solexa_example.fastq"), "fastq-solexa", "fastq-solexa")
        write_read(os.path.join("Quality", "solexa_example.fastq"), "fastq-solexa", "fastq-illumina")
        write_read(os.path.join("Quality", "solexa_example.fastq"), "fastq-solexa", "qual")

    def test_illumina_faked(self) :
        """Write and read back illumina_faked.fastq"""
        write_read(os.path.join("Quality", "illumina_faked.fastq"), "fastq-illumina", "fasta")
        write_read(os.path.join("Quality", "illumina_faked.fastq"), "fastq-illumina", "fastq-sanger")
        write_read(os.path.join("Quality", "illumina_faked.fastq"), "fastq-illumina", "fastq-solexa")
        write_read(os.path.join("Quality", "illumina_faked.fastq"), "fastq-illumina", "fastq-illumina")
        write_read(os.path.join("Quality", "illumina_faked.fastq"), "fastq-illumina", "qual")

class MappingTests(unittest.TestCase) :
    def setUp(self):
        warnings.resetwarnings()

    def test_solexa_quality_from_phred(self):
        """Mapping check for function solexa_quality_from_phred"""
        self.assertEqual(-5, round(QualityIO.solexa_quality_from_phred(0)))
        self.assertEqual(-5, round(QualityIO.solexa_quality_from_phred(1)))
        self.assertEqual(-2, round(QualityIO.solexa_quality_from_phred(2)))
        self.assertEqual(0, round(QualityIO.solexa_quality_from_phred(3)))
        self.assertEqual(2, round(QualityIO.solexa_quality_from_phred(4)))
        self.assertEqual(3, round(QualityIO.solexa_quality_from_phred(5)))
        self.assertEqual(5, round(QualityIO.solexa_quality_from_phred(6)))
        self.assertEqual(6, round(QualityIO.solexa_quality_from_phred(7)))
        self.assertEqual(7, round(QualityIO.solexa_quality_from_phred(8)))
        self.assertEqual(8, round(QualityIO.solexa_quality_from_phred(9)))
        for i in range(10,100) :
            self.assertEqual(i, round(QualityIO.solexa_quality_from_phred(i)))
        
    def test_phred_quality_from_solexa(self):
        """Mapping check for function phred_quality_from_solexa"""
        self.assertEqual(1, round(QualityIO.phred_quality_from_solexa(-5)))
        self.assertEqual(1, round(QualityIO.phred_quality_from_solexa(-4)))
        self.assertEqual(2, round(QualityIO.phred_quality_from_solexa(-3)))
        self.assertEqual(2, round(QualityIO.phred_quality_from_solexa(-2)))
        self.assertEqual(3, round(QualityIO.phred_quality_from_solexa(-1)))
        self.assertEqual(3, round(QualityIO.phred_quality_from_solexa(0)))
        self.assertEqual(4, round(QualityIO.phred_quality_from_solexa(1)))
        self.assertEqual(4, round(QualityIO.phred_quality_from_solexa(2)))
        self.assertEqual(5, round(QualityIO.phred_quality_from_solexa(3)))
        self.assertEqual(5, round(QualityIO.phred_quality_from_solexa(4)))
        self.assertEqual(6, round(QualityIO.phred_quality_from_solexa(5)))
        self.assertEqual(7, round(QualityIO.phred_quality_from_solexa(6)))
        self.assertEqual(8, round(QualityIO.phred_quality_from_solexa(7)))
        self.assertEqual(9, round(QualityIO.phred_quality_from_solexa(8)))
        self.assertEqual(10, round(QualityIO.phred_quality_from_solexa(9)))
        for i in range(10,100) :
            self.assertEqual(i, round(QualityIO.phred_quality_from_solexa(i)))

    def test_sanger_to_solexa(self):
        """Mapping check for FASTQ Sanger (0 to 93) to Solexa (-5 to 62)"""
        #The point of this test is the writing code doesn't actually use the
        #solexa_quality_from_phred function directly. For speed it uses a
        #cached dictionary of the mappings.
        seq = "N"*94
        qual = "".join(chr(33+q) for q in range(0,94))
        expected_sol = [min(62,int(round(QualityIO.solexa_quality_from_phred(q)))) \
                        for q in range(0,94)]
        in_handle = StringIO("@Test\n%s\n+\n%s" % (seq,qual))
        out_handle = StringIO("")
        #Want to ignore the data loss warning
        #(on Python 2.6 we could check for it!)
        warnings.simplefilter('ignore', UserWarning)
        SeqIO.write(SeqIO.parse(in_handle, "fastq-sanger"),
                    out_handle, "fastq-solexa")
        warnings.resetwarnings()
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-solexa")
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.letter_annotations["solexa_quality"],
                         expected_sol)

    def test_solexa_to_sanger(self):
        """Mapping check for FASTQ Solexa (-5 to 62) to Sanger (0 to 62)"""
        #The point of this test is the writing code doesn't actually use the
        #solexa_quality_from_phred function directly. For speed it uses a
        #cached dictionary of the mappings.
        seq = "N"*68
        qual = "".join(chr(64+q) for q in range(-5,63))
        expected_phred = [round(QualityIO.phred_quality_from_solexa(q)) \
                          for q in range(-5,63)]
        in_handle = StringIO("@Test\n%s\n+\n%s" % (seq,qual))
        out_handle = StringIO("")
        #Want to ignore the data loss warning
        #(on Python 2.6 we could check for it!)
        warnings.simplefilter('ignore', UserWarning)
        SeqIO.write(SeqIO.parse(in_handle, "fastq-solexa"),
                    out_handle, "fastq-sanger")
        warnings.resetwarnings()
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-sanger")
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.letter_annotations["phred_quality"],
                         expected_phred)

    def test_sanger_to_illumina(self):
        """Mapping check for FASTQ Sanger (0 to 93) to Illumina (0 to 62)"""
        seq = "N"*94
        qual = "".join(chr(33+q) for q in range(0,94))
        expected_phred = [min(62,q) for q in range(0,94)]
        in_handle = StringIO("@Test\n%s\n+\n%s" % (seq,qual))
        out_handle = StringIO("")
        #Want to ignore the data loss warning
        #(on Python 2.6 we could check for it!)
        warnings.simplefilter('ignore', UserWarning)
        SeqIO.write(SeqIO.parse(in_handle, "fastq-sanger"),
                    out_handle, "fastq-illumina")
        warnings.resetwarnings()
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-illumina")
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.letter_annotations["phred_quality"],
                         expected_phred)

    def test_illumina_to_sanger(self):
        """Mapping check for FASTQ Illumina (0 to 62) to Sanger (0 to 62)"""
        seq = "N"*63
        qual = "".join(chr(64+q) for q in range(0,63))
        expected_phred = range(63)
        in_handle = StringIO("@Test\n%s\n+\n%s" % (seq,qual))
        out_handle = StringIO("")
        SeqIO.write(SeqIO.parse(in_handle, "fastq-illumina"),
                    out_handle, "fastq-sanger")
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-sanger")
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.letter_annotations["phred_quality"],
                         expected_phred)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
