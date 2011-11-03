# Copyright 2009-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Additional unit tests for Bio.SeqIO.QualityIO (covering FASTQ and QUAL)."""
import os
import unittest
import warnings

from StringIO import StringIO
try:
    #This is in Python 2.6+, but we need it on Python 3
    from io import BytesIO
except ImportError:
    BytesIO = StringIO

from Bio import BiopythonWarning
from Bio.Alphabet import generic_dna
from Bio.SeqIO import QualityIO
from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Data.IUPACData import ambiguous_dna_letters, ambiguous_rna_letters

BINARY_FORMATS = ["sff", "sff-trim"]

def truncation_expected(format):
    if format in ["fastq-solexa", "fastq-illumina"] :
        return 62
    elif format in ["fastq", "fastq-sanger"]:
        return 93
    else:
        assert format in ["fasta", "qual", "phd", "sff"]
        return None

#Top level function as this makes it easier to use for debugging:
def write_read(filename, in_format, out_format):
    if in_format in BINARY_FORMATS:
        mode = "rb"
    else:
        mode = "r"
    records = list(SeqIO.parse(open(filename, mode),in_format))
    #Write it out...
    if out_format in BINARY_FORMATS:
        handle = BytesIO()
    else :
        handle = StringIO()
    SeqIO.write(records, handle, out_format)
    handle.seek(0)
    #Now load it back and check it agrees,
    records2 = list(SeqIO.parse(handle,out_format))
    compare_records(records, records2, truncation_expected(out_format))

def compare_record(old, new, truncate=None):
    """Quality aware SeqRecord comparision.

    This will check the mapping between Solexa and PHRED scores.
    It knows to ignore UnknownSeq objects for string matching (i.e. QUAL files).
    """
    if old.id != new.id:
        raise ValueError("'%s' vs '%s' " % (old.id, new.id))
    if old.description != new.description \
    and (old.id+" "+old.description).strip() != new.description:
        raise ValueError("'%s' vs '%s' " % (old.description, new.description))
    if len(old.seq) != len(new.seq):
        raise ValueError("%i vs %i" % (len(old.seq), len(new.seq)))
    if isinstance(old.seq, UnknownSeq) or isinstance(new.seq, UnknownSeq):
        pass
    elif str(old.seq) != str(new.seq):
        if len(old.seq) < 200:
            raise ValueError("'%s' vs '%s'" % (old.seq, new.seq))
        else:
            raise ValueError("'%s...' vs '%s...'" % (old.seq[:100], new.seq[:100]))
    if "phred_quality" in old.letter_annotations \
    and "phred_quality" in new.letter_annotations \
    and old.letter_annotations["phred_quality"] != new.letter_annotations["phred_quality"]:
        if truncate and [min(q,truncate) for q in old.letter_annotations["phred_quality"]] == \
                        [min(q,truncate) for q in new.letter_annotations["phred_quality"]]:
            pass
        else:
            raise ValuerError("Mismatch in phred_quality")
    if "solexa_quality" in old.letter_annotations \
    and "solexa_quality" in new.letter_annotations \
    and old.letter_annotations["solexa_quality"] != new.letter_annotations["solexa_quality"]:
        if truncate and [min(q,truncate) for q in old.letter_annotations["solexa_quality"]] == \
                        [min(q,truncate) for q in new.letter_annotations["solexa_quality"]]:
            pass
        else:
            raise ValueError("Mismatch in phred_quality")
    if "phred_quality" in old.letter_annotations \
    and "solexa_quality" in new.letter_annotations:
        #Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
        #Assume "old" is the original, and "new" has been converted.
        converted = [round(QualityIO.solexa_quality_from_phred(q)) \
                     for q in old.letter_annotations["phred_quality"]]
        if truncate:
            converted = [min(q,truncate) for q in converted]
        if converted != new.letter_annotations["solexa_quality"]:
            print
            print old.letter_annotations["phred_quality"]
            print converted
            print new.letter_annotations["solexa_quality"]
            raise ValueError("Mismatch in phred_quality vs solexa_quality")
    if "solexa_quality" in old.letter_annotations \
    and "phred_quality" in new.letter_annotations:
        #Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
        #Assume "old" is the original, and "new" has been converted.
        converted = [round(QualityIO.phred_quality_from_solexa(q)) \
                     for q in old.letter_annotations["solexa_quality"]]
        if truncate:
            converted = [min(q,truncate) for q in converted]
        if converted != new.letter_annotations["phred_quality"]:
            print old.letter_annotations["solexa_quality"]
            print converted
            print new.letter_annotations["phred_quality"]
            raise ValueError("Mismatch in solexa_quality vs phred_quality")
    return True

def compare_records(old_list, new_list, truncate_qual=None):
    """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list):
        raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list):
        if not compare_record(old,new,truncate_qual):
            return False
    return True


class TestFastqErrors(unittest.TestCase):
    """Test reject invalid FASTQ files."""
    def check_fails(self, filename, good_count, formats=None, raw=True):
        if not formats:
            formats = ["fastq-sanger", "fastq-solexa", "fastq-illumina"]
        for format in formats:
            handle = open(filename, "rU")
            records = SeqIO.parse(handle, format)
            for i in range(good_count):
                record = records.next() #Make sure no errors!
                self.assertTrue(isinstance(record, SeqRecord))
            self.assertRaises(ValueError, records.next)
            handle.close()

    def check_general_fails(self, filename, good_count):
        handle = open(filename, "rU")
        tuples = QualityIO.FastqGeneralIterator(handle)
        for i in range(good_count):
            title, seq, qual = tuples.next() #Make sure no errors!
        self.assertRaises(ValueError, tuples.next)
        handle.close()

    def check_general_passes(self, filename, record_count):
        handle = open(filename, "rU")
        tuples = QualityIO.FastqGeneralIterator(handle)
        #This "raw" parser doesn't check the ASCII characters which means
        #certain invalid FASTQ files will get parsed without errors.
        count = 0
        for title, seq, qual in tuples:
            self.assertEqual(len(seq), len(qual))
            count += 1
        self.assertEqual(count, record_count)
        handle.close()

    def check_all_fail(self, filename, count):
        self.check_fails(filename, count)
        self.check_general_fails(filename, count)

    def check_qual_char(self, filename, good_count, count):
        self.check_fails(filename, good_count)
        self.check_general_passes(filename, count)

#Now add methods at run time... these FASTQ files will be rejected
#by both the low level parser AND the high level SeqRecord parser:
tests = [("diff_ids", 2),
         ("no_qual", 0),
         ("long_qual", 3),
         ("short_qual", 2),
         ("double_seq", 3),
         ("double_qual", 2),
         ("tabs", 0),
         ("spaces", 0),
         ("trunc_in_title", 4),
         ("trunc_in_seq", 4),
         ("trunc_in_plus", 4),
         ("trunc_in_qual", 4),
         ("trunc_at_seq", 4),
         ("trunc_at_plus", 4),
         ("trunc_at_qual", 4)]
for base_name, good_count in tests:
    def funct(name,c):
        f = lambda x : x.check_all_fail("Quality/error_%s.fastq" % name,c)
        f.__doc__ = "Reject FASTQ with %s" % name.replace("_"," ")
        return f
    setattr(TestFastqErrors, "test_%s" % (base_name),
            funct(base_name, good_count))
    del funct        

#Now add methods for FASTQ files which will be rejected by the high
#level SeqRecord parser, but will be accepted by the low level parser:
tests = [("del", 3, 5),
         ("space", 3, 5),
         ("vtab", 0, 5),
         ("escape", 4, 5),
         ("unit_sep", 2, 5),
         ("tab", 4, 5),
         ("null", 0, 5)]
for base_name, good_count, full_count in tests:
    def funct(name,c1,c2):
        f = lambda x : x.check_qual_char("Quality/error_qual_%s.fastq"%name,c1,c2)
        f.__doc__ = "Reject FASTQ with %s in quality" % name.replace("_"," ")
        return f
    setattr(TestFastqErrors, "test_qual_%s" % (base_name),
            funct(base_name, good_count, full_count))
    del funct        


class TestReferenceSffConversions(unittest.TestCase):
    def check(self, sff_name, sff_format, out_name, format) :
        wanted = list(SeqIO.parse(open(out_name), format))
        data = StringIO()
        count = SeqIO.convert(sff_name, sff_format, data, format)
        self.assertEqual(count, len(wanted))
        data.seek(0)
        converted = list(SeqIO.parse(data, format))
        self.assertEqual(len(wanted), len(converted))
        for old, new in zip(wanted, converted) :
            self.assertEqual(old.id, new.id)
            self.assertEqual(old.name, new.name)
            if format!="qual" :
                self.assertEqual(str(old.seq), str(new.seq))
            elif format!="fasta" :
                self.assertEqual(old.letter_annotations["phred_quality"],
                                 new.letter_annotations["phred_quality"])

    def check_sff(self, sff_name):
        self.check(sff_name, "sff", "Roche/E3MFGYR02_random_10_reads_no_trim.fasta", "fasta")
        self.check(sff_name, "sff", "Roche/E3MFGYR02_random_10_reads_no_trim.qual", "qual")
        self.check(sff_name, "sff-trim", "Roche/E3MFGYR02_random_10_reads.fasta", "fasta")
        self.check(sff_name, "sff-trim", "Roche/E3MFGYR02_random_10_reads.qual", "qual")

    def test_original(self) :
        """Test converting E3MFGYR02_random_10_reads.sff into FASTA+QUAL"""
        self.check_sff("Roche/E3MFGYR02_random_10_reads.sff")
        
    def test_no_manifest(self) :
        """Test converting E3MFGYR02_no_manifest.sff into FASTA+QUAL"""
        self.check_sff("Roche/E3MFGYR02_no_manifest.sff")
        
    def test_alt_index_at_start(self) :
        """Test converting E3MFGYR02_alt_index_at_start into FASTA+QUAL"""
        self.check_sff("Roche/E3MFGYR02_alt_index_at_start.sff")

    def test_alt_index_in_middle(self) :
        """Test converting E3MFGYR02_alt_index_in_middle into FASTA+QUAL"""
        self.check_sff("Roche/E3MFGYR02_alt_index_in_middle.sff")

    def test_alt_index_at_end(self) :
        """Test converting E3MFGYR02_alt_index_at_end into FASTA+QUAL"""
        self.check_sff("Roche/E3MFGYR02_alt_index_at_end.sff")

    def test_index_at_start(self) :
        """Test converting E3MFGYR02_index_at_start into FASTA+QUAL"""
        self.check_sff("Roche/E3MFGYR02_index_at_start.sff")

    def test_index_at_end(self) :
        """Test converting E3MFGYR02_index_in_middle into FASTA+QUAL"""
        self.check_sff("Roche/E3MFGYR02_index_in_middle.sff")

class TestReferenceFastqConversions(unittest.TestCase):
    """Tests where we have reference output."""
    def simple_check(self, base_name, in_variant):
        for out_variant in ["sanger", "solexa", "illumina"]:
            if out_variant != "sanger":
                #Ignore data loss warnings from max qualities
                warnings.simplefilter('ignore', BiopythonWarning)
            in_filename = "Quality/%s_original_%s.fastq" \
                          % (base_name, in_variant)
            self.assertTrue(os.path.isfile(in_filename))
            #Load the reference output...  
            expected = open("Quality/%s_as_%s.fastq" \
                            % (base_name, out_variant),
                            "rU").read()
            #Check matches using convert...
            handle = StringIO()
            SeqIO.convert(in_filename, "fastq-"+in_variant,
                          handle, "fastq-"+out_variant)
            self.assertEqual(expected, handle.getvalue())
            #Check matches using parse/write
            handle = StringIO()
            SeqIO.write(SeqIO.parse(open(in_filename), "fastq-"+in_variant),
                        handle, "fastq-"+out_variant)
            self.assertEqual(expected, handle.getvalue())
            if out_variant != "sanger":
                warnings.filters.pop()

#Now add methods at run time...
tests = [("illumina_full_range", "illumina"),
         ("sanger_full_range", "sanger"),
         ("longreads", "sanger"),
         ("solexa_full_range", "solexa"),
         ("misc_dna", "sanger"),
         ("wrapping", "sanger"),
         ("misc_rna", "sanger")]
for base_name, variant in tests:
    assert variant in ["sanger", "solexa", "illumina"]
    def funct(bn,var):
        f = lambda x : x.simple_check(bn,var)
        f.__doc__ = "Reference conversions of %s file %s" % (var, bn)
        return f
    setattr(TestReferenceFastqConversions, "test_%s_%s" % (base_name, variant),
            funct(base_name, variant))
    del funct        

class TestQual(unittest.TestCase):
    """Tests with QUAL files."""
    def test_paired(self):
        """Check FASTQ parsing matches FASTA+QUAL parsing"""
        records1 = list(\
            QualityIO.PairedFastaQualIterator(open("Quality/example.fasta"),
                                              open("Quality/example.qual")))
        records2 = list(SeqIO.parse(open("Quality/example.fastq"),"fastq"))
        self.assertTrue(compare_records(records1, records2))

    def test_qual(self):
        """Check FASTQ parsing matches QUAL parsing"""
        records1 = list(SeqIO.parse(open("Quality/example.qual"),"qual"))
        records2 = list(SeqIO.parse(open("Quality/example.fastq"),"fastq"))
        #Will ignore the unknown sequences :)
        self.assertTrue(compare_records(records1, records2))

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
        self.assertTrue(compare_records(records1, records2))

    def test_fasta_out(self):
        """Check FASTQ to FASTA output"""
        records = SeqIO.parse(open("Quality/example.fastq"),"fastq")
        h = StringIO("")
        SeqIO.write(records, h, "fasta")
        self.assertEqual(h.getvalue(),open("Quality/example.fasta").read())

    def test_qual_negative(self):
        """Check QUAL negative scores mapped to PHRED zero"""
        data = """>1117_10_107_F3
23 31 -1 -1 -1 29 -1 -1 20 32 -1 18 25 7 -1 6 -1 -1 -1 30 -1 20 13 7 -1 -1 21 30 -1 24 -1 22 -1 -1 22 14 -1 12 26 21 -1 5 -1 -1 -1 20 -1 -1 12 28 
>1117_10_146_F3
20 33 -1 -1 -1 29 -1 -1 28 28 -1 7 16 5 -1 30 -1 -1 -1 14 -1 4 13 4 -1 -1 11 13 -1 5 -1 7 -1 -1 10 16 -1 4 12 15 -1 8 -1 -1 -1 16 -1 -1 10 4 
>1117_10_1017_F3
33 33 -1 -1 -1 27 -1 -1 17 16 -1 28 24 11 -1 6 -1 -1 -1 29 -1 8 29 24 -1 -1 8 8 -1 20 -1 13 -1 -1 8 13 -1 28 10 24 -1 10 -1 -1 -1 4 -1 -1 7 6 
>1117_11_136_F3
16 22 -1 -1 -1 33 -1 -1 30 27 -1 27 28 32 -1 29 -1 -1 -1 27 -1 18 9 6 -1 -1 23 16 -1 26 -1 5 7 -1 22 7 -1 18 14 8 -1 8 -1 -1 -1 11 -1 -1 4 24"""
        h = StringIO(data)
        h2 = StringIO()
        self.assertEqual(4, SeqIO.convert(h, "qual", h2, "fastq"))
        self.assertEqual(h2.getvalue(), """@1117_10_107_F3
??????????????????????????????????????????????????
+
8@!!!>!!5A!3:(!'!!!?!5.(!!6?!9!7!!7/!-;6!&!!!5!!-=
@1117_10_146_F3
??????????????????????????????????????????????????
+
5B!!!>!!==!(1&!?!!!/!%.%!!,.!&!(!!+1!%-0!)!!!1!!+%
@1117_10_1017_F3
??????????????????????????????????????????????????
+
BB!!!<!!21!=9,!'!!!>!)>9!!))!5!.!!).!=+9!+!!!%!!('
@1117_11_136_F3
??????????????????????????????????????????????????
+
17!!!B!!?<!<=A!>!!!<!3*'!!81!;!&(!7(!3/)!)!!!,!!%9
""")



class TestReadWrite(unittest.TestCase):
    """Test can read and write back files."""
    def test_fastq_2000(self):
        """Read and write back simple example with upper case 2000bp read"""
        data = "@%s\n%s\n+\n%s\n" \
               % ("id descr goes here", "ACGT"*500, "!@a~"*500)
        handle = StringIO("")
        self.assertEqual(1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq"))
        self.assertEqual(data, handle.getvalue())

    def test_fastq_1000(self):
        """Read and write back simple example with mixed case 1000bp read"""
        data = "@%s\n%s\n+\n%s\n" \
               % ("id descr goes here", "ACGTNncgta"*100, "abcd!!efgh"*100)
        handle = StringIO("")
        self.assertEqual(1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq"))
        self.assertEqual(data, handle.getvalue())

    def test_fastq_dna(self):
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

    def test_fastq_rna(self):
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


class TestWriteRead(unittest.TestCase):
    """Test can write and read back files."""
    def test_generated(self):
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
        warnings.simplefilter('ignore', BiopythonWarning)
        #TODO - Include phd output?
        for format in ["fasta", "fastq", "fastq-solexa", "fastq-illumina", "qual"]:
            handle = StringIO()
            SeqIO.write(records, handle, format)
            handle.seek(0)
            compare_records(records,
                            list(SeqIO.parse(handle, format)),
                            truncation_expected(format))
        warnings.filters.pop()
            
    def check(self, filename, format, out_formats):
        for f in out_formats:
            write_read(filename, format, f)

    def test_tricky(self):
        """Write and read back tricky.fastq"""
        self.check(os.path.join("Quality", "tricky.fastq"), "fastq",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd"])

    def test_sanger_93(self):
        """Write and read back sanger_93.fastq"""
        self.check(os.path.join("Quality", "sanger_93.fastq"), "fastq",
                   ["fastq", "fastq-sanger", "fasta", "qual", "phd"])
        #TODO - Have a Biopython defined "DataLossWarning?"
        #TODO - On Python 2.6+ we can check this warning is really triggered
        warnings.simplefilter('ignore', BiopythonWarning)
        self.check(os.path.join("Quality", "sanger_93.fastq"), "fastq",
                   ["fastq-solexa","fastq-illumina"])
        warnings.filters.pop()

    def test_sanger_faked(self):
        """Write and read back sanger_faked.fastq"""
        self.check(os.path.join("Quality", "sanger_faked.fastq"), "fastq",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd"])

    def test_example_fasta(self):
        """Write and read back example.fasta"""
        write_read(os.path.join("Quality", "example.fasta"), "fasta", "fasta")
        #TODO - tests to check can't write FASTQ or QUAL...

    def test_example_fastq(self):
        """Write and read back example.fastq"""
        self.check(os.path.join("Quality", "example.fastq"), "fastq",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd"])

    def test_example_qual(self):
        """Write and read back example.qual"""
        self.check(os.path.join("Quality", "example.qual"), "qual",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd"])

    def test_solexa_faked(self):
        """Write and read back solexa_faked.fastq"""
        self.check(os.path.join("Quality", "solexa_faked.fastq"), "fastq-solexa",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd"])

    def test_solexa_example(self):
        """Write and read back solexa_example.fastq"""
        self.check(os.path.join("Quality", "solexa_example.fastq"), "fastq-solexa",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd"])

    def test_illumina_faked(self):
        """Write and read back illumina_faked.fastq"""
        self.check(os.path.join("Quality", "illumina_faked.fastq"), "fastq-illumina",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd"])

    def test_greek_sff(self) :
        """Write and read back greek.sff"""
        self.check(os.path.join("Roche", "greek.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_paired_sff(self) :
        """Write and read back paired.sff"""
        self.check(os.path.join("Roche", "paired.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_E3MFGYR02(self) :
        """Write and read back E3MFGYR02_random_10_reads.sff"""
        self.check(os.path.join("Roche", "E3MFGYR02_random_10_reads.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_E3MFGYR02_no_manifest(self) :
        """Write and read back E3MFGYR02_no_manifest.sff"""
        self.check(os.path.join("Roche", "E3MFGYR02_no_manifest.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_E3MFGYR02_index_at_start(self) :
        """Write and read back E3MFGYR02_index_at_start.sff"""
        self.check(os.path.join("Roche", "E3MFGYR02_index_at_start.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_E3MFGYR02_index_in_middle(self) :
        """Write and read back E3MFGYR02_index_in_middle.sff"""
        self.check(os.path.join("Roche", "E3MFGYR02_index_in_middle.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_E3MFGYR02_alt_index_at_start(self) :
        """Write and read back E3MFGYR02_alt_index_at_start.sff"""
        self.check(os.path.join("Roche", "E3MFGYR02_alt_index_at_start.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_E3MFGYR02_alt_index_in_middle(self) :
        """Write and read back E3MFGYR02_alt_index_in_middle.sff"""
        self.check(os.path.join("Roche", "E3MFGYR02_alt_index_in_middle.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_E3MFGYR02_alt_index_at_end(self) :
        """Write and read back E3MFGYR02_alt_index_at_end.sff"""
        self.check(os.path.join("Roche", "E3MFGYR02_alt_index_at_end.sff"), "sff",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd", "sff"])

    def test_E3MFGYR02_trimmed(self) :
        """Write and read back E3MFGYR02_random_10_reads.sff (trimmed)"""
        self.check(os.path.join("Roche", "E3MFGYR02_random_10_reads.sff"), "sff-trim",
                   ["fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa",
                    "fasta", "qual", "phd"]) #not sff as output


class MappingTests(unittest.TestCase):
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
        for i in range(10,100):
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
        for i in range(10,100):
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
        warnings.simplefilter('ignore', BiopythonWarning)
        SeqIO.write(SeqIO.parse(in_handle, "fastq-sanger"),
                    out_handle, "fastq-solexa")
        warnings.filters.pop()
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
        warnings.simplefilter('ignore', BiopythonWarning)
        SeqIO.write(SeqIO.parse(in_handle, "fastq-solexa"),
                    out_handle, "fastq-sanger")
        warnings.filters.pop()
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
        warnings.simplefilter('ignore', BiopythonWarning)
        SeqIO.write(SeqIO.parse(in_handle, "fastq-sanger"),
                    out_handle, "fastq-illumina")
        warnings.filters.pop()
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
