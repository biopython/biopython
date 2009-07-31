# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Additional unit tests for Bio.SeqIO.QualityIO (covering FASTQ and QUAL)."""
import os
import unittest
from Bio.Alphabet import generic_dna
from Bio.SeqIO import QualityIO
from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, MutableSeq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO

#Top level function as this makes it easier to use for debugging:
def write_read(filename, in_format, out_format) :
    records = list(SeqIO.parse(open(filename),in_format))
    #Write it out...
    handle = StringIO()
    SeqIO.write(records, handle, out_format)
    handle.seek(0)
    #Now load it back and check it agrees,
    records2 = list(SeqIO.parse(handle,out_format))
    compare_records(records, records2)

def compare_record(old, new) :
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
        raise ValuerError("Mismatch in phred_quality")
    if "solexa_quality" in old.letter_annotations \
    and "solexa_quality" in new.letter_annotations \
    and old.letter_annotations["solexa_quality"] != new.letter_annotations["solexa_quality"] :
            raise ValuerError("Mismatch in phred_quality")
    if "phred_quality" in old.letter_annotations \
    and "solexa_quality" in new.letter_annotations :
        #Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
        #Assume "old" is the original, and "new" has been converted.
        converted = [round(QualityIO.solexa_quality_from_phred(q)) \
                     for q in old.letter_annotations["phred_quality"]]
        if converted != new.letter_annotations["solexa_quality"] :
            raise ValueError("Mismatch in phred_quality vs solexa_quality")
    if "solexa_quality" in old.letter_annotations \
    and "phred_quality" in new.letter_annotations :
        #Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
        #Assume "old" is the original, and "new" has been converted.
        converted = [round(QualityIO.phred_quality_from_solexa(q)) \
                     for q in old.letter_annotations["solexa_quality"]]
        if converted != new.letter_annotations["phred_quality"] :
            raise ValueError("Mismatch in solexa_quality vs phred_quality")
    return True

def compare_records(old_list, new_list) :
    """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list) :
        raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list) :
        if not compare_record(old,new) :
            return False
    return True

class TestWriteRead(unittest.TestCase) :
    """Test can write and read back files."""

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
        #TODO - Record with no identifier?
        records = [record1, record2, record3, record4, record5, record6, record7]
        for format in ["fasta", "fastq", "fastq-solexa", "fastq-illumina", "qual"] :
            handle = StringIO()
            SeqIO.write(records, handle, format)
            handle.seek(0)
            compare_records(records, list(SeqIO.parse(handle, format)))
            
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
        #TODO - tests with truncation of the high PHRED values...

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


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
