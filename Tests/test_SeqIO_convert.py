# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.SeqIO.convert(...) function."""
import os
import unittest
import warnings
from Bio.Seq import UnknownSeq
from Bio import SeqIO
from Bio.SeqIO import QualityIO
from Bio.SeqIO._convert import _converter as converter_dict
from StringIO import StringIO
from Bio.Alphabet import generic_protein, generic_nucleotide, generic_dna

#TODO - share this with the QualityIO tests...
def truncation_expected(format):
    if format in ["fastq-solexa", "fastq-illumina"]:
        return 62
    elif format in ["fastq", "fastq-sanger"]:
        return 93
    else:
        return None

#Top level function as this makes it easier to use for debugging:
def check_convert(in_filename, in_format, out_format, alphabet=None):
    records = list(SeqIO.parse(in_filename,in_format, alphabet))
    #Write it out...
    handle = StringIO()
    qual_truncate = truncation_expected(out_format)
    if qual_truncate:
        warnings.simplefilter('ignore', UserWarning)
    SeqIO.write(records, handle, out_format)
    if qual_truncate:
        warnings.filters.pop()
    handle.seek(0)
    #Now load it back and check it agrees,
    records2 = list(SeqIO.parse(handle, out_format, alphabet))
    compare_records(records, records2, qual_truncate)
    #Finally, use the convert fuction, and check that agrees:
    handle2 = StringIO()
    if qual_truncate:
        warnings.simplefilter('ignore', UserWarning)
    SeqIO.convert(in_filename, in_format, handle2, out_format, alphabet)
    if qual_truncate:
        warnings.filters.pop()
    #We could re-parse this, but it is simpler and stricter:
    assert handle.getvalue() == handle2.getvalue()

def check_convert_fails(in_filename, in_format, out_format, alphabet=None):
    qual_truncate = truncation_expected(out_format)
    #We want the SAME error message from parse/write as convert!
    err1 = None
    try:
        records = list(SeqIO.parse(in_filename,in_format, alphabet))
        handle = StringIO()
        if qual_truncate:
            warnings.simplefilter('ignore', UserWarning)
        SeqIO.write(records, handle, out_format)
        if qual_truncate:
            warnings.filters.pop()
        handle.seek(0)
        assert False, "Parse or write should have failed!"
    except ValueError, err:
        err1 = err
    #Now do the conversion...
    try:
        handle2 = StringIO()
        if qual_truncate:
            warnings.simplefilter('ignore', UserWarning)
        SeqIO.convert(in_filename, in_format, handle2, out_format, alphabet)
        if qual_truncate:
            warnings.filters.pop()
        assert False, "Convert should have failed!"
    except ValueError, err2:
        assert str(err1) == str(err2), \
               "Different failures, parse/write:\n%s\nconvert:\n%s" \
               % (err1, err2)
    #print err
    
#TODO - move this to a shared test module...
def compare_record(old, new, truncate=None):
    """Quality aware SeqRecord comparision.

    This will check the mapping between Solexa and PHRED scores.
    It knows to ignore UnknownSeq objects for string matching (i.e. QUAL files).
    """
    if old.id != new.id:
        raise ValueError("'%s' vs '%s' " % (old.id, new.id))
    if old.description != new.description \
    and (old.id+" "+old.description).strip() != new.description \
    and new.description != "<unknown description>" \
    and new.description != "" : #e.g. tab format
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

class ConvertTests(unittest.TestCase):
    """Cunning unit test where methods are added at run time."""
    def simple_check(self, filename, in_format, out_format, alphabet):
        check_convert(filename, in_format, out_format, alphabet)
    def failure_check(self, filename, in_format, out_format, alphabet):
        check_convert_fails(filename, in_format, out_format, alphabet)

tests = [
    ("Quality/example.fastq", "fastq", None),
    ("Quality/example.fastq", "fastq-sanger", generic_dna),
    ("Quality/tricky.fastq", "fastq", generic_nucleotide),
    ("Quality/sanger_93.fastq", "fastq-sanger", None),
    ("Quality/sanger_faked.fastq", "fastq-sanger", generic_dna),
    ("Quality/solexa_faked.fastq", "fastq-solexa", generic_dna),
    ("Quality/illumina_faked.fastq", "fastq-illumina", generic_dna),
    ("EMBL/U87107.embl", "embl", None),
    ("EMBL/TRBG361.embl", "embl", None),
    ("GenBank/NC_005816.gb", "gb", None),
    ("GenBank/cor6_6.gb", "genbank", None),
    ]
for filename, format, alphabet in tests:
    for (in_format, out_format) in converter_dict:
        if in_format != format : continue
        def funct(fn,fmt1, fmt2, alpha):
            f = lambda x : x.simple_check(fn, fmt1, fmt2, alpha)
            f.__doc__ = "Convert %s from %s to %s" % (fn, fmt1, fmt2)
            return f
        setattr(ConvertTests, "test_%s_%s_to_%s" \
                % (filename.replace("/","_").replace(".","_"), in_format, out_format),
                funct(filename, in_format, out_format, alphabet))
        del funct

#Fail tests:
tests = [
    ("Quality/error_diff_ids.fastq", "fastq", None),
    ("Quality/error_long_qual.fastq", "fastq", None),
    ("Quality/error_no_qual.fastq", "fastq", None),
    ("Quality/error_qual_del.fastq", "fastq", None),
    ("Quality/error_qual_escape.fastq", "fastq", None),
    ("Quality/error_qual_null.fastq", "fastq", None),
    ("Quality/error_qual_space.fastq", "fastq", None),
    ("Quality/error_qual_tab.fastq", "fastq", None),
    ("Quality/error_qual_unit_sep.fastq", "fastq", None),
    ("Quality/error_qual_vtab.fastq", "fastq", None),
    ("Quality/error_short_qual.fastq", "fastq", None),
    ("Quality/error_spaces.fastq", "fastq", None),
    ("Quality/error_tabs.fastq", "fastq", None),
    ("Quality/error_trunc_at_plus.fastq", "fastq", None),
    ("Quality/error_trunc_at_qual.fastq", "fastq", None),
    ("Quality/error_trunc_at_seq.fastq", "fastq", None),
    ("Quality/error_trunc_in_title.fastq", "fastq", generic_dna),
    ("Quality/error_trunc_in_seq.fastq", "fastq", generic_nucleotide),
    ("Quality/error_trunc_in_plus.fastq", "fastq", None),
    ("Quality/error_trunc_in_qual.fastq", "fastq", generic_dna),
    ("Quality/error_double_seq.fastq", "fastq", generic_dna),
    ("Quality/error_double_qual.fastq", "fastq", generic_dna),
    ]
for filename, format, alphabet in tests:
    for (in_format, out_format) in converter_dict:
        if in_format != format : continue
        if in_format in ["fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina"] \
        and out_format in ["fasta", "tab"] and filename.startswith("Quality/error_qual_"):
            #TODO? These conversions don't check for bad characters in the quality,
            #and in order to pass this strict test they should.
            continue
        def funct(fn,fmt1, fmt2, alpha):
            f = lambda x : x.failure_check(fn, fmt1, fmt2, alpha)
            f.__doc__ = "Convert %s from %s to %s" % (fn, fmt1, fmt2)
            return f
        setattr(ConvertTests, "test_%s_%s_to_%s" \
                % (filename.replace("/","_").replace(".","_"), in_format, out_format),
                funct(filename, in_format, out_format, alphabet))
    del funct


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
