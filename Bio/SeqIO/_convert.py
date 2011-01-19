# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Optimised sequence conversion code (PRIVATE).

You are not expected to access this module, or any of its code, directly. This
is all handled internally by the Bio.SeqIO.convert(...) function which is the
public interface for this.

The idea here is rather while doing this will work:

from Bio import SeqIO
records = SeqIO.parse(in_handle, in_format)
count = SeqIO.write(records, out_handle, out_format)

it is shorter to write:

from Bio import SeqIO
count = SeqIO.convert(in_handle, in_format, out_handle, out_format)

Also, the convert function can take a number of special case optimisations. This
means that using Bio.SeqIO.convert() may be faster, as well as more convenient.
All these file format specific optimisations are handled by this (private) module.
"""

from Bio import SeqIO
#NOTE - Lots of lazy imports further on...

def _genbank_convert_fasta(in_handle, out_handle, alphabet=None):
    """Fast GenBank to FASTA (PRIVATE)."""
    #We don't need to parse the features...
    from Bio.GenBank.Scanner import GenBankScanner
    records = GenBankScanner().parse_records(in_handle, do_features=False)
    #For FASTA output we can ignore the alphabet too
    return SeqIO.write(records, out_handle, "fasta")


def _embl_convert_fasta(in_handle, out_handle, alphabet=None):
    """Fast EMBL to FASTA (PRIVATE)."""
    #We don't need to parse the features...
    from Bio.GenBank.Scanner import EmblScanner
    records = EmblScanner().parse_records(in_handle, do_features=False)
    #For FASTA output we can ignore the alphabet too
    return SeqIO.write(records, out_handle, "fasta")


def _fastq_generic(in_handle, out_handle, mapping):
    """FASTQ helper function where can't have data loss by truncation (PRIVATE)."""
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    #For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    null = chr(0)
    for title, seq, old_qual in FastqGeneralIterator(in_handle):
        count += 1
        #map the qual...
        qual = old_qual.translate(mapping)
        if null in qual:
            raise ValueError("Invalid character in quality string")
        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    return count

    
def _fastq_generic2(in_handle, out_handle, mapping, truncate_char, truncate_msg):
    """FASTQ helper function where there could be data loss by truncation (PRIVATE)."""
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    #For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    null = chr(0)
    for title, seq, old_qual in FastqGeneralIterator(in_handle):
        count += 1
        #map the qual...
        qual = old_qual.translate(mapping)
        if null in qual:
            raise ValueError("Invalid character in quality string")
        if truncate_char in qual:
            qual = qual.replace(truncate_char, chr(126))
            import warnings
            warnings.warn(truncate_msg)
        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    return count


def _fastq_sanger_convert_fastq_sanger(in_handle, out_handle, alphabet=None):
    """Fast Sanger FASTQ to Sanger FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    #Map unexpected chars to null
    mapping = "".join([chr(0) for ascii in range(0, 33)] \
                     +[chr(ascii) for ascii in range(33, 127)] \
                     +[chr(0) for ascii in range(127, 256)])
    assert len(mapping)==256
    return _fastq_generic(in_handle, out_handle, mapping)


def _fastq_solexa_convert_fastq_solexa(in_handle, out_handle, alphabet=None):
    """Fast Solexa FASTQ to Solexa FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.
    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    #Map unexpected chars to null
    mapping = "".join([chr(0) for ascii in range(0, 59)] \
                     +[chr(ascii) for ascii in range(59, 127)] \
                     +[chr(0) for ascii in range(127, 256)])
    assert len(mapping)==256
    return _fastq_generic(in_handle, out_handle, mapping)


def _fastq_illumina_convert_fastq_illumina(in_handle, out_handle, alphabet=None):
    """Fast Illumina 1.3+ FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.
    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    #Map unexpected chars to null
    mapping = "".join([chr(0) for ascii in range(0, 64)] \
                     +[chr(ascii) for ascii in range(64,127)] \
                     +[chr(0) for ascii in range(127,256)])
    assert len(mapping)==256
    return _fastq_generic(in_handle, out_handle, mapping)


def _fastq_illumina_convert_fastq_sanger(in_handle, out_handle, alphabet=None):
    """Fast Illumina 1.3+ FASTQ to Sanger FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    #Map unexpected chars to null
    mapping = "".join([chr(0) for ascii in range(0, 64)] \
                     +[chr(33+q) for q in range(0, 62+1)] \
                     +[chr(0) for ascii in range(127, 256)])
    assert len(mapping)==256
    return _fastq_generic(in_handle, out_handle, mapping)


def _fastq_sanger_convert_fastq_illumina(in_handle, out_handle, alphabet=None):
    """Fast Sanger FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion. Will issue a warning if the scores had to be truncated at 62
    (maximum possible in the Illumina 1.3+ FASTQ format)
    """
    #Map unexpected chars to null
    trunc_char = chr(1)
    mapping = "".join([chr(0) for ascii in range(0, 33)] \
                     +[chr(64+q) for q in range(0, 62+1) ] \
                     +[trunc_char for ascii in range(96,127)] \
                     +[chr(0) for ascii in range(127, 256)])
    assert len(mapping)==256
    return _fastq_generic2(in_handle, out_handle, mapping, trunc_char,
                          "Data loss - max PHRED quality 62 in Illumina 1.3+ FASTQ")


def _fastq_solexa_convert_fastq_sanger(in_handle, out_handle, alphabet=None):
    """Fast Solexa FASTQ to Sanger FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    #Map unexpected chars to null
    from Bio.SeqIO.QualityIO import phred_quality_from_solexa
    mapping = "".join([chr(0) for ascii in range(0, 59)] \
                     +[chr(33+int(round(phred_quality_from_solexa(q)))) \
                       for q in range(-5, 62+1)]\
                      +[chr(0) for ascii in range(127, 256)])
    assert len(mapping)==256
    return _fastq_generic(in_handle, out_handle, mapping)

def _fastq_sanger_convert_fastq_solexa(in_handle, out_handle, alphabet=None):
    """Fast Sanger FASTQ to Solexa FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion. Will issue a warning if the scores had to be truncated at 62
    (maximum possible in the Solexa FASTQ format)
    """
    #Map unexpected chars to null
    from Bio.SeqIO.QualityIO import solexa_quality_from_phred
    trunc_char = chr(1)
    mapping = "".join([chr(0) for ascii in range(0, 33)] \
                     +[chr(64+int(round(solexa_quality_from_phred(q)))) \
                       for q in range(0, 62+1)] \
                     +[trunc_char for ascii in range(96, 127)] \
                     +[chr(0) for ascii in range(127, 256)])
    assert len(mapping)==256
    return _fastq_generic2(in_handle, out_handle, mapping, trunc_char,
                          "Data loss - max Solexa quality 62 in Solexa FASTQ")


def _fastq_solexa_convert_fastq_illumina(in_handle, out_handle, alphabet=None):
    """Fast Solexa FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    #Map unexpected chars to null
    from Bio.SeqIO.QualityIO import phred_quality_from_solexa
    mapping = "".join([chr(0) for ascii in range(0, 59)] \
                     +[chr(64+int(round(phred_quality_from_solexa(q)))) \
                       for q in range(-5, 62+1)]\
                      +[chr(0) for ascii in range(127, 256)])
    assert len(mapping)==256
    return _fastq_generic(in_handle, out_handle, mapping)


def _fastq_illumina_convert_fastq_solexa(in_handle, out_handle, alphabet=None):
    """Fast Illumina 1.3+ FASTQ to Solexa FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    #Map unexpected chars to null
    from Bio.SeqIO.QualityIO import solexa_quality_from_phred
    trunc_char = chr(1)
    mapping = "".join([chr(0) for ascii in range(0, 64)] \
                     +[chr(64+int(round(solexa_quality_from_phred(q)))) \
                       for q in range(0, 62+1)] \
                     +[chr(0) for ascii in range(127, 256)])
    assert len(mapping)==256
    return _fastq_generic(in_handle, out_handle, mapping)


def _fastq_convert_fasta(in_handle, out_handle, alphabet=None):
    """Fast FASTQ to FASTA conversion (PRIVATE).

    Avoids dealing with the FASTQ quality encoding, and creating SeqRecord and
    Seq objects in order to speed up this conversion.

    NOTE - This does NOT check the characters used in the FASTQ quality string
    are valid!
    """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    #For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    for title, seq, qual in FastqGeneralIterator(in_handle):
        count += 1
        out_handle.write(">%s\n" % title)
        #Do line wrapping
        for i in range(0, len(seq), 60):
            out_handle.write(seq[i:i+60] + "\n")
    return count

def _fastq_convert_tab(in_handle, out_handle, alphabet=None):
    """Fast FASTQ to simple tabbed conversion (PRIVATE).

    Avoids dealing with the FASTQ quality encoding, and creating SeqRecord and
    Seq objects in order to speed up this conversion.

    NOTE - This does NOT check the characters used in the FASTQ quality string
    are valid!
    """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    #For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    for title, seq, qual in FastqGeneralIterator(in_handle):
        count += 1
        out_handle.write("%s\t%s\n" % (title.split(None, 1)[0], seq))
    return count

def _fastq_convert_qual(in_handle, out_handle, mapping):
    """FASTQ helper function for QUAL output (PRIVATE).

    Mapping should be a dictionary mapping expected ASCII characters from the
    FASTQ quality string to PHRED quality scores (as strings).
    """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    #For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    for title, seq, qual in FastqGeneralIterator(in_handle):
        count += 1
        out_handle.write(">%s\n" % title)
        #map the qual...
        try:
            qualities_strs = [mapping[ascii] for ascii in qual]
        except KeyError:
            raise ValueError("Invalid character in quality string")
        data = " ".join(qualities_strs)
        while True:
            if len(data) <= 60:
                out_handle.write(data + "\n")
                break
            else:
                #By construction there must be spaces in the first 60 chars
                #(unless we have 60 digit or higher quality scores!)
                i = data.rfind(" ", 0, 60)
                out_handle.write(data[:i] + "\n")
                data = data[i+1:]
    return count

    
def _fastq_sanger_convert_qual(in_handle, out_handle, alphabet=None):
    """Fast Sanger FASTQ to QUAL conversion (PRIVATE)."""
    mapping = dict((chr(q+33), str(q)) for q in range(0,93+1))
    return _fastq_convert_qual(in_handle, out_handle, mapping)


def _fastq_solexa_convert_qual(in_handle, out_handle, alphabet=None):
    """Fast Solexa FASTQ to QUAL conversion (PRIVATE)."""
    from Bio.SeqIO.QualityIO import phred_quality_from_solexa
    mapping = dict((chr(q+64), str(int(round(phred_quality_from_solexa(q))))) \
                   for q in range(-5,62+1))
    return _fastq_convert_qual(in_handle, out_handle, mapping)


def _fastq_illumina_convert_qual(in_handle, out_handle, alphabet=None):
    """Fast Illumina 1.3+ FASTQ to QUAL conversion (PRIVATE)."""
    mapping = dict((chr(q+64), str(q)) for q in range(0,62+1))
    return _fastq_convert_qual(in_handle, out_handle, mapping)


#TODO? - Handling aliases explicitly would let us shorten this list:
_converter = {
    ("genbank", "fasta") : _genbank_convert_fasta,
    ("gb", "fasta") : _genbank_convert_fasta,
    ("embl", "fasta") : _embl_convert_fasta,
    ("fastq", "fasta") : _fastq_convert_fasta,
    ("fastq-sanger", "fasta") : _fastq_convert_fasta,
    ("fastq-solexa", "fasta") : _fastq_convert_fasta,
    ("fastq-illumina", "fasta") : _fastq_convert_fasta,
    ("fastq", "tab") : _fastq_convert_tab,
    ("fastq-sanger", "tab") : _fastq_convert_tab,
    ("fastq-solexa", "tab") : _fastq_convert_tab,
    ("fastq-illumina", "tab") : _fastq_convert_tab,
    ("fastq", "fastq") : _fastq_sanger_convert_fastq_sanger,
    ("fastq-sanger", "fastq") : _fastq_sanger_convert_fastq_sanger,
    ("fastq-solexa", "fastq") : _fastq_solexa_convert_fastq_sanger,
    ("fastq-illumina", "fastq") : _fastq_illumina_convert_fastq_sanger,
    ("fastq", "fastq-sanger") : _fastq_sanger_convert_fastq_sanger,
    ("fastq-sanger", "fastq-sanger") : _fastq_sanger_convert_fastq_sanger,
    ("fastq-solexa", "fastq-sanger") : _fastq_solexa_convert_fastq_sanger,
    ("fastq-illumina", "fastq-sanger") : _fastq_illumina_convert_fastq_sanger,
    ("fastq", "fastq-solexa") : _fastq_sanger_convert_fastq_solexa,
    ("fastq-sanger", "fastq-solexa") : _fastq_sanger_convert_fastq_solexa,
    ("fastq-solexa", "fastq-solexa") : _fastq_solexa_convert_fastq_solexa,
    ("fastq-illumina", "fastq-solexa") : _fastq_illumina_convert_fastq_solexa,
    ("fastq", "fastq-illumina") : _fastq_sanger_convert_fastq_illumina,
    ("fastq-sanger", "fastq-illumina") : _fastq_sanger_convert_fastq_illumina,
    ("fastq-solexa", "fastq-illumina") : _fastq_solexa_convert_fastq_illumina,
    ("fastq-illumina", "fastq-illumina") : _fastq_illumina_convert_fastq_illumina,
    ("fastq", "qual") : _fastq_sanger_convert_qual,
    ("fastq-sanger", "qual") : _fastq_sanger_convert_qual,
    ("fastq-solexa", "qual") : _fastq_solexa_convert_qual,
    ("fastq-illumina", "qual") : _fastq_illumina_convert_qual,
    }

def _handle_convert(in_handle, in_format, out_handle, out_format, alphabet=None):
    """SeqIO conversion function (PRIVATE)."""
    try:
        f = _converter[(in_format, out_format)]
    except KeyError:
        f = None
    if f:
        return f(in_handle, out_handle, alphabet)
    else:
        records = SeqIO.parse(in_handle, in_format, alphabet)
        return SeqIO.write(records, out_handle, out_format)
