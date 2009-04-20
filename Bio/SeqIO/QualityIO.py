# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading and writing FASTQ and QUAL format files as
# SeqRecord objects, and is expected to be used via the Bio.SeqIO API.

"""Bio.SeqIO support for the "fastq" and "qual" file formats.

Note that you are expected to use this code via the Bio.SeqIO interface, as
shown below.

The FASTQ file format is used frequently at the Wellcome Trust Sanger Institute
to bundle a FASTA sequence and its PHRED quality data (integers between 0 and
90).  Rather than using a single FASTQ file, often paired FASTA and QUAL files
are used containing the sequence and the quality information separately.

The PHRED software reads DNA sequencing trace files, calls bases, and
assigns a quality value between 0 and 90 to each called base using a logged
transformation of the error probability, Q = -10 log10( Pe ), for example::

    Pe = 0.0,         Q =  0
    Pe = 0.1,         Q = 10
    Pe = 0.01,        Q = 20
    ...
    Pe = 0.00000001,  Q = 80
    Pe = 0.000000001, Q = 90

In the QUAL format these quality values are held as space separated text in
a FASTA like file format.  In the FASTQ format, each quality values is encoded
with a single ASCI character using chr(Q+33), meaning zero maps to the
character "!" and for example 80 maps to "q".  The sequences and quality are
then stored in pairs in a FASTA like format.

Unfortunately there is no official document describing the FASTQ file format,
and worse, several related but different variants exist.  Reasonable
documentation exists at: http://maq.sourceforge.net/fastq.shtml

Solexa/Illumina quality scores use Q = - 10 log10 ( Pe / (1-Pe) ), which can
be negative or easily exceed 90.  PHRED scores and Solexa scores are NOT
interchangeable (but a reasonable mapping can be achieved between them).
Confusingly Solexa produces a FASTQ like file but using their own score
mapping instead.

Also note that Roche 454 sequencers can output files in the QUAL format, and
thankfully they use PHREP style scores like Sanger.  To extract QUAL files from
a Roche 454 SFF binary file, use the Roche off instrument command line tool
"sffinfo" with the -q or -qual argument.  You can extract a matching FASTA file
using the -s or -seq argument instead.

You are expected to use this module via the Bio.SeqIO functions, with the
following format names:
 - "fastq" means Sanger style FASTQ files using PHRED scores.
 - "fastq-solexa" means Solexa/Illumina style FASTQ files.
 - "qual" means simple quality files using PHRED scores.

For example, consider the following short FASTQ file (extracted from a real
NCBI dataset)::

    @EAS54_6_R1_2_1_413_324
    CCCTTCTTGTCTTCAGCGTTTCTCC
    +
    ;;3;;;;;;;;;;;;7;;;;;;;88
    @EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    +
    ;;;;;;;;;;;7;;;;;-;;;3;83
    @EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG
    +
    ;;;;;;;;;;;9;7;;.7;393333

This contains three reads of length 25.  From the read length these were
probably originally from an early Solexa/Illumina sequencer but NCBI have
followed the Sanger FASTQ convention and this actually uses PHRED style
qualities.  This means we can parse this file using Bio.SeqIO using "fastq"
as the format name:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse(open("Quality/example.fastq"), "fastq") :
    ...     print record.id, record.seq
    EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
    EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
    EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG

The qualities are held as a list of integers in each record's annotation:

    >>> print record
    ID: EAS54_6_R1_2_1_443_348
    Name: EAS54_6_R1_2_1_443_348
    Description: EAS54_6_R1_2_1_443_348
    Number of features: 0
    Per letter annotation for: phred_quality
    Seq('GTTGCTTCTGGCGTGGGTGGGGGGG', SingleLetterAlphabet())
    >>> print record.letter_annotations["phred_quality"]
    [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

You can use the SeqRecord format method you can show this in the QUAL format:

    >>> print record.format("qual")
    >EAS54_6_R1_2_1_443_348
    26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 13 22 26 18
    24 18 18 18 18
    <BLANKLINE>

Or go back to the FASTQ format,

    >>> print record.format("fastq")
    @EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG
    +
    ;;;;;;;;;;;9;7;;.7;393333
    <BLANKLINE>

You can also get Biopython to convert the scores and show a Solexa style
FASTQ file:

    >>> print record.format("fastq-solexa")
    @EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG
    +
    ZZZZZZZZZZZXZVZZMVZRXRRRR
    <BLANKLINE>

If you wanted to trim your sequences (perhaps to remove low quality regions,
or to remove a primer sequence), try slicing the SeqRecord objects.  e.g.

    >>> sub_rec = record[5:15]
    >>> print sub_rec
    ID: EAS54_6_R1_2_1_443_348
    Name: EAS54_6_R1_2_1_443_348
    Description: EAS54_6_R1_2_1_443_348
    Number of features: 0
    Per letter annotation for: phred_quality
    Seq('TTCTGGCGTG', SingleLetterAlphabet())
    >>> print sub_rec.letter_annotations["phred_quality"]
    [26, 26, 26, 26, 26, 26, 24, 26, 22, 26]
    >>> print sub_rec.format("fastq")
    @EAS54_6_R1_2_1_443_348
    TTCTGGCGTG
    +
    ;;;;;;9;7;
    <BLANKLINE>
    
If you wanted to, you could read in this FASTQ file, and save it as a QUAL file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse(open("Quality/example.fastq"), "fastq")
    >>> out_handle = open("Quality/temp.qual", "w")
    >>> SeqIO.write(record_iterator, out_handle, "qual")
    3
    >>> out_handle.close()

You can of course read in a QUAL file, such as the one we just created:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse(open("Quality/temp.qual"), "qual") :
    ...     print record.id, record.seq
    EAS54_6_R1_2_1_413_324 ?????????????????????????
    EAS54_6_R1_2_1_540_792 ?????????????????????????
    EAS54_6_R1_2_1_443_348 ?????????????????????????

Notice that QUAL files don't have a proper sequence present!  But the quality
information is there:

    >>> print record
    ID: EAS54_6_R1_2_1_443_348
    Name: EAS54_6_R1_2_1_443_348
    Description: EAS54_6_R1_2_1_443_348
    Number of features: 0
    Per letter annotation for: phred_quality
    UnknownSeq(25, alphabet = SingleLetterAlphabet(), character = '?')
    >>> print record.letter_annotations["phred_quality"]
    [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

Just to keep things tidy, if you are following this example yourself, you can
delete this temporary file now:

    >>> import os
    >>> os.remove("Quality/temp.qual")

Sometimes you won't have a FASTQ file, but rather just a pair of FASTA and QUAL
files.  Because the Bio.SeqIO system is designed for reading single files, you
would have to read the two in separately and then combine the data.  However,
since this is such a common thing to want to do, there is a helper iterator
defined in this module that does this for you - PairedFastaQualIterator.

Alternatively, if you have enough RAM to hold all the records in memory at once,
then a simple dictionary approach would work:

    >>> from Bio import SeqIO
    >>> reads = SeqIO.to_dict(SeqIO.parse(open("Quality/example.fasta"), "fasta"))
    >>> for rec in SeqIO.parse(open("Quality/example.qual"), "qual") :
    ...     reads[rec.id].letter_annotations["phred_quality"]=rec.letter_annotations["phred_quality"]

You can then access any record by its key, and get both the sequence and the
quality scores.

    >>> print reads["EAS54_6_R1_2_1_540_792"].format("fastq")
    @EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    +
    ;;;;;;;;;;;7;;;;;-;;;3;83
    <BLANKLINE>

It is important that you explicitly tell Bio.SeqIO which FASTQ variant you are
using ("fastq" for the Sanger standard using PHRED values, or "fastq-solexa"
for the Solexa/Illumina variant), as this cannot be detected reliably
automatically.
"""
__docformat__ = "epytext en" #Don't just use plain text in epydoc API pages!

#See also http://blog.malde.org/index.php/2008/09/09/the-fastq-file-format-for-sequences/

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Interfaces import SequentialSequenceWriter
from math import log

# define score offsets. See discussion for differences between Sanger and
# Solexa offsets.
SANGER_SCORE_OFFSET = 33
SOLEXA_SCORE_OFFSET = 64

def solexa_quality_from_phred(phred_quality) :
    """Covert a PHRED quality (range 0 to about 90) to a Solexa quality.

    This will return a floating point number, it is up to you to round this to
    the nearest integer if appropriate.  e.g.

    >>> print "%0.2f" % round(solexa_quality_from_phred(80),2)
    80.00
    >>> print "%0.2f" % round(solexa_quality_from_phred(50),2)
    50.00
    >>> print "%0.2f" % round(solexa_quality_from_phred(20),2)
    19.96
    >>> print "%0.2f" % round(solexa_quality_from_phred(10),2)
    9.54
    >>> print "%0.2f" % round(solexa_quality_from_phred(1),2)
    -5.87
    """
    return 10*log(10**(phred_quality/10.0) - 1, 10)

def phred_quality_from_solexa(solexa_quality) :
    """Convert a Solexa quality (which can be negative) to a PHRED quality.

    This will return a floating point number, it is up to you to round this to
    the nearest integer if appropriate.  e.g.

    >>> print "%0.2f" % round(phred_quality_from_solexa(80),2)
    80.00
    >>> print "%0.2f" % round(phred_quality_from_solexa(20),2)
    20.04
    >>> print "%0.2f" % round(phred_quality_from_solexa(10),2)
    10.41
    >>> print "%0.2f" % round(phred_quality_from_solexa(0),2)
    3.01
    >>> print "%0.2f" % round(phred_quality_from_solexa(-10),2)
    0.41
    """
    return 10*log(10**(solexa_quality/10.0) + 1, 10)

def _get_phred_quality(record) :
    """Extract PHRED qualities from a SeqRecord's letter_annotations (PRIVATE).

    If there are no PHRED qualities, but there are Solexa qualities, those are
    used instead after conversion.
    """
    try :
        return record.letter_annotations["phred_quality"]
    except KeyError :
        pass
    try :
        return [phred_quality_from_solexa(q) for \
                q in record.letter_annotations["solexa_quality"]]
    except KeyError :
        raise ValueError("No suitable quality scores found in letter_annotations "
                         "of SeqRecord (id=%s)." % record.id)

def _get_solexa_quality(record) :
    """Extract Solexa qualities from a SeqRecord's letter_annotations (PRIVATE).

    If there are no Solexa qualities, but there are PHRED qualities, those are
    used instead after conversion.
    """
    try :
        return record.letter_annotations["solexa_quality"]
    except KeyError :
        pass
    try :
        return [solexa_quality_from_phred(q) for \
                q in record.letter_annotations["phred_quality"]]
    except KeyError :
        raise ValueError("No suitable quality scores found in letter_annotation "
                         "of SeqRecord (id=%s)." % record.id)


#TODO - Default to nucleotide or even DNA?
def FastqGeneralIterator(handle) :
    """Iterate over Fastq records as string tuples (not as SeqRecord objects).

    This code does not try to interpret the quality string numerically.  It
    just returns tuples of the title, sequence and quality as strings.  For
    the sequence and quality, any whitespace (such as new lines) is removed.

    Our SeqRecord based FASTQ iterators call this function internally, and then
    turn the strings into a SeqRecord objects, mapping the quality string into
    a list of numerical scores.  If you want to do a custom quality mapping,
    then you might consider calling this function directly.

    For parsing FASTQ files, the title string from the "@" line at the start
    of each record can optionally be omitted on the "+" lines.  If it is
    repeated, it must be identical.

    The sequence string and the quality string can optionally be split over
    multiple lines, although several sources discourage this.  In comparison,
    for the FASTA file format line breaks between 60 and 80 characters are
    the norm.

    WARNING - Because the "@" character can appear in the quality string,
    this can cause problems as this is also the marker for the start of
    a new sequence.  In fact, the "+" sign can also appear as well.  Some
    sources recommended having no line breaks in the  quality to avoid this,
    but even that is not enough, consider this example::

        @071113_EAS56_0053:1:1:998:236
        TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA
        +071113_EAS56_0053:1:1:998:236
        IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
        @071113_EAS56_0053:1:1:182:712
        ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG
        +
        @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
        @071113_EAS56_0053:1:1:153:10
        TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT
        +
        IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
        @071113_EAS56_0053:1:3:990:501
        TGGGAGGTTTTATGTGGA
        AAGCAGCAATGTACAAGA
        +
        IIIIIII.IIIIII1@44
        @-7.%<&+/$/%4(++(%

    This is four PHRED encoded FASTQ entries originally from an NCBI source
    (given the read length of 36, these are probably Solexa Illumna reads where
    the quality has been mapped onto the PHRED values).

    This example has been edited to illustrate some of the nasty things allowed
    in the FASTQ format.  Firstly, on the "+" lines most but not all of the
    (redundant) identifiers are ommited.  In real files it is likely that all or
    none of these extra identifiers will be present.

    Secondly, while the first three sequences have been shown without line
    breaks, the last has been split over multiple lines.  In real files any line
    breaks are likely to be consistent.

    Thirdly, some of the quality string lines start with an "@" character.  For
    the second record this is unavoidable.  However for the fourth sequence this
    only happens because its quality string is split over two lines.  A naive
    parser could wrongly treat any line starting with an "@" as the beginning of
    a new sequence!  This code copes with this possible ambiguity by keeping track
    of the length of the sequence which gives the expected length of the quality
    string.

    Using this tricky example file as input, this short bit of code demonstrates
    what this parsing function would return:

    >>> handle = open("Quality/tricky.fastq", "rU")
    >>> for (title, sequence, quality) in FastqGeneralIterator(handle) :
    ...     print title
    ...     print sequence, quality
    071113_EAS56_0053:1:1:998:236
    TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
    071113_EAS56_0053:1:1:182:712
    ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
    071113_EAS56_0053:1:1:153:10
    TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
    071113_EAS56_0053:1:3:990:501
    TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA IIIIIII.IIIIII1@44@-7.%<&+/$/%4(++(%
    >>> handle.close()

    Finally we note that some sources state that the quality string should
    start with "!" (which using the PHRED mapping means the first letter always
    has a quality score of zero).  This rather restrictive rule is not widely
    observed, so is therefore ignored here.  One plus point about this "!" rule
    is that (provided there are no line breaks in the quality sequence) it
    would prevent the above problem with the "@" character.
    """
    #Skip any text before the first record (e.g. blank lines, comments?)
    while True :
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line[0] == "@" :
            break

    while True :
        if line[0]!="@" :
            raise ValueError("Records in Fastq files should start with '@' character")
        title_line = line[1:].rstrip()

        seq_lines = []
        line = handle.readline()
        while True:
            if not line :
                raise ValueError("End of file without quality information.")
            if line[0] == "+":
                #The title here is optional, but if present must match!
                if line[1:].rstrip() and line[1:].rstrip() != title_line :
                    raise ValueError("Sequence and quality captions differ.")
                break
            seq_lines.extend(line.split()) #removes any whitespace
            line = handle.readline()

        seq_string = "".join(seq_lines)
        del seq_lines

        quality_lines = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == "@":
                #This COULD be the start of a new sequence. However, it MAY just
                #be a line of quality data which starts with a "@" character.  We
                #should be able to check this by looking at the sequence length
                #and the amount of quality data found so far.
                if len("".join(quality_lines)) >= len(seq_string) :
                    #We expect it to be equal if this is the start of a new record.
                    #If the quality data is longer, we'll raise an error below.
                    break
                #Continue - its just some (more) sequence data.
                
            quality_lines.extend(line.split()) #removes any whitespace
            line = handle.readline()

        quality_string = "".join(quality_lines)
        del quality_lines
        
        if len(seq_string) != len(quality_string) :
            raise ValueError("Lengths of sequence and quality values differs "
                             " for %s (%i and %i)." \
                             % (title_line, len(seq_string), len(quality_string)))

        #Return the record and then continue...
        yield (title_line, seq_string, quality_string)
        if not line : return #StopIteration at end of file
    assert False, "Should not reach this line"
        
#This is a generator function!
def FastqPhredIterator(handle, alphabet = single_letter_alphabet, title2ids = None) :
    """Generator function to iterate over FASTQ records (as SeqRecord objects).

     - handle - input file
     - alphabet - optional alphabet
     - title2ids - A function that, when given the title line from the FASTQ
                   file (without the beginning >), will return the id, name and
                   description (in that order) for the record as a tuple of
                   strings.  If this is not given, then the entire title line
                   will be used as the description, and the first word as the
                   id and name.

    Note that use of title2ids matches that of Bio.SeqIO.FastaIO.

    For each sequence in a (Sanger style) FASTQ file there is a matching string
    encoding the PHRED qualities (integers between 0 and about 90) using ASCII
    values with an offset of 33.

    For example, consider a file containing three short reads::

        @EAS54_6_R1_2_1_413_324
        CCCTTCTTGTCTTCAGCGTTTCTCC
        +
        ;;3;;;;;;;;;;;;7;;;;;;;88
        @EAS54_6_R1_2_1_540_792
        TTGGCAGGCCAAGGCCGATGGATCA
        +
        ;;;;;;;;;;;7;;;;;-;;;3;83
        @EAS54_6_R1_2_1_443_348
        GTTGCTTCTGGCGTGGGTGGGGGGG
        +
        ;;;;;;;;;;;9;7;;.7;393333

    For each sequence (e.g. "CCCTTCTTGTCTTCAGCGTTTCTCC") there is a matching
    string encoding the PHRED qualities using a ASCI values with an offset of
    33 (e.g. ";;3;;;;;;;;;;;;7;;;;;;;88").

    Using this module directly you might run:

    >>> handle = open("Quality/example.fastq", "rU")
    >>> for record in FastqPhredIterator(handle) :
    ...     print record.id, record.seq
    EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
    EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
    EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG
    >>> handle.close()

    Typically however, you would call this via Bio.SeqIO instead with "fastq" as
    the format:

    >>> from Bio import SeqIO
    >>> handle = open("Quality/example.fastq", "rU")
    >>> for record in SeqIO.parse(handle, "fastq") :
    ...     print record.id, record.seq
    EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
    EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
    EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG
    >>> handle.close()

    If you want to look at the qualities, they are record in each record's
    per-letter-annotation dictionary as a simple list of integers:

    >>> print record.letter_annotations["phred_quality"]
    [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]
    """
    for title_line, seq_string, quality_string in FastqGeneralIterator(handle) :
        if title2ids :
            id, name, descr = title2ids(title_line)
        else :
            descr = title_line
            id   = descr.split()[0]
            name = id
        record = SeqRecord(Seq(seq_string, alphabet),
                           id=id, name=name, description=descr)
        
        assert SANGER_SCORE_OFFSET == ord("!")
        #According to BioPerl documentation at least, the first character should
        #be an "!" (and therefore quality zero).  This seems crazy - what if the
        #sequence has been trimmed to remove any poor quality sequence?  In any
        #case real examples from the NCBI don't follow this practice, so we
        #won't enforce it here.
        #e.g. ftp://ftp.ncbi.nih.gov/pub/TraceDB/ShortRead/SRA000271/fastq/200x36x36-071113_EAS56_0053-s_1_1.fastq.gz
        #
        #if quality_string[0] != "!" :
        #    raise ValueError("The quality string should always start with a ! character.")
        qualities = [ord(letter)-SANGER_SCORE_OFFSET for letter in quality_string]
        if qualities :
            if min(qualities) < 0 or max(qualities) > 90 :
                raise ValueError("Quality score outside 0 to 90 found - these are perhaps "
                                 "in a Solexa/Illumina format, not the Sanger FASTQ format "
                                 "which uses PHRED scores.")
        record.letter_annotations["phred_quality"] = qualities
        yield record

#This is a generator function!
def FastqSolexaIterator(handle, alphabet = single_letter_alphabet, title2ids = None) :
    """Parsing the Solexa/Illumina FASTQ like files (which differ in the quality mapping).

    The optional arguments are the same as those for the FastqPhredIterator.

    For each sequence in Solexa/Illumina FASTQ files there is a matching string
    encoding the Solexa integer qualities using ASCII values with an offset
    of 64.  Solexa scores are scaled differently to PHRED scores, and Biopython
    will NOT perform any automatic conversion when loading.

    For example, consider a file containing these five records::
        
        @SLXA-B3_649_FC8437_R1_1_1_610_79
        GATGTGCAATACCTTTGTAGAGGAA
        +SLXA-B3_649_FC8437_R1_1_1_610_79
        YYYYYYYYYYYYYYYYYYWYWYYSU
        @SLXA-B3_649_FC8437_R1_1_1_397_389
        GGTTTGAGAAAGAGAAATGAGATAA
        +SLXA-B3_649_FC8437_R1_1_1_397_389
        YYYYYYYYYWYYYYWWYYYWYWYWW
        @SLXA-B3_649_FC8437_R1_1_1_850_123
        GAGGGTGTTGATCATGATGATGGCG
        +SLXA-B3_649_FC8437_R1_1_1_850_123
        YYYYYYYYYYYYYWYYWYYSYYYSY
        @SLXA-B3_649_FC8437_R1_1_1_362_549
        GGAAACAAAGTTTTTCTCAACATAG
        +SLXA-B3_649_FC8437_R1_1_1_362_549
        YYYYYYYYYYYYYYYYYYWWWWYWY
        @SLXA-B3_649_FC8437_R1_1_1_183_714
        GTATTATTTAATGGCATACACTCAA
        +SLXA-B3_649_FC8437_R1_1_1_183_714
        YYYYYYYYYYWYYYYWYWWUWWWQQ
        
    Using this module directly you might run:

    >>> handle = open("Quality/solexa_example.fastq", "rU")
    >>> for record in FastqSolexaIterator(handle) :
    ...     print record.id, record.seq
    SLXA-B3_649_FC8437_R1_1_1_610_79 GATGTGCAATACCTTTGTAGAGGAA
    SLXA-B3_649_FC8437_R1_1_1_397_389 GGTTTGAGAAAGAGAAATGAGATAA
    SLXA-B3_649_FC8437_R1_1_1_850_123 GAGGGTGTTGATCATGATGATGGCG
    SLXA-B3_649_FC8437_R1_1_1_362_549 GGAAACAAAGTTTTTCTCAACATAG
    SLXA-B3_649_FC8437_R1_1_1_183_714 GTATTATTTAATGGCATACACTCAA
    >>> handle.close()

    Typically however, you would call this via Bio.SeqIO instead with "fastq" as
    the format:

    >>> from Bio import SeqIO
    >>> handle = open("Quality/solexa_example.fastq", "rU")
    >>> for record in SeqIO.parse(handle, "fastq-solexa") :
    ...     print record.id, record.seq
    SLXA-B3_649_FC8437_R1_1_1_610_79 GATGTGCAATACCTTTGTAGAGGAA
    SLXA-B3_649_FC8437_R1_1_1_397_389 GGTTTGAGAAAGAGAAATGAGATAA
    SLXA-B3_649_FC8437_R1_1_1_850_123 GAGGGTGTTGATCATGATGATGGCG
    SLXA-B3_649_FC8437_R1_1_1_362_549 GGAAACAAAGTTTTTCTCAACATAG
    SLXA-B3_649_FC8437_R1_1_1_183_714 GTATTATTTAATGGCATACACTCAA
    >>> handle.close()

    If you want to look at the qualities, they are recorded in each record's
    per-letter-annotation dictionary as a simple list of integers:

    >>> print record.letter_annotations["solexa_quality"]
    [25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 23, 25, 25, 25, 25, 23, 25, 23, 23, 21, 23, 23, 23, 17, 17]

    These scores aren't very good, but they are high enough that they map
    almost exactly onto PHRED scores:

    >>> print "%0.2f" % phred_quality_from_solexa(25)
    25.01

    Let's look at another example read which is even worse, where there are
    more noticeable differences between the Solexa and PHRED scores::

         @slxa_0013_1_0001_24
         ACAAAAATCACAAGCATTCTTATACACC
         +slxa_0013_1_0001_24
         ??????????????????:??<?<-6%.

    Again, you would typically use Bio.SeqIO to read this file in (rather than
    calling the Bio.SeqIO.QualtityIO module directly).  Most FASTQ files will
    contain thousands of reads, so you would normally use Bio.SeqIO.parse()
    as shown above.  This example has only as one entry, so instead we can
    use the Bio.SeqIO.read() function:

    >>> from Bio import SeqIO
    >>> handle = open("Quality/solexa.fastq", "rU")
    >>> record = SeqIO.read(handle, "fastq-solexa")
    >>> handle.close()
    >>> print record.id, record.seq
    slxa_0013_1_0001_24 ACAAAAATCACAAGCATTCTTATACACC
    >>> print record.letter_annotations["solexa_quality"]
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6, -1, -1, -4, -1, -4, -19, -10, -27, -18]

    These quality scores are so low that when converted from the Solexa scheme
    into PHRED scores they look quite different:

    >>> print "%0.2f" % phred_quality_from_solexa(-1)
    2.54

    Note you can use the Bio.SeqIO.write() function or the SeqRecord's format
    method to output the record(s):

    >>> print record.format("fastq-solexa")
    @slxa_0013_1_0001_24
    ACAAAAATCACAAGCATTCTTATACACC
    +
    ??????????????????:??<?<-6%.
    <BLANKLINE>

    Note this output is slightly different from the input file as Biopython
    has left out the optional repetition of the sequence identifier on the "+"
    line.  If you want the to use PHRED scores, use "fastq" or "qual" as the
    output format instead, and Biopython will do the conversion for you:

    >>> print record.format("fastq")
    @slxa_0013_1_0001_24
    ACAAAAATCACAAGCATTCTTATACACC
    +
    $$$$$$$$$$$$$$$$$$"$$"$"!!!!
    <BLANKLINE>

    >>> print record.format("qual")
    >slxa_0013_1_0001_24
    3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 3 3 1 3 1 0 0 0 0
    <BLANKLINE>
    """
    for title_line, seq_string, quality_string in FastqGeneralIterator(handle) :
        if title2ids :
            id, name, descr = title_line
        else :
            descr = title_line
            id   = descr.split()[0]
            name = id
        record = SeqRecord(Seq(seq_string, alphabet),
                           id=id, name=name, description=descr)
        qualities = [ord(letter)-SOLEXA_SCORE_OFFSET for letter in quality_string]
        #DO NOT convert these into PHRED qualities automatically!
        record.letter_annotations["solexa_quality"] = qualities
        yield record

def QualPhredIterator(handle, alphabet = single_letter_alphabet, title2ids = None) :
    """For QUAL files which include PHRED quality scores, but no sequence.

    For example, consider this short QUAL file::

        >EAS54_6_R1_2_1_413_324
        26 26 18 26 26 26 26 26 26 26 26 26 26 26 26 22 26 26 26 26
        26 26 26 23 23
        >EAS54_6_R1_2_1_540_792
        26 26 26 26 26 26 26 26 26 26 26 22 26 26 26 26 26 12 26 26
        26 18 26 23 18
        >EAS54_6_R1_2_1_443_348
        26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 13 22 26 18
        24 18 18 18 18
    
    Using this module directly you might run:

    >>> handle = open("Quality/example.qual", "rU")
    >>> for record in QualPhredIterator(handle) :
    ...     print record.id, record.seq
    EAS54_6_R1_2_1_413_324 ?????????????????????????
    EAS54_6_R1_2_1_540_792 ?????????????????????????
    EAS54_6_R1_2_1_443_348 ?????????????????????????
    >>> handle.close()

    Typically however, you would call this via Bio.SeqIO instead with "qual"
    as the format:

    >>> from Bio import SeqIO
    >>> handle = open("Quality/example.qual", "rU")
    >>> for record in SeqIO.parse(handle, "qual") :
    ...     print record.id, record.seq
    EAS54_6_R1_2_1_413_324 ?????????????????????????
    EAS54_6_R1_2_1_540_792 ?????????????????????????
    EAS54_6_R1_2_1_443_348 ?????????????????????????
    >>> handle.close()

    Becase QUAL files don't contain the sequence string itself, the seq
    property is set to an UnknownSeq object.  As no alphabet was given, this
    has defaulted to a generic single letter alphabet and the character "?"
    used.

    By specifying a nucleotide alphabet, "N" is used instead:

    >>> from Bio import SeqIO
    >>> from Bio.Alphabet import generic_dna
    >>> handle = open("Quality/example.qual", "rU")
    >>> for record in SeqIO.parse(handle, "qual", alphabet=generic_dna) :
    ...     print record.id, record.seq
    EAS54_6_R1_2_1_413_324 NNNNNNNNNNNNNNNNNNNNNNNNN
    EAS54_6_R1_2_1_540_792 NNNNNNNNNNNNNNNNNNNNNNNNN
    EAS54_6_R1_2_1_443_348 NNNNNNNNNNNNNNNNNNNNNNNNN
    >>> handle.close()

    However, the quality scores themselves are available as a list of integers
    in each record's per-letter-annotation:

    >>> print record.letter_annotations["phred_quality"]
    [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

    You can still slice one of these SeqRecord objects with an UnknownSeq:

    >>> sub_record = record[5:10]
    >>> print sub_record.id, sub_record.letter_annotations["phred_quality"]
    EAS54_6_R1_2_1_443_348 [26, 26, 26, 26, 26]
    """
    #Skip any text before the first record (e.g. blank lines, comments)
    while True :
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line[0] == ">" :
            break

    while True :
        if line[0]!=">" :
            raise ValueError("Records in Fasta files should start with '>' character")
        if title2ids :
            id, name, descr = title2ids(line[1:].rstrip())
        else :
            descr = line[1:].rstrip()
            id   = descr.split()[0]
            name = id

        qualities = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">": break
            qualities.extend([int(word) for word in line.split()])
            line = handle.readline()

        if qualities :
            if min(qualities) < 0 or max(qualities) > 90 :
                raise ValueError(("Quality score range for %s is %i to %i, outside the " \
                                 +"expected 0 to 90.  Perhaps these are Solexa/Illumina " \
                                 +"scores, and not PHRED scores?") \
                                 % (id, min(qualities), max(qualities)))
        
        #Return the record and then continue...
        record = SeqRecord(UnknownSeq(len(qualities), alphabet),
                           id = id, name = name, description = descr)
        record.letter_annotations["phred_quality"] = qualities
        yield record

        if not line : return #StopIteration
    assert False, "Should not reach this line"

class FastqPhredWriter(SequentialSequenceWriter):
    """Class to write FASTQ format files (using PHRED quality scores).

    Although you can use this class directly, you are strongly encouraged
    to use the Bio.SeqIO.write() function instead.  For example, this code
    reads in a FASTQ (PHRED) file and re-saves it as another FASTQ (PHRED)
    file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse(open("Quality/example.fastq"), "fastq")
    >>> out_handle = open("Quality/temp.fastq", "w")
    >>> SeqIO.write(record_iterator, out_handle, "fastq")
    3
    >>> out_handle.close()

    You might want to do this if the original file included extra line breaks,
    which while valid may not be supported by all tools.  The output file from
    Biopython will have each sequence on a single line, and each quality
    string on a single line (which is considered desirable for maximum
    compatibility).

    In this next example, a Solexa FASTQ file is converted into a standard
    Sanger style FASTQ file using PHRED qualities:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse(open("Quality/solexa.fastq"), "fastq-solexa")
    >>> out_handle = open("Quality/temp.fastq", "w")
    >>> SeqIO.write(record_iterator, out_handle, "fastq")
    1
    >>> out_handle.close()

    This code is also called if you use the .format("fastq") method of a
    SeqRecord.

    P.S. To avoid cluttering up your working directory, you can delete this
    temporary file now:

    >>> import os
    >>> os.remove("Quality/temp.fastq")

    """
    def write_record(self, record):
        """Write a single FASTQ record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        #TODO - Is an empty sequence allowed in FASTQ format?
        assert SANGER_SCORE_OFFSET == ord("!")
        #This rounds to the nearest integer:
        qualities = "".join([chr(int(round(q+SANGER_SCORE_OFFSET,0))) for q \
                             in _get_phred_quality(record)])
        if record.seq is None:
            raise ValueError("No sequence for record %s" % record.id)
        if len(qualities) != len(record) :
            raise ValueError("Record %s has sequence length %i but %i quality scores" \
                             % (record.id, len(record), len(qualities)))

        title = self.clean(record.id) #TODO - add the description too? cf Fasta output
        self.handle.write("@%s\n%s\n+\n%s\n" % (title, record.seq, qualities))

class QualPhredWriter(SequentialSequenceWriter):
    """Class to write QUAL format files (using PHRED quality scores).

    Although you can use this class directly, you are strongly encouraged
    to use the Bio.SeqIO.write() function instead.  For example, this code
    reads in a FASTQ file and saves the quality scores into a QUAL file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse(open("Quality/example.fastq"), "fastq")
    >>> out_handle = open("Quality/temp.qual", "w")
    >>> SeqIO.write(record_iterator, out_handle, "qual")
    3
    >>> out_handle.close()

    This code is also called if you use the .format("qual") method of a
    SeqRecord.

    P.S. Don't forget to clean up the temp file if you don't need it anymore:

    >>> import os
    >>> os.remove("Quality/temp.qual")
    """
    def __init__(self, handle, wrap=60, record2title=None):
        """Create a QUAL writer.

        Arguments:
         - handle - Handle to an output file, e.g. as returned
                    by open(filename, "w")
         - wrap   - Optional line length used to wrap sequence lines.
                    Defaults to wrapping the sequence at 60 characters
                    Use zero (or None) for no wrapping, giving a single
                    long line for the sequence.
         - record2title - Optional function to return the text to be
                    used for the title line of each record.  By default
                    a combination of the record.id and record.description
                    is used.  If the record.description starts with the
                    record.id, then just the record.description is used.

        The record2title argument is present for consistency with the
        Bio.SeqIO.FastaIO writer class.
        """
        SequentialSequenceWriter.__init__(self, handle)
        #self.handle = handle
        self.wrap = None
        if wrap :
            if wrap < 1 :
                raise ValueError
        self.wrap = wrap
        self.record2title = record2title

    def write_record(self, record):
        """Write a single QUAL record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        if self.record2title :
            title=self.clean(self.record2title(record))
        else :
            id = self.clean(record.id)
            description = self.clean(record.description)

            #if description[:len(id)]==id :
            if description and description.split(None,1)[0]==id :
                #The description includes the id at the start
                title = description
            else :
                title = "%s %s" % (id, description)

        assert "\n" not in title
        assert "\r" not in title
        self.handle.write(">%s\n" % title)

        #This rounds to the nearest integer.
        #TODO - can we put a float in a qual file?
        qualities = [("%i" % round(q,0)) for q in _get_phred_quality(record)]

        if self.wrap :
            while qualities :
                line=qualities.pop(0)
                while qualities \
                and len(line) + 1 + len(qualities[0]) < self.wrap :
                    line += " " + qualities.pop(0)
                self.handle.write(line + "\n")
        else :
            data = " ".join(qualities)
            self.handle.write(data + "\n")

class FastqSolexaWriter(SequentialSequenceWriter):
    """Class to write FASTQ format files (using Solexa quality scores).

    Although you can use this class directly, you are strongly encouraged
    to use the Bio.SeqIO.write() function instead.  For example, this code
    reads in a FASTQ file and re-saves it as another FASTQ file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse(open("Quality/solexa.fastq"), "fastq-solexa")
    >>> out_handle = open("Quality/temp.fastq", "w")
    >>> SeqIO.write(record_iterator, out_handle, "fastq-solexa")
    1
    >>> out_handle.close()

    You might want to do this if the original file included extra line
    breaks, which (while valid) may not be supported by all tools.  The
    output file from Biopython will have each sequence on a single line, and
    each quality string on a single line (which is considered desirable for
    maximum compatibility).

    This code is also called if you use the .format("fastq-solexa") method of
    a SeqRecord.

    P.S. Don't forget to delete the temp file if you don't need it anymore:

    >>> import os
    >>> os.remove("Quality/temp.fastq")
    """
    def write_record(self, record):
        """Write a single FASTQ record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        #TODO - Is an empty sequence allowed in FASTQ format?
        qualities = "".join([chr(int(round(q+SOLEXA_SCORE_OFFSET,0))) for q \
                             in _get_solexa_quality(record)])
        if record.seq is None:
            raise ValueError("No sequence for record %s" % record.id)
        if len(qualities) != len(record) :
            raise ValueError("Record %s has sequence length %i but %i quality scores" \
                             % (record.id, len(record), len(qualities)))

        title = self.clean(record.id) #TODO - add the description too? cf Fasta output
        self.handle.write("@%s\n%s\n+\n%s\n" % (title, record.seq, qualities))
        
def PairedFastaQualIterator(fasta_handle, qual_handle, alphabet = single_letter_alphabet, title2ids = None) :
    """Iterate over matched FASTA and QUAL files as SeqRecord objects.

    For example, consider this short QUAL file::

        >EAS54_6_R1_2_1_413_324
        26 26 18 26 26 26 26 26 26 26 26 26 26 26 26 22 26 26 26 26
        26 26 26 23 23
        >EAS54_6_R1_2_1_540_792
        26 26 26 26 26 26 26 26 26 26 26 22 26 26 26 26 26 12 26 26
        26 18 26 23 18
        >EAS54_6_R1_2_1_443_348
        26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 13 22 26 18
        24 18 18 18 18
    
    And a matching FASTA file::

        >EAS54_6_R1_2_1_413_324
        CCCTTCTTGTCTTCAGCGTTTCTCC
        >EAS54_6_R1_2_1_540_792
        TTGGCAGGCCAAGGCCGATGGATCA
        >EAS54_6_R1_2_1_443_348
        GTTGCTTCTGGCGTGGGTGGGGGGG

    You can parse these separately using Bio.SeqIO with the "qual" and
    "fasta" formats, but then you'll get a group of SeqRecord objects with
    no sequence, and a matching group with the sequence but not the
    qualities.  Because it only deals with one input file handle, Bio.SeqIO
    can't be used to read the two files together - but this function can!
    For example,
    
    >>> rec_iter = PairedFastaQualIterator(open("Quality/example.fasta", "rU"),
    ...                                    open("Quality/example.qual", "rU"))
    >>> for record in rec_iter :
    ...     print record.id, record.seq
    EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
    EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
    EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG

    As with the FASTQ or QUAL parsers, if you want to look at the qualities,
    they are in each record's per-letter-annotation dictionary as a simple
    list of integers:

    >>> print record.letter_annotations["phred_quality"]
    [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

    If you have access to data as a FASTQ format file, using that directly
    would be simpler and more straight forward.  Note that you can easily use
    this function to convert paired FASTA and QUAL files into FASTQ files:

    >>> from Bio import SeqIO
    >>> rec_iter = PairedFastaQualIterator(open("Quality/example.fasta", "rU"),
    ...                                    open("Quality/example.qual", "rU"))
    >>> out_handle = open("Quality/temp.fastq", "w")
    >>> SeqIO.write(rec_iter, out_handle, "fastq")
    3
    >>> out_handle.close()

    And don't forget to clean up the temp file if you don't need it anymore:

    >>> import os
    >>> os.remove("Quality/temp.fastq")    
    """
    from Bio.SeqIO.FastaIO import FastaIterator    
    fasta_iter = FastaIterator(fasta_handle, alphabet=alphabet, \
                               title2ids=title2ids)
    qual_iter = QualPhredIterator(qual_handle, alphabet=alphabet, \
                                  title2ids=title2ids)

    #Using zip(...) would create a list loading everything into memory!
    #It would also not catch any extra records found in only one file.
    while True :
        try :
            f_rec = fasta_iter.next()
        except StopIteration :
            f_rec = None
        try :
            q_rec = qual_iter.next()
        except StopIteration :
            q_rec = None
        if f_rec is None and q_rec is None :
            #End of both files
            break
        if f_rec is None :
            raise ValueError("FASTA file has more entries than the QUAL file.")
        if q_rec is None :
            raise ValueError("QUAL file has more entries than the FASTA file.")
        if f_rec.id != q_rec.id :
            raise ValueError("FASTA and QUAL entries do not match (%s vs %s)." \
                             % (f_rec.id, q_rec.id))
        if len(f_rec) != len(q_rec.letter_annotations["phred_quality"]) :
            raise ValueError("Sequence length and number of quality scores disagree for %s" \
                             % f_rec.id)
        #Merge the data....
        f_rec.letter_annotations["phred_quality"] = q_rec.letter_annotations["phred_quality"]
        yield f_rec
    #Done
    

def _test():
    """Run the Bio.SeqIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","Tests")) :
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","Tests"))
        assert os.path.isfile("Quality/example.fastq")
        assert os.path.isfile("Quality/example.fasta")
        assert os.path.isfile("Quality/example.qual")
        assert os.path.isfile("Quality/tricky.fastq")
        assert os.path.isfile("Quality/solexa.fastq")
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
        
if __name__ == "__main__" :
    _test()

