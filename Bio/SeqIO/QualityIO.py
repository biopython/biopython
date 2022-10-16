# Copyright 2009-2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the FASTQ and QUAL file formats.

Note that you are expected to use this code via the Bio.SeqIO interface, as
shown below.

The FASTQ file format is used frequently at the Wellcome Trust Sanger Institute
to bundle a FASTA sequence and its PHRED quality data (integers between 0 and
90).  Rather than using a single FASTQ file, often paired FASTA and QUAL files
are used containing the sequence and the quality information separately.

The PHRED software reads DNA sequencing trace files, calls bases, and
assigns a non-negative quality value to each called base using a logged
transformation of the error probability, Q = -10 log10( Pe ), for example::

    Pe = 1.0,         Q =  0
    Pe = 0.1,         Q = 10
    Pe = 0.01,        Q = 20
    ...
    Pe = 0.00000001,  Q = 80
    Pe = 0.000000001, Q = 90

In typical raw sequence reads, the PHRED quality valuea will be from 0 to 40.
In the QUAL format these quality values are held as space separated text in
a FASTA like file format.  In the FASTQ format, each quality values is encoded
with a single ASCI character using chr(Q+33), meaning zero maps to the
character "!" and for example 80 maps to "q".  For the Sanger FASTQ standard
the allowed range of PHRED scores is 0 to 93 inclusive. The sequences and
quality are then stored in pairs in a FASTA like format.

Unfortunately there is no official document describing the FASTQ file format,
and worse, several related but different variants exist. For more details,
please read this open access publication::

    The Sanger FASTQ file format for sequences with quality scores, and the
    Solexa/Illumina FASTQ variants.
    P.J.A.Cock (Biopython), C.J.Fields (BioPerl), N.Goto (BioRuby),
    M.L.Heuer (BioJava) and P.M. Rice (EMBOSS).
    Nucleic Acids Research 2010 38(6):1767-1771
    https://doi.org/10.1093/nar/gkp1137

The good news is that Roche 454 sequencers can output files in the QUAL format,
and sensibly they use PHREP style scores like Sanger.  Converting a pair of
FASTA and QUAL files into a Sanger style FASTQ file is easy. To extract QUAL
files from a Roche 454 SFF binary file, use the Roche off instrument command
line tool "sffinfo" with the -q or -qual argument.  You can extract a matching
FASTA file using the -s or -seq argument instead.

The bad news is that Solexa/Illumina did things differently - they have their
own scoring system AND their own incompatible versions of the FASTQ format.
Solexa/Illumina quality scores use Q = - 10 log10 ( Pe / (1-Pe) ), which can
be negative.  PHRED scores and Solexa scores are NOT interchangeable (but a
reasonable mapping can be achieved between them, and they are approximately
equal for higher quality reads).

Confusingly early Solexa pipelines produced a FASTQ like file but using their
own score mapping and an ASCII offset of 64. To make things worse, for the
Solexa/Illumina pipeline 1.3 onwards, they introduced a third variant of the
FASTQ file format, this time using PHRED scores (which is more consistent) but
with an ASCII offset of 64.

i.e. There are at least THREE different and INCOMPATIBLE variants of the FASTQ
file format: The original Sanger PHRED standard, and two from Solexa/Illumina.

The good news is that as of CASAVA version 1.8, Illumina sequencers will
produce FASTQ files using the standard Sanger encoding.

You are expected to use this module via the Bio.SeqIO functions, with the
following format names:

    - "qual" means simple quality files using PHRED scores (e.g. from Roche 454)
    - "fastq" means Sanger style FASTQ files using PHRED scores and an ASCII
      offset of 33 (e.g. from the NCBI Short Read Archive and Illumina 1.8+).
      These can potentially hold PHRED scores from 0 to 93.
    - "fastq-sanger" is an alias for "fastq".
    - "fastq-solexa" means old Solexa (and also very early Illumina) style FASTQ
      files, using Solexa scores with an ASCII offset 64. These can hold Solexa
      scores from -5 to 62.
    - "fastq-illumina" means newer Illumina 1.3 to 1.7 style FASTQ files, using
      PHRED scores but with an ASCII offset 64, allowing PHRED scores from 0
      to 62.

We could potentially add support for "qual-solexa" meaning QUAL files which
contain Solexa scores, but thus far there isn't any reason to use such files.

For example, consider the following short FASTQ file::

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
probably originally from an early Solexa/Illumina sequencer but this file
follows the Sanger FASTQ convention (PHRED style qualities with an ASCII
offset of 33).  This means we can parse this file using Bio.SeqIO using
"fastq" as the format name:

>>> from Bio import SeqIO
>>> for record in SeqIO.parse("Quality/example.fastq", "fastq"):
...     print("%s %s" % (record.id, record.seq))
EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG

The qualities are held as a list of integers in each record's annotation:

>>> print(record)
ID: EAS54_6_R1_2_1_443_348
Name: EAS54_6_R1_2_1_443_348
Description: EAS54_6_R1_2_1_443_348
Number of features: 0
Per letter annotation for: phred_quality
Seq('GTTGCTTCTGGCGTGGGTGGGGGGG')
>>> print(record.letter_annotations["phred_quality"])
[26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

You can use the SeqRecord format method to show this in the QUAL format:

>>> print(record.format("qual"))
>EAS54_6_R1_2_1_443_348
26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 13 22 26 18
24 18 18 18 18
<BLANKLINE>

Or go back to the FASTQ format, use "fastq" (or "fastq-sanger"):

>>> print(record.format("fastq"))
@EAS54_6_R1_2_1_443_348
GTTGCTTCTGGCGTGGGTGGGGGGG
+
;;;;;;;;;;;9;7;;.7;393333
<BLANKLINE>

Or, using the Illumina 1.3+ FASTQ encoding (PHRED values with an ASCII offset
of 64):

>>> print(record.format("fastq-illumina"))
@EAS54_6_R1_2_1_443_348
GTTGCTTCTGGCGTGGGTGGGGGGG
+
ZZZZZZZZZZZXZVZZMVZRXRRRR
<BLANKLINE>

You can also get Biopython to convert the scores and show a Solexa style
FASTQ file:

>>> print(record.format("fastq-solexa"))
@EAS54_6_R1_2_1_443_348
GTTGCTTCTGGCGTGGGTGGGGGGG
+
ZZZZZZZZZZZXZVZZMVZRXRRRR
<BLANKLINE>

Notice that this is actually the same output as above using "fastq-illumina"
as the format! The reason for this is all these scores are high enough that
the PHRED and Solexa scores are almost equal. The differences become apparent
for poor quality reads. See the functions solexa_quality_from_phred and
phred_quality_from_solexa for more details.

If you wanted to trim your sequences (perhaps to remove low quality regions,
or to remove a primer sequence), try slicing the SeqRecord objects.  e.g.

>>> sub_rec = record[5:15]
>>> print(sub_rec)
ID: EAS54_6_R1_2_1_443_348
Name: EAS54_6_R1_2_1_443_348
Description: EAS54_6_R1_2_1_443_348
Number of features: 0
Per letter annotation for: phred_quality
Seq('TTCTGGCGTG')
>>> print(sub_rec.letter_annotations["phred_quality"])
[26, 26, 26, 26, 26, 26, 24, 26, 22, 26]
>>> print(sub_rec.format("fastq"))
@EAS54_6_R1_2_1_443_348
TTCTGGCGTG
+
;;;;;;9;7;
<BLANKLINE>

If you wanted to, you could read in this FASTQ file, and save it as a QUAL file:

>>> from Bio import SeqIO
>>> record_iterator = SeqIO.parse("Quality/example.fastq", "fastq")
>>> with open("Quality/temp.qual", "w") as out_handle:
...     SeqIO.write(record_iterator, out_handle, "qual")
3

You can of course read in a QUAL file, such as the one we just created:

>>> from Bio import SeqIO
>>> for record in SeqIO.parse("Quality/temp.qual", "qual"):
...     print("%s read of length %d" % (record.id, len(record.seq)))
EAS54_6_R1_2_1_413_324 read of length 25
EAS54_6_R1_2_1_540_792 read of length 25
EAS54_6_R1_2_1_443_348 read of length 25

Notice that QUAL files don't have a proper sequence present!  But the quality
information is there:

>>> print(record)
ID: EAS54_6_R1_2_1_443_348
Name: EAS54_6_R1_2_1_443_348
Description: EAS54_6_R1_2_1_443_348
Number of features: 0
Per letter annotation for: phred_quality
Undefined sequence of length 25
>>> print(record.letter_annotations["phred_quality"])
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
>>> reads = SeqIO.to_dict(SeqIO.parse("Quality/example.fasta", "fasta"))
>>> for rec in SeqIO.parse("Quality/example.qual", "qual"):
...     reads[rec.id].letter_annotations["phred_quality"]=rec.letter_annotations["phred_quality"]

You can then access any record by its key, and get both the sequence and the
quality scores.

>>> print(reads["EAS54_6_R1_2_1_540_792"].format("fastq"))
@EAS54_6_R1_2_1_540_792
TTGGCAGGCCAAGGCCGATGGATCA
+
;;;;;;;;;;;7;;;;;-;;;3;83
<BLANKLINE>

It is important that you explicitly tell Bio.SeqIO which FASTQ variant you are
using ("fastq" or "fastq-sanger" for the Sanger standard using PHRED values,
"fastq-solexa" for the original Solexa/Illumina variant, or "fastq-illumina"
for the more recent variant), as this cannot be detected reliably
automatically.

To illustrate this problem, let's consider an artificial example:

>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> test = SeqRecord(Seq("NACGTACGTA"), id="Test", description="Made up!")
>>> print(test.format("fasta"))
>Test Made up!
NACGTACGTA
<BLANKLINE>
>>> print(test.format("fastq"))
Traceback (most recent call last):
 ...
ValueError: No suitable quality scores found in letter_annotations of SeqRecord (id=Test).

We created a sample SeqRecord, and can show it in FASTA format - but for QUAL
or FASTQ format we need to provide some quality scores. These are held as a
list of integers (one for each base) in the letter_annotations dictionary:

>>> test.letter_annotations["phred_quality"] = [0, 1, 2, 3, 4, 5, 10, 20, 30, 40]
>>> print(test.format("qual"))
>Test Made up!
0 1 2 3 4 5 10 20 30 40
<BLANKLINE>
>>> print(test.format("fastq"))
@Test Made up!
NACGTACGTA
+
!"#$%&+5?I
<BLANKLINE>

We can check this FASTQ encoding - the first PHRED quality was zero, and this
mapped to a exclamation mark, while the final score was 40 and this mapped to
the letter "I":

>>> ord('!') - 33
0
>>> ord('I') - 33
40
>>> [ord(letter)-33 for letter in '!"#$%&+5?I']
[0, 1, 2, 3, 4, 5, 10, 20, 30, 40]

Similarly, we could produce an Illumina 1.3 to 1.7 style FASTQ file using PHRED
scores with an offset of 64:

>>> print(test.format("fastq-illumina"))
@Test Made up!
NACGTACGTA
+
@ABCDEJT^h
<BLANKLINE>

And we can check this too - the first PHRED score was zero, and this mapped to
"@", while the final score was 40 and this mapped to "h":

>>> ord("@") - 64
0
>>> ord("h") - 64
40
>>> [ord(letter)-64 for letter in "@ABCDEJT^h"]
[0, 1, 2, 3, 4, 5, 10, 20, 30, 40]

Notice how different the standard Sanger FASTQ and the Illumina 1.3 to 1.7 style
FASTQ files look for the same data! Then we have the older Solexa/Illumina
format to consider which encodes Solexa scores instead of PHRED scores.

First let's see what Biopython says if we convert the PHRED scores into Solexa
scores (rounding to one decimal place):

>>> for q in [0, 1, 2, 3, 4, 5, 10, 20, 30, 40]:
...     print("PHRED %i maps to Solexa %0.1f" % (q, solexa_quality_from_phred(q)))
PHRED 0 maps to Solexa -5.0
PHRED 1 maps to Solexa -5.0
PHRED 2 maps to Solexa -2.3
PHRED 3 maps to Solexa -0.0
PHRED 4 maps to Solexa 1.8
PHRED 5 maps to Solexa 3.3
PHRED 10 maps to Solexa 9.5
PHRED 20 maps to Solexa 20.0
PHRED 30 maps to Solexa 30.0
PHRED 40 maps to Solexa 40.0

Now here is the record using the old Solexa style FASTQ file:

>>> print(test.format("fastq-solexa"))
@Test Made up!
NACGTACGTA
+
;;>@BCJT^h
<BLANKLINE>

Again, this is using an ASCII offset of 64, so we can check the Solexa scores:

>>> [ord(letter)-64 for letter in ";;>@BCJT^h"]
[-5, -5, -2, 0, 2, 3, 10, 20, 30, 40]

This explains why the last few letters of this FASTQ output matched that using
the Illumina 1.3 to 1.7 format - high quality PHRED scores and Solexa scores
are approximately equal.

"""
import warnings

from math import log

from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
from Bio import BiopythonDeprecationWarning
from Bio import StreamModeError
from Bio.File import as_handle
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import _clean
from .Interfaces import _get_seq_string
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter


# define score offsets. See discussion for differences between Sanger and
# Solexa offsets.
SANGER_SCORE_OFFSET = 33
SOLEXA_SCORE_OFFSET = 64


def solexa_quality_from_phred(phred_quality):
    """Convert a PHRED quality (range 0 to about 90) to a Solexa quality.

    PHRED and Solexa quality scores are both log transformations of a
    probality of error (high score = low probability of error). This function
    takes a PHRED score, transforms it back to a probability of error, and
    then re-expresses it as a Solexa score. This assumes the error estimates
    are equivalent.

    How does this work exactly? Well the PHRED quality is minus ten times the
    base ten logarithm of the probability of error::

        phred_quality = -10*log(error,10)

    Therefore, turning this round::

        error = 10 ** (- phred_quality / 10)

    Now, Solexa qualities use a different log transformation::

        solexa_quality = -10*log(error/(1-error),10)

    After substitution and a little manipulation we get::

         solexa_quality = 10*log(10**(phred_quality/10.0) - 1, 10)

    However, real Solexa files use a minimum quality of -5. This does have a
    good reason - a random base call would be correct 25% of the time,
    and thus have a probability of error of 0.75, which gives 1.25 as the PHRED
    quality, or -4.77 as the Solexa quality. Thus (after rounding), a random
    nucleotide read would have a PHRED quality of 1, or a Solexa quality of -5.

    Taken literally, this logarithic formula would map a PHRED quality of zero
    to a Solexa quality of minus infinity. Of course, taken literally, a PHRED
    score of zero means a probability of error of one (i.e. the base call is
    definitely wrong), which is worse than random! In practice, a PHRED quality
    of zero usually means a default value, or perhaps random - and therefore
    mapping it to the minimum Solexa score of -5 is reasonable.

    In conclusion, we follow EMBOSS, and take this logarithmic formula but also
    apply a minimum value of -5.0 for the Solexa quality, and also map a PHRED
    quality of zero to -5.0 as well.

    Note this function will return a floating point number, it is up to you to
    round this to the nearest integer if appropriate.  e.g.

    >>> print("%0.2f" % round(solexa_quality_from_phred(80), 2))
    80.00
    >>> print("%0.2f" % round(solexa_quality_from_phred(50), 2))
    50.00
    >>> print("%0.2f" % round(solexa_quality_from_phred(20), 2))
    19.96
    >>> print("%0.2f" % round(solexa_quality_from_phred(10), 2))
    9.54
    >>> print("%0.2f" % round(solexa_quality_from_phred(5), 2))
    3.35
    >>> print("%0.2f" % round(solexa_quality_from_phred(4), 2))
    1.80
    >>> print("%0.2f" % round(solexa_quality_from_phred(3), 2))
    -0.02
    >>> print("%0.2f" % round(solexa_quality_from_phred(2), 2))
    -2.33
    >>> print("%0.2f" % round(solexa_quality_from_phred(1), 2))
    -5.00
    >>> print("%0.2f" % round(solexa_quality_from_phred(0), 2))
    -5.00

    Notice that for high quality reads PHRED and Solexa scores are numerically
    equal. The differences are important for poor quality reads, where PHRED
    has a minimum of zero but Solexa scores can be negative.

    Finally, as a special case where None is used for a "missing value", None
    is returned:

    >>> print(solexa_quality_from_phred(None))
    None
    """
    if phred_quality is None:
        # Assume None is used as some kind of NULL or NA value; return None
        # e.g. Bio.SeqIO gives Ace contig gaps a quality of None.
        return None
    elif phred_quality > 0:
        # Solexa uses a minimum value of -5, which after rounding matches a
        # random nucleotide base call.
        return max(-5.0, 10 * log(10 ** (phred_quality / 10.0) - 1, 10))
    elif phred_quality == 0:
        # Special case, map to -5 as discussed in the docstring
        return -5.0
    else:
        raise ValueError(
            f"PHRED qualities must be positive (or zero), not {phred_quality!r}"
        )


def phred_quality_from_solexa(solexa_quality):
    """Convert a Solexa quality (which can be negative) to a PHRED quality.

    PHRED and Solexa quality scores are both log transformations of a
    probality of error (high score = low probability of error). This function
    takes a Solexa score, transforms it back to a probability of error, and
    then re-expresses it as a PHRED score. This assumes the error estimates
    are equivalent.

    The underlying formulas are given in the documentation for the sister
    function solexa_quality_from_phred, in this case the operation is::

        phred_quality = 10*log(10**(solexa_quality/10.0) + 1, 10)

    This will return a floating point number, it is up to you to round this to
    the nearest integer if appropriate.  e.g.

    >>> print("%0.2f" % round(phred_quality_from_solexa(80), 2))
    80.00
    >>> print("%0.2f" % round(phred_quality_from_solexa(20), 2))
    20.04
    >>> print("%0.2f" % round(phred_quality_from_solexa(10), 2))
    10.41
    >>> print("%0.2f" % round(phred_quality_from_solexa(0), 2))
    3.01
    >>> print("%0.2f" % round(phred_quality_from_solexa(-5), 2))
    1.19

    Note that a solexa_quality less then -5 is not expected, will trigger a
    warning, but will still be converted as per the logarithmic mapping
    (giving a number between 0 and 1.19 back).

    As a special case where None is used for a "missing value", None is
    returned:

    >>> print(phred_quality_from_solexa(None))
    None
    """
    if solexa_quality is None:
        # Assume None is used as some kind of NULL or NA value; return None
        return None
    if solexa_quality < -5:
        warnings.warn(
            f"Solexa quality less than -5 passed, {solexa_quality!r}", BiopythonWarning
        )
    return 10 * log(10 ** (solexa_quality / 10.0) + 1, 10)


def _get_phred_quality(record):
    """Extract PHRED qualities from a SeqRecord's letter_annotations (PRIVATE).

    If there are no PHRED qualities, but there are Solexa qualities, those are
    used instead after conversion.
    """
    try:
        return record.letter_annotations["phred_quality"]
    except KeyError:
        pass
    try:
        return [
            phred_quality_from_solexa(q)
            for q in record.letter_annotations["solexa_quality"]
        ]
    except KeyError:
        raise ValueError(
            "No suitable quality scores found in "
            "letter_annotations of SeqRecord (id=%s)." % record.id
        ) from None


# Only map 0 to 93, we need to give a warning on truncating at 93
_phred_to_sanger_quality_str = {
    qp: chr(min(126, qp + SANGER_SCORE_OFFSET)) for qp in range(0, 93 + 1)
}
# Only map -5 to 93, we need to give a warning on truncating at 93
_solexa_to_sanger_quality_str = {
    qs: chr(min(126, int(round(phred_quality_from_solexa(qs)) + SANGER_SCORE_OFFSET)))
    for qs in range(-5, 93 + 1)
}


def _get_sanger_quality_str(record):
    """Return a Sanger FASTQ encoded quality string (PRIVATE).

    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> r = SeqRecord(Seq("ACGTAN"), id="Test",
    ...               letter_annotations = {"phred_quality":[50, 40, 30, 20, 10, 0]})
    >>> _get_sanger_quality_str(r)
    'SI?5+!'

    If as in the above example (or indeed a SeqRecord parser with Bio.SeqIO),
    the PHRED qualities are integers, this function is able to use a very fast
    pre-cached mapping. However, if they are floats which differ slightly, then
    it has to do the appropriate rounding - which is slower:

    >>> r2 = SeqRecord(Seq("ACGTAN"), id="Test2",
    ...      letter_annotations = {"phred_quality":[50.0, 40.05, 29.99, 20, 9.55, 0.01]})
    >>> _get_sanger_quality_str(r2)
    'SI?5+!'

    If your scores include a None value, this raises an exception:

    >>> r3 = SeqRecord(Seq("ACGTAN"), id="Test3",
    ...               letter_annotations = {"phred_quality":[50, 40, 30, 20, 10, None]})
    >>> _get_sanger_quality_str(r3)
    Traceback (most recent call last):
       ...
    TypeError: A quality value of None was found

    If (strangely) your record has both PHRED and Solexa scores, then the PHRED
    scores are used in preference:

    >>> r4 = SeqRecord(Seq("ACGTAN"), id="Test4",
    ...               letter_annotations = {"phred_quality":[50, 40, 30, 20, 10, 0],
    ...                                     "solexa_quality":[-5, -4, 0, None, 0, 40]})
    >>> _get_sanger_quality_str(r4)
    'SI?5+!'

    If there are no PHRED scores, but there are Solexa scores, these are used
    instead (after the appropriate conversion):

    >>> r5 = SeqRecord(Seq("ACGTAN"), id="Test5",
    ...      letter_annotations = {"solexa_quality":[40, 30, 20, 10, 0, -5]})
    >>> _get_sanger_quality_str(r5)
    'I?5+$"'

    Again, integer Solexa scores can be looked up in a pre-cached mapping making
    this very fast. You can still use approximate floating point scores:

    >>> r6 = SeqRecord(Seq("ACGTAN"), id="Test6",
    ...      letter_annotations = {"solexa_quality":[40.1, 29.7, 20.01, 10, 0.0, -4.9]})
    >>> _get_sanger_quality_str(r6)
    'I?5+$"'

    Notice that due to the limited range of printable ASCII characters, a
    PHRED quality of 93 is the maximum that can be held in an Illumina FASTQ
    file (using ASCII 126, the tilde). This function will issue a warning
    in this situation.
    """
    # TODO - This functions works and is fast, but it is also ugly
    # and there is considerable repetition of code for the other
    # two FASTQ variants.
    try:
        # These take priority (in case both Solexa and PHRED scores found)
        qualities = record.letter_annotations["phred_quality"]
    except KeyError:
        # Fall back on solexa scores...
        pass
    else:
        # Try and use the precomputed mapping:
        try:
            return "".join(_phred_to_sanger_quality_str[qp] for qp in qualities)
        except KeyError:
            # Could be a float, or a None in the list, or a high value.
            pass
        if None in qualities:
            raise TypeError("A quality value of None was found")
        if max(qualities) >= 93.5:
            warnings.warn(
                "Data loss - max PHRED quality 93 in Sanger FASTQ", BiopythonWarning
            )
        # This will apply the truncation at 93, giving max ASCII 126
        return "".join(
            chr(min(126, int(round(qp)) + SANGER_SCORE_OFFSET)) for qp in qualities
        )
    # Fall back on the Solexa scores...
    try:
        qualities = record.letter_annotations["solexa_quality"]
    except KeyError:
        raise ValueError(
            "No suitable quality scores found in "
            "letter_annotations of SeqRecord (id=%s)." % record.id
        ) from None
    # Try and use the precomputed mapping:
    try:
        return "".join(_solexa_to_sanger_quality_str[qs] for qs in qualities)
    except KeyError:
        # Either no PHRED scores, or something odd like a float or None
        pass
    if None in qualities:
        raise TypeError("A quality value of None was found")
    # Must do this the slow way, first converting the PHRED scores into
    # Solexa scores:
    if max(qualities) >= 93.5:
        warnings.warn(
            "Data loss - max PHRED quality 93 in Sanger FASTQ", BiopythonWarning
        )
    # This will apply the truncation at 93, giving max ASCII 126
    return "".join(
        chr(min(126, int(round(phred_quality_from_solexa(qs))) + SANGER_SCORE_OFFSET))
        for qs in qualities
    )


# Only map 0 to 62, we need to give a warning on truncating at 62
assert 62 + SOLEXA_SCORE_OFFSET == 126
_phred_to_illumina_quality_str = {
    qp: chr(qp + SOLEXA_SCORE_OFFSET) for qp in range(0, 62 + 1)
}
# Only map -5 to 62, we need to give a warning on truncating at 62
_solexa_to_illumina_quality_str = {
    qs: chr(int(round(phred_quality_from_solexa(qs))) + SOLEXA_SCORE_OFFSET)
    for qs in range(-5, 62 + 1)
}


def _get_illumina_quality_str(record):
    """Return an Illumina 1.3 to 1.7 FASTQ encoded quality string (PRIVATE).

    Notice that due to the limited range of printable ASCII characters, a
    PHRED quality of 62 is the maximum that can be held in an Illumina FASTQ
    file (using ASCII 126, the tilde). This function will issue a warning
    in this situation.
    """
    # TODO - This functions works and is fast, but it is also ugly
    # and there is considerable repetition of code for the other
    # two FASTQ variants.
    try:
        # These take priority (in case both Solexa and PHRED scores found)
        qualities = record.letter_annotations["phred_quality"]
    except KeyError:
        # Fall back on solexa scores...
        pass
    else:
        # Try and use the precomputed mapping:
        try:
            return "".join(_phred_to_illumina_quality_str[qp] for qp in qualities)
        except KeyError:
            # Could be a float, or a None in the list, or a high value.
            pass
        if None in qualities:
            raise TypeError("A quality value of None was found")
        if max(qualities) >= 62.5:
            warnings.warn(
                "Data loss - max PHRED quality 62 in Illumina FASTQ", BiopythonWarning
            )
        # This will apply the truncation at 62, giving max ASCII 126
        return "".join(
            chr(min(126, int(round(qp)) + SOLEXA_SCORE_OFFSET)) for qp in qualities
        )
    # Fall back on the Solexa scores...
    try:
        qualities = record.letter_annotations["solexa_quality"]
    except KeyError:
        raise ValueError(
            "No suitable quality scores found in "
            "letter_annotations of SeqRecord (id=%s)." % record.id
        ) from None
    # Try and use the precomputed mapping:
    try:
        return "".join(_solexa_to_illumina_quality_str[qs] for qs in qualities)
    except KeyError:
        # Either no PHRED scores, or something odd like a float or None
        pass
    if None in qualities:
        raise TypeError("A quality value of None was found")
    # Must do this the slow way, first converting the PHRED scores into
    # Solexa scores:
    if max(qualities) >= 62.5:
        warnings.warn(
            "Data loss - max PHRED quality 62 in Illumina FASTQ", BiopythonWarning
        )
    # This will apply the truncation at 62, giving max ASCII 126
    return "".join(
        chr(min(126, int(round(phred_quality_from_solexa(qs))) + SOLEXA_SCORE_OFFSET))
        for qs in qualities
    )


# Only map 0 to 62, we need to give a warning on truncating at 62
assert 62 + SOLEXA_SCORE_OFFSET == 126
_solexa_to_solexa_quality_str = {
    qs: chr(min(126, qs + SOLEXA_SCORE_OFFSET)) for qs in range(-5, 62 + 1)
}
# Only map -5 to 62, we need to give a warning on truncating at 62
_phred_to_solexa_quality_str = {
    qp: chr(min(126, int(round(solexa_quality_from_phred(qp))) + SOLEXA_SCORE_OFFSET))
    for qp in range(0, 62 + 1)
}


def _get_solexa_quality_str(record):
    """Return a Solexa FASTQ encoded quality string (PRIVATE).

    Notice that due to the limited range of printable ASCII characters, a
    Solexa quality of 62 is the maximum that can be held in a Solexa FASTQ
    file (using ASCII 126, the tilde). This function will issue a warning
    in this situation.
    """
    # TODO - This functions works and is fast, but it is also ugly
    # and there is considerable repetition of code for the other
    # two FASTQ variants.
    try:
        # These take priority (in case both Solexa and PHRED scores found)
        qualities = record.letter_annotations["solexa_quality"]
    except KeyError:
        # Fall back on PHRED scores...
        pass
    else:
        # Try and use the precomputed mapping:
        try:
            return "".join(_solexa_to_solexa_quality_str[qs] for qs in qualities)
        except KeyError:
            # Could be a float, or a None in the list, or a high value.
            pass
        if None in qualities:
            raise TypeError("A quality value of None was found")
        if max(qualities) >= 62.5:
            warnings.warn(
                "Data loss - max Solexa quality 62 in Solexa FASTQ", BiopythonWarning
            )
        # This will apply the truncation at 62, giving max ASCII 126
        return "".join(
            chr(min(126, int(round(qs)) + SOLEXA_SCORE_OFFSET)) for qs in qualities
        )
    # Fall back on the PHRED scores...
    try:
        qualities = record.letter_annotations["phred_quality"]
    except KeyError:
        raise ValueError(
            "No suitable quality scores found in "
            "letter_annotations of SeqRecord (id=%s)." % record.id
        ) from None
    # Try and use the precomputed mapping:
    try:
        return "".join(_phred_to_solexa_quality_str[qp] for qp in qualities)
    except KeyError:
        # Either no PHRED scores, or something odd like a float or None
        # or too big to be in the cache
        pass
    if None in qualities:
        raise TypeError("A quality value of None was found")
    # Must do this the slow way, first converting the PHRED scores into
    # Solexa scores:
    if max(qualities) >= 62.5:
        warnings.warn(
            "Data loss - max Solexa quality 62 in Solexa FASTQ", BiopythonWarning
        )
    return "".join(
        chr(min(126, int(round(solexa_quality_from_phred(qp))) + SOLEXA_SCORE_OFFSET))
        for qp in qualities
    )


# TODO - Default to nucleotide or even DNA?
def FastqGeneralIterator(source):
    """Iterate over Fastq records as string tuples (not as SeqRecord objects).

    Arguments:
     - source - input stream opened in text mode, or a path to a file

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

    **WARNING** - Because the "@" character can appear in the quality string,
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
    (given the read length of 36, these are probably Solexa Illumina reads where
    the quality has been mapped onto the PHRED values).

    This example has been edited to illustrate some of the nasty things allowed
    in the FASTQ format.  Firstly, on the "+" lines most but not all of the
    (redundant) identifiers are omitted.  In real files it is likely that all or
    none of these extra identifiers will be present.

    Secondly, while the first three sequences have been shown without line
    breaks, the last has been split over multiple lines.  In real files any line
    breaks are likely to be consistent.

    Thirdly, some of the quality string lines start with an "@" character.  For
    the second record this is unavoidable.  However for the fourth sequence this
    only happens because its quality string is split over two lines.  A naive
    parser could wrongly treat any line starting with an "@" as the beginning of
    a new sequence!  This code copes with this possible ambiguity by keeping
    track of the length of the sequence which gives the expected length of the
    quality string.

    Using this tricky example file as input, this short bit of code demonstrates
    what this parsing function would return:

    >>> with open("Quality/tricky.fastq") as handle:
    ...     for (title, sequence, quality) in FastqGeneralIterator(handle):
    ...         print(title)
    ...         print("%s %s" % (sequence, quality))
    ...
    071113_EAS56_0053:1:1:998:236
    TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
    071113_EAS56_0053:1:1:182:712
    ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
    071113_EAS56_0053:1:1:153:10
    TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
    071113_EAS56_0053:1:3:990:501
    TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA IIIIIII.IIIIII1@44@-7.%<&+/$/%4(++(%

    Finally we note that some sources state that the quality string should
    start with "!" (which using the PHRED mapping means the first letter always
    has a quality score of zero).  This rather restrictive rule is not widely
    observed, so is therefore ignored here.  One plus point about this "!" rule
    is that (provided there are no line breaks in the quality sequence) it
    would prevent the above problem with the "@" character.
    """
    try:
        handle = open(source)
    except TypeError:
        handle = source
        if handle.read(0) != "":
            raise StreamModeError("Fastq files must be opened in text mode") from None
    try:
        try:
            line = next(handle)
        except StopIteration:
            return  # Premature end of file, or just empty?

        while True:
            if line[0] != "@":
                raise ValueError(
                    "Records in Fastq files should start with '@' character"
                )
            title_line = line[1:].rstrip()
            seq_string = ""
            # There will now be one or more sequence lines; keep going until we
            # find the "+" marking the quality line:
            for line in handle:
                if line[0] == "+":
                    break
                seq_string += line.rstrip()
            else:
                if seq_string:
                    raise ValueError("End of file without quality information.")
                else:
                    raise ValueError("Unexpected end of file")
            # The title here is optional, but if present must match!
            second_title = line[1:].rstrip()
            if second_title and second_title != title_line:
                raise ValueError("Sequence and quality captions differ.")
            # This is going to slow things down a little, but assuming
            # this isn't allowed we should try and catch it here:
            if " " in seq_string or "\t" in seq_string:
                raise ValueError("Whitespace is not allowed in the sequence.")
            seq_len = len(seq_string)

            # There will now be at least one line of quality data, followed by
            # another sequence, or EOF
            line = None
            quality_string = ""
            for line in handle:
                if line[0] == "@":
                    # This COULD be the start of a new sequence. However, it MAY just
                    # be a line of quality data which starts with a "@" character.  We
                    # should be able to check this by looking at the sequence length
                    # and the amount of quality data found so far.
                    if len(quality_string) >= seq_len:
                        # We expect it to be equal if this is the start of a new record.
                        # If the quality data is longer, we'll raise an error below.
                        break
                    # Continue - its just some (more) quality data.
                quality_string += line.rstrip()
            else:
                if line is None:
                    raise ValueError("Unexpected end of file")
                line = None

            if seq_len != len(quality_string):
                raise ValueError(
                    "Lengths of sequence and quality values differs for %s (%i and %i)."
                    % (title_line, seq_len, len(quality_string))
                )

            # Return the record and then continue...
            yield (title_line, seq_string, quality_string)

            if line is None:
                break
    finally:
        if handle is not source:
            handle.close()


class FastqPhredIterator(SequenceIterator):
    """Parser for FASTQ files."""

    def __init__(self, source, alphabet=None, title2ids=None):
        """Iterate over FASTQ records as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
         - alphabet - optional alphabet, no longer used. Leave as None.
         - title2ids (DEPRECATED) - A function that, when given the title line
           from the FASTQ file (without the beginning >), will return the id,
           name and description (in that order) for the record as a tuple of
           strings.  If this is not given, then the entire title line will be
           used as the description, and the first word as the id and name.

        The use of title2ids matches that of Bio.SeqIO.FastaIO.

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
        string encoding the PHRED qualities using a ASCII values with an offset of
        33 (e.g. ";;3;;;;;;;;;;;;7;;;;;;;88").

        Using this module directly you might run:

        >>> with open("Quality/example.fastq") as handle:
        ...     for record in FastqPhredIterator(handle):
        ...         print("%s %s" % (record.id, record.seq))
        EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
        EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
        EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG

        Typically however, you would call this via Bio.SeqIO instead with "fastq"
        (or "fastq-sanger") as the format:

        >>> from Bio import SeqIO
        >>> with open("Quality/example.fastq") as handle:
        ...     for record in SeqIO.parse(handle, "fastq"):
        ...         print("%s %s" % (record.id, record.seq))
        EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
        EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
        EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG

        If you want to look at the qualities, they are record in each record's
        per-letter-annotation dictionary as a simple list of integers:

        >>> print(record.letter_annotations["phred_quality"])
        [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

        The title2ids argument is deprecated. Instead, please use a generator
        function to modify the records returned by the parser. For example, to
        store the mean PHRED quality in the record description, use

        >>> from statistics import mean
        >>> def modify_records(records):
        ...     for record in records:
        ...         record.description = mean(record.letter_annotations['phred_quality'])
        ...         yield record
        ...
        >>> with open('Quality/example.fastq') as handle:
        ...     for record in modify_records(FastqPhredIterator(handle)):
        ...         print(record.id, record.description)
        ...
        EAS54_6_R1_2_1_413_324 25.28
        EAS54_6_R1_2_1_540_792 24.52
        EAS54_6_R1_2_1_443_348 23.4

        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        if title2ids is not None:
            warnings.warn(
                "The title2ids argument is deprecated. Instead, please use a "
                "generator function to modify records returned by the parser. "
                "For example, to change the record description to a counter, "
                "use\n"
                "\n"
                ">>> from statistics import mean\n"
                ">>> def modify_records(records):\n"
                "...     for record in records:\n"
                "...         record.description = mean(record.letter_annotations['phred_quality'])\n"
                "...         yield record\n"
                "...\n"
                ">>> with open('Quality/example.fastq') as handle:\n"
                "...     for record in modify_records(FastqPhredIterator(handle)):\n"
                "...         print(record.id, record.description)\n"
                "\n",
                BiopythonDeprecationWarning,
            )
        self.title2ids = title2ids
        super().__init__(source, mode="t", fmt="Fastq")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = self.iterate(handle)
        return records

    def iterate(self, handle):
        """Parse the file and generate SeqRecord objects."""
        title2ids = self.title2ids
        assert SANGER_SCORE_OFFSET == ord("!")
        # Originally, I used a list expression for each record:
        #
        # qualities = [ord(letter)-SANGER_SCORE_OFFSET for letter in quality_string]
        #
        # Precomputing is faster, perhaps partly by avoiding the subtractions.
        q_mapping = {
            chr(letter): letter - SANGER_SCORE_OFFSET
            for letter in range(SANGER_SCORE_OFFSET, 94 + SANGER_SCORE_OFFSET)
        }

        for title_line, seq_string, quality_string in FastqGeneralIterator(handle):
            if title2ids:
                id, name, descr = title2ids(title_line)
            else:
                descr = title_line
                id = descr.split()[0]
                name = id
            record = SeqRecord(Seq(seq_string), id=id, name=name, description=descr)
            try:
                qualities = [q_mapping[letter] for letter in quality_string]
            except KeyError:
                raise ValueError("Invalid character in quality string") from None
            # For speed, will now use a dirty trick to speed up assigning the
            # qualities. We do this to bypass the length check imposed by the
            # per-letter-annotations restricted dict (as this has already been
            # checked by FastqGeneralIterator). This is equivalent to:
            # record.letter_annotations["phred_quality"] = qualities
            dict.__setitem__(record._per_letter_annotations, "phred_quality", qualities)
            yield record


def FastqSolexaIterator(source, alphabet=None, title2ids=None):
    r"""Parse old Solexa/Illumina FASTQ like files (which differ in the quality mapping).

    The optional arguments are the same as those for the FastqPhredIterator.

    For each sequence in Solexa/Illumina FASTQ files there is a matching string
    encoding the Solexa integer qualities using ASCII values with an offset
    of 64.  Solexa scores are scaled differently to PHRED scores, and Biopython
    will NOT perform any automatic conversion when loading.

    NOTE - This file format is used by the OLD versions of the Solexa/Illumina
    pipeline. See also the FastqIlluminaIterator function for the NEW version.

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

    >>> with open("Quality/solexa_example.fastq") as handle:
    ...     for record in FastqSolexaIterator(handle):
    ...         print("%s %s" % (record.id, record.seq))
    SLXA-B3_649_FC8437_R1_1_1_610_79 GATGTGCAATACCTTTGTAGAGGAA
    SLXA-B3_649_FC8437_R1_1_1_397_389 GGTTTGAGAAAGAGAAATGAGATAA
    SLXA-B3_649_FC8437_R1_1_1_850_123 GAGGGTGTTGATCATGATGATGGCG
    SLXA-B3_649_FC8437_R1_1_1_362_549 GGAAACAAAGTTTTTCTCAACATAG
    SLXA-B3_649_FC8437_R1_1_1_183_714 GTATTATTTAATGGCATACACTCAA

    Typically however, you would call this via Bio.SeqIO instead with
    "fastq-solexa" as the format:

    >>> from Bio import SeqIO
    >>> with open("Quality/solexa_example.fastq") as handle:
    ...     for record in SeqIO.parse(handle, "fastq-solexa"):
    ...         print("%s %s" % (record.id, record.seq))
    SLXA-B3_649_FC8437_R1_1_1_610_79 GATGTGCAATACCTTTGTAGAGGAA
    SLXA-B3_649_FC8437_R1_1_1_397_389 GGTTTGAGAAAGAGAAATGAGATAA
    SLXA-B3_649_FC8437_R1_1_1_850_123 GAGGGTGTTGATCATGATGATGGCG
    SLXA-B3_649_FC8437_R1_1_1_362_549 GGAAACAAAGTTTTTCTCAACATAG
    SLXA-B3_649_FC8437_R1_1_1_183_714 GTATTATTTAATGGCATACACTCAA

    If you want to look at the qualities, they are recorded in each record's
    per-letter-annotation dictionary as a simple list of integers:

    >>> print(record.letter_annotations["solexa_quality"])
    [25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 23, 25, 25, 25, 25, 23, 25, 23, 23, 21, 23, 23, 23, 17, 17]

    These scores aren't very good, but they are high enough that they map
    almost exactly onto PHRED scores:

    >>> print("%0.2f" % phred_quality_from_solexa(25))
    25.01

    Let's look at faked example read which is even worse, where there are
    more noticeable differences between the Solexa and PHRED scores::

         @slxa_0001_1_0001_01
         ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
         +slxa_0001_1_0001_01
         hgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@?>=<;

    Again, you would typically use Bio.SeqIO to read this file in (rather than
    calling the Bio.SeqIO.QualtityIO module directly).  Most FASTQ files will
    contain thousands of reads, so you would normally use Bio.SeqIO.parse()
    as shown above.  This example has only as one entry, so instead we can
    use the Bio.SeqIO.read() function:

    >>> from Bio import SeqIO
    >>> with open("Quality/solexa_faked.fastq") as handle:
    ...     record = SeqIO.read(handle, "fastq-solexa")
    >>> print("%s %s" % (record.id, record.seq))
    slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
    >>> print(record.letter_annotations["solexa_quality"])
    [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]

    These quality scores are so low that when converted from the Solexa scheme
    into PHRED scores they look quite different:

    >>> print("%0.2f" % phred_quality_from_solexa(-1))
    2.54
    >>> print("%0.2f" % phred_quality_from_solexa(-5))
    1.19

    Note you can use the Bio.SeqIO.write() function or the SeqRecord's format
    method to output the record(s):

    >>> print(record.format("fastq-solexa"))
    @slxa_0001_1_0001_01
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
    +
    hgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@?>=<;
    <BLANKLINE>

    Note this output is slightly different from the input file as Biopython
    has left out the optional repetition of the sequence identifier on the "+"
    line.  If you want the to use PHRED scores, use "fastq" or "qual" as the
    output format instead, and Biopython will do the conversion for you:

    >>> print(record.format("fastq"))
    @slxa_0001_1_0001_01
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
    +
    IHGFEDCBA@?>=<;:9876543210/.-,++*)('&&%%$$##""
    <BLANKLINE>

    >>> print(record.format("qual"))
    >slxa_0001_1_0001_01
    40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21
    20 19 18 17 16 15 14 13 12 11 10 10 9 8 7 6 5 5 4 4 3 3 2 2
    1 1
    <BLANKLINE>

    As shown above, the poor quality Solexa reads have been mapped to the
    equivalent PHRED score (e.g. -5 to 1 as shown earlier).
    """
    if alphabet is not None:
        raise ValueError("The alphabet argument is no longer supported")

    q_mapping = {
        chr(letter): letter - SOLEXA_SCORE_OFFSET
        for letter in range(SOLEXA_SCORE_OFFSET - 5, 63 + SOLEXA_SCORE_OFFSET)
    }

    for title_line, seq_string, quality_string in FastqGeneralIterator(source):
        if title2ids:
            id, name, descr = title2ids(title_line)
        else:
            descr = title_line
            id = descr.split()[0]
            name = id
        record = SeqRecord(Seq(seq_string), id=id, name=name, description=descr)
        try:
            qualities = [q_mapping[letter] for letter in quality_string]
        # DO NOT convert these into PHRED qualities automatically!
        except KeyError:
            raise ValueError("Invalid character in quality string") from None
        # Dirty trick to speed up this line:
        # record.letter_annotations["solexa_quality"] = qualities
        dict.__setitem__(record._per_letter_annotations, "solexa_quality", qualities)
        yield record


def FastqIlluminaIterator(source, alphabet=None, title2ids=None):
    """Parse Illumina 1.3 to 1.7 FASTQ like files (which differ in the quality mapping).

    The optional arguments are the same as those for the FastqPhredIterator.

    For each sequence in Illumina 1.3+ FASTQ files there is a matching string
    encoding PHRED integer qualities using ASCII values with an offset of 64.

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("Quality/illumina_faked.fastq", "fastq-illumina")
    >>> print("%s %s" % (record.id, record.seq))
    Test ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN
    >>> max(record.letter_annotations["phred_quality"])
    40
    >>> min(record.letter_annotations["phred_quality"])
    0

    NOTE - Older versions of the Solexa/Illumina pipeline encoded Solexa scores
    with an ASCII offset of 64. They are approximately equal but only for high
    quality reads. If you have an old Solexa/Illumina file with negative
    Solexa scores, and try and read this as an Illumina 1.3+ file it will fail:

    >>> record2 = SeqIO.read("Quality/solexa_faked.fastq", "fastq-illumina")
    Traceback (most recent call last):
       ...
    ValueError: Invalid character in quality string

    NOTE - True Sanger style FASTQ files use PHRED scores with an offset of 33.
    """
    if alphabet is not None:
        raise ValueError("The alphabet argument is no longer supported")

    q_mapping = {
        chr(letter): letter - SOLEXA_SCORE_OFFSET
        for letter in range(SOLEXA_SCORE_OFFSET, 63 + SOLEXA_SCORE_OFFSET)
    }

    for title_line, seq_string, quality_string in FastqGeneralIterator(source):
        if title2ids:
            id, name, descr = title2ids(title_line)
        else:
            descr = title_line
            id = descr.split()[0]
            name = id
        record = SeqRecord(Seq(seq_string), id=id, name=name, description=descr)
        try:
            qualities = [q_mapping[letter] for letter in quality_string]
        except KeyError:
            raise ValueError("Invalid character in quality string") from None
        # Dirty trick to speed up this line:
        # record.letter_annotations["phred_quality"] = qualities
        dict.__setitem__(record._per_letter_annotations, "phred_quality", qualities)
        yield record


class QualPhredIterator(SequenceIterator):
    """Parser for QUAL files with PHRED quality scores but no sequence."""

    def __init__(self, source, alphabet=None, title2ids=None):
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

        >>> with open("Quality/example.qual") as handle:
        ...     for record in QualPhredIterator(handle):
        ...         print("%s read of length %d" % (record.id, len(record.seq)))
        EAS54_6_R1_2_1_413_324 read of length 25
        EAS54_6_R1_2_1_540_792 read of length 25
        EAS54_6_R1_2_1_443_348 read of length 25

        Typically however, you would call this via Bio.SeqIO instead with "qual"
        as the format:

        >>> from Bio import SeqIO
        >>> with open("Quality/example.qual") as handle:
        ...     for record in SeqIO.parse(handle, "qual"):
        ...         print("%s read of length %d" % (record.id, len(record.seq)))
        EAS54_6_R1_2_1_413_324 read of length 25
        EAS54_6_R1_2_1_540_792 read of length 25
        EAS54_6_R1_2_1_443_348 read of length 25

        Only the sequence length is known, as the QUAL file does not contain
        the sequence string itself.

        The quality scores themselves are available as a list of integers
        in each record's per-letter-annotation:

        >>> print(record.letter_annotations["phred_quality"])
        [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

        You can still slice one of these SeqRecord objects:

        >>> sub_record = record[5:10]
        >>> print("%s %s" % (sub_record.id, sub_record.letter_annotations["phred_quality"]))
        EAS54_6_R1_2_1_443_348 [26, 26, 26, 26, 26]

        As of Biopython 1.59, this parser will accept files with negatives quality
        scores but will replace them with the lowest possible PHRED score of zero.
        This will trigger a warning, previously it raised a ValueError exception.
        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        self.title2ids = title2ids
        super().__init__(source, mode="t", fmt="QUAL")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = self.iterate(handle)
        return records

    def iterate(self, handle):
        """Parse the file and generate SeqRecord objects."""
        title2ids = self.title2ids
        # Skip any text before the first record (e.g. blank lines, comments)
        for line in handle:
            if line[0] == ">":
                break
        else:
            return

        while True:
            if line[0] != ">":
                raise ValueError(
                    "Records in Fasta files should start with '>' character"
                )
            if title2ids:
                id, name, descr = title2ids(line[1:].rstrip())
            else:
                descr = line[1:].rstrip()
                id = descr.split()[0]
                name = id

            qualities = []
            for line in handle:
                if line[0] == ">":
                    break
                qualities.extend(int(word) for word in line.split())
            else:
                line = None

            if qualities and min(qualities) < 0:
                warnings.warn(
                    "Negative quality score %i found, substituting PHRED zero instead."
                    % min(qualities),
                    BiopythonParserWarning,
                )
                qualities = [max(0, q) for q in qualities]

            # Return the record and then continue...
            sequence = Seq(None, length=len(qualities))
            record = SeqRecord(sequence, id=id, name=name, description=descr)
            # Dirty trick to speed up this line:
            # record.letter_annotations["phred_quality"] = qualities
            dict.__setitem__(record._per_letter_annotations, "phred_quality", qualities)
            yield record

            if line is None:
                return  # StopIteration
        raise ValueError("Unrecognised QUAL record format.")


class FastqPhredWriter(SequenceWriter):
    """Class to write standard FASTQ format files (using PHRED quality scores) (OBSOLETE).

    Although you can use this class directly, you are strongly encouraged
    to use the ``as_fastq`` function, or top level ``Bio.SeqIO.write()``
    function instead via the format name "fastq" or the alias "fastq-sanger".

    For example, this code reads in a standard Sanger style FASTQ file
    (using PHRED scores) and re-saves it as another Sanger style FASTQ file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("Quality/example.fastq", "fastq")
    >>> with open("Quality/temp.fastq", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "fastq")
    3

    You might want to do this if the original file included extra line breaks,
    which while valid may not be supported by all tools.  The output file from
    Biopython will have each sequence on a single line, and each quality
    string on a single line (which is considered desirable for maximum
    compatibility).

    In this next example, an old style Solexa/Illumina FASTQ file (using Solexa
    quality scores) is converted into a standard Sanger style FASTQ file using
    PHRED qualities:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("Quality/solexa_example.fastq", "fastq-solexa")
    >>> with open("Quality/temp.fastq", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "fastq")
    5

    This code is also called if you use the .format("fastq") method of a
    SeqRecord, or .format("fastq-sanger") if you prefer that alias.

    Note that Sanger FASTQ files have an upper limit of PHRED quality 93, which is
    encoded as ASCII 126, the tilde. If your quality scores are truncated to fit, a
    warning is issued.

    P.S. To avoid cluttering up your working directory, you can delete this
    temporary file now:

    >>> import os
    >>> os.remove("Quality/temp.fastq")
    """

    assert SANGER_SCORE_OFFSET == ord("!")

    def write_record(self, record):
        """Write a single FASTQ record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True
        # TODO - Is an empty sequence allowed in FASTQ format?
        seq = record.seq
        if seq is None:
            raise ValueError(f"No sequence for record {record.id}")
        qualities_str = _get_sanger_quality_str(record)
        if len(qualities_str) != len(seq):
            raise ValueError(
                "Record %s has sequence length %i but %i quality scores"
                % (record.id, len(seq), len(qualities_str))
            )

        # FASTQ files can include a description, just like FASTA files
        # (at least, this is what the NCBI Short Read Archive does)
        id = self.clean(record.id)
        description = self.clean(record.description)
        if description and description.split(None, 1)[0] == id:
            # The description includes the id at the start
            title = description
        elif description:
            title = f"{id} {description}"
        else:
            title = id

        self.handle.write(f"@{title}\n{seq}\n+\n{qualities_str}\n")


def as_fastq(record):
    """Turn a SeqRecord into a Sanger FASTQ formatted string.

    This is used internally by the SeqRecord's .format("fastq")
    method and by the SeqIO.write(..., ..., "fastq") function,
    and under the format alias "fastq-sanger" as well.
    """
    seq_str = _get_seq_string(record)
    qualities_str = _get_sanger_quality_str(record)
    if len(qualities_str) != len(seq_str):
        raise ValueError(
            "Record %s has sequence length %i but %i quality scores"
            % (record.id, len(seq_str), len(qualities_str))
        )
    id = _clean(record.id)
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id:
        title = description
    elif description:
        title = f"{id} {description}"
    else:
        title = id
    return f"@{title}\n{seq_str}\n+\n{qualities_str}\n"


class QualPhredWriter(SequenceWriter):
    """Class to write QUAL format files (using PHRED quality scores) (OBSOLETE).

    Although you can use this class directly, you are strongly encouraged
    to use the ``as_qual`` function, or top level ``Bio.SeqIO.write()``
    function instead.

    For example, this code reads in a FASTQ file and saves the quality scores
    into a QUAL file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("Quality/example.fastq", "fastq")
    >>> with open("Quality/temp.qual", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "qual")
    3

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
           Defaults to wrapping the sequence at 60 characters. Use
           zero (or None) for no wrapping, giving a single long line
           for the sequence.
         - record2title - Optional function to return the text to be
           used for the title line of each record.  By default a
           combination of the record.id and record.description is
           used.  If the record.description starts with the record.id,
           then just the record.description is used.

        The record2title argument is present for consistency with the
        Bio.SeqIO.FastaIO writer class.
        """
        super().__init__(handle)
        # self.handle = handle
        self.wrap = None
        if wrap:
            if wrap < 1:
                raise ValueError
        self.wrap = wrap
        self.record2title = record2title

    def write_record(self, record):
        """Write a single QUAL record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        handle = self.handle
        wrap = self.wrap

        if self.record2title:
            title = self.clean(self.record2title(record))
        else:
            id = self.clean(record.id)
            description = self.clean(record.description)
            if description and description.split(None, 1)[0] == id:
                # The description includes the id at the start
                title = description
            elif description:
                title = f"{id} {description}"
            else:
                title = id
        handle.write(f">{title}\n")

        qualities = _get_phred_quality(record)
        try:
            # This rounds to the nearest integer.
            # TODO - can we record a float in a qual file?
            qualities_strs = [("%i" % round(q, 0)) for q in qualities]
        except TypeError:
            if None in qualities:
                raise TypeError("A quality value of None was found") from None
            else:
                raise

        if wrap > 5:
            # Fast wrapping
            data = " ".join(qualities_strs)
            while True:
                if len(data) <= wrap:
                    self.handle.write(data + "\n")
                    break
                else:
                    # By construction there must be spaces in the first X chars
                    # (unless we have X digit or higher quality scores!)
                    i = data.rfind(" ", 0, wrap)
                    handle.write(data[:i] + "\n")
                    data = data[i + 1 :]
        elif wrap:
            # Safe wrapping
            while qualities_strs:
                line = qualities_strs.pop(0)
                while qualities_strs and len(line) + 1 + len(qualities_strs[0]) < wrap:
                    line += " " + qualities_strs.pop(0)
                handle.write(line + "\n")
        else:
            # No wrapping
            data = " ".join(qualities_strs)
            handle.write(data + "\n")


def as_qual(record):
    """Turn a SeqRecord into a QUAL formatted string.

    This is used internally by the SeqRecord's .format("qual")
    method and by the SeqIO.write(..., ..., "qual") function.
    """
    id = _clean(record.id)
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id:
        title = description
    elif description:
        title = f"{id} {description}"
    else:
        title = id
    lines = [f">{title}\n"]

    qualities = _get_phred_quality(record)
    try:
        # This rounds to the nearest integer.
        # TODO - can we record a float in a qual file?
        qualities_strs = [("%i" % round(q, 0)) for q in qualities]
    except TypeError:
        if None in qualities:
            raise TypeError("A quality value of None was found") from None
        else:
            raise

    # Safe wrapping
    while qualities_strs:
        line = qualities_strs.pop(0)
        while qualities_strs and len(line) + 1 + len(qualities_strs[0]) < 60:
            line += " " + qualities_strs.pop(0)
        lines.append(line + "\n")
    return "".join(lines)


class FastqSolexaWriter(SequenceWriter):
    r"""Write old style Solexa/Illumina FASTQ format files (with Solexa qualities) (OBSOLETE).

    This outputs FASTQ files like those from the early Solexa/Illumina
    pipeline, using Solexa scores and an ASCII offset of 64. These are
    NOT compatible with the standard Sanger style PHRED FASTQ files.

    If your records contain a "solexa_quality" entry under letter_annotations,
    this is used, otherwise any "phred_quality" entry will be used after
    conversion using the solexa_quality_from_phred function. If neither style
    of quality scores are present, an exception is raised.

    Although you can use this class directly, you are strongly encouraged
    to use the ``as_fastq_solexa`` function, or top-level ``Bio.SeqIO.write()``
    function instead.  For example, this code reads in a FASTQ file and re-saves
    it as another FASTQ file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("Quality/solexa_example.fastq", "fastq-solexa")
    >>> with open("Quality/temp.fastq", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "fastq-solexa")
    5

    You might want to do this if the original file included extra line breaks,
    which (while valid) may not be supported by all tools.  The output file
    from Biopython will have each sequence on a single line, and each quality
    string on a single line (which is considered desirable for maximum
    compatibility).

    This code is also called if you use the .format("fastq-solexa") method of
    a SeqRecord. For example,

    >>> record = SeqIO.read("Quality/sanger_faked.fastq", "fastq-sanger")
    >>> print(record.format("fastq-solexa"))
    @Test PHRED qualities from 40 to 0 inclusive
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN
    +
    hgfedcba`_^]\[ZYXWVUTSRQPONMLKJHGFECB@>;;
    <BLANKLINE>

    Note that Solexa FASTQ files have an upper limit of Solexa quality 62, which is
    encoded as ASCII 126, the tilde.  If your quality scores must be truncated to fit,
    a warning is issued.

    P.S. Don't forget to delete the temp file if you don't need it anymore:

    >>> import os
    >>> os.remove("Quality/temp.fastq")
    """

    def write_record(self, record):
        """Write a single FASTQ record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        # TODO - Is an empty sequence allowed in FASTQ format?
        seq = record.seq
        if seq is None:
            raise ValueError(f"No sequence for record {record.id}")
        qualities_str = _get_solexa_quality_str(record)
        if len(qualities_str) != len(seq):
            raise ValueError(
                "Record %s has sequence length %i but %i quality scores"
                % (record.id, len(seq), len(qualities_str))
            )

        # FASTQ files can include a description, just like FASTA files
        # (at least, this is what the NCBI Short Read Archive does)
        id = self.clean(record.id)
        description = self.clean(record.description)
        if description and description.split(None, 1)[0] == id:
            # The description includes the id at the start
            title = description
        elif description:
            title = f"{id} {description}"
        else:
            title = id

        self.handle.write(f"@{title}\n{seq}\n+\n{qualities_str}\n")


def as_fastq_solexa(record):
    """Turn a SeqRecord into a Solexa FASTQ formatted string.

    This is used internally by the SeqRecord's .format("fastq-solexa")
    method and by the SeqIO.write(..., ..., "fastq-solexa") function.
    """
    seq_str = _get_seq_string(record)
    qualities_str = _get_solexa_quality_str(record)
    if len(qualities_str) != len(seq_str):
        raise ValueError(
            "Record %s has sequence length %i but %i quality scores"
            % (record.id, len(seq_str), len(qualities_str))
        )
    id = _clean(record.id)
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id:
        # The description includes the id at the start
        title = description
    elif description:
        title = f"{id} {description}"
    else:
        title = id
    return f"@{title}\n{seq_str}\n+\n{qualities_str}\n"


class FastqIlluminaWriter(SequenceWriter):
    r"""Write Illumina 1.3+ FASTQ format files (with PHRED quality scores) (OBSOLETE).

    This outputs FASTQ files like those from the Solexa/Illumina 1.3+ pipeline,
    using PHRED scores and an ASCII offset of 64. Note these files are NOT
    compatible with the standard Sanger style PHRED FASTQ files which use an
    ASCII offset of 32.

    Although you can use this class directly, you are strongly encouraged to
    use the ``as_fastq_illumina`` or top-level ``Bio.SeqIO.write()`` function
    with format name "fastq-illumina" instead. This code is also called if you
    use the .format("fastq-illumina") method of a SeqRecord. For example,

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("Quality/sanger_faked.fastq", "fastq-sanger")
    >>> print(record.format("fastq-illumina"))
    @Test PHRED qualities from 40 to 0 inclusive
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN
    +
    hgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@
    <BLANKLINE>

    Note that Illumina FASTQ files have an upper limit of PHRED quality 62, which is
    encoded as ASCII 126, the tilde. If your quality scores are truncated to fit, a
    warning is issued.
    """

    def write_record(self, record):
        """Write a single FASTQ record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        # TODO - Is an empty sequence allowed in FASTQ format?
        seq = record.seq
        if seq is None:
            raise ValueError(f"No sequence for record {record.id}")
        qualities_str = _get_illumina_quality_str(record)
        if len(qualities_str) != len(seq):
            raise ValueError(
                "Record %s has sequence length %i but %i quality scores"
                % (record.id, len(seq), len(qualities_str))
            )

        # FASTQ files can include a description, just like FASTA files
        # (at least, this is what the NCBI Short Read Archive does)
        id = self.clean(record.id)
        description = self.clean(record.description)
        if description and description.split(None, 1)[0] == id:
            # The description includes the id at the start
            title = description
        elif description:
            title = f"{id} {description}"
        else:
            title = id

        self.handle.write(f"@{title}\n{seq}\n+\n{qualities_str}\n")


def as_fastq_illumina(record):
    """Turn a SeqRecord into an Illumina FASTQ formatted string.

    This is used internally by the SeqRecord's .format("fastq-illumina")
    method and by the SeqIO.write(..., ..., "fastq-illumina") function.
    """
    seq_str = _get_seq_string(record)
    qualities_str = _get_illumina_quality_str(record)
    if len(qualities_str) != len(seq_str):
        raise ValueError(
            "Record %s has sequence length %i but %i quality scores"
            % (record.id, len(seq_str), len(qualities_str))
        )
    id = _clean(record.id)
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id:
        title = description
    elif description:
        title = f"{id} {description}"
    else:
        title = id
    return f"@{title}\n{seq_str}\n+\n{qualities_str}\n"


def PairedFastaQualIterator(fasta_source, qual_source, alphabet=None, title2ids=None):
    """Iterate over matched FASTA and QUAL files as SeqRecord objects.

    For example, consider this short QUAL file with PHRED quality scores::

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

    >>> with open("Quality/example.fasta") as f:
    ...     with open("Quality/example.qual") as q:
    ...         for record in PairedFastaQualIterator(f, q):
    ...             print("%s %s" % (record.id, record.seq))
    ...
    EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
    EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
    EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG

    As with the FASTQ or QUAL parsers, if you want to look at the qualities,
    they are in each record's per-letter-annotation dictionary as a simple
    list of integers:

    >>> print(record.letter_annotations["phred_quality"])
    [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

    If you have access to data as a FASTQ format file, using that directly
    would be simpler and more straight forward.  Note that you can easily use
    this function to convert paired FASTA and QUAL files into FASTQ files:

    >>> from Bio import SeqIO
    >>> with open("Quality/example.fasta") as f:
    ...     with open("Quality/example.qual") as q:
    ...         SeqIO.write(PairedFastaQualIterator(f, q), "Quality/temp.fastq", "fastq")
    ...
    3

    And don't forget to clean up the temp file if you don't need it anymore:

    >>> import os
    >>> os.remove("Quality/temp.fastq")
    """
    if alphabet is not None:
        raise ValueError("The alphabet argument is no longer supported")

    from Bio.SeqIO.FastaIO import FastaIterator

    fasta_iter = FastaIterator(fasta_source, title2ids=title2ids)
    qual_iter = QualPhredIterator(qual_source, title2ids=title2ids)

    # Using zip wouldn't load everything into memory, but also would not catch
    # any extra records found in only one file.
    while True:
        try:
            f_rec = next(fasta_iter)
        except StopIteration:
            f_rec = None
        try:
            q_rec = next(qual_iter)
        except StopIteration:
            q_rec = None
        if f_rec is None and q_rec is None:
            # End of both files
            break
        if f_rec is None:
            raise ValueError("FASTA file has more entries than the QUAL file.")
        if q_rec is None:
            raise ValueError("QUAL file has more entries than the FASTA file.")
        if f_rec.id != q_rec.id:
            raise ValueError(
                f"FASTA and QUAL entries do not match ({f_rec.id} vs {q_rec.id})."
            )
        if len(f_rec) != len(q_rec.letter_annotations["phred_quality"]):
            raise ValueError(
                f"Sequence length and number of quality scores disagree for {f_rec.id}"
            )
        # Merge the data....
        f_rec.letter_annotations["phred_quality"] = q_rec.letter_annotations[
            "phred_quality"
        ]
        yield f_rec
    # Done


def _fastq_generic(in_file, out_file, mapping):
    """FASTQ helper function where can't have data loss by truncation (PRIVATE)."""
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    null = chr(0)
    with as_handle(out_file, "w") as out_handle:
        for title, seq, old_qual in FastqGeneralIterator(in_file):
            count += 1
            # map the qual...
            qual = old_qual.translate(mapping)
            if null in qual:
                raise ValueError("Invalid character in quality string")
            out_handle.write(f"@{title}\n{seq}\n+\n{qual}\n")
    return count


def _fastq_generic2(in_file, out_file, mapping, truncate_char, truncate_msg):
    """FASTQ helper function where there could be data loss by truncation (PRIVATE)."""
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    null = chr(0)
    with as_handle(out_file, "w") as out_handle:
        for title, seq, old_qual in FastqGeneralIterator(in_file):
            count += 1
            # map the qual...
            qual = old_qual.translate(mapping)
            if null in qual:
                raise ValueError("Invalid character in quality string")
            if truncate_char in qual:
                qual = qual.replace(truncate_char, chr(126))
                warnings.warn(truncate_msg, BiopythonWarning)
            out_handle.write(f"@{title}\n{seq}\n+\n{qual}\n")
    return count


def _fastq_sanger_convert_fastq_sanger(in_file, out_file):
    """Fast Sanger FASTQ to Sanger FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(0, 33)]
        + [chr(ascii) for ascii in range(33, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_solexa_convert_fastq_solexa(in_file, out_file):
    """Fast Solexa FASTQ to Solexa FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.
    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(0, 59)]
        + [chr(ascii) for ascii in range(59, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_illumina_convert_fastq_illumina(in_file, out_file):
    """Fast Illumina 1.3+ FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.
    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(0, 64)]
        + [chr(ascii) for ascii in range(64, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_illumina_convert_fastq_sanger(in_file, out_file):
    """Fast Illumina 1.3+ FASTQ to Sanger FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(0, 64)]
        + [chr(33 + q) for q in range(0, 62 + 1)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_sanger_convert_fastq_illumina(in_file, out_file):
    """Fast Sanger FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion. Will issue a warning if the scores had to be truncated at 62
    (maximum possible in the Illumina 1.3+ FASTQ format)
    """
    # Map unexpected chars to null
    trunc_char = chr(1)
    mapping = "".join(
        [chr(0) for ascii in range(0, 33)]
        + [chr(64 + q) for q in range(0, 62 + 1)]
        + [trunc_char for ascii in range(96, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic2(
        in_file,
        out_file,
        mapping,
        trunc_char,
        "Data loss - max PHRED quality 62 in Illumina 1.3+ FASTQ",
    )


def _fastq_solexa_convert_fastq_sanger(in_file, out_file):
    """Fast Solexa FASTQ to Sanger FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(0, 59)]
        + [
            chr(33 + int(round(phred_quality_from_solexa(q))))
            for q in range(-5, 62 + 1)
        ]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_sanger_convert_fastq_solexa(in_file, out_file):
    """Fast Sanger FASTQ to Solexa FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion. Will issue a warning if the scores had to be truncated at 62
    (maximum possible in the Solexa FASTQ format)
    """
    # Map unexpected chars to null
    trunc_char = chr(1)
    mapping = "".join(
        [chr(0) for ascii in range(0, 33)]
        + [chr(64 + int(round(solexa_quality_from_phred(q)))) for q in range(0, 62 + 1)]
        + [trunc_char for ascii in range(96, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic2(
        in_file,
        out_file,
        mapping,
        trunc_char,
        "Data loss - max Solexa quality 62 in Solexa FASTQ",
    )


def _fastq_solexa_convert_fastq_illumina(in_file, out_file):
    """Fast Solexa FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(0, 59)]
        + [
            chr(64 + int(round(phred_quality_from_solexa(q))))
            for q in range(-5, 62 + 1)
        ]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_illumina_convert_fastq_solexa(in_file, out_file):
    """Fast Illumina 1.3+ FASTQ to Solexa FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(0, 64)]
        + [chr(64 + int(round(solexa_quality_from_phred(q)))) for q in range(0, 62 + 1)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_convert_fasta(in_file, out_file):
    """Fast FASTQ to FASTA conversion (PRIVATE).

    Avoids dealing with the FASTQ quality encoding, and creating SeqRecord and
    Seq objects in order to speed up this conversion.

    NOTE - This does NOT check the characters used in the FASTQ quality string
    are valid!
    """
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    with as_handle(out_file, "w") as out_handle:
        for title, seq, qual in FastqGeneralIterator(in_file):
            count += 1
            out_handle.write(f">{title}\n")
            # Do line wrapping
            for i in range(0, len(seq), 60):
                out_handle.write(seq[i : i + 60] + "\n")
    return count


def _fastq_convert_tab(in_file, out_file):
    """Fast FASTQ to simple tabbed conversion (PRIVATE).

    Avoids dealing with the FASTQ quality encoding, and creating SeqRecord and
    Seq objects in order to speed up this conversion.

    NOTE - This does NOT check the characters used in the FASTQ quality string
    are valid!
    """
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    with as_handle(out_file, "w") as out_handle:
        for title, seq, qual in FastqGeneralIterator(in_file):
            count += 1
            out_handle.write(f"{title.split(None, 1)[0]}\t{seq}\n")
    return count


def _fastq_convert_qual(in_file, out_file, mapping):
    """FASTQ helper function for QUAL output (PRIVATE).

    Mapping should be a dictionary mapping expected ASCII characters from the
    FASTQ quality string to PHRED quality scores (as strings).
    """
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    with as_handle(out_file, "w") as out_handle:
        for title, seq, qual in FastqGeneralIterator(in_file):
            count += 1
            out_handle.write(f">{title}\n")
            # map the qual... note even with Sanger encoding max 2 digits
            try:
                qualities_strs = [mapping[ascii] for ascii in qual]
            except KeyError:
                raise ValueError("Invalid character in quality string") from None
            data = " ".join(qualities_strs)
            while len(data) > 60:
                # Know quality scores are either 1 or 2 digits, so there
                # must be a space in any three consecutive characters.
                if data[60] == " ":
                    out_handle.write(data[:60] + "\n")
                    data = data[61:]
                elif data[59] == " ":
                    out_handle.write(data[:59] + "\n")
                    data = data[60:]
                else:
                    assert data[58] == " ", "Internal logic failure in wrapping"
                    out_handle.write(data[:58] + "\n")
                    data = data[59:]
            out_handle.write(data + "\n")
    return count


def _fastq_sanger_convert_qual(in_file, out_file):
    """Fast Sanger FASTQ to QUAL conversion (PRIVATE)."""
    mapping = {chr(q + 33): str(q) for q in range(0, 93 + 1)}
    return _fastq_convert_qual(in_file, out_file, mapping)


def _fastq_solexa_convert_qual(in_file, out_file):
    """Fast Solexa FASTQ to QUAL conversion (PRIVATE)."""
    mapping = {
        chr(q + 64): str(int(round(phred_quality_from_solexa(q))))
        for q in range(-5, 62 + 1)
    }
    return _fastq_convert_qual(in_file, out_file, mapping)


def _fastq_illumina_convert_qual(in_file, out_file):
    """Fast Illumina 1.3+ FASTQ to QUAL conversion (PRIVATE)."""
    mapping = {chr(q + 64): str(q) for q in range(0, 62 + 1)}
    return _fastq_convert_qual(in_file, out_file, mapping)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
