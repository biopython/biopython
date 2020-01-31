# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Command line wrapper for samtools."""
# Last Checked with samtools [0.1.20 and 1.2]
# TODO samtools 1.x has additional options over 0.x which
# are missing from this wrapper


from Bio.Application import _Option, _Argument, _Switch
from Bio.Application import AbstractCommandline, _ArgumentList
from Bio.Application import _StaticArgument


class SamtoolsViewCommandline(AbstractCommandline):
    """Command line wrapper for samtools view.

    Extract/print all or sub alignments in SAM or BAM format, equivalent to::

        $ samtools view [-bchuHS] [-t in.refList] [-o output] [-f reqFlag]
                        [-F skipFlag] [-q minMapQ] [-l library] [-r readGroup]
                        [-R rgFile] <in.bam>|<in.sam> [region1 [...]]

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsViewCommandline
    >>> input_file = "/path/to/sam_or_bam_file"
    >>> samtools_view_cmd = SamtoolsViewCommandline(input_file=input_file)
    >>> print(samtools_view_cmd)
    samtools view /path/to/sam_or_bam_file

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("view"),
            _Switch(["-b", "b"], "Output in the BAM format"),
            _Switch(
                ["-c", "c"],
                """Instead of printing the alignments, only count them and
                    print the total number.

                    All filter options, such as '-f', '-F' and '-q',
                    are taken into account""",
            ),
            _Switch(["-h", "h"], "Include the header in the output"),
            _Switch(
                ["-u", "u"],
                """Output uncompressed BAM.

                    This option saves time spent on compression/decompression
                    and is thus preferred when the output is piped to
                    another samtools command""",
            ),
            _Switch(["-H", "H"], "Output the header only"),
            _Switch(
                ["-S", "S"],
                """Input is in SAM.
                    If @SQ header lines are absent,
                    the '-t' option is required.""",
            ),
            _Option(
                ["-t", "t"],
                """This file is TAB-delimited.
                    Each line must contain the reference name and the
                    length of the reference, one line for each
                    distinct reference; additional fields are ignored.

                    This file also defines the order of the reference
                    sequences in sorting.
                    If you run   'samtools faidx <ref.fa>',
                    the resultant index file <ref.fa>.fai can be used
                    as this <in.ref_list> file.""",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-o", "o"],
                "Output file",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-f", "f"],
                """Only output alignments with all bits in
                    INT present in the FLAG field""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-F", "F"],
                "Skip alignments with bits present in INT",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-q", "q"],
                "Skip alignments with MAPQ smaller than INT",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-r", "r"],
                "Only output reads in read group STR",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-R", "R"],
                "Output reads in read groups listed in FILE",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-l", "l"],
                "Only output reads in library STR",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Switch(
                ["-1", "fast_bam"],
                "Use zlib compression level 1 to compress the output",
            ),
            _Argument(
                ["input", "input_file"],
                "Input File Name",
                filename=True,
                is_required=True,
            ),
            _Argument(["region"], "Region", is_required=False),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsMpileupCommandline(AbstractCommandline):
    """Command line wrapper for samtools mpileup.

    Generate BCF or pileup for one or multiple BAM files, equivalent to::

        $ samtools mpileup [-EBug] [-C capQcoef] [-r reg] [-f in.fa]
                           [-l list] [-M capMapQ] [-Q minBaseQ]
                           [-q minMapQ] in.bam [in2.bam [...]]

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsMpileupCommandline
    >>> input = ["/path/to/sam_or_bam_file"]
    >>> samtools_mpileup_cmd = SamtoolsMpileupCommandline(input_file=input)
    >>> print(samtools_mpileup_cmd)
    samtools mpileup /path/to/sam_or_bam_file

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("mpileup"),
            _Switch(
                ["-E", "E"],
                """Extended BAQ computation.
                    This option helps sensitivity especially
                    for MNPs, but may hurt specificity a little bit""",
            ),
            _Switch(
                ["-B", "B"],
                """Disable probabilistic realignment for the
                    computation of base alignment quality (BAQ).

                    BAQ is the Phred-scaled probability of a read base being
                    misaligned.
                    Applying this option greatly helps to reduce false SNPs
                    caused by misalignments""",
            ),
            _Switch(
                ["-g", "g"],
                """Compute genotype likelihoods and output them in the
                    binary call format (BCF)""",
            ),
            _Switch(
                ["-u", "u"],
                """Similar to -g except that the output is
                    uncompressed BCF, which is preferred for piping""",
            ),
            _Option(
                ["-C", "C"],
                """Coefficient for downgrading mapping quality for
                    reads containing excessive mismatches.

                    Given a read with a phred-scaled probability q of
                    being generated from the mapped position,
                    the new mapping quality is about sqrt((INT-q)/INT)*INT.
                    A zero value disables this functionality;
                    if enabled, the recommended value for BWA is 50""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-r", "r"],
                "Only generate pileup in region STR",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-f", "f"],
                """The faidx-indexed reference file in the FASTA format.

                    The file can be optionally compressed by razip""",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-l", "l"],
                """BED or position list file containing a list of regions
                    or sites where pileup or BCF should be generated""",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-M", "M"],
                "Cap Mapping Quality at M",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-q", "q"],
                "Minimum mapping quality for an alignment to be used",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-Q", "Q"],
                "Minimum base quality for a base to be considered",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Switch(
                ["-6", "illumina_13"],
                "Assume the quality is in the Illumina 1.3+ encoding",
            ),
            _Switch(
                ["-A", "A"], "Do not skip anomalous read pairs in variant calling."
            ),
            _Option(
                ["-b", "b"],
                "List of input BAM files, one file per line",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-d", "d"],
                "At a position, read maximally INT reads per input BAM",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Switch(["-D", "D"], "Output per-sample read depth"),
            _Switch(
                ["-S", "S"],
                """Output per-sample Phred-scaled
                                strand bias P-value""",
            ),
            _Option(
                ["-e", "e"],
                """Phred-scaled gap extension sequencing error probability.

                    Reducing INT leads to longer indels""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-h", "h"],
                """Coefficient for modeling homopolymer errors.

                    Given an l-long homopolymer run, the sequencing error
                    of an indel of size s is modeled as INT*s/l""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Switch(["-I", "I"], "Do not perform INDEL calling"),
            _Option(
                ["-L", "L"],
                """Skip INDEL calling if the average per-sample
                    depth is above INT""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-o", "o"],
                """Phred-scaled gap open sequencing error probability.

                    Reducing INT leads to more indel calls.""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-p", "p"],
                """Comma delimited list of platforms (determined by @RG-PL)
                    from which indel candidates are obtained.

                    It is recommended to collect indel candidates from
                    sequencing technologies that have low indel error rate
                    such as ILLUMINA""",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _ArgumentList(
                ["input_file"],
                "Input File for generating mpileup",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsReheaderCommandline(AbstractCommandline):
    """Command line wrapper for samtools reheader.

    Replace the header in in.bam with the header
    in in.header.sam, equivalent to::

    $ samtools reheader <in.header.sam> <in.bam>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsReheaderCommandline
    >>> input_header = "/path/to/header_sam_file"
    >>> input_bam = "/path/to/input_bam_file"
    >>> reheader_cmd = SamtoolsReheaderCommandline(input_header=input_header,
    ...                                            input_bam=input_bam)
    >>> print(reheader_cmd)
    samtools reheader /path/to/header_sam_file /path/to/input_bam_file

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("reheader"),
            _Argument(
                ["input_header", "header_sam", "sam_file"],
                "Sam file with header",
                filename=True,
                is_required=True,
            ),
            _Argument(
                ["input_bam", "input_file", "bam_file"],
                "BAM file for writing header to",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsCatCommandline(AbstractCommandline):
    """Command line wrapper for samtools cat.

    Concatenate BAMs, equivalent to::

        $ samtools cat [-h header.sam] [-o out.bam] <in1.bam> <in2.bam> [ ... ]

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsCatCommandline
    >>> input_bam1 = "/path/to/input_bam1"
    >>> input_bam2 = "/path/to/input_bam2"
    >>> input_bams = [input_bam1, input_bam2]
    >>> samtools_cat_cmd = SamtoolsCatCommandline(input_bam=input_bams)
    >>> print(samtools_cat_cmd)
    samtools cat /path/to/input_bam1 /path/to/input_bam2

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("cat"),
            _Option(
                ["-h", "h"],
                "Header SAM file",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-o", "o"],
                "Output SAM file",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _ArgumentList(
                ["input", "input_bam", "bams"],
                "Input BAM files",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsVersion0xSortCommandline(AbstractCommandline):
    """Command line wrapper for samtools version 0.1.x sort.

    Concatenate BAMs, equivalent to::

    $ samtools sort [-no] [-m maxMem] <in.bam> <out.prefix>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsVersion0xSortCommandline
    >>> input_bam = "/path/to/input_bam"
    >>> out_prefix = "/path/to/out_prefix"
    >>> samtools_sort_cmd = SamtoolsVersion0xSortCommandline(input=input_bam, out_prefix=out_prefix)
    >>> print(samtools_sort_cmd)
    samtools sort /path/to/input_bam /path/to/out_prefix

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd

        # options for version samtools 0.0.19
        self.parameters = [
            _StaticArgument("sort"),
            _Switch(
                ["-o", "o"],
                """Output the final alignment
                                    to the standard output""",
            ),
            _Switch(
                ["-n", "n"],
                """Sort by read names rather
                                    than by chromosomal coordinates""",
            ),
            _Option(
                ["-m", "m"],
                "Approximately the maximum required memory",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Argument(["input"], "Input BAM file", filename=True, is_required=True),
            _Argument(["out_prefix"], "Output prefix", filename=True, is_required=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsVersion1xSortCommandline(AbstractCommandline):
    """Command line wrapper for samtools version 1.3.x sort.

    Concatenate BAMs, equivalent to::

    $ samtools sort [-n] [-T FREFIX] [-o file] [-I INT] [-m maxMem] <in.bam>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsVersion1xSortCommandline
    >>> input_bam = "/path/to/input_bam"
    >>> FREFIX = "/path/to/out_prefix"
    >>> file_name = "/path/to/out_file"
    >>> samtools_sort_cmd = SamtoolsVersion1xSortCommandline(input=input_bam, T=FREFIX, o=file_name)
    >>> print(samtools_sort_cmd)
    samtools sort -o /path/to/out_file -T /path/to/out_prefix /path/to/input_bam

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd

        # options for version samtools 1.3.1
        self.parameters = [
            _StaticArgument("sort"),
            _Switch(
                ["-n", "n"],
                """Sort by read names rather
                                    than by chromosomal coordinates""",
            ),
            _Option(
                ["-o", "o"],
                """(file) Write the final sorted output to FILE,
                    rather than to standard output""",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-O", "O"],
                """(FORMAT) Write the final output as sam, bam, or cram""",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-T", "T"],
                """(PREFIX) Write temporary files to PREFIX.nnnn.bam, or if the specified PREFIX
                    is an existing directory, to PREFIX/samtools.mmm.mmm.tmp.nnnn.bam,
                    where mmm is unique to this invocation of the sort command""",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-I", "I"],
                """(INT) Set the desired compression level for the final output file,
                    ranging from 0 (uncompressed) or 1 (fastest but minimal compression)
                    to 9 (best compression but slowest to write), similarly to gzip(1)'s compression level setting.""",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-m", "m"],
                "Approximately the maximum required memory",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Argument(
                ["input"], "Input SAM/BAM/CRAM file", filename=True, is_required=True
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsMergeCommandline(AbstractCommandline):
    """Command line wrapper for samtools merge.

    Merge multiple sorted alignments, equivalent to::

        $ samtools merge [-nur1f] [-h inh.sam] [-R reg]
                         <out.bam> <in1.bam> <in2.bam> [...]

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsMergeCommandline
    >>> out_bam = "/path/to/out_bam"
    >>> in_bam = ["/path/to/input_bam1", "/path/to/input_bam2"]
    >>> merge_cmd = SamtoolsMergeCommandline(out_bam=out_bam,
    ...                                      input_bam=in_bam)
    >>> print(merge_cmd)
    samtools merge /path/to/out_bam /path/to/input_bam1 /path/to/input_bam2

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("merge"),
            _Switch(
                ["-n", "n"],
                """The input alignments are sorted by read names
                    rather than by chromosomal coordinates""",
            ),
            _Switch(
                ["-r", "r"],
                """Attach an RG tag to each alignment.
                    The tag value is inferred from file names""",
            ),
            _Switch(["-u", "u"], "Uncompressed BAM output"),
            _Switch(
                ["-1", "fast_bam"],
                """Use zlib compression level 1
                                           to compress the output""",
            ),
            _Switch(
                ["-f", "f"],
                """Force to overwrite the
                                    output file if present""",
            ),
            _Option(
                ["-h", "h"],
                """Use the lines of FILE as '@'
                                    headers to be copied to out.bam""",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-R", "R"],
                "Merge files in the specified region indicated by STR",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Argument(
                ["output_bam", "out_bam", "out", "output"],
                "Output BAM file",
                filename=True,
                is_required=True,
            ),
            _ArgumentList(
                ["input_bam", "in_bam", "input", "bam"],
                "Input BAM",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsIndexCommandline(AbstractCommandline):
    """Command line wrapper for samtools index.

    Index sorted alignment for fast random access, equivalent to::

    $ samtools index <aln.bam>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsIndexCommandline
    >>> input = "/path/to/aln_bam"
    >>> samtools_index_cmd = SamtoolsIndexCommandline(input_bam=input)
    >>> print(samtools_index_cmd)
    samtools index /path/to/aln_bam

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("index"),
            _Argument(["input", "in_bam", "input_bam"], "BAM file to be indexed"),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsIdxstatsCommandline(AbstractCommandline):
    """Command line wrapper for samtools idxstats.

    Retrieve and print stats in the index file, equivalent to::

    $ samtools idxstats <aln.bam>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsIdxstatsCommandline
    >>> input = "/path/to/aln_bam"
    >>> samtools_idxstats_cmd = SamtoolsIdxstatsCommandline(input_bam=input)
    >>> print(samtools_idxstats_cmd)
    samtools idxstats /path/to/aln_bam

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("idxstats"),
            _Argument(["input", "in_bam", "input_bam"], "BAM file to be indexed"),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsFaidxCommandline(AbstractCommandline):
    """Command line wrapper for samtools faidx.

    Retrieve and print stats in the index file, equivalent to::

    $ samtools faidx <ref.fasta> [region1 [...]]

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsFaidxCommandline
    >>> reference = "/path/to/reference.fasta"
    >>> samtools_faidx_cmd = SamtoolsFaidxCommandline(reference=reference)
    >>> print(samtools_faidx_cmd)
    samtools faidx /path/to/reference.fasta

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("faidx"),
            _Argument(
                ["reference", "reference_fasta", "ref"],
                "Reference FASTA to be indexed",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsFixmateCommandline(AbstractCommandline):
    """Command line wrapper for samtools fixmate.

    Fill in mate coordinates, ISIZE and mate related
    flags from a name-sorted alignment, equivalent to::

    $ samtools fixmate <in.nameSrt.bam> <out.bam>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsFixmateCommandline
    >>> in_bam = "/path/to/in.nameSrt.bam"
    >>> out_bam = "/path/to/out.bam"
    >>> fixmate_cmd = SamtoolsFixmateCommandline(input_bam=in_bam,
    ...                                          out_bam=out_bam)
    >>> print(fixmate_cmd)
    samtools fixmate /path/to/in.nameSrt.bam /path/to/out.bam

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("fixmate"),
            _Argument(
                ["in_bam", "sorted_bam", "input_bam", "input", "input_file"],
                "Name Sorted Alignment File ",
                filename=True,
                is_required=True,
            ),
            _Argument(
                ["out_bam", "output_bam", "output", "output_file"],
                "Output file",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsRmdupCommandline(AbstractCommandline):
    """Command line wrapper for samtools rmdup.

    Remove potential PCR duplicates, equivalent to::

    $ samtools rmdup [-sS] <input.srt.bam> <out.bam>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsRmdupCommandline
    >>> input_sorted_bam = "/path/to/input.srt.bam"
    >>> out_bam = "/path/to/out.bam"
    >>> rmdup_cmd = SamtoolsRmdupCommandline(input_bam=input_sorted_bam,
    ...                                      out_bam=out_bam)
    >>> print(rmdup_cmd)
    samtools rmdup /path/to/input.srt.bam /path/to/out.bam

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("rmdup"),
            _Switch(
                ["-s", "s"],
                """Remove duplicates for single-end reads.

                    By default, the command works for paired-end
                    reads only""",
            ),
            _Switch(
                ["-S", "S"],
                """Treat paired-end reads
                                    as single-end reads""",
            ),
            _Argument(
                ["in_bam", "sorted_bam", "input_bam", "input", "input_file"],
                "Name Sorted Alignment File ",
                filename=True,
                is_required=True,
            ),
            _Argument(
                ["out_bam", "output_bam", "output", "output_file"],
                "Output file",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsCalmdCommandline(AbstractCommandline):
    """Command line wrapper for samtools calmd.

    Generate the MD tag, equivalent to::

    $ samtools calmd [-EeubSr] [-C capQcoef] <aln.bam> <ref.fasta>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsCalmdCommandline
    >>> input_bam = "/path/to/aln.bam"
    >>> reference_fasta = "/path/to/reference.fasta"
    >>> calmd_cmd = SamtoolsCalmdCommandline(input_bam=input_bam,
    ...                                      reference=reference_fasta)
    >>> print(calmd_cmd)
    samtools calmd /path/to/aln.bam /path/to/reference.fasta

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("calmd"),
            _Switch(
                ["-E", "E"],
                """Extended BAQ calculation.
                    This option trades specificity for sensitivity,
                    though the effect is minor.""",
            ),
            _Switch(
                ["-e", "e"],
                """Convert the read base to = if it is
                    identical to the aligned reference base.

                    Indel caller does not support the = bases
                    at the moment.""",
            ),
            _Switch(["-u", "u"], "Output uncompressed BAM"),
            _Switch(["-b", "b"], "Output compressed BAM "),
            _Switch(["-S", "S"], "The input is SAM with header lines "),
            _Switch(
                ["-r", "r"],
                """Compute the BQ tag (without -A)
                    or cap base quality by BAQ (with -A).""",
            ),
            _Switch(
                ["-A", "A"],
                """When used jointly with -r this option overwrites
                    the original base quality""",
            ),
            _Option(
                ["-C", "C"],
                """Coefficient to cap mapping quality
                    of poorly mapped reads.

                    See the pileup command for details.""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Argument(
                ["input", "input_file", "in_bam", "infile", "input_bam"],
                "Input BAM",
                filename=True,
                is_required=True,
            ),
            _Argument(
                ["reference", "reference_fasta", "ref"],
                "Reference FASTA to be indexed",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsTargetcutCommandline(AbstractCommandline):
    """Command line wrapper for samtools targetcut.

    This command identifies target regions by examining the continuity
    of read depth, computes haploid consensus sequences of targets
    and outputs a SAM with each sequence corresponding to a target,
    equivalent to::

        $ samtools targetcut [-Q minBaseQ] [-i inPenalty] [-0 em0]
                             [-1 em1] [-2 em2] [-f ref] <in.bam>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsTargetcutCommandline
    >>> input_bam = "/path/to/aln.bam"
    >>> samtools_targetcut_cmd = SamtoolsTargetcutCommandline(input_bam=input_bam)
    >>> print(samtools_targetcut_cmd)
    samtools targetcut /path/to/aln.bam

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("targetcut"),
            _Option(
                ["-Q", "Q"],
                "Minimum Base Quality ",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-i", "i"],
                "Insertion Penalty",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-f", "f"],
                "Reference Filename",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-0", "em0"],
                "em0",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-1", "em1"],
                "em1",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Option(
                ["-2", "em2"],
                "em2",
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Argument(
                ["input", "input_bam", "in_bam"],
                "Input file",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SamtoolsPhaseCommandline(AbstractCommandline):
    """Command line wrapper for samtools phase.

    Call and phase heterozygous SNPs, equivalent to::

        $ samtools phase [-AF] [-k len] [-b prefix]
                         [-q minLOD] [-Q minBaseQ] <in.bam>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Examples
    --------
    >>> from Bio.Sequencing.Applications import SamtoolsPhaseCommandline
    >>> input_bam = "/path/to/in.bam"
    >>> samtools_phase_cmd = SamtoolsPhaseCommandline(input_bam=input_bam)
    >>> print(samtools_phase_cmd)
    samtools phase /path/to/in.bam

    """

    def __init__(self, cmd="samtools", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("phase"),
            _Argument(
                ["input", "input_bam", "in_bam"],
                "Input file",
                filename=True,
                is_required=True,
            ),
            _Switch(["-A", "A"], "Drop reads with ambiguous phase"),
            _Option(
                ["-b", "b"],
                "Prefix of BAM output",
                filename=True,
                equate=False,
                checker_function=lambda x: isinstance(x, str),
            ),
            _Switch(["-F", "F"], "Do not attempt to fix chimeric reads"),
            _Option(
                ["-k", "k"],
                "Maximum length for local phasing",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-q", "q"],
                """Minimum Phred-scaled LOD to
                    call a heterozygote""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-Q", "Q"],
                """Minimum base quality to be
                    used in het calling""",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
