"""Command line wrapper for bwa
"""

__docformat__ = "epytext en"

from Bio.Application import _Option, _Argument, _Switch, AbstractCommandline

class BwaIndexCommandline(AbstractCommandline):
    """Command line wrapper for Burrows Wheeler Aligner

    http://bio-bwa.sourceforge.net/

    Command:
		bwa index [-p prefix] [-a algoType] [-c] <in.db.fasta>

    Description:
        Index database sequences in the FASTA format.

    Options:
		-c	 Build color-space index. The input fast should be in nucleotide space.
		-p STR	 Prefix of the output database [same as db filename]
		-a STR	 Algorithm for constructing BWT index. Available options are:
				   is	 IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.
				   bwtsw	 Algorithm implemented in BWT-SW. This method works with the whole human genome, but it does not work with database smaller than 10MB and it is usually slower than IS.

    Example:
        >>> from Bio.Sequencing.Applications import BwaIndexCommandline
        >>> reference_genome = "/path/to/reference_genome.fasta"
        >>> index = BwaIndexCommandline(infile=reference_genome, algorithm="bwtsw")
        >>> index()




    """
    def __init__(self, cmd="bwa index", **kwargs):
        self.program_name = cmd
        self.parameters = \
                [
                    _Option(["-a","a","algorithm"],"Algorithm for constructing BWT index. Available options are:\n \
                            is:   IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.\n \
                            bwtsw:   Algorithm implemented in BWT-SW. This method works with the whole human genome, but it does not work with database smaller than 10MB and it is usually slower than IS.", \
                            filename=False,equate=False,checker_function=lambda x:x in ["is","bwtsw"],is_required=True ),
                    _Option(["-p","p","prefix"],"Prefix of the output database [same as db filename]",filename=False,equate=False,is_required=False),
                    _Argument(["infile"],"Input file name", filename=True, is_required=True),
                    _Switch(["-c","c"],"Build color-space index. The input fast should be in nucleotide space.")
                ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class BwaAlignCommandline(AbstractCommandline):
    """Command line wrapper for Burrows Wheeler Aligner

    http://bio-bwa.sourceforge.net/

	Command:
		bwa aln [-n maxDiff] [-o maxGapO] [-e maxGapE] [-d nDelTail] [-i nIndelEnd]
		[-k maxSeedDiff] [-l seedLen] [-t nThrds] [-cRN] [-M misMsc] [-O gapOsc]
		[-E gapEsc] [-q trimQual] <in.db.fasta> <in.query.fq> > <out.sai>

    Description:
        Find the SA coordinates of the input reads. Maximum maxSeedDiff differences are allowed
        in the first seedLen subsequence and maximum maxDiff differences are allowed in the whole sequence.

	Options:

		-n NUM	 Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]
		-o INT	 Maximum number of gap opens [1]
		-e INT	 Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]
		-d INT	 Disallow a long deletion within INT bp towards the 3-end [16]
		-i INT	 Disallow an indel within INT bp towards the ends [5]
		-l INT	 Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for '-k 2'. [inf]
		-k INT	 Maximum edit distance in the seed [2]
		-t INT	 Number of threads (multi-threading mode) [1]
		-M INT	 Mismatch penalty. BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc). [3]
		-O INT	 Gap open penalty [11]
		-E INT	 Gap extension penalty [4]
		-R INT	 Proceed with suboptimal alignments if there are no more than INT equally best hits. This option only affects paired-end mapping. Increasing this threshold helps to improve the pairing accuracy at the cost of speed, especially for short reads (~32bp).
		-c	 Reverse query but not complement it, which is required for alignment in the color space.
		-N	 Disable iterative search. All hits with no more than maxDiff differences will be found. This mode is much slower than the default.
		-q INT	 Parameter for read trimming. BWA trims a read down to argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT where l is the original read length. [0]
		-I	 The input is in the Illumina 1.3+ read format (quality equals ASCII-64).
		-B INT	 Length of barcode starting from the 5'-end. When INT is positive, the barcode of each read will be trimmed before mapping and will be written at the BC SAM tag. For paired-end reads, the barcode from both ends are concatenated. [0]
		-b	 Specify the input read sequence file is the BAM format. For paired-end data, two ends in a pair must be grouped together and options -1 or -2 are usually applied to specify which end should be mapped. Typical command lines for mapping pair-end data in the BAM format are:
			bwa aln ref.fa -b1 reads.bam > 1.sai
			bwa aln ref.fa -b2 reads.bam > 2.sai
			bwa sampe ref.fa 1.sai 2.sai reads.bam reads.bam > aln.sam

			-0	 When -b is specified, only use single-end reads in mapping.
			-1	 When -b is specified, only use the first read in a read pair in mapping (skip single-end reads and the second reads).
			-2	 When -b is specified, only use the second read in a read pair in mapping.

	 Example:
        >>> from Bio.Sequencing.Applications import BwaAlignCommandline
        >>> reference_genome = "/path/to/reference_genome.fasta"
        >>> read_file = "/path/to/read_1.fq"
        >>> output_sai_file = "/path/to/read_1.sai"
        >>> read_group="@RG\tID:foo\tSM:bar"
        >>> align = BwaAlignCommandline(reference=reference_genome, read_file=read_file)
        >>> output = align(stdout=output_sai_file)



    """
    def __init__(self, cmd="bwa aln", **kwargs):
        self.program_name = cmd
        self.parameters = \
                [
                    _Argument(["reference"],"Reference file name", filename=True, is_required=True),
                    _Argument(["read_file"],"Read  file name", filename=True, is_required=True),
                    _Option(["-n","n"],"Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,(int,float))),
                    _Option(["-o","o"],"Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,(int,float))),
                    _Option(["-e","e"],"Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-d","d"],"Disallow a long deletion within INT bp towards the 3-end [16]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-i","i"],"Disallow an indel within INT bp towards the ends [5]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-l","l"],"Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for -k 2. [inf]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-k","k"],"Maximum edit distance in the seed [2]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-t","t"],"Number of threads (multi-threading mode) [1]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-M","M"],"Mismatch penalty. BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc). [3]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-O","O"],"Gap open penalty [11]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-E","E"],"Gap extension penalty [4]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-R","R"],"Proceed with suboptimal alignments if there are no more than INT equally best hits. This option only affects paired-end mapping. Increasing this threshold helps to improve the pairing accuracy at the cost of speed, especially for short reads (~32bp).",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-q","q"],"Parameter for read trimming. BWA trims a read down to argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT where l is the original read length. [0]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Option(["-B","B"],"Length of barcode starting from the 5-end. When INT is positive, the barcode of each read will be trimmed before mapping and will be written at the BC SAM tag. For paired-end reads, the barcode from both ends are concatenated. [0]",filename=False, equate=False,checker_function= lambda x :  isinstance(x,int)),
                    _Switch(["-c","c"],"Reverse query but not complement it, which is required for alignment in the color space."),
                    _Switch(["-N","N"],"Disable iterative search. All hits with no more than maxDiff differences will be found. This mode is much slower than the default."),
                    _Switch(["-I","I"],"The input is in the Illumina 1.3+ read format (quality equals ASCII-64)."),
                    _Switch(["-b","b"],"Specify the input read sequence file is the BAM format"),
                    _Switch(["-b1","b1"],"When -b is specified, only use the first read in a read pair in mapping (skip single-end reads and the second reads)."),
                    _Switch(["-b2","b2"],"When -b is specified, only use the second read in a read pair in mapping.")
                  ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
class BwaSamseCommandline(AbstractCommandline):
    """Command line wrapper for Burrows Wheeler Aligner

    http://bio-bwa.sourceforge.net/

    Command:
		bwa samse [-n maxOcc] <in.db.fasta> <in.sai> <in.fq> > <out.sam>

    Description:
        Generate alignments in the SAM format given single-end reads.
        Repetitive hits will be randomly chosen.

	Options:
		-n INT	 Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]
		-r STR	 Specify the read group in a format like '@RG\tID:foo\tSM:bar'. [null]

	Example:
        >>> from Bio.Sequencing.Applications import BwaSamseCommandline
        >>> reference_genome = "/path/to/reference_genome.fasta"
        >>> read_file = "/path/to/read_1.fq"
        >>> sai_file = "/path/to/read_1.sai"
        >>> output_sam_file = "/path/to/read_1.sam"
        >>> samse = BwaSamseCommandline(reference=reference_genome, read_file=read_file, sai_file=sai_file)
        >>> output = samse(stdout=output_sam_file)





    """
    def __init__(self, cmd="bwa samse", **kwargs):
        self.program_name = cmd
        self.parameters = \
                [
                    _Argument(["reference"],"Reference file name", filename=True, is_required=True),
                    _Argument(["sai_file"],"Sai file name", filename=True, is_required=True),
                    _Argument(["read_file"],"Read  file name", filename=True, is_required=True),
                    _Option(["-n","n"],"Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-r","r"],"Specify the read group in a format like '@RG\tID:foo\tSM:bar'. [null]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,basestring))

                  ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class BwaSampeCommandline(AbstractCommandline):
    """Command line wrapper for Burrows Wheeler Aligner

    http://bio-bwa.sourceforge.net/

    Command:
		bwa sampe [-a maxInsSize] [-o maxOcc] [-n maxHitPaired] [-N maxHitDis] [-P] <in.db.fasta> <in1.sai> <in2.sai> <in1.fq> <in2.fq> > <out.sam>

    Description:
        Generate alignments in the SAM format given paired-end
        reads. Repetitive read pairs will be placed randomly.

	Options:
		-a INT	 Maximum insert size for a read pair to be considered being mapped properly. Since 0.4.5, this option is only used when there are not enough good alignment to infer the distribution of insert sizes. [500]
        -o INT	 Maximum occurrences of a read for pairing. A read with more occurrneces will be treated as a single-end read. Reducing this parameter helps faster pairing. [100000]
        -P	     Load the entire FM-index into memory to reduce disk operations (base-space reads only). With this option, at least 1.25N bytes of memory are required, where N is the length of the genome.
        -n INT	 Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]
        -N INT	 Maximum number of alignments to output in the XA tag for disconcordant read pairs (excluding singletons). If a read has more than INT hits, the XA tag will not be written. [10]
        -r STR	 Specify the read group in a format like '@RG\tID:foo\tSM:bar'. [null]

    Example:
        >>> from Bio.Sequencing.Applications import BwaSampeCommandline
        >>> reference_genome = "/path/to/reference_genome.fasta"
        >>> read_file1 = "/path/to/read_1.fq"
        >>> read_file2 = "/path/to/read_2.fq"
        >>> sai_file1 = "/path/to/read_1.sai"
        >>> sai_file2 = "/path/to/read_2.sai"
        >>> output_sam_file = "/path/to/output.sam"
        >>> read_group="@RG\tID:foo\tSM:bar"
        >>> sampe = BwaSampeCommandline(reference=reference_genome, sai_file1=sai_file1, sai_file2=sai_file2, read_file1=read_file1, read_file2=read_file2, r=read_group)
        >>> output = sampe(stdout=output_sam_file)




    """
    def __init__(self, cmd="bwa sampe", **kwargs):
        self.program_name = cmd
        self.parameters = \
                [
                    _Argument(["reference"],"Reference file name", filename=True, is_required=True),
                    _Argument(["sai_file1"],"Sai file 1", filename=True, is_required=True),
                    _Argument(["sai_file2"],"Sai file 2", filename=True, is_required=True),
                    _Argument(["read_file1"],"Read  file 1", filename=True, is_required=True),
                    _Argument(["read_file2"],"Read  file 2", filename=True, is_required=True),
                    _Option(["-a","a"],"Maximum insert size for a read pair to be considered being mapped properly. Since 0.4.5, this option is only used when there are not enough good alignment to infer the distribution of insert sizes. [500]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-o","o"],"Maximum occurrences of a read for pairing. A read with more occurrneces will be treated as a single-end read. Reducing this parameter helps faster pairing. [100000]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-n","n"],"Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-N","N"],"Maximum number of alignments to output in the XA tag for disconcordant read pairs (excluding singletons). If a read has more than INT hits, the XA tag will not be written. [10]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-r","r"],"Specify the read group in a format like '@RG\tID:foo\tSM:bar'. [null]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,basestring)),


                  ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class BwaBwaswCommandline(AbstractCommandline):
    """Command line wrapper for Burrows Wheeler Aligner

    http://bio-bwa.sourceforge.net/

	Command:
		bwa bwasw [-a matchScore] [-b mmPen] [-q gapOpenPen] [-r gapExtPen] [-t nThreads] [-w bandWidth] [-T thres] [-s hspIntv] [-z zBest] [-N nHspRev] [-c thresCoef] <in.db.fasta> <in.fq>

    Description:
        Align query sequences in the in.fq file. When mate.fq is present, perform paired-end alignment.
        The paired-end mode only works for reads Illumina short-insert libraries. In the paired-end mode,
        BWA-SW may still output split alignments but they are all marked as not properly paired;
        the mate positions will not be written if the mate has multiple local hits


	Options:
	    -a INT	 Score of a match [1]
	    -b INT	 Mismatch penalty [3]
	    -q INT	 Gap open penalty [5]
	    -r INT	 Gap extension penalty. The penalty for a contiguous gap of size k is q+k*r. [2]
	    -t INT	 Number of threads in the multi-threading mode [1]
	    -w INT	 Band width in the banded alignment [33]
	    -T INT	 Minimum score threshold divided by a [37]
	    -c FLOAT	 Coefficient for threshold adjustment according to query length. Given an l-long query, the threshold for a hit to be retained is a*max{T,c*log(l)}. [5.5]
	    -z INT	 Z-best heuristics. Higher -z increases accuracy at the cost of speed. [1]
	    -s INT	 Maximum SA interval size for initiating a seed. Higher -s increases accuracy at the cost of speed. [3]
	    -N INT	 Minimum number of seeds supporting the resultant alignment to skip reverse alignment. [5]

	Example:
        >>> from Bio.Sequencing.Applications import BwaBwaswCommandline
        >>> reference_genome = "/path/to/reference_genome.fasta"
        >>> read_file = "/path/to/read_1.fq"
        >>> bwasw = BwaBwaswCommandline(reference=reference_genome, read_file=read_file)
        >>> output = bwasw()




    """
    def __init__(self, cmd="bwa bwasw", **kwargs):
        self.program_name = cmd
        self.parameters = \
                [
                    _Argument(["reference"],"Reference file name", filename=True, is_required=True),
                    _Argument(["read_file"],"Read file", filename=True, is_required=True),
                    _Argument(["mate_file"],"Mate file", filename=True, is_required=False),
                    _Option(["-a","a"],"Score of a match [1]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-b","b"],"Mismatch penalty [3]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-q","q"],"Gap open penalty [5]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-r","r"],"Gap extension penalty. The penalty for a contiguous gap of size k is q+k*r. [2]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-t","t"],"Number of threads in the multi-threading mode [1]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-w","w"],"Band width in the banded alignment [33]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-T","T"],"Minimum score threshold divided by a [37]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-c","c"],"Coefficient for threshold adjustment according to query length. Given an l-long query, the threshold for a hit to be retained is a*max{T,c*log(l)}. [5.5]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,float)),
                    _Option(["-z","z"],"Z-best heuristics. Higher -z increases accuracy at the cost of speed. [1]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-s","s"],"Maximum SA interval size for initiating a seed. Higher -s increases accuracy at the cost of speed. [3]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),
                    _Option(["-N","N"],"Minimum number of seeds supporting the resultant alignment to skip reverse alignment. [5]",filename=False, equate=False,checker_function=lambda x :  isinstance(x,int)),

                  ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


def _test():
        """Run the module's doctests (PRIVATE)."""
        print "Running modules doctests..."
        import doctest
        doctest.testmod()
        print "Done"
if __name__ == "__main__":
    _test()







