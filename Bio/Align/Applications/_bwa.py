"""Command line wrapper for bwa
"""

__docformat__ = "epytext en"

from Bio.Application import _Option, _Argument, _Switch, AbstractCommandline

class BwaCommandline(AbstractCommandline):
    """Command line wrapper for Burrows Wheeler Aligner

    http://bio-bwa.sourceforge.net/

    Example to reference a fastq file :
        >>> from Bio.Align.Applications import BwaCommandline
        >>> reference_genome = "reference_genome.fasta"
        >>> index_gnome =
        >>> BwaCommandline(infile=reference_genome, algorithm="bwt,prefix=None, color_space=False)
        bwa index -a bwtsw reference_genome.fastadsadsa
    """
    def __init__(self, cmd="bwa", **kwargs):
        self.program_name = cmd
        self.parameters = \
                [
                    _Argument(["action"],"Perform on eof the following actions :\
                              1.index \n \
                              2.aln\n \
                              3.samse \n \
                              4.sampe \n \
                              5.bwasw",checker_function= lambda x : x in ["index","aln","samse","sampe","bwasw"],is_filename=False,is_requred=True,equate=False),
                    _Option(["-a","a"],"Algorithm for constructing BWT index. Available options are:\n \
                            is:   IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.\n \
                            bwtsw:   Algorithm implemented in BWT-SW. This method works with the whole human genome, but it does not work with database smaller than 10MB and it is usually slower than IS.", \
                            filename=False,equate=False,checker_function=lambda x:x in ["is","bwtsw"]

                            ),
                    _Option(["-p","p"],"Prefix of the output database [same as db filename]",filename=False,equate=False,is_required=False),
                    _Argument(["input"],"Input file name", filename=True, is_required=True),
                    _Switch(["-c","c"],"     Build color-space index. The input fast should be in nucleotide space.",is_required=False)
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







