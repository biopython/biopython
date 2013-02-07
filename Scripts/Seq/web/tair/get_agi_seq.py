"""Demonstration use-case of the Bio.Seq.web.tair module.
See comments in source code for pointers on the use of this module.
"""

from Bio.Seq.web import tair
from Bio import SeqIO
from optparse import OptionParser


# Get the commandline options, this section is unrelated to the function of the
# tair module.
parser = OptionParser()
parser.add_option('-f', '--file', dest='filename', default=None)
parser.add_option('-d', '--dataset', dest='dataset', default="transcript",
        help="one of: upstream_500, downstream_3000, intergenic, 5prime_utr, \
        upstream_1000, intron, downstream_500, cds, 3prime_utr, genomic, \
        protein, gene, transcript, upstream_3000, downstream_1000")
parser.add_option('-a', '--agis', dest='agis',
        default="AT5G63980.1,AT5G63280.1,AT5G23980.1",
        help='Comma seperated list of AGIs to fetch')
(options, args) = parser.parse_args()


# The tair module get functions require a AGIs to be given as a python list
# of strings, containing AGIs. This gets a list of strings from the CSV
# string given.
agis = options.agis.split(",")

# Get the sequences
seqs = tair.get(agis, options.dataset, "representative")

if options.filename is None:
    # If a filename is not specified, print a fasta to stdout
    for seq in seqs:
        print seq.format("fasta")
else:
    # If we have a filename, write a fasta to the file specified
    file_handle = open(options.filename, "wb")
    SeqIO.write(seqs, file_handle, "fasta")
    file_handle.close()
