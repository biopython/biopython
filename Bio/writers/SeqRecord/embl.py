"""Part of an old unused and undocumented sequence writing framework (DEPRECATED)."""
# Not clear on the distinction, if any, between 'embl' and 'embl/65'.  This
# code might apply to either or both.

# See 'http://www.ebi.ac.uk/embl/Documentation/User_manual/usrman.html' for a
# definition of this file format.

# This code only makes a best effort--the output may not be strictly valid.
# So, for example, the EMBL ID is supposed to be alphanumeric, starting with a
# letter, but we don't check for this, etc.


# Example:
# ID   AA03518    standard; DNA; FUN; 237 BP.
# XX
# AC   U03518;
# XX
# DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
# DE   rRNA and 5.8S rRNA genes, partial sequence.
# XX
# SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
#      aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc        60
#      tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg       120
#      ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc       180
#      tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc          237
# //


import textwrap

from Bio import Alphabet
from Bio import Writer

class WriteEmbl(Writer.Writer):
    def __init__(self, outfile):
        Writer.Writer.__init__(self, outfile)
        
    def write(self, record):
        seq = record.seq
        assert seq.alphabet.size == 1, "cannot handle alphabet of size %d" % \
               seq.alphabet.size
        data = seq.data
        upperdata = data.upper()

# It'd be nice if the alphabet was usefully set, but for many interesting
# cases (e.g., reading from FASTA files), it's not.

        if isinstance(seq.alphabet, Alphabet.RNAAlphabet):
            molecule = 'mRNA'
            letters = ['A', 'C', 'G', 'U']
        else:
            molecule = 'DNA'
            letters = ['A', 'C', 'G', 'T']

        division = 'UNC'                # unknown

        self.outfile.write("ID   %s  standard; %s; %s; %d BP.\n"
                           % (record.id, molecule, division, len(data)))

        desclist = textwrap.wrap(record.description, 74)
        for l in desclist:
            self.outfile.write("DE   %s\n" % l)

        counts = [ upperdata.count(l) for l in letters ]
        othercount = len(upperdata) - sum(counts)

        countstring = ''.join([ " %d %s;" % p for p in zip(counts, letters) ])

        self.outfile.write("SQ   Sequence %s BP;%s %d other;\n"
                           % (len(data), countstring, othercount))

        rowlength = 60
        blocklength = 10
        for i in xrange(0, len(data), rowlength):
            self.outfile.write(" " * 5)
            row = data[i:i+rowlength]
            for b in xrange(0, rowlength, blocklength):
                block = row[b:b+blocklength]
                self.outfile.write("%-*s" % (blocklength+1, block))
            self.outfile.write("%9d\n" % min(i+rowlength, len(data)))

        self.outfile.write("//\n")


make_writer = WriteEmbl
