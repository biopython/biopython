from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC

from Bio.CodonAlign import build

from Bio.CodonAlign.CodonAlignment import mktest

pro_aln = AlignIO.read('adh.aln', 'clustal', alphabet=IUPAC.protein)
p = SeqIO.index('drosophilla.fasta', 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA())

codon_aln = build(pro_aln, p)
print codon_aln
print mktest([codon_aln[1:12], codon_aln[12:16], codon_aln[16:]])
