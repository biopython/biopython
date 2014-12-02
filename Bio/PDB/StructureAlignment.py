# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Map the residues of two structures to each other based on a FASTA alignment
file.
"""

from __future__ import print_function

from Bio.Data import SCOPData

from Bio.PDB import Selection
from Bio.PDB.Polypeptide import is_aa

__docformat__ = "restructuredtext en"

class StructureAlignment(object):
    """
    This class aligns two structures based on an alignment of their
    sequences.
    """
    def __init__(self, fasta_align, m1, m2, si=0, sj=1):
        """
        Attributes:
    
            - fasta_align --- Alignment object
            - m1, m2 --- two models
            - si, sj --- the sequences in the Alignment object that
              correspond to the structures
        """
        l=fasta_align.get_alignment_length()
        # Get the residues in the models
        rl1=Selection.unfold_entities(m1, 'R')
        rl2=Selection.unfold_entities(m2, 'R')
        # Residue positions
        p1=0
        p2=0
        # Map equivalent residues to each other
        map12={}
        map21={}
        # List of residue pairs (None if -)
        duos=[]
        for i in range(0, l):
            column=fasta_align.get_column(i)
            aa1=column[si]
            aa2=column[sj]
            if aa1!="-":
                # Position in seq1 is not -
                while True:
                    # Loop until an aa is found
                    r1=rl1[p1]
                    p1=p1+1
                    if is_aa(r1):
                        break
                self._test_equivalence(r1, aa1)
            else:
                r1=None
            if aa2!="-":
                # Position in seq2 is not -
                while True:
                    # Loop until an aa is found
                    r2=rl2[p2]
                    p2=p2+1
                    if is_aa(r2):
                        break
                self._test_equivalence(r2, aa2)
            else:
                r2=None
            if r1:
                # Map residue in seq1 to its equivalent in seq2
                map12[r1]=r2
            if r2:
                # Map residue in seq2 to its equivalent in seq1
                map21[r2]=r1
            # Append aligned pair (r is None if gap)
            duos.append((r1, r2))
        self.map12=map12
        self.map21=map21
        self.duos=duos

    def _test_equivalence(self, r1, aa1):
        "Test if aa in sequence fits aa in structure."
        resname=r1.get_resname()
        resname=SCOPData.protein_letters_3to1[resname]
        assert(aa1==resname)

    def get_maps(self):
        """
        Return two dictionaries that map a residue in one structure to
        the equivealent residue in the other structure.
        """
        return self.map12, self.map21

    def get_iterator(self):
        """
        Iterator over all residue pairs.
        """
        for i in range(0, len(self.duos)):
            yield self.duos[i]


if __name__=="__main__":
    import sys
    from Bio.Alphabet import generic_protein
    from Bio import AlignIO
    from Bio.PDB import PDBParser

    if len(sys.argv) != 4:
        print("Expects three arguments,")
        print(" - FASTA alignment filename (expect two sequences)")
        print(" - PDB file one")
        print(" - PDB file two")
        sys.exit()

    # The alignment
    fa=AlignIO.read(open(sys.argv[1]), "fasta", generic_protein)

    pdb_file1=sys.argv[2]
    pdb_file2=sys.argv[3]

    # The structures
    p=PDBParser()
    s1=p.get_structure('1', pdb_file1)
    p=PDBParser()
    s2=p.get_structure('2', pdb_file2)

    # Get the models
    m1=s1[0]
    m2=s2[0]

    al=StructureAlignment(fa, m1, m2)

    # Print aligned pairs (r is None if gap)
    for (r1, r2) in al.get_iterator():
        print("%s %s" % (r1, r2))
