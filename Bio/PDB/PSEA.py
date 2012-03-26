# Copyright (C) 2006, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Wrappers for PSEA, a program for secondary structure assignment.

See this citation for P-SEA, PMID: 9183534

Labesse G, Colloc'h N, Pothier J, Mornon J-P:  P-SEA: a new efficient
assignment of secondary structure from C_alpha.
Comput Appl Biosci 1997 , 13:291-295

ftp://ftp.lmcp.jussieu.fr/pub/sincris/software/protein/p-sea/
"""

import os

from Bio.PDB.Polypeptide import is_aa


def run_psea(fname):
    """Run PSEA and return output filename.
    
    Note that this assumes the P-SEA binary is called "psea" and that it is
    on the path.
    
    Note that P-SEA will write an output file in the current directory using
    the input filename with extension ".sea".
    
    Note that P-SEA will write output to the terminal while run.
    """
    os.system("psea "+fname)
    last=fname.split("/")[-1]
    base=last.split(".")[0]
    return base+".sea"

def psea(pname):
    """Parse PSEA output file."""
    fname=run_psea(pname)
    start=0
    ss=""
    fp=open(fname, 'r')
    for l in fp.readlines():
        if l[0:6]==">p-sea":
            start=1
            continue
        if not start:
            continue
        if l[0]=="\n":
            break
        ss=ss+l[0:-1]
    fp.close()
    return ss

def psea2HEC(pseq):
    """Translate PSEA secondary structure string into HEC."""
    seq=[]
    for ss in pseq:
        if ss=="a":
            n="H"
        elif ss=="b":
            n="E"
        elif ss=="c":
            n="C"
        seq.append(n)
    return seq

def annotate(m, ss_seq):
    """Apply seconardary structure information to residues in model."""
    c=m.get_list()[0]
    all=c.get_list()
    residues=[]
    # Now remove HOH etc.
    for res in all:
        if is_aa(res):
            residues.append(res)
    L=len(residues)
    if not (L==len(ss_seq)):
        raise ValueError("Length mismatch %i %i" % (L, len(ss_seq)))
    for i in range(0, L):
        residues[i].xtra["SS_PSEA"]=ss_seq[i]
    #os.system("rm "+fname)

class PSEA(object):
    def __init__(self, model, filename):
        ss_seq=psea(filename)
        ss_seq=psea2HEC(ss_seq)
        annotate(model, ss_seq)
        self.ss_seq=ss_seq

    def get_seq(self):
        """
        Return secondary structure string.
        """
        return self.ss_seq
        

if __name__=="__main__":

    import sys
    from Bio.PDB import PDBParser

    # Parse PDB file
    p=PDBParser()
    s=p.get_structure('X', sys.argv[1])

    # Annotate structure with PSEA sceondary structure info
    PSEA(s[0], sys.argv[1])
