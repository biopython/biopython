# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import tempfile
from Bio.PDB import *
from PDBExceptions import PDBException
from AbstractPropertyMap import AbstractResiduePropertyMap
import re


__doc__="""
Use the DSSP program to calculate secondary structure and accessibility.
You need to have a working version of DSSP (and a license, free for 
academic use) in order to use this. For DSSP, see U{http://www.cmbi.kun.nl/gv/dssp/}.

The DSSP codes for secondary structure used here are:

    - H        Alpha helix (4-12)
    - B        Isolated beta-bridge residue
    - E        Strand
    - G        3-10 helix
    - I        pi helix
    - T        Turn
    - S        Bend
    - -        None
"""

# Match C in DSSP
_dssp_cys=re.compile('[a-z]')

# Maximal ASA of amino acids
# Values from Sander & Rost, (1994), Proteins, 20:216-226
# Used for relative accessibility
MAX_ACC={}
MAX_ACC["ALA"]=106.0
MAX_ACC["CYS"]=135.0
MAX_ACC["ASP"]=163.0
MAX_ACC["GLU"]=194.0
MAX_ACC["PHE"]=197.0
MAX_ACC["GLY"]=84.0
MAX_ACC["HIS"]=184.0
MAX_ACC["ILE"]=169.0
MAX_ACC["LYS"]=205.0
MAX_ACC["LEU"]=164.0
MAX_ACC["MET"]=188.0
MAX_ACC["ASN"]=157.0
MAX_ACC["PRO"]=136.0
MAX_ACC["GLN"]=198.0
MAX_ACC["ARG"]=248.0
MAX_ACC["SER"]=130.0
MAX_ACC["THR"]=142.0
MAX_ACC["VAL"]=142.0
MAX_ACC["TRP"]=227.0
MAX_ACC["TYR"]=222.0


def ss_to_index(ss):
    """
    Secondary structure symbol to index.
    H=0
    E=1
    C=2
    """
    if ss=='H':
        return 0
    if ss=='E':
        return 1
    if ss=='C':
        return 2
    assert(0)

def dssp_dict_from_pdb_file(in_file, DSSP="dssp"):
    """
    Create a DSSP dictionary from a PDB file.

    Example:
        >>> dssp_dict=dssp_dict_from_pdb_file("1fat.pdb")
        >>> aa, ss, acc=dssp_dict[('A', 1)]

    @param in_file: pdb file
    @type in_file: string

    @param DSSP: DSSP executable (argument to os.system)
    @type DSSP: string

    @return: a dictionary that maps (chainid, resid) to 
        amino acid type, secondary structure code and 
        accessibility.
    @rtype: {}
    """
    out_file=tempfile.mktemp()
    os.system(DSSP+" %s > %s" % (in_file, out_file))
    dict, keys=make_dssp_dict(out_file)
    # This can be dangerous...
    #os.system("rm "+out_file)
    return dict, keys

def make_dssp_dict(filename):
    """
    Return a DSSP dictionary that maps (chainid, resid) to
    aa, ss and accessibility, from a DSSP file.

    @param filename: the DSSP output file
    @type filename: string
    """
    dssp={}
    fp=open(filename, "r")
    start=0
    keys=[]
    for l in fp.readlines():
        sl=l.split()
        if sl[1]=="RESIDUE":
            # start
            start=1
            continue
        if not start:
            continue
        if l[9]==" ":
            # skip -- missing residue
            continue
        resseq=int(l[5:10])
        icode=l[10]
        chainid=l[11]
        aa=l[13]
        ss=l[16]
        if ss==" ":
            ss="-"
        acc=int(l[34:38])
        res_id=(" ", resseq, icode)
        dssp[(chainid, res_id)]=(aa, ss, acc)
        keys.append((chainid, res_id))
    fp.close()
    return dssp, keys


class DSSP(AbstractResiduePropertyMap):
    """
    Run DSSP on a pdb file, and provide a handle to the 
    DSSP secondary structure and accessibility.

    Note that DSSP can only handle one model.

    Example:
        >>> p=PDBParser()
        >>> structure=parser.get_structure("1fat.pdb")
        >>> model=structure[0]
        >>> dssp=DSSP(model, "1fat.pdb")
        >>> # print dssp data for a residue
        >>> secondary_structure, accessibility=dssp[(chain_id, res_id)]
    """
    def __init__(self, model, pdb_file, dssp="dssp"):
        """
        @param model: the first model of the structure
        @type model: L{Model}

        @param pdb_file: a PDB file
        @type pdb_file: string

        @param dssp: the dssp executable (ie. the argument to os.system)
        @type dssp: string
        """
        # create DSSP dictionary
        dssp_dict, dssp_keys=dssp_dict_from_pdb_file(pdb_file, dssp)
        dssp_map={}
        dssp_list=[]
        # Now create a dictionary that maps Residue objects to 
        # secondary structure and accessibility, and a list of 
        # (residue, (secondary structure, accessibility)) tuples
        for key in dssp_keys:
            chain_id, res_id=key
            chain=model[chain_id]
            res=chain[res_id]
            aa, ss, acc=dssp_dict[key]
            res.xtra["SS_DSSP"]=ss
            res.xtra["EXP_DSSP_ASA"]=acc
            # relative accessibility
            resname=res.get_resname()
            rel_acc=acc/MAX_ACC[resname]
            if rel_acc>1.0:
                rel_acc=1.0
            res.xtra["EXP_DSSP_RASA"]=rel_acc
            # Verify if AA in DSSP == AA in Structure
            # Something went wrong if this is not true!
            resname=to_one_letter_code[resname]
            if resname=="C":
                # DSSP renames C in C-bridges to a,b,c,d,...
                # - we rename it back to 'C'
                if _dssp_cys.match(aa):
                    aa='C'
            if not (resname==aa):
                raise PDBException("Structure/DSSP mismatch at "+str(res))
            dssp_map[key]=((res, ss, acc, rel_acc))
            dssp_list.append((res, ss, acc, rel_acc))
        AbstractResiduePropertyMap.__init__(self, dssp_map, dssp_keys, dssp_list)


if __name__=="__main__":

    import sys

    p=PDBParser()
    s=p.get_structure('X', sys.argv[1])

    model=s[0]

    d=DSSP(model, sys.argv[1])

    for r in d:
        print r

    print d.keys()

    print len(d)

    print d.has_key(('A', 1))

    print d[('A', 1)]

    print s[0]['A'][1].xtra



