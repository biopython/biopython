# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from types import StringType

from Bio.Alphabet import ProteinAlphabet
from Bio.Seq import Seq
from Bio.SCOP.Raf import to_one_letter_code
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Residue import Residue, DisorderedResidue
from Vector import calc_dihedral, calc_angle

__doc__="""
Polypeptide related classes (construction and representation).

Example:

    >>> ppb=PPBuilder()
    >>> for pp in ppb.build_peptides(structure):
    >>>     print pp.get_sequence()
"""

standard_aa_names=["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", 
                   "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",
                   "TRP", "TYR"]


aa1="ACDEFGHIKLMNPQRSTVWY"
aa3=standard_aa_names

d1_to_index={}
dindex_to_1={}
d3_to_index={}
dindex_to_3={}

# Create some lookup tables
for i in range(0, 20):
    n1=aa1[i]
    n3=aa3[i]
    d1_to_index[n1]=i
    dindex_to_1[i]=n1
    d3_to_index[n3]=i
    dindex_to_3[i]=n3

def index_to_one(index):
    """
    Index to corresponding one letter amino acid name.
    For example: 0 to A.
    """
    return dindex_to_1[index]

def one_to_index(s):
    """
    One letter code to index.
    For example: A to 0.
    """
    return d1_to_index[s]

def index_to_three(i):
    """
    Index to corresponding three letter amino acid name.
    For example: 0 to ALA.
    """
    return dindex_to_3[i]

def three_to_index(s):
    """
    Three letter code to index.
    For example: ALA to 0.
    """
    return d3_to_index[s]

def three_to_one(s):
    """
    Three letter code to one letter code.
    For example: ALA to A.
    """
    i=d3_to_index[s]
    return dindex_to_1[i]

def one_to_three(s):
    """
    One letter code to three letter code.
    For example: A to ALA.
    """
    i=d1_to_index[s]
    return dindex_to_3[i]

def is_aa(residue, standard=0):
    """
    Return 1 if residue object/string is an amino acid.

    @param residue: a L{Residue} object OR a three letter amino acid code
    @type residue: L{Residue} or string

    @param standard: flag to check for the 20 AA (default false) 
    @type standard: boolean
    """
    if not type(residue)==StringType:
        residue=residue.get_resname()
    residue=residue.upper()
    if standard:
        return d3_to_index.has_key(residue)
    else:
        return to_one_letter_code.has_key(residue)


class Polypeptide(list):
    """
    A polypeptide is simply a list of L{Residue} objects.
    """
    def get_ca_list(self):
        """
        @return: the list of C-alpha atoms
        @rtype: [L{Atom}, L{Atom}, ...]
        """
        ca_list=[]
        for res in self:
            ca=res["CA"]
            ca_list.append(ca)
        return ca_list

    def get_phi_psi_list(self):
        """
        Return the list of phi/psi dihedral angles
        """
        ppl=[]
        lng=len(self)
        for i in range(0, lng):
            res=self[i]
            try:
                n=res['N'].get_vector()
                ca=res['CA'].get_vector()
                c=res['C'].get_vector()
            except:
                # Some atoms are missing
                # Phi/Psi cannot be calculated for this residue
                ppl.append((None, None))
                res.xtra["PHI"]=None
                res.xtra["PSI"]=None
                continue
            # Phi
            if i>0:
                rp=self[i-1]
                try:
                    cp=rp['C'].get_vector()
                    phi=calc_dihedral(cp, n, ca, c)
                except:
                    phi=None
            else:
                # No phi for residue 0!
                phi=None
            # Psi
            if i<(lng-1):
                rn=self[i+1]
                try:
                    nn=rn['N'].get_vector()
                    psi=calc_dihedral(n, ca, c, nn)
                except:
                    psi=None
            else:
                # No psi for last residue!
                psi=None
            ppl.append((phi, psi))
            # Add Phi/Psi to xtra dict of residue
            res.xtra["PHI"]=phi
            res.xtra["PSI"]=psi
        return ppl

    def get_tau_list(self):
        """
        Return list of tau torsions angles for all 4 consecutive
        Calpha atoms.
        """
        ca_list=self.get_ca_list()
        tau_list=[]
        for i in range(0, len(ca_list)-3):
            atom_list=[ca_list[i], ca_list[i+1], ca_list[i+2], ca_list[i+3]]
            vector_list=map(lambda a: a.get_vector(), atom_list)
            v1, v2, v3, v4=vector_list
            tau=calc_dihedral(v1, v2, v3, v4)
            tau_list.append(tau)
            # Put tau in xtra dict of residue
            res=ca_list[i+2].get_parent()
            res.xtra["TAU"]=tau
        return tau_list

    def get_theta_list(self):
        """
        Return list of theta angles for all 3 consecutive
        Calpha atoms.
        """
        theta_list=[]
        ca_list=self.get_ca_list()
        for i in range(0, len(ca_list)-2):
            atom_list=[ca_list[i], ca_list[i+1], ca_list[i+2]]
            vector_list=map(lambda a: a.get_vector(), atom_list)
            v1, v2, v3=vector_list
            theta=calc_angle(v1, v2, v3)
            theta_list.append(theta)
            # Put tau in xtra dict of residue
            res=ca_list[i+1].get_parent()
            res.xtra["THETA"]=theta
        return theta_list

    def get_sequence(self):
        """
        Return the AA sequence.

        @return: polypeptide sequence 
        @rtype: L{Seq}
        """
        s=""
        for res in self:
            resname=res.get_resname()
            if to_one_letter_code.has_key(resname):
                resname=to_one_letter_code[resname]
            else:
                resname='X'
            s=s+resname
        seq=Seq(s, ProteinAlphabet)
        return seq

    def __repr__(self):
        """
        Return <Polypeptide start=START end=END>, where START
        and END are sequence identifiers of the outer residues.
        """
        start=self[0].get_id()[1]
        end=self[-1].get_id()[1]
        s="<Polypeptide start=%s end=%s>" % (start, end)
        return s

class _PPBuilder:
    """
    Base class to extract polypeptides.
    It checks if two consecutive residues in a chain 
    are connected. The connectivity test is implemented by a 
    subclass.
    """
    def __init__(self, radius):
        """
        @param radius: distance
        @type radius: float
        """
        self.radius=radius

    def _accept(self, residue):
        "Check if the residue is an amino acid."
        if is_aa(residue):
            return 1
        else:
            if "CA" in residue.child_dict :
                #It has an alpha carbon...
                #We probably need to update the hard coded list of
                #non-standard residues, see function is_aa for details.
                import warnings
                warnings.warn("Assuming residue %s is an unknown modified "
                              "amino acid" % residue.get_resname())
                return 1
            # not a standard AA so skip
            return 0
    
    def build_peptides(self, entity, aa_only=1):
        """
        Build and return a list of Polypeptide objects.

        @param entity: polypeptides are searched for in this object
        @type entity: L{Structure}, L{Model} or L{Chain}

        @param aa_only: if 1, the residue needs to be a standard AA
        @type aa_only: int
        """
        is_connected=self._is_connected
        accept=self._accept
        level=entity.get_level()
        # Decide wich entity we are dealing with
        if level=="S":
            model=entity[0]
            chain_list=model.get_list()
        elif level=="M":
            chain_list=entity.get_list()
        elif level=="C":
            chain_list=[entity]
        else:
            raise PDBException("Entity should be Structure, Model or Chain.")
        pp_list=[]
        for chain in chain_list:
            chain_it=iter(chain)
            prev=chain_it.next()
            pp=None
            for next in chain_it:
                if aa_only and not accept(prev):
                    prev=next
                    continue
                if is_connected(prev, next):
                    if pp is None:
                        pp=Polypeptide()
                        pp.append(prev)
                        pp_list.append(pp)
                    pp.append(next)
                else:
                    pp=None
                prev=next
        return pp_list


class CaPPBuilder(_PPBuilder):
    """
    Use CA--CA distance to find polypeptides.
    """
    def __init__(self, radius=4.3):
        _PPBuilder.__init__(self, radius)

    def _is_connected(self, prev, next):
        for r in [prev, next]:
            if not r.has_id("CA"):
                return 0
        n=next["CA"]
        p=prev["CA"]
        # Unpack disordered
        if n.is_disordered():
            nlist=n.disordered_get_list()
        else:
            nlist=[n]
        if p.is_disordered():
            plist=p.disordered_get_list()
        else:
            plist=[p]
        for nn in nlist:
            for pp in plist:
                if (nn-pp)<self.radius:
                    return 1
        return 0


class PPBuilder(_PPBuilder):
    """
    Use C--N distance to find polypeptides.
    """
    def __init__(self, radius=1.8):
        _PPBuilder.__init__(self, radius)

    def _is_connected(self, prev, next):
        if not prev.has_id("C"):
            return 0
        if not next.has_id("N"):
            return 0
        test_dist=self._test_dist
        c=prev["C"]
        n=next["N"]
        # Test all disordered atom positions!
        if c.is_disordered():
            clist=c.disordered_get_list()
        else:
            clist=[c]
        if n.is_disordered():
            nlist=n.disordered_get_list()
        else:
            nlist=[n]
        for nn in nlist:
            for cc in clist:
                # To form a peptide bond, N and C must be 
                # within radius and have the same altloc
                # identifier or one altloc blank
                n_altloc=nn.get_altloc()
                c_altloc=cc.get_altloc()
                if n_altloc==c_altloc or n_altloc==" " or c_altloc==" ": 
                    if test_dist(nn, cc):
                        # Select the disordered atoms that
                        # are indeed bonded
                        if c.is_disordered():
                            c.disordered_select(c_altloc)
                        if n.is_disordered():
                            n.disordered_select(n_altloc)
                        return 1
        return 0

    def _test_dist(self, c, n):
        "Return 1 if distance between atoms<radius"
        if (c-n)<self.radius:
            return 1
        else:
            return 0
    

if __name__=="__main__":

    import sys

    from Bio.PDB.PDBParser import PDBParser

    p=PDBParser(PERMISSIVE=1)

    s=p.get_structure("scr", sys.argv[1])

    ppb=PPBuilder()

    print "C-N"
    for pp in ppb.build_peptides(s):
        print pp.get_sequence()
    for pp in ppb.build_peptides(s[0]):
        print pp.get_sequence()
    for pp in ppb.build_peptides(s[0]["A"]):
        print pp.get_sequence()

    for pp in ppb.build_peptides(s):
        for phi, psi in pp.get_phi_psi_list():
            print phi, psi

    ppb=CaPPBuilder()

    print "CA-CA"
    for pp in ppb.build_peptides(s):
        print pp.get_sequence()
    for pp in ppb.build_peptides(s[0]):
        print pp.get_sequence()
    for pp in ppb.build_peptides(s[0]["A"]):
        print pp.get_sequence()



