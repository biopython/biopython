from Numeric import sum
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

def is_aa(residue):
    """
    Return 1 if residue object/string is an amino acid.

    @param residue: a L{Residue} object OR a three letter amino acid code
    @type residue: L{Residue} or string
    """
    if not type(residue)==StringType:
        residue=residue.get_resname()
    residue=residue.upper()
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
        return tau_list

    def get_omega_list(self):
        """
        Return list of omega angles for all 3 consecutive
        Calpha atoms.
        """
        omega_list=[]
        ca_list=self.get_ca_list()
        for i in range(0, len(ca_list)-2):
            atom_list=[ca_list[i], ca_list[i+1], ca_list[i+2]]
            vector_list=map(lambda a: a.get_vector(), atom_list)
            v1, v2, v3=vector_list
            omega=calc_angle(v1, v2, v3)
            omega_list.append(omega)
        return omega_list

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
            # not a standard AA so skip
            return 1
        else:
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
            raise PDBException, "Entity should be Structure, Model or Chain."
        pp_list=[]
        for chain in chain_list:
            prev=None
            pp=None
            for next in chain.get_list():
                if prev:
                    if is_connected(prev, next):
                        if pp is None:
                            pp=Polypeptide()
                            pp.append(prev)
                            pp_list.append(pp)
                        pp.append(next)
                    else:
                        pp=None
                if aa_only:
                    if accept(next):
                        prev=next
                    else:
                        prev=None
                        pp=None
                else:
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
        n_ca=next["CA"]
        p_ca=prev["CA"]
        for a in [n_ca, p_ca]:
            # if a CA is disordered we consider
            # it as not part of a polypeptide 
            # chain
            if a.is_disordered():
                return 0
        if (n_ca-p_ca)<self.radius:
            return 1
        else:
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
                # identifier
                n_altloc=nn.get_altloc()
                c_altloc=cc.get_altloc()
                if n_altloc==c_altloc: 
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

    ppb=CaPPBuilder()

    print "CA-CA"
    for pp in ppb.build_peptides(s):
        print pp.get_sequence()
    for pp in ppb.build_peptides(s[0]):
        print pp.get_sequence()
    for pp in ppb.build_peptides(s[0]["A"]):
        print pp.get_sequence()



