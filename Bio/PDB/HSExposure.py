from Numeric import matrixmultiply, transpose, array, sqrt
from math import pi
from MLab import eye

from Bio.PDB import *


class HSExposure:
    """
    Calculates the Half-Sphere Exposure (HSE).
    """
    # unit vector along z
    # CA-CB vectors or the CA-CA-CA approximation will
    # be rotated onto this unit vector
    unit_z=Vector(0.0, 0.0, 1.0)

    def __init__(self):
        # List of CA-CB direction calculated from CA-CA-CA
        # Used for the PyMol script
        self.ca_cb_list=[]

    def write_pymol_script(self, filename="pymol.py"):
        """
        Write a PyMol script that visualizes the pseudo CB-CA directions 
        at the CA coordinates.
        """
        if len(self.ca_cb_list)==0:
            print "Nothing to draw."
            return
        fp=open(filename, "w")
        fp.write("from pymol.cgo import *\n")
        fp.write("from pymol import cmd\n")
        fp.write("obj=[\n")
        fp.write("BEGIN, LINES,\n")
        fp.write("COLOR, %.2f, %.2f, %.2f,\n" % (1.0, 1.0, 1.0))
        for (ca, cb) in self.ca_cb_list:
            x,y,z=ca.get_array()
            fp.write("VERTEX, %.2f, %.2f, %.2f,\n" % (x,y,z))
            x,y,z=cb.get_array()
            fp.write("VERTEX, %.2f, %.2f, %.2f,\n" % (x,y,z))
        fp.write("END]\n")
        fp.write("cmd.load_cgo(obj, 'HS')\n")
        fp.close()

    def _get_gly_cb_coord(self, residue):
        """
        Return a pseudo CB coord for a Gly residue.
        
        CB coord=N coord rotated over -120 degrees 
        along the CA-C axis.
        """
        n=residue["N"].get_vector()
        c=residue["C"].get_vector()
        ca=residue["CA"].get_vector()
        # center at origin
        n=n-ca
        c=c-ca
        # rotation around c-ca over -120 deg
        rot=rotaxis(-pi*120.0/180.0, c)
        cb_at_origin=n.left_multiply(rot)
        # move back to ca position
        cb=cb_at_origin+ca
        return cb.get_array()

    def _get_cb_from_ca(self, ca1, ca2, ca3):
        """
        Calculate the approximate CA-CB direction for a central
        CA atom based on the two flanking CA positions. 
        
        The CA-CB vector is centered at the origin.
        """
        # center
        ca1=ca1-ca2
        ca3=ca3-ca2
        ca1=ca1.normalize()
        ca3=ca3.normalize()
        # bisection
        b=(ca1+ca3).normalize()
        b=-b
        # Add to ca_cb_list for drawing
        self.ca_cb_list.append((ca2, b+ca2))
        # b is centered at the origin!
        return b

    def _get_data_list_from_cb(self, residue_list):
        """
        Return a list of (translation, rotation, ca, residue)
        tuples (using CB). 
        
        The translation and rotation are calculated 
        using the CA-CB coordinates.

        This list is used in the following way:

        A certain atom is found within radius of a given position.
        The coordinates are translated by -translation.
        The coordinates are rotated by rotation.
        If the z coordinate is >0 it lies in the side-chain sphere.
        Otherwise it lies in the backbone sphere.
        """
        # list of (translation, rotation, ca, residue) tuples
        data_list=[]
        for residue in residue_list:
            if is_aa(residue):
                ca=residue["CA"]
                ca_coord=ca.get_coord()
                if residue.has_id("CB"):
                    cb=residue["CB"]
                    cb_coord=cb.get_coord()
                else:
                    # GLY has no CB - calculate pseudo CB position
                    cb_coord=self._get_gly_cb_coord(residue)
                x,y,z=cb_coord-ca_coord
                # CB-CA vector
                cb_ca=Vector(x,y,z)
                # Rotate CB-CA vector to unit vector along Z
                rotation=rotmat(cb_ca, self.unit_z)
                translation=ca_coord
                data_list.append((translation, rotation, ca, residue))
        return data_list

    def _get_data_list_from_ca(self, structure):
        """
        Return a list of (translation, rotation, ca, residue)
        tuples (using CA). 
        
        The translation and rotation are calculated 
        using the CA-CA-CA coordinates (ie. the side chain 
        direction is approximated using the CA-CA-CA vectors).

        This list is used in the following way:

        A certain atom is found within radius of a given position.
        The coordinates are translated by -translation.
        The coordinates are rotated by rotation.
        If the z coordinate is >0 it lies in the side-chain sphere.
        Otherwise it lies in the backbone sphere.
        """
        # list of (translation, rotation, ca, residue) tuples
        data_list=[]
        ppb=PPBuilder()
        for pp in ppb.build_peptides(structure):
            ca_list=[]
            ca_list=pp.get_ca_list()
            for i in range(1, len(ca_list)-1):
                r=pp[i]
                ca1=ca_list[i-1]
                ca2=ca_list[i]
                ca3=ca_list[i+1]
                ca_v1=ca1.get_vector()
                ca_v2=ca2.get_vector()
                ca_v3=ca3.get_vector()
                # pseudo cb centred at origin
                cb_v=self._get_cb_from_ca(ca_v1, ca_v2, ca_v3)
                # rotate cb to unit vector along z
                rotation=rotmat(cb_v, self.unit_z)
                translation=ca2.get_coord()
                data_list.append((translation, rotation, ca2, r))
        return data_list

    def _calc_hs_exposure(self, data_list, radius):
        d={}
        for tran1, rot1, ca1, r1 in data_list:
            hs_sidechain=0
            hs_backbone=0
            for tran2, rot2, ca2, r2 in data_list:
                if r1 is r2:
                    continue
                if (ca1-ca2)<radius:
                    ca_coord2=ca2.get_coord()
                    rot_coord=matrixmultiply(rot1, ca_coord2-tran1)
                    if rot_coord[2]>0:
                        hs_sidechain+=1
                    else:
                        hs_backbone+=1
            d[r1]=(hs_sidechain, hs_backbone)
        return d
    
    def calc_ca_exposure(self, structure, radius=12.0):
        data_list=self._get_data_list_from_ca(structure)
        return self._calc_hs_exposure(data_list, radius)
                
    def calc_cb_exposure(self, structure, radius=12.0):
        residue_list=Selection.unfold_entities(structure, 'R')
        data_list=self._get_data_list_from_cb(residue_list)
        return self._calc_hs_exposure(data_list, radius)

def calc_simple_exposure(residue_list, radius):
    """
    A residue's exposure is defined as the number of CA atoms around 
    that residues CA atom. A dictionary is returned that uses a Residue
    object as key, and the residue exposure as corresponding value.
    """
    ca_list=[]
    # Extract the CA coordinates
    for r in residue_list:
        if is_aa(r) and r.has_id("CA"):
            ca=r["CA"]
            ca_list.append((ca, r))
    # Calculate exposure
    d={}
    for ca1, res1 in ca_list:
        for ca2, res2 in ca_list:
            if ca1 is ca2:
                continue
            if (ca1-ca2)<=radius:
                if d.has_key(res1):
                    d[res1]=d[res1]+1
                else:
                    d[res1]=1
    return d


if __name__=="__main__":

    import sys
    import os

    p=PDBParser(PERMISSIVE=1)
    s=p.get_structure('X', sys.argv[1])

    #residue_list=Selection.unfold_entities(s, 'R')
    #exp_simple=calc_simple_exposure(residue_list, 12.0)

    hse=HSExposure()
    exp=hse.calc_ca_exposure(s, 12.0)
    #exp=hse.calc_cb_exposure(s, 12.0)

    keys=exp.keys()
    keys.sort()

    os.system("touch exp_ca3.txt")

    fp=open("exp_ca3.txt", "a")
    for key in keys:
        if exp.has_key(key):
            name=key.get_resname()
            up, down=exp[key][0], exp[key][1]
            fp.write("%s %.2f %.2f\n" % (name, up, down))
    fp.close()

    print "%s done."  % sys.argv[1]
    sys.stdout.flush()

    #hse.write_pymol_script()



