from Numeric import matrixmultiply, transpose, array, sqrt
from math import pi
from MLab import eye

from Bio.PDB import *


class HSExposure:
    """
    Calculates the Half-Sphere Exposure (HSE).
    """
    # unit vector along z
    unit=Vector(0.0, 0.0, 1.0)

    def __init__(self):
        self.ca_cb_list=[]

    def write_pymol_script(self, filename="pymol.py"):
        if len(self.ca_cb_list)==0:
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
        "CB=N rotated over -120 degrees along CA-C axis"
        n=residue["N"].get_vector()
        c=residue["C"].get_vector()
        ca=residue["CA"].get_vector()
        n=n-ca
        c=c-ca
        rot=rotaxis(-pi*120.0/180.0, c)
        cb=n.left_multiply(rot)
        cb=cb+ca
        return cb.get_array()

    def _get_cb_from_ca(self, ca1, ca2, ca3):
        """
        Calculate the approximate CA-CB direction based
        on two flanking CA positions.
        """
        # center
        ca1=ca1-ca2
        ca3=ca3-ca2
        ca1=ca1.normalize()
        ca3=ca3.normalize()
        # bisection
        b=(ca1+ca3).normalize()
        b=-b
        # for drawing
        self.ca_cb_list.append((ca2, b+ca2))
        return b

    def _get_data_list_from_cb(self, residue_list):
        # Unit vector along Z
        data_list=[]
        for r in residue_list:
            if is_aa(r):
                ca=r["CA"]
                ca_coord=ca.get_coord()
                if r.has_id("CB"):
                    cb=r["CB"]
                    cb_coord=cb.get_coord()
                else:
                    # GLY
                    cb_coord=self._get_gly_cb_coord(r)
                x,y,z=cb_coord-ca_coord
                cb_ca=Vector(x,y,z)
                m=rotmat(cb_ca, self.unit)
                data_list.append((ca_coord, m, ca, r))
        return data_list

    def _get_data_list_from_ca(self, structure):
        ppb=PPBuilder()
        data_list=[]
        for pp in ppb.build_peptides(structure):
            ca_list=[]
            for r in pp:
                ca_list.append(r["CA"])
            for i in range(1, len(ca_list)-1):
                r=pp[i]
                ca1=ca_list[i-1]
                ca2=ca_list[i]
                ca3=ca_list[i+1]
                ca_v1=ca1.get_vector()
                ca_v2=ca2.get_vector()
                ca_v3=ca3.get_vector()
                cb_v=self._get_cb_from_ca(ca_v1, ca_v2, ca_v3)
                m=rotmat(cb_v, self.unit)
                ca_coord2=ca2.get_coord()
                data_list.append((ca_coord2, m, ca2, r))
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



