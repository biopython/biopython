from Numeric import matrixmultiply, transpose, array, sqrt
from math import pi
from MLab import eye
import sys

from Bio.PDB import *

__doc__="Half sphere exposure and coordination number calculation."


class HSExposure:
    """
    Calculates the Half-Sphere Exposure (HSE).
    """
    # unit vector along z
    # CA-CB vectors or the CA-CA-CA approximation will
    # be rotated onto this unit vector
    _unit_z=Vector(0.0, 0.0, 1.0)

    def __init__(self, OFFSET=0.0):
        # List of CA-CB direction calculated from CA-CA-CA
        # Used for the PyMol script
        self.ca_cb_list=[]
        # Dummy 
        self.angles={}
        self.OFFSET=OFFSET

    def write_pymol_script(self, filename="hs_exp.py"):
        """
        Write a PyMol script that visualizes the pseudo CB-CA directions 
        at the CA coordinates.
        """
        if len(self.ca_cb_list)==0:
            sys.stderr.write("Nothing to draw.\n")
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

    def _get_gly_cb_vector(self, residue):
        """
        Return a pseudo CB coord for a Gly residue.
        
        The pseudoCB vector is centered at CA (ie. it's the
        atom position itself).

        CB coord=N coord rotated over -120 degrees 
        along the CA-C axis.
        """
        n_v=residue["N"].get_vector()
        c_v=residue["C"].get_vector()
        ca_v=residue["CA"].get_vector()
        # center at origin
        n_v=n_v-ca_v
        c_v=c_v-ca_v
        # rotation around c-ca over -120 deg
        rot=rotaxis(-pi*120.0/180.0, c_v)
        cb_at_origin_v=n_v.left_multiply(rot)
        # move back to ca position
        cb_v=cb_at_origin_v+ca_v
        # This is for PyMol visualization
        self.ca_cb_list.append((ca_v, cb_v))
        return cb_v

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
        # vector b is centered at the origin!
        return b

    def _get_rotran_list_from_cb(self, residue_list):
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
        rotran_list=[]
        for residue in residue_list:
            if is_aa(residue):
                ca=residue["CA"]
                ca_v=ca.get_vector()
                if residue.has_id("CB"):
                    cb=residue["CB"]
                    cb_v=cb.get_vector()
                    # Call _get_gly_cb_vector here for PyMol output
                    # self._get_gly_cb_vector(residue)
                else:
                    # GLY has no CB - calculate pseudo CB position
                    # based on N/CA/C positions
                    cb_v=self._get_gly_cb_vector(residue)
                # CB-CA vector
                cb_ca_v=cb_v-ca_v
                # Rotate CB-CA vector to unit vector along Z
                rotation=rotmat(cb_ca_v, self._unit_z)
                translation=ca.get_coord()
                rotran_list.append((translation, rotation, ca, residue))
        return rotran_list

    def _get_rotran_list_from_ca(self, model):
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
        rotran_list=[]
        ppb=PPBuilder()
        # angles between pseudo-CB-CA and CA-CB
        angles={}
        for pp in ppb.build_peptides(model):
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
                rotation=rotmat(cb_v, self._unit_z)
                translation=ca2.get_coord()
                rotran_list.append((translation, rotation, ca2, r))
                # Calculate angle between pseudo-CB-CA and CA-CB
                angles[r]=self._calc_delta_angle(r, cb_v)
        self.angles=angles
        return rotran_list

    def _calc_delta_angle(self, res, pseudo_cb_v):
        """
        Caculate delta angle between CB-CA and pseudoCB-CA
        """
        if res.get_resname()=="GLY":
            # Calculate pseudo CA for GLY
            cb_v=self._get_gly_cb_vector(res)
            ca_v=res["CA"].get_vector()
        else:
            # Calculate CB-CA vector
            cb_v=res["CB"].get_vector()
            ca_v=res["CA"].get_vector()
        real_cb_v=cb_v-ca_v
        angle=360*real_cb_v.angle(pseudo_cb_v)/(2*pi)
        return angle

    def _calc_hs_exposure(self, rotran_list, residue_list, radius):
        """
        Calculate for each residue how many CA atoms are present in 
        the half sphere in the side chain direction, and in the half
        sphere on the opposite side.'
        """
        d={}
        for tran, rot, ca1, r1 in rotran_list:
            hs_sidechain=0
            hs_mainchain=0
            for r2 in residue_list:
                if not is_aa(r2):
                    continue
                if r1 is r2:
                    continue
                ca2=r2["CA"]
                if (ca1-ca2)<radius:
                    neighbor_coord=ca2.get_coord()
                    # Rotate neighbor to the CB-CA direction
                    rot_coord=matrixmultiply(rot, neighbor_coord-tran)
                    if rot_coord[2]>self.OFFSET:
                        # in side chain half sphere
                        hs_sidechain+=1
                    else:
                        # in main chain half sphere
                        hs_mainchain+=1
            d[r1]=(hs_sidechain, hs_mainchain)
        return d
    
    def calc_hs_exposure(self, model, radius=12.0, option='CB'):
        residue_list=Selection.unfold_entities(model, 'R')
        if option=='CA3':
            rotran_list=self._get_rotran_list_from_ca(model)
        elif option=='CB':
            rotran_list=self._get_rotran_list_from_cb(residue_list)
        else:
            raise "Options: CA3 or CB"
        return self._calc_hs_exposure(rotran_list, residue_list, radius)

    def calc_fs_exposure(self, model, radius=12.0):
        """
        A residue's exposure is defined as the number of CA atoms around 
        that residues CA atom. A dictionary is returned that uses a Residue
        object as key, and the residue exposure as corresponding value.
        """
        residue_list=Selection.unfold_entities(model, 'R')
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

    def get_angles(self):
        """
        Return delta angle between CA-CB and pseudoCB-CA
        """
        return self.angles


if __name__=="__main__":

    import sys

    p=PDBParser()
    s=p.get_structure('X', sys.argv[1])

    model=s[0]

    # Neighbor sphere radius
    RADIUS=13.0

    hse=HSExposure()

    # Calculate CA-CA-CA based HS-exposure
    exp_ca=hse.calc_hs_exposure(model, RADIUS, option='CA3')

    # Calculate CB based HS-exposure
    exp_cb=hse.calc_hs_exposure(model, RADIUS, option='CB')

    # Calculate classical coordination number
    exp_fs=hse.calc_fs_exposure(model, RADIUS)

    angles=hse.get_angles()

    # All residues in the model
    residue_list=Selection.unfold_entities(model, 'R')

    # Print accessibilities for each residue
    for r in residue_list:
        print r

        # CA hs-exposure
        if exp_ca.has_key(r):
            print "CA3  ", exp_ca[r]

        # CB hs-exposure
        if exp_cb.has_key(r):
            print "CB   ", exp_cb[r]

        # Classical sphere coordination number
        if exp_fs.has_key(r):
            print "SPHE ", exp_fs[r]

        if angles.has_key(r):
            print "DELTA %.2f" % angles[r]

        print "--------------------"

    hse.write_pymol_script()




