# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Half-sphere exposure and coordination number calculation."""

import warnings
from math import pi

from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import CaPPBuilder, is_aa
from Bio.PDB.Vector import rotaxis


class _AbstractHSExposure(AbstractPropertyMap):
    """
    Abstract class to calculate Half-Sphere Exposure (HSE).

    The HSE can be calculated based on the CA-CB vector, or the pseudo CB-CA
    vector based on three consecutive CA atoms. This is done by two separate 
    subclasses. 
    """
    def __init__(self, model, radius, offset, hse_up_key, hse_down_key, 
            angle_key=None):
        """
        @param model: model
        @type model: L{Model}

        @param radius: HSE radius
        @type radius: float

        @param offset: number of flanking residues that are ignored in the calculation
        of the number of neighbors
        @type offset: int

        @param hse_up_key: key used to store HSEup in the entity.xtra attribute
        @type hse_up_key: string

        @param hse_down_key: key used to store HSEdown in the entity.xtra attribute
        @type hse_down_key: string

        @param angle_key: key used to store the angle between CA-CB and CA-pCB in 
            the entity.xtra attribute
        @type angle_key: string
        """
        assert(offset>=0)
        # For PyMOL visualization
        self.ca_cb_list=[]
        ppb=CaPPBuilder()
        ppl=ppb.build_peptides(model)
        hse_map={}
        hse_list=[]
        hse_keys=[]
        for pp1 in ppl:
            for i in range(0, len(pp1)):
                if i==0:
                    r1=None
                else:
                    r1=pp1[i-1]
                r2=pp1[i]
                if i==len(pp1)-1:
                    r3=None
                else:
                    r3=pp1[i+1]
                # This method is provided by the subclasses to calculate HSE
                result=self._get_cb(r1, r2, r3)
                if result is None:
                    # Missing atoms, or i==0, or i==len(pp1)-1
                    continue
                pcb, angle=result
                hse_u=0
                hse_d=0
                ca2=r2['CA'].get_vector()
                for pp2 in ppl:
                    for j in range(0, len(pp2)):
                        if pp1 is pp2 and abs(i-j)<=offset:
                            # neighboring residues in the chain are ignored 
                            continue
                        ro=pp2[j]
                        if not is_aa(ro) or not ro.has_id('CA'):
                            continue
                        cao=ro['CA'].get_vector()
                        d=(cao-ca2)
                        if d.norm()<radius:
                            if d.angle(pcb)<(pi/2):
                                hse_u+=1
                            else:
                                hse_d+=1
                res_id=r2.get_id()
                chain_id=r2.get_parent().get_id()
                # Fill the 3 data structures
                hse_map[(chain_id, res_id)]=(hse_u, hse_d, angle)
                hse_list.append((r2, (hse_u, hse_d, angle)))
                hse_keys.append((chain_id, res_id))
                # Add to xtra
                r2.xtra[hse_up_key]=hse_u
                r2.xtra[hse_down_key]=hse_d
                if angle_key:
                    r2.xtra[angle_key]=angle
        AbstractPropertyMap.__init__(self, hse_map, hse_keys, hse_list)

    def _get_cb(self, r1, r2, r3):
        """This method is provided by the subclasses to calculate HSE."""
        return NotImplemented

    def _get_gly_cb_vector(self, residue):
        """
        Return a pseudo CB vector for a Gly residue.
        The pseudoCB vector is centered at the origin.

        CB coord=N coord rotated over -120 degrees 
        along the CA-C axis.
        """
        try:
            n_v=residue["N"].get_vector()
            c_v=residue["C"].get_vector()
            ca_v=residue["CA"].get_vector()
        except:
            return None
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
        return cb_at_origin_v



class HSExposureCA(_AbstractHSExposure):
    """
    Class to calculate HSE based on the approximate CA-CB vectors,
    using three consecutive CA positions.
    """
    def __init__(self, model, radius=12, offset=0):
        """
        @param model: the model that contains the residues
        @type model: L{Model}

        @param radius: radius of the sphere (centred at the CA atom)
        @type radius: float

        @param offset: number of flanking residues that are ignored in the calculation            of the number of neighbors
        @type offset: int
        """
        _AbstractHSExposure.__init__(self, model, radius, offset, 
                'EXP_HSE_A_U', 'EXP_HSE_A_D', 'EXP_CB_PCB_ANGLE')

    def _get_cb(self, r1, r2, r3):
        """
        Calculate the approximate CA-CB direction for a central
        CA atom based on the two flanking CA positions, and the angle
        with the real CA-CB vector. 
        
        The CA-CB vector is centered at the origin.

        @param r1, r2, r3: three consecutive residues
        @type r1, r2, r3: L{Residue}
        """
        if r1 is None or r3 is None:
            return None
        try:
            ca1=r1['CA'].get_vector()
            ca2=r2['CA'].get_vector()
            ca3=r3['CA'].get_vector()
        except:
            return None
        # center
        d1=ca2-ca1
        d3=ca2-ca3
        d1.normalize()
        d3.normalize()
        # bisection
        b=(d1+d3)
        b.normalize()
        # Add to ca_cb_list for drawing
        self.ca_cb_list.append((ca2, b+ca2))
        if r2.has_id('CB'):
            cb=r2['CB'].get_vector()
            cb_ca=cb-ca2
            cb_ca.normalize()
            angle=cb_ca.angle(b)
        elif r2.get_resname()=='GLY':
            cb_ca=self._get_gly_cb_vector(r2)
            if cb_ca is None:
                angle=None
            else:
                angle=cb_ca.angle(b)
        else:
            angle=None
        # vector b is centered at the origin!
        return b, angle

    def pcb_vectors_pymol(self, filename="hs_exp.py"):
        """
        Write a PyMol script that visualizes the pseudo CB-CA directions 
        at the CA coordinates.

        @param filename: the name of the pymol script file
        @type filename: string
        """
        if len(self.ca_cb_list)==0:
            warnings.warn("Nothing to draw.", RuntimeWarning)
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


class HSExposureCB(_AbstractHSExposure):
    """
    Class to calculate HSE based on the real CA-CB vectors.
    """
    def __init__(self, model, radius=12, offset=0):
        """
        @param model: the model that contains the residues
        @type model: L{Model}

        @param radius: radius of the sphere (centred at the CA atom)
        @type radius: float

        @param offset: number of flanking residues that are ignored in the calculation            of the number of neighbors
        @type offset: int
        """
        _AbstractHSExposure.__init__(self, model, radius, offset,
                'EXP_HSE_B_U', 'EXP_HSE_B_D')

    def _get_cb(self, r1, r2, r3):
        """
        Method to calculate CB-CA vector.

        @param r1, r2, r3: three consecutive residues (only r2 is used)
        @type r1, r2, r3: L{Residue}
        """
        if r2.get_resname()=='GLY':
            return self._get_gly_cb_vector(r2), 0.0
        else:
            if r2.has_id('CB') and r2.has_id('CA'):
                vcb=r2['CB'].get_vector()
                vca=r2['CA'].get_vector()
                return (vcb-vca), 0.0
        return None


class ExposureCN(AbstractPropertyMap):
    def __init__(self, model, radius=12.0, offset=0):
        """
        A residue's exposure is defined as the number of CA atoms around 
        that residues CA atom. A dictionary is returned that uses a L{Residue}
        object as key, and the residue exposure as corresponding value.

        @param model: the model that contains the residues
        @type model: L{Model}

        @param radius: radius of the sphere (centred at the CA atom)
        @type radius: float

        @param offset: number of flanking residues that are ignored in the calculation            of the number of neighbors
        @type offset: int

        """
        assert(offset>=0)
        ppb=CaPPBuilder()
        ppl=ppb.build_peptides(model)
        fs_map={}
        fs_list=[]
        fs_keys=[]
        for pp1 in ppl:
            for i in range(0, len(pp1)):
                fs=0
                r1=pp1[i]
                if not is_aa(r1) or not r1.has_id('CA'):
                    continue
                ca1=r1['CA']
                for pp2 in ppl:
                    for j in range(0, len(pp2)):
                        if pp1 is pp2 and abs(i-j)<=offset:
                            continue
                        r2=pp2[j]
                        if not is_aa(r2) or not r2.has_id('CA'):
                            continue
                        ca2=r2['CA']
                        d=(ca2-ca1)
                        if d<radius:
                            fs+=1
                res_id=r1.get_id()
                chain_id=r1.get_parent().get_id()
                # Fill the 3 data structures
                fs_map[(chain_id, res_id)]=fs
                fs_list.append((r1, fs))
                fs_keys.append((chain_id, res_id))
                # Add to xtra
                r1.xtra['EXP_CN']=fs
        AbstractPropertyMap.__init__(self, fs_map, fs_keys, fs_list)


if __name__=="__main__":

    import sys

    p=PDBParser()
    s=p.get_structure('X', sys.argv[1])
    model=s[0]

    # Neighbor sphere radius
    RADIUS=13.0
    OFFSET=0

    hse=HSExposureCA(model, radius=RADIUS, offset=OFFSET)
    for l in hse:
        print l
    print

    hse=HSExposureCB(model, radius=RADIUS, offset=OFFSET)
    for l in hse:
        print l
    print

    hse=ExposureCN(model, radius=RADIUS, offset=OFFSET)
    for l in hse:
        print l
    print

    for c in model:
        for r in c:
            try:
                print r.xtra['PCB_CB_ANGLE']
            except:
                pass


