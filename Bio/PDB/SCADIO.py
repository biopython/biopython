# Copyright 2019 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""SCADIO: write OpenSCAD program to create protein structure."""
# import re

from Bio.File import as_handle
from Bio.PDB.PDBExceptions import PDBException

from Bio.PDB.ic_classes import IC_Residue, IC_Chain
from Bio.PDB.internal_coords import internal_to_atom_coordinates
from Bio.PDB.internal_coords import atom_to_internal_coordinates
from Bio.PDB.vectors import homog_scale_mtx


def write_SCAD(entity, file, scale=None, handle='protein', pdbid=None,
               backboneOnly=False, includeCode=True):
    """Write hedron assembly to file as OpenSCAD matrices.

    Output data format is primarily:
        matrix for each hedron: len1, angle2, len3, atom covalent bond class,
            flags toindicate atom/bond represented in previous hedron (OpenSCAD
            very slow for overlapping elements)
        matrices for each residue to assemble hedra into dihedrons
        matrices to transform each residue set of dihedra to position in chain

    :param entity: Biopython PDB structure entity
        structure data to export
    :param file: Bipoython as_handle filename or open file pointer
        file to write data to
    :param scale: float
        units (usually mm) per angstrom for STL output, written in output
    :param handle: str, default 'protein'
        name for top level of generated OpenSCAD matrix structure,
    :param pdbid: str
        PDB idcode, written in output. Defaults to '0PDB' if not supplied
        and no 'idcode' set in entity
    :param backboneOnly: bool default False
        Do not output side chain data if True
    :param includeCode: bool default True
        Include Bio/PDB/peptide.scad so output file can be loaded into
        OpenSCAD; if False, output data matrices only

    """
    # step one need IC_Residue atom_coords loaded in order to scale
    have_PIC_Atoms = False
    if 'S' == entity.level or 'M' == entity.level:
        for chn in entity.get_chains():
            if not hasattr(chn, 'internal_coord'):
                chn.internal_coord = IC_Chain(chn)
                have_PIC_Atoms = True
    elif 'C' == entity.level:
        if not hasattr(entity, 'internal_coord'):
            entity.internal_coord = IC_Chain(entity)
            have_PIC_Atoms = True
    elif 'R' == entity.level:
        if not hasattr(entity, 'internal_coord'):
            entity.internal_coord = IC_Residue(entity)
            have_PIC_Atoms = True
    else:
        raise PDBException("level not S, M. C or R: " + str(entity.level))

    if not have_PIC_Atoms and scale is not None:
        # if loaded pic file and need to scale, generate atom coords
        internal_to_atom_coordinates(entity)

    if scale is not None:
        scaleMtx = homog_scale_mtx(scale)
        for res in entity.get_residues():
            if hasattr(res, 'internal_coord'):
                res.internal_coord.applyMtx(scaleMtx)
                if res.internal_coord.gly_Cbeta:
                    res.internal_coord.scale = scale

    # generate internal coords for scaled entity
    # -- hedron bond lengths have changed
    # if not scaling, still need to generate internal coordinate
    # bonds for ring sidechains
    atom_to_internal_coordinates(entity, allBonds=True)

    # clear initNCaC - want at origin, not match PDB file
    for chn in entity.get_chains():
        chn.internal_coord.initNCaC = {}

    internal_to_atom_coordinates(entity)

    with as_handle(file, 'w') as fp:

        fp.write('protein_scale=' + str(scale) + ';\n')
        if includeCode:
            fp.write('$fn=20;\nchain(protein);\n')
            fp.write(peptide_scad)
            # codeFile = re.sub(r"pic.py\Z", "peptide.scad", __file__)
            # with as_handle(codeFile, 'r') as cf:
            #     for line in cf.readlines():
            #        fp.write(line)

        if not pdbid and hasattr(entity, 'header'):
            pdbid = entity.header.get('idcode', None)
        if pdbid is None or '' == pdbid:
            pdbid = '0PDB'
        fp.write('protein = [ "' + pdbid + '", protein_scale,\n')

        if 'S' == entity.level or 'M' == entity.level:
            for chn in entity.get_chains():
                fp.write(' [\n')
                chn.internal_coord.write_SCAD(fp, scale, backboneOnly)
                fp.write(' ]\n')
        elif 'C' == entity.level:
            fp.write(' [\n')
            entity.internal_coord.write_SCAD(fp, scale, backboneOnly)
            fp.write(' ]\n')
        elif 'R' == entity.level:
            raise NotImplementedError(
                'writescad single residue not yet implemented.')

        fp.write('\n];\n')


peptide_scad = """
//
// peptide.scad
// Copyright(c) 2019 Robert T. Miller.  All rights reserved.
// This file is part of the Biopython distribution and governed by your
// choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
// Please see the LICENSE file that should have been included as part of this
// package.

// This is the support file to build an OpenSCAD (http://www.openscad.org/)
// model of a protein from internal coordinates.
//
// Lines similar to:
//    protein_scale=2;
//    $fn=20;
//    chain(protein);
//
// should be pre-pended to the beginning of this file, and data matrices
// should be appended below to form a program ready to load into the OpenSCAD
// application.
//
//  "protein_scale" is the value supplied when generating the data for build
//  units per PDB angstrom.  You may wish to modify it here to adjust the
//  appearance of the model in terms of atom sphere or bond cylinder diameter.
//
//  "$fn" (fragment number) is an OpenSCAD parameter controlling the smoothness
//  of the model surface.  Smaller values will render faster, but yield more
//  'blocky' models.
//
//  "chain(protein);" above is the command to call the chain() function below
//  to render the data in the "protein" array at the end of this file.
//
//  This is intended to be a working example, you are encouraged to modify the
//  OpenSCAD subroutines below to generate a model to your liking.  For more
//  information, start with http://www.openscad.org/cheatsheet/index.html
//

// output parameters
bondRadius=0.4;
clearance=0.25;
extrusionWidth=0.4;
minThickness=2*extrusionWidth;
atomScale=0.8;

//
// Generate a sphere to represent an atom.
// Colour and size determined for the atom covalent radius specified by the
//  parameter 'a' by lookup in the atomData table below, then scaled by the
//  supplied parameter 'scal'.
module atom(a,scal)
{
    ad = atomData[search([a],atomData)[0]];
    color(ad[1]) {
        sphere(r=((ad[2]*atomScale)*scal));
    }
}


// Generate a 'hedron', one plane of 3 points, consisting of 3 atoms joined by
//  two bonds.
// In some cases the sequence of atoms in the h[] array is reversed, as
//  detailed in the comments.
module hedron(h,rev=0,scal)
{
    newh = (rev ?
        // reversed
        //    0     1     2     3     4     5     6      7
            [ h[2], h[0], h[5], h[3], h[8], h[6], h[10], h[9] ]
            :
        // not reversed
        //    0     1     2     3     4     5      6     7
            [ h[0], h[2], h[3], h[5], h[6], h[8],  h[9], h[10] ]
            );

    scal=scal;

    if (h[7]) {
        // central atom at 0,0,0
        atom(h[4],scal);
    }

    if (newh[5]) {
        // comments for non-reversed case
        // atom 3 is len3 up on +z
        translate([0,0,newh[1]]) atom(newh[3],scal);
    }

    if (newh[7]) {
        // atom 2 - atom 3 bond from origin up +z distance len3
        cylinder(h=newh[1],r=bondRadius*scal,center=false);
    }

    // atom 1 - atom 2 bond rotated about Y
    rotate([0, h[1], 0]) {   // rotate following elements by angle2 about Y

        if (newh[6]) {
            // a1-a2 bond is len1
            cylinder(h=newh[0],r=bondRadius*scal,center=false);
        }

        if (newh[4]) {
            // put atom1 sphere len1 away on Z
            translate([0,0,newh[0]]) atom(newh[2],scal);
        }
    }
}


// Generate a hedron rotated to specific angle d
module d2(d,hedra,scal)
{
    // get h2 len1 depending on reversed state
    tz = (d[d_reversed] ? hedra[d[d_h2ndx]][2] : hedra[d[d_h2ndx]][0]);

    // 4. rotate h2 to specified dihedral angle
    rotate(d[d_dangle1]) {
        // 3. translate h2 h2:len1 up +z
        translate([0,0,tz]) {
            // 2. rotate h2r about X so h2:a3 in +z and h2:a1 in -z
            rotate([180, 0, 0]) {
                // 1. reverse hedron 2 orientation = h2r
                hedron(hedra[d[d_h2ndx]],(!d[d_reversed]),scal);
            }
        }
    }
}

// Generate two hedra at specified dihedral angle d
module dihedron(d,hedra,scal)
{
    // reverse h1 if dihedral reversed
    hedron(hedra[d[d_h1ndx]],d[d_reversed],scal);
    d2(d,hedra,scal);
}

// Generate a residue consisting of the set of dihedra in the parameter 'r'.
//   r references hedra in the table speicified in the parameter 'hedra'.
module residue(r,hedra, scal)
{
    for (d = r) {
        multmatrix(d[d_dihedralTransform]) {
            dihedron(d, hedra, scal);
        }
    }
}

// Generate a chain of residues, each positioned by a supplied
//  rotation/translation matrix.
module chain(protein)
{
    chnD = protein[p_chainData];
    c = chnD[c_residues];
    dihedra = chnD[c_dihedra];
    hedra = chnD[c_hedra];
    for (r = c) {
        multmatrix(r[r_resTransform]) {
            residue(dihedra[r[r_resNdx]],hedra, protein[p_proteinScale]);
        }
    }
}

//
// OpenSCAD array indices (offsets) to reference protein data
//

// protein base level
p_pdbid = 0;
p_proteinScale = 1;
p_chainData = 2;

// chain level data
c_chainID = 0;
c_dihedra = 1;
c_hedra = 2;
c_residues = 3;

// hedra definitions
h_len1 = 0;
h_angle2 = 1;
h_len3 = 2;
h_atom1 = 3;
h_atom2 = 4;
h_atom3 = 5;

// dihedra specifications for each residue in sequence, dihedral array
d_dangle1 = 0;
d_h1ndx = 1;
d_h2ndx = 2;
d_reversed = 3;
d_dihedralTransform = 4;

// residueSet: world transform for each residue in sequence array
r_resNdx = 0;
r_resID = 1;
r_resTransform = 2;

// covalent radii from Heyrovska, Raji : 'Atomic Structures of all the Twenty
//   Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of
//   Atomic Covalent Radii'
// https://arxiv.org/pdf/0804.2488.pdf

atomData = [
    ["Csb","green" , 0.77],
    ["Cres","green" , 0.72],
    ["Cdb","green" , 0.67],
    ["Osb","red" , 0.67],
    ["Ores","red" , 0.635],
    ["Odb","red" , 0.60],
    ["Nsb","blue" , 0.70],
    ["Nres","blue" , 0.66],
    ["Ndb","blue" , 0.62],
    ["Hsb","gray" , 0.37],
    ["Ssb","yellow" , 1.04]
];

// Protein specific array data below.
"""
