# Copyright 2019-21 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""SCADIO: write OpenSCAD program to create protein structure 3D model.

3D printing a protein structure is a non-trivial exercise due to the
overall complexity and the general requirement for supporting overhang regions
while printing.  This software is a path to generating a model for printing
(e.g. an STL file), and does not address the issues around converting the
model to a physical product.  OpenSCAD <http://www.openscad.org/> can create
a printable model from the script this software produces.  MeshMixer
<http://www.meshmixer.com/>, various slicer software, and the 3D printer
technology available to you provide options for addressing the problems around
physically rendering the model.

The model generated here consists of OpenSCAD primitives, e.g. spheres and
cylinders, representing individual atoms and bonds in an explicit model of a
protein structure.  The benefit is that individual atoms/bonds may be selected
for specific print customizations relevant to 3D printing (such as rotatable
bond mechanisms or hydrogen bond magnets).  Alternatively, use e.g. Chimera to
render a structure as ribbons or similar for printing as a single object.

I suggest generating your initial model using the OpenSCAD script provided
here, then modifying that script according to your needs.  Changing the
atomScale and bondRadius values can simplify the model by removing gaps and
the corresponding need for supports, or you may wish to modify the
hedronDispatch() routine to select residues or chain sections for printing
separately and subsequently joining with rotatable bonds.  During this
development phase you will likely have your version include only the data
matrices generated here, by using the `includeCode=False` option to
write_SCAD().  An example project using rotatable backbone and magnetic
hydrogen bonds is at <https://www.thingiverse.com/thing:3957471>.
"""
# import re

from Bio.File import as_handle
from Bio.PDB.PDBExceptions import PDBException

from Bio.PDB.internal_coords import IC_Residue, IC_Chain

# from Bio.PDB.Structure import Structure
# from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import homog_scale_mtx

import numpy as np  # type: ignore


def _scale_residue(res, scale, scaleMtx):
    if res.internal_coord:
        res.internal_coord.applyMtx(scaleMtx)
        if res.internal_coord.gly_Cbeta:
            res.internal_coord.scale = scale


def write_SCAD(
    entity,
    file,
    scale=None,
    pdbid=None,
    backboneOnly=False,
    includeCode=True,
    maxPeptideBond=None,
    start=None,
    fin=None,
    handle="protein",
):
    """Write hedron assembly to file as OpenSCAD matrices.

    This routine calls both :meth:`.IC_Chain.internal_to_atom_coordinates` and
    :meth:`.IC_Chain.atom_to_internal_coordinates` due to requirements for
    scaling, explicit bonds around rings, and setting the coordinate space of
    the output model.

    Output data format is primarily:

    - matrix for each hedron:
        len1, angle2, len3, atom covalent bond class, flags to indicate
        atom/bond represented in previous hedron (OpenSCAD very slow with
        redundant overlapping elements), flags for bond features
    - transform matrices to assemble each hedron into residue dihedra sets
    - transform matrices for each residue to position in chain

    OpenSCAD software is included in this Python file to process these
    matrices into a model suitable for a 3D printing project.

    :param entity: Biopython PDB :class:`.Structure` entity
        structure data to export
    :param file: Bipoython :func:`.as_handle` filename or open file pointer
        file to write data to
    :param float scale:
        units (usually mm) per angstrom for STL output, written in output
    :param str pdbid:
        PDB idcode, written in output. Defaults to '0PDB' if not supplied
        and no 'idcode' set in entity
    :param bool backboneOnly: default False.
        Do not output side chain data past Cbeta if True
    :param bool includeCode: default True.
        Include OpenSCAD software (inline below) so output file can be loaded
        into OpenSCAD; if False, output data matrices only
    :param float maxPeptideBond: Optional default None.
        Override the cut-off in IC_Chain class (default 1.4) for detecting
        chain breaks.  If your target has chain breaks, pass a large number
        here to create a very long 'bond' spanning the break.
    :param int start,fin: default None
        Parameters for internal_to_atom_coords() to limit chain segment.
    :param str handle: default 'protein'
        name for top level of generated OpenSCAD matrix structure

    See :meth:`.IC_Residue.set_flexible` to set flags for specific residues to
    have rotatable bonds, and :meth:`.IC_Residue.set_hbond` to include cavities
    for small magnets to work as hydrogen bonds.
    See <https://www.thingiverse.com/thing:3957471> for implementation example.

    The OpenSCAD code explicitly creates spheres and cylinders to
    represent atoms and bonds in a 3D model.  Options are available
    to support rotatable bonds and magnetic hydrogen bonds.

    Matrices are written to link, enumerate and describe residues,
    dihedra, hedra, and chains, mirroring contents of the relevant IC_*
    data structures.

    The OpenSCAD matrix of hedra has additional information as follows:

    * the atom and bond state (single, double, resonance) are logged
        so that covalent radii may be used for atom spheres in the 3D models

    * bonds and atoms are tracked so that each is only created once

    * bond options for rotation and magnet holders for hydrogen bonds
        may be specified (see :meth:`.IC_Residue.set_flexible` and
        :meth:`.IC_Residue.set_hbond` )

    Note the application of :data:`Bio.PDB.internal_coords.IC_Chain.MaxPeptideBond`
    :  missing residues may be linked (joining chain segments with arbitrarily
    long bonds) by setting this to a large value.

    Note this uses the serial assembly per residue, placing each residue at
    the origin and supplying the coordinate space transform to OpenaSCAD

    All ALTLOC (disordered) residues and atoms are written to the output
    model.  (see :data:`Bio.PDB.internal_coords.IC_Residue.no_altloc`)
    """
    if maxPeptideBond is not None:
        mpbStash = IC_Chain.MaxPeptideBond
        IC_Chain.MaxPeptideBond = float(maxPeptideBond)

    # step one need IC_Residue atom_coords loaded in order to scale
    # so if no internal_coords, initialise from Atom coordinates
    added_IC_Atoms = False
    if "S" == entity.level or "M" == entity.level:
        for chn in entity.get_chains():
            if not chn.internal_coord:
                chn.internal_coord = IC_Chain(chn)
                added_IC_Atoms = True
    elif "C" == entity.level:
        if not entity.internal_coord:  # entity.internal_coord:
            entity.internal_coord = IC_Chain(entity)
            added_IC_Atoms = True
    else:
        raise PDBException("level not S, M or C: " + str(entity.level))

    if added_IC_Atoms:
        # if loaded pdb, need to scale, and asm, gen atomArray
        entity.atom_to_internal_coordinates()
    else:
        # if loaded pic file and need to scale, generate atom coords
        entity.internal_to_atom_coordinates(None)

    if scale is not None:
        scaleMtx = homog_scale_mtx(scale)

        if "C" == entity.level:
            entity.internal_coord.atomArray = np.dot(
                entity.internal_coord.atomArray[:], scaleMtx
            )
            entity.internal_coord.hAtoms_needs_update[:] = True
            entity.internal_coord.scale = scale
        else:
            for chn in entity.get_chains():
                if hasattr(chn.internal_coord, "atomArray"):
                    chn.internal_coord.atomArray = np.dot(
                        chn.internal_coord.atomArray[:], scaleMtx
                    )
                    chn.internal_coord.hAtoms_needs_update[:] = True
                    chn.internal_coord.scale = scale

    # generate internal coords for scaled entity
    # (hedron bond lengths have changed if scaled)
    # if not scaling, still need to generate internal coordinate
    # bonds for ring sidechains
    # AllBonds is a class attribute for IC_Residue.atom_to_internal_coordinates
    # to generate explicit hedra covering all bonds

    allBondsStash = IC_Residue._AllBonds
    IC_Residue._AllBonds = True
    # trigger rebuild of hedra for AllBonds
    if "C" == entity.level:
        entity.internal_coord.ordered_aa_ic_list[0].hedra = {}
        delattr(entity.internal_coord, "hAtoms_needs_update")
        delattr(entity.internal_coord, "hedraLen")
    else:
        for chn in entity.get_chains():
            chn.internal_coord.ordered_aa_ic_list[0].hedra = {}
            delattr(chn.internal_coord, "hAtoms_needs_update")
            delattr(chn.internal_coord, "hedraLen")
    entity.atom_to_internal_coordinates()
    IC_Residue._AllBonds = allBondsStash

    # rebuild atom coordinates now with chain starting at origin: in OpenSCAD
    # code, each residue model is transformed to N-Ca-C start position instead
    # of updating transform matrix along chain
    entity.internal_to_atom_coordinates()

    with as_handle(file, "w") as fp:
        if includeCode:
            fp.write(peptide_scad)

        if not pdbid and hasattr(entity, "header"):
            pdbid = entity.header.get("idcode", None)
        if pdbid is None or "" == pdbid:
            pdbid = "0PDB"
        fp.write(
            'protein = [ "' + pdbid + '", ' + str(scale) + ",  // ID, protein_scale\n"
        )

        if "S" == entity.level or "M" == entity.level:
            for chn in entity.get_chains():
                fp.write(" [\n")
                chn.internal_coord._write_SCAD(
                    fp, backboneOnly=backboneOnly, start=start, fin=fin
                )
                fp.write(" ]\n")
        elif "C" == entity.level:
            fp.write(" [\n")
            entity.internal_coord._write_SCAD(
                fp, backboneOnly=backboneOnly, start=start, fin=fin
            )
            fp.write(" ]\n")
        elif "R" == entity.level:
            raise NotImplementedError("writescad single residue not yet implemented.")

        fp.write("\n];\n")

    if maxPeptideBond is not None:
        IC_Chain.MaxPeptideBond = mpbStash


peptide_scad = """
/*
//
// peptide.scad
// Copyright (c) 2019 Robert T. Miller.  All rights reserved.
// This file is part of the Biopython distribution and governed by your
// choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
// Please see the LICENSE file that should have been included as part of this
// package.
//
// This is the support file to build an OpenSCAD (http://www.openscad.org/) model
// of a protein from internal coordinates.  The resulting model may be constructed
// on a 3D printer.
//
// data matrices should be appended below to form a program ready to load into
// the OpenSCAD application.
//
//  The protein_scale value used throughout is the second element of the
//    protein[] array appended below.
//    This is the value supplied when generating the data for build units per
//    PDB angstrom.
//    You may wish to modify it here to adjust the appearance of the model in
//    terms of atom sphere or bond cylinder diameter, however the bond lengths
//    are fixed with the supplied value when the data matrices are generated.
//    Atom sphere and bond cylinder radii may be individually adjusted below as
//    well.
//
//  $fn (fragment number) is an OpenSCAD parameter controlling the smoothness
//    of the model surface.  Smaller values will render faster, but yield more
//    'blocky' models.
//
//  This is intended to be a working example, you are encouraged to modify the
//    OpenSCAD subroutines below to generate a model to your liking.  For more
//    information, start with http://www.openscad.org/cheatsheet/index.html
//
//  Note especially the hedronDispatch() subroutine below: here you may select
//    hedra based on residue, sequence position, and class (hedron atoms) for
//    special handling.  Also see the per hedron render options in the hedra[]
//    array.
//
//  If you modify this file, you may find it useful to generate the data
//    matrices without this OpenSCAD code by calling write_SCAD() with the
//    includeCode=False option, then use the OpenSCAD 'include <>' facility at
//    the end of your modified OpenSCAD program.
*/

rotate([-90,0,0])  // convenient for default location (no N-Ca-C start coordinates)
    chain(protein);   // this is the main subroutine call to  build the structure

// top-level OpenSCAD $fn for visible surfaces.  Rotatable bonds use $fn=8
// inside, regardless of this setting.
$fn = 0;  // 0 yields OpenSCAD default of 30.  $n=8 should print with minimal support

tubes=false;     // style: render atoms and bonds as constant diameter cylinders, preferred for rotatable bonds / h-bonds
support=false;   // enable print-in-place internal support for rotatable bonds
// N.B. rotatable bonds must be parallel to build plate for internal support
// structures to be generated correctly by slicer

// output parameters
atomScale=1.0;  // 0.8 better for rotatable bonds
defaultAtomRadius = 0.77;  // used if tubes = true

bondRadius = (tubes ? defaultAtomRadius * atomScale : 0.4);
jBondRadius = defaultAtomRadius * atomScale;  // radius for rotatable bonds

// general printer, slicer, print settings
layerHeight=0.15;  // must match slicer setting for print-in-place support
clearance=0.3;     // sliding clearance - can be smaller (0.2) if not doing print-in-place
pClearance=0.2;    // press-fit clearance (magnets for h-bonds)
shim=0.05;         // extra to make OpenSCAD surfaces distinct in difference()
nozzleDiameter=0.4;

// need one magnet for each side of hydrogen bond, suggest 3mm x 5mm e.g. from eBay
// use compass to identify poles if you care, North pointing (red) repels compass North pointing
magR=3/2;    // magnet radius
magL=5;      // magnet length

// for $fn=8 which works nice on fdm printer
oRot = 22.5;              // 45/2, rotate to make fn=8 spheres and cylinders flat on build plate
apmFac = cos(180/8);      // apothem factor - multiply by radius for center to octagon side distance
octSide = 2* tan(180/8);  // multiply by radius to get length of octagon side
// for values of $fn:
fnRot = ($fn ? 90-(180/$fn) : 90-(180/30));

bondLenFac = 0.6;         // fraction of bond length to extend from atom for each arm of hedron in join

hblen = 1.97;             // hydrogen bond length

wall = 3*nozzleDiameter;
joinerStep = 1;           // radius difference between rotatable bond axle and end knob inside bond cylinder

caTop = false;     // only make top of N_C-alpha_C hedron plus C-beta (see hedron() and hedron_dispatch() examples)

/*
//
// Generate a sphere to represent an atom.
// Colour and size determined for the atom covalent radius specified by the
//   parameter 'a' by lookup in the atomData table below, then scaled by the
//   supplied parameter 'scal'.
//
// scal : protein_scale
// clr : additional radius if used to create clearance for rotatable bonds
//
*/
module atom(a,scal,clr=0)
{
    ad = atomData[search([a],atomData)[0]];
    color(ad[1]) {
        rotate([0,0,fnRot]) sphere(r=((ad[2]*atomScale)*scal)+clr);
    }
}

/*
//
// a hedron (below) may be 'reversed' in terms of the order of its two bonds;
// this function fixes the ordering
//
*/
function hFlip(h,rev) =
        //   yes reversed                                     :  not reversed
        //    0    1     2     3     4     5     6      7     :     0     1     2     3    4     5      6      7
        //  len1  len3  atom1 atom3  a1    a2   a1-a2  a2-a3      len1  len3  atom1 atom3   a1    a3  a1-a2  a2-a3
    (rev ? [ h[2], h[0], h[5], h[3], h[8], h[6], h[10], h[9] ] : [ h[0], h[2], h[3], h[5], h[6], h[8],  h[9], h[10] ]);
    // h[1] = angle2 for both cases


/*
//
// generate the male or female interior cylinders of a rotating bond
//
*/
module joinUnit(cOuterLen, cOuterRad, cInnerLen, cInnerRad, male=false) {
    if (male) {
        rotate([0,0,oRot]) {
            cylinder(h=cInnerLen, r=cInnerRad, center=false, $fn=8);
            cylinder(h=cOuterLen, r=cOuterRad, center=false, $fn=8);
        }
    } else {
        rotate([0,0,fnRot]) {
            cylinder(h=cInnerLen, r=cInnerRad, center=false, $fn=30);
            cylinder(h=cOuterLen, r=cOuterRad, center=false, $fn=30);
        }
    }
}

/*
//
// create a rotatable bond
//
// supportSel : 0 for no support, 1 or 2 for support on top or bottom (needed
// for reversed hedra)
//
*/
module joiner(bondlen, scal, male=0, ver=0, supportSelect=0) {  // ver = differentiate joiner part lengths to guide assembly, but not used
    lenfac = bondLenFac;
    jClr = clearance+0.05;

    cOuterRad = (jBondRadius * scal) - (2*wall + (male ? jClr/2 : -jClr/2));
    cInnerRad = cOuterRad - joinerStep;  // m/f jClr already in cOuterRad;  - (male ? 0 : -0*jClr/2);

    hArmLen = (bondlen * lenfac);
    lenClr = 0.6*jClr;  // length clearance applied to male and female both, so effective clearance is 2x this value
    cOuterLen = hArmLen * lenfac + (ver ? 0.5 : - 0.5) - (wall+ (male ? lenClr*2 : -lenClr*2  ));

    joinerOffset = (hArmLen * (1 - lenfac)) + (male ? lenClr : -lenClr) - (ver ? 1 : 0);

    i=supportSelect-1;
    oside = cOuterRad*octSide;
    wid = oside+2*wall+4*jClr+1;

    if (male) {
        rotate([0,180,0])
        translate([0,0,-(bondlen-joinerOffset)]) {
            difference() {
                joinUnit(cOuterLen, cOuterRad, bondlen, cInnerRad, male=true);
                if (supportSelect) {
                    rotate([0,0,i*180]) {
                        translate([0,(cOuterRad*apmFac)-0.5*layerHeight,cOuterLen/2]) {
                                cube([oside+2*shim,layerHeight+shim,cOuterLen+2*shim],center=true);
                        }
                    }
                }
            }
            if (supportSelect) {
                rotate([0,0,i*180]) {
                    translate([0,(cOuterRad*apmFac)-0.5*layerHeight,cOuterLen/2]) {
                        for (j=[0:1]) {
                            rotate([0,(j?60:-60),0])
                                cube([wid,layerHeight,2*nozzleDiameter],center=true);
                        }
                    }
                }
            }
        }
    } else {
        translate([0,0,joinerOffset]) {
            joinUnit(cOuterLen, cOuterRad, bondlen, cInnerRad);
            if (supportSelect) {  // extra gap top and bottom because filament sags
                supHeight = max(5*layerHeight,2*(cOuterRad-cOuterRad*apmFac));  // double because center=true below
                for(j=[0:1]) {
                    rotate([0,0,j*180]) {
                        translate([0,(cOuterRad*apmFac),cOuterLen/2]) {
                            cube([oside+2*shim,supHeight+shim,cOuterLen+2*shim],center=true);
                        }
                    }
                }
            }
        }
    }
}


/*
//
// create bond with different options (regular, skinny, h-bond atom, rotatable
// male or female
//
//  parameters:
//  bl : bond length
//  br : bond radius
//  scal : protein_scale
//  key : option symbols defined below
//  atm : atomic element symbol, used for color and radius by atom() routine above
//  ver : make rotatable bonds slightly different based on value; currently unused
//  supporSel : enable print-in-place support for rotatable bonds
//
*/

// option symbols - these names generated in BioPython code so avoid changing without thought
StdBond = 1;
FemaleJoinBond = 2;
MaleJoinBond = 3;
SkinnyBond = 4;        // Calpha - Cbeta bond cylinder needs to be skinny for clearance with rotating bonds
HBond = 5;             // make room inside atom/bond to insert magnet to appropriate depth

module bond(bl, br, scal, key, atm, ver, supportSel=0) {

    br = (key == FemaleJoinBond ? jBondRadius * scal : br)  * (key == SkinnyBond ? 0.65 : 1);   // bond radius smaller for skinnyBond
    bl = (key == FemaleJoinBond ? bl * bondLenFac : bl);  // make female joiner shorter
    if (key == MaleJoinBond) { // male join is direct solid, others need difference()
        joiner(bl, scal, male = true, ver = ver, supportSelect=supportSel);
    } else {  // regular bond / skinny / h-bond / female join
        bhblen = bl +(hblen/2 * scal);
        rotate([0,0,fnRot]) {
            difference() {
                union() {
                    cylinder(h=bl,r=br,center=false);
                    if (key == HBond) {  // make extension collar for h-bond magnet
                        rotate([0,0,oRot-fnRot]) cylinder(h=bhblen-1,r=(magR + clearance +wall),center=false, $fn=8);
                    }
                }
                atom(atm,scal,-clearance);  // remove overlap with atom to clear area for female join
                if (key == HBond) {     // make space to insert magnet inside bond cylinder
                    translate([0,0,(bhblen-magL)-pClearance])
                        cylinder(h=magL+pClearance+shim, r=magR+pClearance, center=false, $fn=8);
                }
            }
        }
    }
}

/*
//
// Generate a 'hedron', one plane of 3 points, consisting of 3 atoms joined by
//   two bonds.
//   Defined as bond length - bond angle - bond length
//
// In some cases the sequence of atoms in the h[] array is reversed (rev flag),
// as detailed in the comments.
//
// other parameters:
//
// h = hedron array data according to rev flag:
//   yes reversed                                     :  not reversed
//    0    1     2     3     4     5     6      7     :     0     1     2     3    4     5      6      7
//  len1  len3  atom1 atom3  a1    a2   a1-a2  a2-a3      len1  len3  atom1 atom3   a1    a3  a1-a2  a2-a3
//
// split: chop half of the hedron - to selectively print parts of a rotating
//   bond to be glued together.  top or bottom half selected by global caTop
//   (C-alpha top) variable, undef by default so bottom half.
//
// supporSel: enable support structure inside rotatable bond to print in place.
//  Please note the bond needs to be exactly parallel to the buildplate and the
//  layerHeight global variable above needs to be set correctly for the
//  structure to be correctly created by your slicer software.
//
 */

module hedron(h,rev=0,scal,split=0, supportSel) {

    newh = hFlip(h, rev);  // make a consistent hedron array regardless of rev flag

    bondRad = bondRadius * scal;
    difference() {
        union(){
            if (h[7]) {
                // central atom at 0,0,0
                atom(h[4],scal);
            }

            if (newh[5] && newh[7] != FemaleJoinBond) {  // not female join
                // comments for non-reversed case
                // atom 3 is len3 up on +z
                translate([0,0,newh[1]])
                    difference() {
                        atom(newh[3],scal * (newh[7] == SkinnyBond ? 0.7 : 1));  // if skinny bond make atom (C-beta) same diameter as bond
                        if (newh[7] == HBond) {  // make room for hbond magnet through atom - this branch not used for backbone N,O
                            translate([0,0,scal*hblen/2-magL-pClearance])
                                cylinder(h=magL+pClearance,r=magR+pClearance,$fn=8);
                        }
                    }
            }

            if (newh[7]) {
                // atom 2 - atom 3 bond from origin up +z distance len3
                bond(newh[1], bondRad, scal, newh[7], h[4], ver=1, supportSel=supportSel);
            }
            rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
                if (newh[6]) {
                    bond(newh[0], bondRad, scal, newh[6], h[4], ver=1, supportSel=supportSel);  // h[4] is center atom (atom 2)
                }
                if (newh[4] && newh[6] != FemaleJoinBond) {   // if draw atom 2 and atom1-atom2 not joiner
                    translate([0,0,newh[0]]) {
                        difference() {
                            atom(newh[2],scal * (newh[6] == SkinnyBond ? 0.7 : 1));  // put atom1 sphere len1 away on Z
                            if (newh[6] == HBond) {  // make room for hbond magnet through atom
                                translate([0,0,scal*hblen/2-magL-pClearance])
                                    cylinder(h=magL+pClearance,r=magR+pClearance,$fn=8);
                            }
                        }
                    }
                }
            }
        }

        if (split) {
            // top / bottom half cutter
            thick = 2*bondRadius * scal;
            Zdim = newh[0];
            Xdim = newh[1];

            cside = 7* defaultAtomRadius * atomScale * scal / 12 + (caTop ? pClearance : -pClearance);
            difference() {
                translate([-Xdim,((rev || caTop) ? 0 : -thick),-Zdim]) {
                    cube([2*Xdim,thick,2*Zdim]);
                }
                if (!caTop) {
                    rotate([0,(rev ? h[1] : 0),0])
                    rotate([45,0,0])
                    cube([cside, cside, cside],center=true);
                }
            }
            if (caTop) {
                //translate([tx+cside,0,tx+cside])
                    rotate([0,(rev ? h[1] : 0),0])
                        rotate([45,0,0])
                        cube([cside, cside, cside], center=true);
            }
        }

        if (newh[7] == FemaleJoinBond) {  // female join
            joiner(newh[1], scal, male=false, ver=1, supportSelect=supportSel);
        }

        if (newh[6] == FemaleJoinBond) {  // female join
            rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
            joiner(newh[0], scal, male=false, ver=1, supportSelect=supportSel);
            translate([0,0,newh[0]])
                atom(newh[2],scal+0.5,clearance);  // clearance for atom against join outer cylinder
            }
        }

        if (newh[7] == FemaleJoinBond || newh[6] == FemaleJoinBond) {  // female join both hedron arms
            translate([0,0,newh[1]]) atom(newh[3],scal+0.5,clearance);  // clearance for atom against join outer cylinder
        }
    }
}

/*
//
// Hook to call custom routines for specific hedra.
//
// Residue is h[h_residue]
// Sequence position is h[h_seqpos]
//
*/
module hedronDispatch(h,rev=0,scal) {
    // default action is just to pass to hedron()

    hedron(h, rev, scal, 0, (support ? 1 : 0));

    /*
    // Some examples for special handling for specific hedra below:
    // note use of h_seqpos, h_residue, h_class for selecting hedra

    // bool flag caTop (for rotatable bond part) needs to be a global variable
    // so hedron() above can see it.

caBase1 = false;   // only make bottom of N_C-alpha_C hedron
caBase2 = false;   // same as caBase1 but for case of reversed hedron (for testing, should be identical to caBase1 result)
amideOnly = false; // make only the first amide

    if (caTop) {
        // these examples select a specific sequence position (h[h_seqpos] == n)
        if (h[h_seqpos] == 1) {
            if (h[h_class] == "NCAC") {
                hedron(h, rev, scal, 1);
            } else if (h[h_class] == "CBCAC") {
                color("yellow") {  // ca-cb
                    hedron(h, rev, scal);
                }
            }
        }
    } else if (caBase1) {
        if (h[h_seqpos] == 1 && (h[h_class] == "NCAC")) {
            hedron(h, rev, scal, true, (support ? 1 : 0));
        }
    } else if (caBase2) {
        if (h[h_seqpos] == 5 && (h[h_class] == "NCAC")) {
            hedron(h, rev, scal, true, (support ? 1 : 0));
        }
    } else if (amideOnly) {
        if (h[h_seqpos] == 1) {
            if (h[h_class] == "CACN") {
                color("darkgray") {
                    hedron(h, rev, scal);
                }
            }  else if (h[h_class] == "CACO") {
                color("red") {   // c=o
                    hedron(h, rev, scal);
                }
            }  else if (h[h_class] == "CNCA") {
                color("cyan") {  // h=n
                    hedron(h, rev, scal);
                }
            }
        } else if ((h[h_seqpos] == 2) && (h[h_class] == "HNCA")) {
            color("cyan") {  // h=n
                hedron(h, rev, scal);
            }
        }
       // actions above select out only a single hedron
    } else {
        // actions below will process hedra all but handle selected ones differently

        if (h[h_class] == "NCAC") {
            if (h[h_seqpos] == 1) {
                if (! CCap && NCap) {  // make split rotatable bond for terminal NH3
                    hedron(h, rev, scal, true, (support ? 1 : 0));
                }
            } else if (h[h_seqpos] == 5) {  // make split rotatable bond for terminal COOH
                hedron(h, rev, scal, true, (support ? 2 : 0));  // note supportSel = 2
            } else {
                hedron(h, rev, scal, 0, (support ? 2 : 0));
            }
        } else if (h[h_class] == "CBCAC") {
            color("yellow") {                     // ca-cb -- color yellow in OpenSCAD renderer
                if (h[h_seqpos] == 1 ) {         // don't make here for N-term
                } else if (h[h_seqpos] == 5 ) {  // don't make here for C-term
                } else {
                    hedron(h, rev, scal);       // otherwise do make here
                }
            }
        } else if (h[h_class] == "HNCA") {
            color("cyan") { // color h-n in OenSCAD renderer
                if (h[h_seqpos] == 1) {
                    if (NCap) {                      // only make at N term if variable NCap is true
                        hedron(h, rev, scal, 0, (support ? 1 : 0));
                    }
                } else {
                    hedron(h, rev, scal, 0, (support ? 1 : 0));
                }
            }
        } else if (h[h_residue] == "P") {
            color("darkgray")   // highlight Prolines in OpenSCAD renderer
                hedron(h, rev, scal);
        } else {
            echo("unrecognised hedron", h[h_class]);
            color("pink")
                hedron(h, rev, scal, 0, (support ? 1 : 0));
        }
    }
    */
}

/*
//
// Generate a hedron rotated to specific angle d
//
*/
module d2(d,hedra,scal)
{
    tz = (d[d_reversed] ? hedra[d[d_h2ndx]][2] : hedra[d[d_h2ndx]][0]);      // get h2 len1 depending on reversed
    rotate(d[d_dangle1]) {                                                   // 4. rotate h2 to specified dihedral angle
        translate([0,0,tz]) {                                               // 3. translate h2 h2:len1 up +z
            rotate([180, 0, 0]) {                                          // 2. rotate h2r about X so h2:a3 in +z and h2:a1 in -z
                hedronDispatch(hedra[d[d_h2ndx]],(!d[d_reversed]),scal);  // 1. reverse hedron 2 orientation = h2r
            }
        }
    }
}

/*
//
// Generate two hedra at specified dihedral angle d
//
*/
module dihedron(d,hedra,scal)
{
    if (d[d_h1new])
        hedronDispatch(hedra[d[d_h1ndx]],d[d_reversed],scal);                // reverse h1 if dihedral reversed
    if (d[d_h2new])
        d2(d,hedra,scal);
}

/*
//
// Generate a residue consisting of the set of dihedra in the parameter 'r',
//   referring to hedra the table specified in the parameter 'hedra'.
//
*/
module residue(r,hedra, scal)
{
    for (d = r) {
        multmatrix(d[d_dihedralTransform]) {
            dihedron(d, hedra, scal);
        }
    }
}

/*
//
// Generate a chain of residues, each positioned by a supplied
// rotation/translation matrix.
//
*/
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

/*
//
// OpenSCAD array indices to reference protein data - tied to BioPython code
//
*/

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
h_atom1class = 3;
h_atom2class = 4;
h_atom3class = 5;
h_atom1state = 6;
h_atom2state = 7;
h_atom3state = 8;
h_bond1state = 9;
h_bond2state = 10;
h_residue = 11;
h_seqpos = 12;  // residue sequence position for first atom in hedra
h_class = 13;

// dihedra specifications for each residue in sequence, dihedral array
d_dangle1 = 0;
d_h1ndx = 1;
d_h2ndx = 2;
d_reversed = 3;
d_h1new = 4;
d_h2new = 5;
d_dihedralTransform = 6;

// residueSet: world transform for each residue in sequence array
r_resNdx = 0;
r_resID = 1;
r_resTransform = 2;


// use single default atom radius for all atoms if tubes = true, else use
// covalent radii from literature
atomData = ( tubes ?
            [   ["Csb","green" , defaultAtomRadius], ["Cres","green" , defaultAtomRadius], ["Cdb","green" , defaultAtomRadius],
                ["Osb","red" , defaultAtomRadius], ["Ores","red" , defaultAtomRadius], ["Odb","red" , defaultAtomRadius],
                ["Nsb","blue" , defaultAtomRadius], ["Nres","blue" , defaultAtomRadius], ["Ndb","blue" , defaultAtomRadius],
                ["Hsb","gray" , defaultAtomRadius],
                ["Ssb","yellow" , defaultAtomRadius] ]
            :

// covalent radii from Heyrovska, Raji : 'Atomic Structures of all the Twenty
// Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic
// Covalent Radii'  https://arxiv.org/pdf/0804.2488.pdf

            [   ["Csb","green" , 0.77], ["Cres","green" , 0.72], ["Cdb","green" , 0.67],
                ["Osb","red" , 0.67], ["Ores","red" , 0.635], ["Odb","red" , 0.60],
                ["Nsb","blue" , 0.70], ["Nres","blue" , 0.66], ["Ndb","blue" , 0.62],
                ["Hsb","gray" , 0.37],
                ["Ssb","yellow" , 1.04] ]
    );


// optionally include protein array data here [ write_SCAD(includeCode=False) ], e.g.:
// include <1rtm.scad>;
// or paste below

"""  # noqa
