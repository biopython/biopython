# Copyright 2019-22 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes to support internal coordinates for protein structures.

Internal coordinates comprise Psi, Omega and Phi dihedral angles along the
protein backbone, Chi angles along the sidechains, and all 3-atom angles and
bond lengths defining a protein chain.  These routines can compute internal
coordinates from atom XYZ coordinates, and compute atom XYZ coordinates from
internal coordinates.

Secondary benefits include the ability to align and compare residue
environments in 3D structures, support for 2D atom distance plots, converting a
distance plot plus chirality information to a structure, generating an OpenSCAD
description of a structure for 3D printing, and reading/writing structures as
internal coordinate data files.

**Usage:**
::

    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.Chain import Chain
    from Bio.PDB.internal_coords import *
    from Bio.PDB.PICIO import write_PIC, read_PIC, read_PIC_seq
    from Bio.PDB.ic_rebuild import write_PDB, IC_duplicate, structure_rebuild_test
    from Bio.PDB.SCADIO import write_SCAD
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.PDB.PDBIO import PDBIO
    import numpy as np

    # load a structure as normal, get first chain
    parser = PDBParser()
    myProtein = parser.get_structure("7rsa", "pdb7rsa.ent")
    myChain = myProtein[0]["A"]

    # compute bond lengths, angles, dihedral angles
    myChain.atom_to_internal_coordinates(verbose=True)

    # check myChain makes sense (can get angles and rebuild same structure)
    resultDict = structure_rebuild_test(myChain)
    assert resultDict['pass'] == True

    # get residue 1 chi2 angle
    r1 = next(myChain.get_residues())
    r1chi2 = r1.internal_coord.get_angle("chi2")

    # rotate residue 1 chi2 angle by 120 degrees (loops w/in +/-180)
    r1.internal_coord.set_angle("chi2", r1chi2 + 120.0)
    # or
    r1.internal_coord.bond_rotate("chi2", 120.0)
    # update myChain XYZ coordinates with chi2 changed
    myChain.internal_to_atom_coordinates()
    # write new conformation with PDBIO
    write_PDB(myProtein, "myChain.pdb")
    # or just the ATOM records without headers:
    io = PDBIO()
    io.set_structure(myProtein)
    io.save("myChain2.pdb")

    # write chain as 'protein internal coordinates' (.pic) file
    write_PIC(myProtein, "myChain.pic")
    # read .pic file
    myProtein2 = read_PIC("myChain.pic")

    # create default structure for random sequence by reading as .pic file
    myProtein3 = read_PIC_seq(
        SeqRecord(
            Seq("GAVLIMFPSTCNQYWDEHKR"),
            id="1RND",
            description="my random sequence",
        )
    )
    myProtein3.internal_to_atom_coordinates()
    write_PDB(myProtein3, "myRandom.pdb")

    # access the all-dihedrals array for the chain, e.g. residue 1 chi2 angle:
    r1chi2_obj = r1.internal_coord.pick_angle("chi2")
    # or same thing: r1chi2_obj = r1.internal_coord.pick_angle("CA:CB:CG:CD")
    r1chi2_key = r1chi2_obj.atomkeys
    # r1chi2_key is tuple of AtomKeys (1_K_CA, 1_K_CB, 1_K_CG, 1_K_CD)
    r1chi2_index = myChain.internal_coord.dihedraNdx[r1chi2_key]
    # or same thing: r1chi2_index = r1chi2_obj.ndx
    r1chi2_value = myChain.internal_coord.dihedraAngle[r1chi2_index]
    # also true: r1chi2_obj == myChain.internal_coord.dihedra[r1chi2_index]

    # access the array of all atoms for the chain, e.g. residue 1 C-beta
    r1_cBeta_index = myChain.internal_coord.atomArrayIndex[AtomKey("1_K_CB")]
    r1_cBeta_coords = myChain.internal_coord.atomArray[r1_cBeta_index]
    # r1_cBeta_coords = [ x, y, z, 1.0 ]

    # the Biopython Atom coord array is now a view into atomArray, so
    assert r1_cBeta_coords[1] == r1["CB"].coord[1]
    r1_cBeta_coords[1] += 1.0  # change the Y coord 1 angstrom
    assert r1_cBeta_coords[1] == r1["CB"].coord[1]
    # they are always the same (they share the same memory)
    r1_cBeta_coords[1] -= 1.0  # restore

    # create a selector to filter just the C-alpha atoms from the all atom array
    atmNameNdx = AtomKey.fields.atm
    atomArrayIndex = myChain.internal_coord.atomArrayIndex
    CaSelect = [
        atomArrayIndex.get(k) for k in atomArrayIndex.keys() if k.akl[atmNameNdx] == "CA"
    ]
    # now the ordered array of C-alpha atom coordinates is:
    CA_coords = myChain.internal_coord.atomArray[CaSelect]
    # note this uses Numpy fancy indexing, so CA_coords is a new copy

    # create a C-alpha distance plot
    caDistances = myChain.internal_coord.distance_plot(CaSelect)
    # display with e.g. MatPlotLib:
    # import matplotlib.pyplot as plt
    # plt.imshow(caDistances, cmap="hot", interpolation="nearest")
    # plt.show()

    # build structure from distance plot:
    ## create the all-atom distance plot
    distances = myChain.internal_coord.distance_plot()
    ## get the sign of the dihedral angles
    chirality = myChain.internal_coord.dihedral_signs()
    ## get new, empty data structure : copy data structure from myChain
    myChain2 = IC_duplicate(myChain)[0]["A"]
    cic2 = myChain2.internal_coord
    ## clear the new atomArray and di/hedra value arrays, just for proof
    cic2.atomArray = np.zeros((cic2.AAsiz, 4), dtype=np.float64)
    cic2.dihedraAngle[:] = 0.0
    cic2.hedraAngle[:] = 0.0
    cic2.hedraL12[:] = 0.0
    cic2.hedraL23[:] = 0.0
    ## copy just the first N-Ca-C coords so structures will superimpose:
    cic2.copy_initNCaCs(myChain.internal_coord)
    ## copy distances to chain arrays:
    cic2.distplot_to_dh_arrays(distances, chirality)
    ## compute angles and dihedral angles from distances:
    cic2.distance_to_internal_coordinates()
    ## generate XYZ coordinates from internal coordinates:
    myChain2.internal_to_atom_coordinates()
    ## confirm result atomArray matches original structure:
    assert np.allclose(cic2.atomArray, myChain.internal_coord.atomArray)

    # superimpose all phe-phe pairs - quick hack just to demonstrate concept
    # for analyzing pairwise residue interactions.  Generates PDB ATOM records
    # placing each PHE at origin and showing all other PHEs in environment
    ## shorthand for key variables:
    cic = myChain.internal_coord
    resNameNdx = AtomKey.fields.resname
    aaNdx = cic.atomArrayIndex
    ## select just PHE atoms:
    pheAtomSelect = [aaNdx.get(k) for k in aaNdx.keys() if k.akl[resNameNdx] == "F"]
    aaF = cic.atomArray[ pheAtomSelect ]  # numpy fancy indexing makes COPY not view

    for ric in cic.ordered_aa_ic_list:  # internal_coords version of get_residues()
        if ric.rbase[2] == "F":  # if PHE, get transform matrices for chi1 dihedral
            chi1 = ric.pick_angle("N:CA:CB:CG")  # chi1 space has C-alpha at origin
            cst = np.transpose(chi1.cst)  # transform TO chi1 space
            # rcst = np.transpose(chi1.rcst)  # transform FROM chi1 space
            cic.atomArray[pheAtomSelect] = aaF.dot(cst)  # transform just the PHEs
            for res in myChain.get_residues():  # print PHEs in new coordinate space
                if res.resname in ["PHE"]:
                    print(res.internal_coord.pdb_residue_string())
            cic.atomArray[pheAtomSelect] = aaF  # restore coordinate space from copy

    # write OpenSCAD program of spheres and cylinders to 3d print myChain backbone
    ## set atom load filter to accept backbone only:
    IC_Residue.accept_atoms = IC_Residue.accept_backbone
    ## delete existing data to force re-read of all atoms:
    myChain.internal_coord = None
    write_SCAD(myChain, "myChain.scad", scale=10.0)

See the `''Internal coordinates module''` section of the `Biopython Tutorial
and Cookbook` for further discussion.

**Terms and key data structures:**
Internal coordinates are defined on sequences of atoms which span
residues or follow accepted nomenclature along sidechains.  To manage these
sequences and support Biopython's disorder mechanisms, :class:`AtomKey`
specifiers are implemented to capture residue, atom and variant identification
in a single object.  A :class:`Hedron` object is specified as three sequential
AtomKeys, comprising two bond lengths and the bond angle between them.  A
:class:`Dihedron` consists of four sequential AtomKeys, linking two Hedra with
a dihedral angle between them.

**Algorithmic overview:**
The Internal Coordinates module combines a specification of connected atoms as
hedra and dihedra in the :mod:`.ic_data` file with routines here to transform
XYZ coordinates of these atom sets between a local coordinate system and the
world coordinates supplied in e.g. a PDB or mmCif data file.  The local
coordinate system places the center atom of a hedron at the origin (0,0,0), one
leg on the +Z axis, and the other leg on the XZ plane (see :class:`Hedron`).
Measurement and creation or manipulation of hedra and dihedra in the local
coordinate space is straightforward, and the calculated transformation matrices
enable assembling these subunits into a protein chain starting from supplied
(PDB) coordinates for the initial N-Ca-C atoms.

Psi and Phi angles are defined on atoms from adjacent residues in a protein
chain, see e.g. :meth:`.pick_angle` and :mod:`.ic_data` for the relevant
mapping between residues and backbone dihedral angles.

Transforms to and from the dihedron local coordinate space described above are
accessible via :data:`IC_Chain.dCoordSpace` and :class:`Dihedron` attributes
.cst and .rcst, and may be applied in the alignment and comparison of residues
and their environments with code along the lines of::

    chi1 = ric0.pick_angle("chi1") # chi1 space defined with CA at origin
    cst = np.transpose(chi1.cst) # transform TO chi1 local space
    newAtomCoords = oldAtomCoords.dot(cst)

The core algorithms were developed independently during 1993-4 for
`''Development and Application of a Three-dimensional Description of Amino Acid
Environments in Protein,'' Miller, Douthart, and Dunker, Advances in Molecular
Bioinformatics, IOS Press, 1994, ISBN 90 5199 172 x, pp. 9-30.
<https://www.google.com/books/edition/Advances_in_Molecular_Bioinformatics/VmFSNNm7k6cC?gbpv=1>`_

A Protein Internal Coordinate (.pic) file format is defined to capture
sufficient detail to reproduce a PDB file from chain starting coordinates
(first residue N, Ca, C XYZ coordinates) and remaining internal coordinates.
These files are used internally to verify that a given structure can be
regenerated from its internal coordinates.  See :mod:`.PICIO` for reading and
writing .pic files and :func:`.structure_rebuild_test` to determine if a
specific PDB or mmCif datafile has sufficient information to interconvert
between cartesian and internal coordinates.

Internal coordinates may also be exported as `OpenSCAD <https://www.openscad.org>`_
data arrays for generating 3D printed protein models.  OpenSCAD software is
provided as a starting point and proof-of-concept for generating such models.
See :mod:`.SCADIO` and this `Thingiverse project <https://www.thingiverse.com/thing:3957471>`_
for a more advanced example.

Refer to :meth:`.distance_plot` and :meth:`.distance_to_internal_coordinates`
for converting structure data to/from 2D distance plots.

The following classes comprise the core functionality for processing internal
coordinates and are sufficiently related and coupled to place them together in
this module:

:class:`IC_Chain`: Extends Biopython Chain on .internal_coord attribute.
    Manages connected sequence of residues and chain breaks; holds numpy arrays
    for all atom coordinates and bond geometries. For 'parallel' processing
    IC_Chain methods operate on these arrays with single numpy commands.

:class:`IC_Residue`: Extends Biopython Residue on .internal_coord attribute.
    Access for per residue views on internal coordinates and methods for serial
    (residue by residue) assembly.

:class:`Dihedron`: four joined atoms forming a dihedral angle.
    Dihedral angle, homogeneous atom coordinates in local coordinate space,
    references to relevant Hedra and IC_Residue.  Getter methods for
    residue dihedral angles, bond angles and bond lengths.

:class:`Hedron`: three joined atoms forming a plane.
    Contains homogeneous atom coordinates in local coordinate space as well as
    bond lengths and angle between them.

:class:`Edron`: base class for Hedron and Dihedron classes.
    Tuple of AtomKeys comprising child, string ID, mainchain membership boolean
    and other routines common for both Hedra and Dihedra.  Implements rich
    comparison.

:class:`AtomKey`: keys (dictionary and string) for referencing atom sequences.
    Capture residue and disorder/occupancy information, provides a
    no-whitespace key for .pic files, and implements rich comparison.

Custom exception classes: :class:`HedronMatchError` and
:class:`MissingAtomError`
"""  # noqa

import re
from collections import deque, namedtuple
import copy

# from numpy import floor, ndarray
from numbers import Integral

import numpy as np  # type: ignore

from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.Data.PDBData import protein_letters_3to1

from Bio.PDB.vectors import multi_coord_space, multi_rot_Z
from Bio.PDB.vectors import coord_space

from Bio.PDB.ic_data import ic_data_backbone, ic_data_sidechains
from Bio.PDB.ic_data import primary_angles
from Bio.PDB.ic_data import ic_data_sidechain_extras, residue_atom_bond_state
from Bio.PDB.ic_data import dihedra_primary_defaults, hedra_defaults

# for type checking only
from typing import (
    List,
    Dict,
    Set,
    TextIO,
    Union,
    Tuple,
    cast,
    TYPE_CHECKING,
    Optional,
)

if TYPE_CHECKING:
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Chain import Chain

HKT = Tuple["AtomKey", "AtomKey", "AtomKey"]  # Hedron key tuple
DKT = Tuple["AtomKey", "AtomKey", "AtomKey", "AtomKey"]  # Dihedron Key Tuple
EKT = Union[HKT, DKT]  # Edron Key Tuple
BKT = Tuple["AtomKey", "AtomKey"]  # Bond Key Tuple

# HACS = Tuple[np.array, np.array, np.array]  # Hedron Atom Coord Set
HACS = np.array  # Hedron Atom Coord Set
DACS = Tuple[np.array, np.array, np.array, np.array]  # Dihedron Atom Coord Set


class IC_Chain:
    """Class to extend Biopython Chain with internal coordinate data.

    Attributes
    ----------
    chain: object reference
        The Biopython :class:`Bio.PDB.Chain.Chain` object this extends

    MaxPeptideBond: float
        **Class** attribute to detect chain breaks.
        Override for fully contiguous chains with some very long bonds - e.g.
        for 3D printing (OpenSCAD output) a structure with missing residues.
        :data:`MaxPeptideBond`

    ParallelAssembleResidues: bool
        **Class** attribute affecting internal_to_atom_coords.
        Short (50 residue and less) chains are faster to assemble without the
        overhead of creating numpy arrays, and the algorithm is easier to
        understand and trace processing a single residue at a time.  Clearing
        (set to False) this flag will switch to the serial algorithm

    ordered_aa_ic_list: list
        IC_Residue objects internal_coords algorithms can process (e.g. no
        waters)

    initNCaC: List of N, Ca, C AtomKey tuples (NCaCKeys).
        NCaCKeys start chain segments (first residue or after chain break).
        These 3 atoms define the coordinate space for a contiguous chain
        segment, as initially specified by PDB or mmCIF file.

    AAsiz = int
        AtomArray size, number of atoms in this chain

    atomArray: numpy array
        homogeneous atom coords ([x,, y, z, 1.0]) for every atom in chain

    atomArrayIndex: dict
        maps AtomKeys to atomArray indexes

    hedra: dict
        Hedra forming residues in this chain; indexed by 3-tuples of AtomKeys.

    hedraLen: int
        length of hedra dict

    hedraNdx: dict
        maps hedra AtomKeys to numeric index into hedra data arrays e.g.
        hedraL12 below

    a2ha_map: [hedraLen x 3]
        atom indexes in hedraNdx order

    dihedra: dict
        Dihedra forming residues in this chain; indexed by 4-tuples of AtomKeys.

    dihedraLen: int
        length of dihedra dict

    dihedraNdx: dict
        maps dihedra AtomKeys to dihedra data arrays e.g. dihedraAngle

    a2da_map : [dihedraLen x 4]
        AtomNdx's in dihedraNdx order

    d2a_map : [dihedraLen x [4]]
        AtomNdx's for each dihedron (reshaped a2da_map)

    Numpy arrays for vector processing of chain di/hedra:

    hedraL12: numpy array
        bond length between hedron 1st and 2nd atom
    hedraAngle: numpy array
        bond angle for each hedron, in degrees
    hedraL23: numpy array
        bond length between hedron 2nd and 3rd atom

    id3_dh_index: dict
        maps hedron key to list of dihedra starting with hedron, used by
        assemble and bond_rotate to find dihedra with h1 key

    id32_dh_index: dict
        like id3_dh_index, find dihedra from h2 key

    hAtoms: numpy array
        homogeneous atom coordinates (3x4) of hedra, central atom at origin

    hAtomsR: numpy array
        hAtoms in reverse orientation

    hAtoms_needs_update: numpy array of bool
        indicates whether hAtoms represent hedraL12/A/L23

    dihedraAngle: numpy array
        dihedral angles (degrees) for each dihedron

    dAtoms: numpy array
        homogeneous atom coordinates (4x4) of dihedra, second atom at origin

    dAtoms_needs_update: numpy array of bool
        indicates whether dAtoms represent dihedraAngle

    dCoordSpace: numpy array
        forward and reverse transform matrices standardising positions of first
        hedron.  See :data:`dCoordSpace`.

    dcsValid: bool
        indicates dCoordSpace up to date

    See also attributes generated by :meth:`build_edraArrays` for indexing
    di/hedra data elements.

    Methods
    -------
    internal_to_atom_coordinates:
        Process ic data to Residue/Atom coordinates; calls assemble_residues()
    assemble_residues:
        Generate IC_Chain atom coords from internal coordinates (parallel)
    assemble_residues_ser:
        Generate IC_Residue atom coords from internal coordinates (serial)
    atom_to_internal_coordinates:
        Calculate dihedrals, angles, bond lengths (internal coordinates) for
        Atom data
    write_SCAD:
        Write OpenSCAD matrices for internal coordinate data comprising chain;
        this is a support routine, see :func:`.SCADIO.write_SCAD` to generate
        OpenSCAD description of a protein chain.
    distance_plot:
        Generate 2D plot of interatomic distances with optional filter
    distance_to_internal_coordinates:
        Compute internal coordinates from distance plot and array of dihedral
        angle signs.
    make_extended:
        Arbitrarily sets all psi and phi backbone angles to 123 and -104 degrees.

    """

    # Class globals
    MaxPeptideBond = 1.4
    """Larger C-N distance than this will be chain break"""

    ParallelAssembleResidues = True
    """Enable parallel internal_to_atom algorithm, is slower for short chains"""

    AAsiz = 0
    """Number of atoms in this chain (size of atomArray)"""

    atomArray: np.array = None
    """AAsiz x [4] of float np.float64 homogeneous atom coordinates, all atoms
    in chain."""

    dCoordSpace = None
    """[2][dihedraLen][4][4] : 2 arrays of 4x4 coordinate space transforms for
    each dihedron.  The first [0] converts TO standard space with first atom on
    the XZ plane, the second atom at the origin, the third on the +Z axis, and
    the fourth placed according to the dihedral angle.  The second [1] transform
    returns FROM the standard space to world coordinates (PDB file input or
    whatever is current).  Also accessible as .cst (forward
    transform) and .rcst (reverse transform) in :class:`Dihedron`."""

    dcsValid = None
    """True if dCoordSpace is up to date.  Use :meth:`.update_dCoordSpace`
    if needed."""

    # for assemble_residues
    _dihedraSelect = np.array([True, True, True, False])
    _dihedraOK = np.array([True, True, True, True])

    def __init__(self, parent: "Chain", verbose: bool = False) -> None:
        """Initialize IC_Chain object, with or without residue/Atom data.

        :param Bio.PDB.Chain parent: Biopython Chain object
            Chain object this extends
        """
        # type hinting parent as Chain leads to import cycle
        self.chain = parent
        self.ordered_aa_ic_list: List[IC_Residue] = []
        # self.initNCaC: Dict[Tuple[str], Dict["AtomKey", np.array]] = {}
        self.initNCaCs = []

        self.sqMaxPeptideBond = np.square(IC_Chain.MaxPeptideBond)
        # need init here for _gen_edra():
        self.hedra = {}
        # self.hedraNdx = {}
        self.dihedra = {}
        # self.dihedraNdx = {}

        # cache of AtomKey results for cak()
        # self.akc: Dict[Tuple(IC_Residue, str), AtomKey] = {}

        self.atomArrayIndex: Dict["AtomKey", int] = {}

        self.bpAtomArray: List["Atom"] = []  # rtm

        self._set_residues(verbose)  # no effect if no residues loaded

    def __deepcopy__(self, memo) -> "IC_Chain":
        """Implement deepcopy for IC_Chain."""
        existing = memo.get(id(self), False)
        if existing:
            return existing
        dup = type(self).__new__(self.__class__)
        memo[id(self)] = dup
        dup.chain = memo[id(self.chain)]
        dup.chain.child_dict = copy.deepcopy(self.chain.child_dict, memo)
        # now have all res and ic_res but ic_res not complete
        dup.chain.child_list = copy.deepcopy(self.chain.child_list, memo)
        dup.akset = copy.deepcopy(self.akset, memo)
        dup.aktuple = copy.deepcopy(self.aktuple, memo)
        # now have all ak w/.ric
        dup.ordered_aa_ic_list = copy.deepcopy(self.ordered_aa_ic_list, memo)

        dup.atomArrayIndex = self.atomArrayIndex.copy()
        dup.atomArrayValid = self.atomArrayValid.copy()
        dup.atomArray = self.atomArray.copy()

        dup.hedra = copy.deepcopy(self.hedra, memo)
        dup.dihedra = copy.deepcopy(self.dihedra, memo)

        dup.id3_dh_index = copy.deepcopy(self.id3_dh_index, memo)
        dup.id32_dh_index = copy.deepcopy(self.id32_dh_index, memo)

        # update missing items in ic_residues and
        # set all bp residue atom coords to be views on dup.atomArray
        # [similar in build_AtomArray() but does not copy from bpAtoms
        # or modify atomArrayValid, and accesses dup]
        dup.AAsiz = self.AAsiz

        dup.bpAtomArray = [None] * dup.AAsiz  # rtm

        def setAtomVw(res, atm):
            ak = AtomKey(res.internal_coord, atm)
            ndx = dup.atomArrayIndex[ak]
            atm.coord = dup.atomArray[ndx, 0:3]  # make view on atomArray

            dup.bpAtomArray[ndx] = atm  # rtm

        def setResAtmVws(res):
            for atm in res.get_atoms():
                # copy not filter so ignore no_altloc
                if atm.is_disordered():
                    for altAtom in atm.child_dict.values():
                        setAtomVw(res, altAtom)
                else:
                    setAtomVw(res, atm)

        for ric in dup.ordered_aa_ic_list:
            setResAtmVws(ric.residue)
            ric.rprev = copy.deepcopy(ric.rprev, memo)
            ric.rnext = copy.deepcopy(ric.rnext, memo)
            ric.ak_set = copy.deepcopy(ric.ak_set, memo)
            ric.akc = copy.deepcopy(ric.akc, memo)
            ric.dihedra = copy.deepcopy(ric.dihedra, memo)
            ric.hedra = copy.deepcopy(ric.hedra, memo)

        dup.sqMaxPeptideBond = self.sqMaxPeptideBond
        dup.initNCaCs = copy.deepcopy(self.initNCaCs, memo)

        dup.hedraLen = self.hedraLen
        dup.hedraL12 = self.hedraL12.copy()
        dup.hedraAngle = self.hedraAngle.copy()
        dup.hedraL23 = self.hedraL23.copy()
        dup.hedraNdx = copy.deepcopy(self.hedraNdx, memo)

        dup.dihedraLen = self.dihedraLen
        dup.dihedraAngle = self.dihedraAngle.copy()
        dup.dihedraAngleRads = self.dihedraAngleRads.copy()
        dup.dihedraNdx = copy.deepcopy(self.dihedraNdx, memo)

        dup.a2da_map = self.a2da_map.copy()
        dup.a2d_map = self.a2d_map.copy()
        dup.d2a_map = self.d2a_map.copy()

        dup.dH1ndx = self.dH1ndx.copy()
        dup.dH2ndx = self.dH2ndx.copy()

        dup.hAtoms = self.hAtoms.copy()
        dup.hAtomsR = self.hAtomsR.copy()
        dup.hAtoms_needs_update = self.hAtoms_needs_update.copy()

        dup.dRev = self.dRev.copy()
        dup.dFwd = self.dFwd.copy()
        dup.dAtoms_needs_update = self.dAtoms_needs_update.copy()

        dup.dAtoms = self.dAtoms.copy()
        dup.a4_pre_rotation = self.a4_pre_rotation.copy()

        dup.dCoordSpace = self.dCoordSpace.copy()
        dup.dcsValid = self.dcsValid.copy()

        for d in dup.dihedra.values():
            d.cst = dup.dCoordSpace[0][d.ndx]
            d.rcst = dup.dCoordSpace[1][d.ndx]

        return dup

    # return True if a0, a1 within supplied cutoff
    def _atm_dist_chk(self, a0: Atom, a1: Atom, cutoff: float, sqCutoff: float) -> bool:
        return sqCutoff > np.sum(np.square(a0.coord - a1.coord))

    # return a string describing issue, or None if OK
    def _peptide_check(self, prev: "Residue", curr: "Residue") -> Optional[str]:
        if 0 == len(curr.child_dict):
            # curr residue with no atoms => reading pic file, no break
            return None
        if (0 != len(curr.child_dict)) and (0 == len(prev.child_dict)):
            # prev residue with no atoms, curr has atoms => reading pic file,
            # have break
            return "PIC data missing atoms"

        # handle non-standard AA not marked as HETATM (1KQF, 1NTH)
        if not prev.internal_coord.isAccept:
            return "previous residue not standard/accepted amino acid"

        # both biopython Residues have Atoms, so check distance
        Natom = curr.child_dict.get("N", None)
        pCatom = prev.child_dict.get("C", None)
        if Natom is None or pCatom is None:
            return f"missing {'previous C' if pCatom is None else 'N'} atom"

        # confirm previous residue has all backbone atoms
        pCAatom = prev.child_dict.get("CA", None)
        pNatom = prev.child_dict.get("N", None)
        if pNatom is None or pCAatom is None:
            return "previous residue missing N or Ca"

        if IC_Residue.no_altloc:
            if Natom.is_disordered():
                Natom = Natom.selected_child
            if pCatom.is_disordered():
                pCatom = pCatom.selected_child

        if IC_Residue.no_altloc or (
            not Natom.is_disordered() and not pCatom.is_disordered()
        ):
            dc = self._atm_dist_chk(
                Natom, pCatom, IC_Chain.MaxPeptideBond, self.sqMaxPeptideBond
            )
            if dc:
                return None
            else:
                return f"MaxPeptideBond ({IC_Chain.MaxPeptideBond} angstroms) exceeded"

        # drop through for else Natom or pCatom is disordered:

        Nlist: List[Atom] = []
        pClist: List[Atom] = []
        if Natom.is_disordered():
            Nlist.extend(Natom.child_dict.values())
        else:
            Nlist = [Natom]
        if pCatom.is_disordered():
            pClist.extend(pCatom.child_dict.values())
        else:
            pClist = [pCatom]

        for n in Nlist:
            for c in pClist:
                if self._atm_dist_chk(
                    n, c, IC_Chain.MaxPeptideBond, self.sqMaxPeptideBond
                ):
                    return None
        return f"MaxPeptideBond ({IC_Chain.MaxPeptideBond} angstroms) exceeded"

    def clear_ic(self):
        """Clear residue internal_coord settings for this chain."""
        for res in self.chain.get_residues():
            res.internal_coord = None

    def _add_residue(
        self,
        res: "Residue",
        last_res: List,
        last_ord_res: List,
        verbose: bool = False,
    ) -> bool:
        """Set rprev, rnext, manage chain break.

        Returns True for no chain break or residue has sufficient data to
        restart at this position after a chain break (sets initNCaC AtomKeys
        in this case).  False return means insufficient data to extend chain
        with this residue.
        """
        # overwrite any existing .internal_coord in case re-initialising chain
        # expected state here is res.internal_coord = None
        res.internal_coord = IC_Residue(res)
        res.internal_coord.cic = self
        ric = res.internal_coord
        if (
            0 < len(last_res)
            and last_ord_res == last_res
            and self._peptide_check(last_ord_res[0].residue, res) is None
        ):
            # no chain break
            for prev in last_ord_res:
                prev.rnext.append(res.internal_coord)
                ric.rprev.append(prev)
            return True
        elif all(atm in res.child_dict for atm in ("N", "CA", "C")):
            # chain break, save coords for restart
            if verbose and len(last_res) != 0:  # not first residue
                if last_ord_res != last_res:
                    reason = f"disordered residues after {last_ord_res.pretty_str()}"
                else:
                    reason = cast(
                        str, self._peptide_check(last_ord_res[0].residue, res)
                    )
                print(f"chain break at {ric.pretty_str()} due to {reason}")

            iNCaC = ric.split_akl(
                (AtomKey(ric, "N"), AtomKey(ric, "CA"), AtomKey(ric, "C"))
            )
            self.initNCaCs.extend(iNCaC)
            return True

        # chain break but do not have N, Ca, C coords to restart from
        return False

    def _set_residues(self, verbose: bool = False) -> None:
        """Initialize .internal_coord for loaded Biopython Residue objects.

        Add IC_Residue as .internal_coord attribute for each :class:`.Residue`
        in parent :class:`Bio.PDB.Chain.Chain`; populate ordered_aa_ic_list with
        :class:`IC_Residue` references for residues which can be built (amino
        acids and some hetatms); set rprev and rnext on each sequential
        IC_Residue, populate initNCaC at start and after chain breaks.

        Generates:
            self.akset : set of :class:`.AtomKey` s in this chain
        """
        # ndx = 0
        last_res: List["IC_Residue"] = []
        last_ord_res: List["IC_Residue"] = []

        # atomCoordDict = {}
        akset = set()
        for res in self.chain.get_residues():
            # select only not hetero or accepted hetero
            if res.id[0] == " " or res.id[0] in IC_Residue.accept_resnames:
                this_res: List["IC_Residue"] = []
                if 2 == res.is_disordered() and not IC_Residue.no_altloc:
                    # print('disordered res:', res.is_disordered(), res)
                    for r in res.child_dict.values():
                        if self._add_residue(r, last_res, last_ord_res, verbose):
                            this_res.append(r.internal_coord)
                            akset.update(r.internal_coord.ak_set)
                else:
                    if self._add_residue(res, last_res, last_ord_res, verbose):
                        this_res.append(res.internal_coord)
                        akset.update(res.internal_coord.ak_set)

                if 0 < len(this_res):
                    self.ordered_aa_ic_list.extend(this_res)
                    last_ord_res = this_res

                last_res = this_res

        self.akset = akset
        self.initNCaCs = sorted(self.initNCaCs)

    def build_atomArray(self) -> None:
        """Build :class:`IC_Chain` numpy coordinate array from biopython atoms.

        See also :meth:`.init_edra` for more complete initialization of IC_Chain.

        Inputs:
            self.akset : set
                :class:`AtomKey` s in this chain

        Generates:
            self.AAsiz : int
                number of atoms in chain (len(akset))
            self.aktuple : AAsiz x AtomKeys
                sorted akset AtomKeys
            self.atomArrayIndex : [AAsiz] of int
                numerical index for each AtomKey in aktuple
            self.atomArrayValid : AAsiz x bool
                atomArray coordinates current with internal coordinates if True
            self.atomArray : AAsiz x np.float64[4]
                homogeneous atom coordinates; Biopython :class:`.Atom`
                coordinates are view into this array after execution
            rak_cache : dict
                lookup cache for AtomKeys for each residue

        """

        def setAtom(res, atm):
            ak = AtomKey(res.internal_coord, atm)
            try:
                ndx = self.atomArrayIndex[ak]
            except KeyError:
                return
            self.atomArray[ndx, 0:3] = atm.coord
            atm.coord = self.atomArray[ndx, 0:3]  # make view on atomArray
            self.atomArrayValid[ndx] = True
            self.bpAtomArray[ndx] = atm  # rtm

        def setResAtms(res):
            for atm in res.get_atoms():
                if atm.is_disordered():
                    if IC_Residue.no_altloc:
                        setAtom(res, atm.selected_child)
                    else:
                        for altAtom in atm.child_dict.values():
                            setAtom(res, altAtom)
                else:
                    setAtom(res, atm)

        self.AAsiz = len(self.akset)
        # sorted(akset) needed here for pdb atom serial number and to maintain
        # consistency between a2ic and i2ac
        self.aktuple = tuple(sorted(self.akset))
        self.atomArrayIndex = dict(zip(self.aktuple, range(self.AAsiz)))
        self.atomArrayValid = np.zeros(self.AAsiz, dtype=bool)
        self.atomArray = np.zeros((self.AAsiz, 4), dtype=np.float64)
        self.atomArray[:, 3] = 1.0
        self.bpAtomArray = [None] * self.AAsiz  # rtm

        for ric in self.ordered_aa_ic_list:
            setResAtms(ric.residue)
            if ric.akc == {}:  # pic file read
                ric._build_rak_cache()

    def build_edraArrays(self) -> None:
        """Build chain level hedra and dihedra arrays.

        Used by :meth:`init_edra` and :meth:`_hedraDict2chain`.  Should be
        private method but exposed for documentation.

        Inputs:
            self.dihedraLen : int
                number of dihedra needed
            self.hedraLen : int
                number of hedra needed
            self.AAsiz : int
                length of atomArray
            self.hedraNdx : dict
                maps hedron keys to range(hedraLen)
            self.dihedraNdx : dict
                maps dihedron keys to range(dihedraLen)
            self.hedra : dict
                maps Hedra keys to Hedra for chain
            self.atomArray : AAsiz x np.float64[4]
                homogeneous atom coordinates for chain
            self.atomArrayIndex : dict
                maps AtomKeys to atomArray
            self.atomArrayValid : AAsiz x bool
                indicates coord is up-to-date

        Generates:
            self.dCoordSpace : [2][dihedraLen][4][4]
                transforms to/from dihedron coordinate space
            self.dcsValid : dihedraLen x bool
                indicates dCoordSpace is current
            self.hAtoms : hedraLen x 3 x np.float64[4]
                atom coordinates in hCoordSpace
            self.hAtomsR : hedraLen x 3 x np.float64[4]
                hAtoms in reverse order (trading space for time)
            self.hAtoms_needs_update : hedraLen x bool
                indicates hAtoms, hAtoms current
            self.a2h_map : AAsiz x [int ...]
                maps atomArrayIndex to hedraNdx's with that atom
            self.a2ha_map : [hedraLen x 3]
                AtomNdx's in hedraNdx order
            self.h2aa : hedraLen x [int ...]
                maps hedraNdx to atomNdx's in hedron (reshaped later)
            Hedron.ndx : int
                self.hedraNdx value stored inside Hedron object
            self.dRev : dihedraLen x bool
                dihedron reversed if true
            self.dH1ndx, dH2ndx : [dihedraLen]
                hedraNdx's for 1st and 2nd hedra
            self.h1d_map : hedraLen x []
                hedraNdx -> [dihedra using hedron]
            Dihedron.h1key, h2key : [AtomKey ...]
                hedron keys for dihedron, reversed as needed
            Dihedron.hedron1, hedron2 : Hedron
                references inside dihedron to hedra
            Dihedron.ndx : int
                self.dihedraNdx info inside Dihedron object
            Dihedron.cst, rcst : np.float64p4][4]
                dCoordSpace references inside Dihedron
            self.a2da_map : [dihedraLen x 4]
                AtomNdx's in dihedraNdx order
            self.d2a_map : [dihedraLen x [4]]
                AtomNdx's for each dihedron (reshaped a2da_map)
            self.dFwd : bool
                dihedron is not Reversed if True
            self.a2d_map : AAsiz x [[dihedraNdx]
                [atom ndx 0-3 of atom in dihedron]], maps atom indexes to
                dihedra and atoms in them
            self.dAtoms_needs_update : dihedraLen x bool
                atoms in h1, h2 are current if False

        """
        # dihedra coord space
        self.dCoordSpace: np.ndarray = np.empty(
            (2, self.dihedraLen, 4, 4), dtype=np.float64
        )
        self.dcsValid: np.ndarray = np.zeros((self.dihedraLen), dtype=bool)

        # hedra atoms
        self.hAtoms: np.ndarray = np.zeros((self.hedraLen, 3, 4), dtype=np.float64)
        self.hAtoms[:, :, 3] = 1.0  # homogeneous
        self.hAtomsR: np.ndarray = np.copy(self.hAtoms)
        self.hAtoms_needs_update = np.full(self.hedraLen, True)

        # maps between hAtoms and atomArray
        a2ha_map = {}
        self.a2h_map = [[] for _ in range(self.AAsiz)]

        h2aa = [[] for _ in range(self.hedraLen)]
        for hk, hndx in self.hedraNdx.items():
            hstep = hndx * 3
            for i in range(3):
                ndx = self.atomArrayIndex[hk[i]]
                a2ha_map[hstep + i] = ndx
            self.hedra[hk].ndx = hndx
            for ak in self.hedra[hk].atomkeys:
                akndx = self.atomArrayIndex[ak]
                h2aa[hndx].append(akndx)
                self.a2h_map[akndx].append(hndx)

        self.a2ha_map = np.array(tuple(a2ha_map.values()))
        self.h2aa = np.array(h2aa)

        # dihedra atoms
        self.dAtoms: np.ndarray = np.empty((self.dihedraLen, 4, 4), dtype=np.float64)
        self.dAtoms[:, :, 3] = 1.0  # homogeneous
        self.a4_pre_rotation = np.empty((self.dihedraLen, 4))

        # maps between dAtoms and atomArray
        # hedra and dihedra
        # dihedra forward/reverse data
        a2da_map = {}
        a2d_map = [[[], []] for _ in range(self.AAsiz)]
        self.dRev: np.ndarray = np.zeros((self.dihedraLen), dtype=bool)

        self.dH1ndx = np.empty(self.dihedraLen, dtype=np.int64)
        self.dH2ndx = np.empty(self.dihedraLen, dtype=np.int64)
        self.h1d_map = [[] for _ in range(self.hedraLen)]

        self.id3_dh_index = {k[0:3]: [] for k in self.dihedraNdx.keys()}
        self.id32_dh_index = {k[1:4]: [] for k in self.dihedraNdx.keys()}

        for dk, dndx in self.dihedraNdx.items():
            # build map between atomArray and dAtoms
            dstep = dndx * 4
            did3 = dk[0:3]
            did32 = dk[1:4]
            d = self.dihedra[dk]
            for i in range(4):
                ndx = self.atomArrayIndex[dk[i]]
                a2da_map[dstep + i] = ndx
                a2d_map[ndx][0].append(dndx)
                a2d_map[ndx][1].append(i)

            try:
                d.h1key = did3
                d.h2key = did32
                h1ndx = self.hedraNdx[d.h1key]
            except KeyError:
                d.h1key = dk[2::-1]
                d.h2key = dk[3:0:-1]
                h1ndx = self.hedraNdx[d.h1key]
                self.dRev[dndx] = True
                d.reverse = True

            h2ndx = self.hedraNdx[d.h2key]
            d.hedron1 = self.hedra[d.h1key]
            d.hedron2 = self.hedra[d.h2key]
            self.dH1ndx[dndx] = h1ndx
            self.dH2ndx[dndx] = h2ndx
            self.h1d_map[h1ndx].append(dndx)

            d.ndx = dndx
            d.cst = self.dCoordSpace[0][dndx]
            d.rcst = self.dCoordSpace[1][dndx]
            self.id3_dh_index[did3].append(dk)
            self.id32_dh_index[did32].append(dk)

        self.a2da_map = np.array(tuple(a2da_map.values()))
        self.d2a_map = self.a2da_map.reshape(-1, 4)
        self.dFwd = self.dRev != True  # noqa: E712

        # manually create np.where(atom in dihedral)
        self.a2d_map = [(np.array(xi[0]), np.array(xi[1])) for xi in a2d_map]
        self.dAtoms_needs_update = np.full(self.dihedraLen, True)

    def _hedraDict2chain(
        self,
        hl12: Dict[str, float],
        ha: Dict[str, float],
        hl23: Dict[str, float],
        da: Dict[str, float],
        bfacs: Dict[str, float],
    ) -> None:
        """Generate chain numpy arrays from :func:`.read_PIC` dicts.

        On entry:
            * chain internal_coord has ordered_aa_ic_list built, akset;
            * residues have rnext, rprev, ak_set and di/hedra dicts initialised
            * Chain, residues do NOT have NCaC info, id3_dh_index
            * Di/hedra have cic, atomkeys set
            * Dihedra do NOT have valid reverse flag, h1/2 info

        """
        for ric in self.ordered_aa_ic_list:
            # log chain starts - beginning and after breaks
            # chain starts are only atom coords in pic files
            # assume valid pic files with all 3 of N, Ca, C coords
            initNCaC = []
            for atm in ric.residue.get_atoms():  # n.b. only few PIC spec atoms
                if 2 == atm.is_disordered():
                    if IC_Residue.no_altloc:
                        initNCaC.append(AtomKey(ric, atm.selected_child))
                    else:
                        for altAtom in atm.child_dict.values():
                            if altAtom.coord is not None:
                                initNCaC.append(AtomKey(ric, altAtom))
                elif atm.coord is not None:
                    initNCaC.append(AtomKey(ric, atm))
            if initNCaC != []:
                self.initNCaCs.append(tuple(initNCaC))

            # next residue NCaCKeys so can do per-residue assemble()
            ric.NCaCKey = []
            ric.NCaCKey.extend(
                ric.split_akl(
                    (AtomKey(ric, "N"), AtomKey(ric, "CA"), AtomKey(ric, "C"))
                )
            )
            ric._link_dihedra()

        # if STILL have no self.initNCacs, assume pic file w/o atoms and grab
        # from first residue
        if self.initNCaCs == []:
            ric = self.ordered_aa_ic_list[0]
            iNCaC = ric.split_akl(
                (AtomKey(ric, "N"), AtomKey(ric, "CA"), AtomKey(ric, "C"))
            )
            self.initNCaCs.extend(iNCaC)

        # set any supplied coordinates from biopython atoms
        # just loaded pic file so only start/chain break residues
        # will have atoms
        self.build_atomArray()

        self.initNCaCs = sorted(self.initNCaCs)
        # now create all biopython atoms for parent chain, setting coords to be
        # view on atomArray entry
        spNdx, icNdx, resnNdx, atmNdx, altlocNdx, occNdx = AtomKey.fields
        sn = None
        for ak, ndx in self.atomArrayIndex.items():
            res = ak.ric.residue  # read_PIC inits with IC_Residue
            atm, altloc = ak.akl[atmNdx], ak.akl[altlocNdx]
            occ = 1.00 if ak.akl[occNdx] is None else float(ak.akl[occNdx])
            bfac = bfacs.get(ak.id, 0.0)
            sn = sn + 1 if sn is not None else ndx + 1
            bpAtm = None
            if res.has_id(atm):
                bpAtm = res[atm]
            if bpAtm is None or (
                2 == bpAtm.is_disordered() and not bpAtm.disordered_has_id(altloc)
            ):
                newAtom = Atom(
                    atm,
                    self.atomArray[ndx][0:3],  # init as view on atomArray
                    bfac,
                    occ,
                    (" " if altloc is None else altloc),
                    atm,
                    sn,
                    atm[0],
                )
                if bpAtm is None:
                    if altloc is None:
                        res.add(newAtom)
                    else:
                        disordered_atom = DisorderedAtom(atm)
                        res.add(disordered_atom)
                        disordered_atom.disordered_add(newAtom)
                        res.flag_disordered()
                else:
                    bpAtm.disordered_add(newAtom)

            else:
                if 2 == bpAtm.is_disordered() and bpAtm.disordered_has_id(altloc):
                    bpAtm.disordered_select(altloc)

                bpAtm.set_bfactor(bfac)
                bpAtm.set_occupancy(occ)
                sn = bpAtm.get_serial_number()

        # hedra
        # dicts sorted on creation by init_edra and maintained by write_PIC
        # python 3.7 minimum for Biopython as of 6 sept 2021 PR #3714
        self.hedraLen = len(ha)
        self.hedraL12 = np.fromiter(hl12.values(), dtype=np.float64)
        self.hedraAngle = np.fromiter(ha.values(), dtype=np.float64)
        self.hedraL23 = np.fromiter(hl23.values(), dtype=np.float64)
        self.hedraNdx = dict(zip(sorted(ha.keys()), range(self.hedraLen)))

        # dihedra
        self.dihedraLen = len(da)
        self.dihedraAngle = np.fromiter(da.values(), dtype=np.float64)
        self.dihedraAngleRads = np.deg2rad(self.dihedraAngle)
        self.dihedraNdx = dict(zip(sorted(da.keys()), range(self.dihedraLen)))

        self.build_edraArrays()

    # @profile
    def assemble_residues(self, verbose: bool = False) -> None:
        """Generate atom coords from internal coords (vectorised).

        This is the 'Numpy parallel' version of :meth:`.assemble_residues_ser`.

        Starting with dihedra already formed by :meth:`.init_atom_coords`, transform
        each from dihedron local coordinate space into protein chain coordinate
        space.  Iterate until all dependencies satisfied.

        Does not update :data:`dCoordSpace` as :meth:`assemble_residues_ser`
        does.  Call :meth:`.update_dCoordSpace` if needed.  Faster to do in
        single operation once all atom coordinates finished.

        :param bool verbose: default False.
            Report number of iterations to compute changed dihedra

        generates:
            self.dSet: AAsiz x dihedraLen x 4
                maps atoms in dihedra to atomArray
            self.dSetValid : [dihedraLen][4] of bool
                map of valid atoms into dihedra to detect 3 or 4 atoms valid

        Output coordinates written to :data:`atomArray`.  Biopython
        :class:`Bio.PDB.Atom` coordinates are a view on this data.
        """
        # dihedron atom positions of chain atom ndxs, maps atomArray to dihedra
        a2da_map = self.a2da_map  # 8468 x int
        # each chain atom to list of [dihedron], [dihedron_position]
        a2d_map = self.a2d_map  # 2000 x ([int], [int])
        # every dihedron atom to chain atoms
        d2a_map = self.d2a_map  # 2117 x [4] ints
        # all chain atoms
        atomArray = self.atomArray  # 2000
        # bool markers for chain atoms with valid coordinates
        atomArrayValid = self.atomArrayValid  # 2000

        # complete array of dihedra atoms
        dAtoms = self.dAtoms  # 2117 x [4][4] float
        # coordinate space transformations optionally supplied
        dCoordSpace1 = self.dCoordSpace[1]
        dcsValid = self.dcsValid

        # dSet is 4-atom arrays for every dihedral, multiple copies of
        # many atoms as the dihedra overlap
        self.dSet = atomArray[a2da_map].reshape(-1, 4, 4)
        dSet = self.dSet
        # dSetValid indicates accurate atom positions in each dSet dihedral
        self.dSetValid = atomArrayValid[a2da_map].reshape(-1, 4)
        dSetValid = self.dSetValid

        # clear any transforms for dihedrals with outdated atoms
        workSelector = (dSetValid == self._dihedraOK).all(axis=1)

        self.dcsValid[np.logical_not(workSelector)] = False

        dihedraWrk = None
        if verbose:
            dihedraWrk = workSelector.size - workSelector.sum()

        # mask for dihedral with 3 valid atoms in dSet, ready to be processed:
        targ = IC_Chain._dihedraSelect
        # select the dihedrals ready for processing
        workSelector = (dSetValid == targ).all(axis=1)

        loopCount = 0
        while np.any(workSelector):
            # indexes of dihedra in dset to update
            workNdxs = np.where(workSelector)
            # subset of dihedra to update
            workSet = dSet[workSelector]
            # will update coordinates of 4th atom in each workSet dihedron
            updateMap = d2a_map[workNdxs, 3][0]

            # get all coordSpace transforms
            if np.all(dcsValid[workSelector]):
                cspace = dCoordSpace1[workSelector]
            else:
                cspace = multi_coord_space(workSet, np.sum(workSelector), True)[1]

            # generate new coords for 4th atoms in workSet dihedra
            initCoords = dAtoms[workSelector].reshape(-1, 4, 4)
            atomArray[updateMap] = np.einsum("ijk,ik->ij", cspace, initCoords[:, 3])

            # mark new computed atom positions as valid
            atomArrayValid[updateMap] = True

            # prep for next iteration
            workSelector[:] = False
            for a in updateMap:
                # copy new atom positions into dihedra atom array
                dSet[a2d_map[a]] = atomArray[a]
                # build new workSelector from only updated dihedra
                adlist = a2d_map[a]
                for d in adlist[0]:
                    dvalid = atomArrayValid[d2a_map[d]]
                    workSelector[d] = (dvalid == targ).all()

            loopCount += 1

        if verbose:
            cid = self.chain.full_id
            print(
                f"{cid[0]} {cid[2]} coordinates for {dihedraWrk} dihedra"
                f" updated in {loopCount} iterations"
            )

    def assemble_residues_ser(
        self,
        verbose: bool = False,
        start: Optional[int] = None,
        fin: Optional[int] = None,
    ) -> None:
        """Generate IC_Residue atom coords from internal coordinates (serial).

        See :meth:`.assemble_residues` for 'numpy parallel' version.

        Filter positions between start and fin if set, find appropriate start
        coordinates for each residue and pass to :meth:`.assemble`

        :param bool verbose: default False.
            Describe runtime problems
        :param int start,fin: default None.
            Sequence position for begin, end of subregion to generate coords
            for.
        """
        self.dcsValid[:] = False

        for ric in self.ordered_aa_ic_list:
            # :  # clear and skip if outside start ... fin
            if (fin and fin < ric.residue.id[1]) or (
                start and start > ric.residue.id[1]
            ):
                ric.ak_set = None
                ric.akc = None
                ric.residue.child_dict = {}
                ric.residue.child_list = []
                continue

            atom_coords = ric.assemble(verbose=verbose)
            if atom_coords:
                ric.ak_set = set(atom_coords.keys())

    def init_edra(self, verbose: bool = False) -> None:
        """Create chain and residue di/hedra structures, arrays, atomArray.

        Inputs:
            self.ordered_aa_ic_list : list of IC_Residue
        Generates:
            * edra objects, self.di/hedra (executes :meth:`._create_edra`)
            * atomArray and support (executes :meth:`.build_atomArray`)
            * self.hedraLen : number of hedra in structure
            * hedraL12 : numpy arrays for lengths, angles (empty)
            * hedraAngle ..
            * hedraL23 ..
            * self.hedraNdx : dict mapping hedrakeys to hedraL12 etc
            * self.dihedraLen : number of dihedra in structure
            * dihedraAngle ..
            * dihedraAngleRads : np arrays for angles (empty)
            * self.dihedraNdx : dict mapping dihedrakeys to dihedraAngle
        """
        if self.ordered_aa_ic_list[0].hedra == {}:
            for ric in self.ordered_aa_ic_list:
                # build di/hedra objects in chain arrays
                ric._create_edra(verbose=verbose)

        if not hasattr(self, "atomArrayValid"):
            self.build_atomArray()  # ric.a2ic added gly CBs to akset

        if not hasattr(self, "hedraLen"):
            # hedra
            self.hedraLen = len(self.hedra)
            self.hedraL12 = np.empty((self.hedraLen), dtype=np.float64)
            self.hedraAngle = np.empty((self.hedraLen), dtype=np.float64)
            self.hedraL23 = np.empty((self.hedraLen), dtype=np.float64)

            # python3.7 sorted dicts
            self.hedraNdx = dict(zip(sorted(self.hedra.keys()), range(len(self.hedra))))

            # dihedra
            self.dihedraLen = len(self.dihedra)
            self.dihedraAngle = np.empty(self.dihedraLen)
            self.dihedraAngleRads = np.empty(self.dihedraLen)

            self.dihedraNdx = dict(
                zip(sorted(self.dihedra.keys()), range(self.dihedraLen))
            )

        if not hasattr(self, "hAtoms_needs_update"):
            self.build_edraArrays()

    # @profile
    def init_atom_coords(self) -> None:
        """Set chain level di/hedra initial coords from angles and distances.

        Initializes atom coordinates in local coordinate space for hedra and
        dihedra, will be transformed appropriately later by :data:`dCoordSpace`
        matrices for assembly.
        """
        # dbg = True
        if not np.all(self.dAtoms_needs_update):
            self.dAtoms_needs_update |= (self.hAtoms_needs_update[self.dH1ndx]) | (
                self.hAtoms_needs_update[self.dH2ndx]
            )
            self.dcsValid &= np.logical_not(self.dAtoms_needs_update)

        # dihedra full size masks:
        mdFwd = self.dFwd & self.dAtoms_needs_update
        mdRev = self.dRev & self.dAtoms_needs_update

        # update size masks
        udFwd = self.dFwd[self.dAtoms_needs_update]
        udRev = self.dRev[self.dAtoms_needs_update]
        """
        if dbg:
            print("mdFwd", mdFwd[0:10])
            print("mdRev", mdRev[0:10])
            print("udFwd", udFwd[0:10])
            print("udRev", udRev[0:10])
        """

        if np.any(self.hAtoms_needs_update):
            # hedra inital coords

            # sar = supplementary angle radian: angles which add to 180
            sar = np.deg2rad(180.0 - self.hedraAngle[self.hAtoms_needs_update])  # angle
            sinSar = np.sin(sar)
            cosSarN = np.cos(sar) * -1
            """
            if dbg:
                print("sar", sar[0:10])
            """
            # a2 is len3 up from a2 on Z axis, X=Y=0
            self.hAtoms[:, 2, 2][self.hAtoms_needs_update] = self.hedraL23[
                self.hAtoms_needs_update
            ]

            # a0 X is sin( sar ) * len12
            self.hAtoms[:, 0, 0][self.hAtoms_needs_update] = (
                sinSar * self.hedraL12[self.hAtoms_needs_update]
            )

            # a0 Z is -(cos( sar ) * len12)
            # (assume angle always obtuse, so a0 is in -Z)
            self.hAtoms[:, 0, 2][self.hAtoms_needs_update] = (
                cosSarN * self.hedraL12[self.hAtoms_needs_update]
            )
            """
            if dbg:
                print("hAtoms_needs_update", self.hAtoms_needs_update[0:10])
                print("self.hAtoms", self.hAtoms[0:10])
            """
            # same again but 'reversed' : a0 on Z axis, a1 at origin, a2 in -Z

            # a0r is len12 up from a1 on Z axis, X=Y=0
            self.hAtomsR[:, 0, 2][self.hAtoms_needs_update] = self.hedraL12[
                self.hAtoms_needs_update
            ]
            # a2r X is sin( sar ) * len23
            self.hAtomsR[:, 2, 0][self.hAtoms_needs_update] = (
                sinSar * self.hedraL23[self.hAtoms_needs_update]
            )
            # a2r Z is -(cos( sar ) * len23)
            self.hAtomsR[:, 2, 2][self.hAtoms_needs_update] = (
                cosSarN * self.hedraL23[self.hAtoms_needs_update]
            )
            """
            if dbg:
                print("self.hAtomsR", self.hAtomsR[0:10])
            """
            self.hAtoms_needs_update[...] = False

            # dihedra parts other than dihedral angle

            dhlen = np.sum(self.dAtoms_needs_update)  # self.dihedraLen

            # only 4th atom takes work:
            # pick 4th atom based on rev flag
            self.a4_pre_rotation[mdRev] = self.hAtoms[self.dH2ndx, 0][mdRev]
            self.a4_pre_rotation[mdFwd] = self.hAtomsR[self.dH2ndx, 2][mdFwd]

            # numpy multiply, add operations below intermediate array but out=
            # not working with masking:
            self.a4_pre_rotation[:, 2][self.dAtoms_needs_update] = np.multiply(
                self.a4_pre_rotation[:, 2][self.dAtoms_needs_update], -1
            )  # a4 to +Z

            a4shift = np.empty(dhlen)
            a4shift[udRev] = self.hedraL23[self.dH2ndx][mdRev]  # len23
            a4shift[udFwd] = self.hedraL12[self.dH2ndx][mdFwd]  # len12

            self.a4_pre_rotation[:, 2][self.dAtoms_needs_update] = np.add(
                self.a4_pre_rotation[:, 2][self.dAtoms_needs_update],
                a4shift,
            )  # so a2 at origin
            """
            if dbg:
                print("dhlen", dhlen)
                print("a4shift", a4shift[0:10])
                print("a4_pre_rotation", self.a4_pre_rotation[0:10])
            """
            # now build dihedra initial coords

            dH1atoms = self.hAtoms[self.dH1ndx]  # fancy indexing so
            dH1atomsR = self.hAtomsR[self.dH1ndx]  # these copy not view

            self.dAtoms[:, :3][mdFwd] = dH1atoms[mdFwd]
            self.dAtoms[:, :3][mdRev] = dH1atomsR[:, 2::-1][mdRev]

            """
            if dbg:
                print("dH1atoms", dH1atoms[0:10])
                print("dH1atosR", dH1atomsR[0:10])
                print("dAtoms", self.dAtoms[0:10])
            """

        # build rz rotation matrix for dihedral angle
        """
        if dbg:
            print("dangle-rads", self.dihedraAngleRads[0:10])
        """
        rz = multi_rot_Z(self.dihedraAngleRads[self.dAtoms_needs_update])

        a4rot = np.matmul(
            rz,
            self.a4_pre_rotation[self.dAtoms_needs_update][:].reshape(-1, 4, 1),
        ).reshape(-1, 4)

        self.dAtoms[:, 3][mdFwd] = a4rot[udFwd]  # [self.dFwd]
        self.dAtoms[:, 3][mdRev] = a4rot[udRev]  # [self.dRev]

        """
        if dbg:
            print("rz", rz[0:3])
            print("dAtoms", self.dAtoms[0:10])
        """

        self.dAtoms_needs_update[...] = False

        # can't start assembly if initial NCaC is not valid, so copy from
        # hAtoms if needed

        """
        if dbg:
            print("initNCaCs", self.initNCaCs)
        """
        for iNCaC in self.initNCaCs:
            invalid = True
            if np.all(self.atomArrayValid[[self.atomArrayIndex[ak] for ak in iNCaC]]):
                invalid = False

            if invalid:
                hatoms = self.hAtoms[self.hedraNdx[iNCaC]]
                for i in range(3):
                    andx = self.atomArrayIndex[iNCaC[i]]
                    self.atomArray[andx] = hatoms[i]
                    self.atomArrayValid[andx] = True
            """
            if dbg:
                hatoms = self.hAtoms[self.hedraNdx[iNCaC]]
                print("hedraNdx iNCaC", self.hedraNdx[iNCaC])
                print("hatoms", hatoms)
            """

    def update_dCoordSpace(self, workSelector: Optional[np.ndarray] = None) -> None:
        """Compute/update coordinate space transforms for chain dihedra.

        Requires all atoms updated so calls :meth:`.assemble_residues`
        (returns immediately if all atoms already assembled).

        :param [bool] workSelector:
            Optional mask to select dihedra for update
        """
        if workSelector is None:
            self.assemble_residues()  # update atoms, fast if nothing to do
            workSelector = np.logical_not(self.dcsValid)
        workSet = self.dSet[workSelector]
        self.dCoordSpace[:, workSelector] = multi_coord_space(
            workSet, np.sum(workSelector), True
        )

        self.dcsValid[workSelector] = True

    def propagate_changes(self) -> None:
        """Track through di/hedra to invalidate dependent atoms."""
        # cs : chain segment
        # csStart, csNext : AtomArray indexes for chain segment
        # process each chain segment
        csNdx = 0
        csLen = len(self.initNCaCs)
        atmNdx = AtomKey.fields.atm
        posNdx = AtomKey.fields.respos
        done = set()

        while csNdx < csLen:  # iterate over chain starts
            startAK = self.initNCaCs[csNdx][0]
            csStart = self.atomArrayIndex[startAK]
            csnTry = csNdx + 1
            # set csNext to be atomArray index of segment end
            if csLen == csnTry:
                csNext = self.AAsiz  # last segment to end of atomArray
            else:  # this segment to next chain start
                finAK = self.initNCaCs[csnTry][0]
                csNext = self.atomArrayIndex[finAK]

            for andx in range(csStart, csNext):
                if not self.atomArrayValid[andx]:
                    ak = self.aktuple[andx]
                    atm = ak.akl[atmNdx]
                    pos = ak.akl[posNdx]  # sequence position = residue number
                    if atm in ("N", "CA", "C"):
                        # backbone moved so all to next start moved
                        self.atomArrayValid[andx:csNext] = False
                        # and done with this invalid_atom_ndxs segment
                        break
                    elif pos not in done and atm != "H":
                        # H is terminal so ignore, not effect subsequent atoms
                        # O is terminal but used to locate CB
                        # atomArray is sorted, sidechain atoms follow backbone
                        for i in range(andx, csNext):
                            if self.aktuple[i].akl[posNdx] == pos:
                                self.atomArrayValid[i] = False
                            else:
                                # done with residue sidechain when find next
                                # seq pos so need not go to fin
                                break
                        done.add(pos)
            csNdx += 1

    # @profile
    def internal_to_atom_coordinates(
        self,
        verbose: bool = False,
        start: Optional[int] = None,
        fin: Optional[int] = None,
    ) -> None:
        """Process IC data to Residue/Atom coords.

        :param bool verbose: default False.
            Describe runtime problems
        :param int start,fin:
            Optional sequence positions for begin, end of subregion
            to process.

        .. note::
            Setting start or fin activates serial :meth:`.assemble_residues_ser`
            instead of (Numpy parallel) :meth:`.assemble_residues`.
            Start C-alpha will be at origin.

        .. seealso::
            :data:`ParallelAssembleResidues`

        """
        if not hasattr(self, "dAtoms_needs_update"):
            return  # escape on no data to process

        # if verbose:
        #    for ric in self.ordered_aa_ic_list:
        #        if not hasattr(ric, "NCaCKey"):
        #            print(
        #                f"no assembly for {ric} due to missing N, Ca"
        #                " and/or C atoms"
        #            )

        if IC_Chain.ParallelAssembleResidues and not (start or fin):
            self.propagate_changes()
            self.init_atom_coords()  # compute initial di/hedra coords
            # transform init di/hedra to chain coord space
            self.assemble_residues(verbose=verbose)

            if verbose and not np.all(self.atomArrayValid):
                dSetValid = self.atomArrayValid[self.a2da_map].reshape(-1, 4)
                for ric in self.ordered_aa_ic_list:
                    for d in ric.dihedra.values():
                        if not dSetValid[d.ndx].all():
                            print(
                                "missing coordinates for chain "
                                f"{ric.cic.chain.id} {ric.pretty_str()} "
                                f"dihedral: {d.id}"
                            )
        else:
            if start:  # set initNCaC tag to build from
                for ric in self.ordered_aa_ic_list:
                    if start != ric.residue.id[1]:
                        continue
                    iNCaC = ric.split_akl(
                        (
                            AtomKey(ric, "N"),
                            AtomKey(ric, "CA"),
                            AtomKey(ric, "C"),
                        )
                    )
                    self.initNCaCs.extend(iNCaC)

            self.init_atom_coords()  # compute initial di/hedra coords
            self.assemble_residues_ser(
                verbose=verbose, start=start, fin=fin
            )  # internal to XYZ coordinates

    # @profile
    def atom_to_internal_coordinates(self, verbose: bool = False) -> None:
        """Calculate dihedrals, angles, bond lengths for Atom data.

        Generates atomArray (through init_edra), value arrays for hedra and
        dihedra, and coordinate space transforms for dihedra.

        Generates Gly C-beta if specified, see :data:`IC_Residue.gly_Cbeta`

        :param bool verbose: default False.
            describe runtime problems
        """
        if self.ordered_aa_ic_list == []:
            return  # escape on no data to process

        self.init_edra(verbose=verbose)

        if self.dihedra == {}:
            return  # escape if no hedra loaded for this chain

        # compute all hedra parameters with law of cosines on 3 atom coords
        ha = self.atomArray[self.a2ha_map].reshape(-1, 3, 4)
        self.hedraL12 = np.linalg.norm(ha[:, 0] - ha[:, 1], axis=1)
        self.hedraL23 = np.linalg.norm(ha[:, 1] - ha[:, 2], axis=1)
        h_a0a2 = np.linalg.norm(ha[:, 0] - ha[:, 2], axis=1)
        np.rad2deg(
            np.arccos(
                (
                    np.square(self.hedraL12)
                    + np.square(self.hedraL23)
                    - np.square(h_a0a2)
                )
                / (2 * self.hedraL12 * self.hedraL23)
            ),
            out=self.hedraAngle,
        )

        # now process dihedra
        dha = self.atomArray[self.a2da_map].reshape(-1, 4, 4)

        # develop coord_space matrix for 1st 3 atoms of dihedra:
        # note use of [...] to modify in place, dihedra cst, rcst remain valid
        self.dCoordSpace[...] = multi_coord_space(dha, self.dihedraLen, True)
        self.dcsValid[:] = True

        # now put atom 4 into that coordinate space
        do4 = np.matmul(self.dCoordSpace[0], dha[:, 3].reshape(-1, 4, 1)).reshape(-1, 4)
        # and read dihedral as azimuth
        np.arctan2(do4[:, 1], do4[:, 0], out=self.dihedraAngleRads)
        np.rad2deg(self.dihedraAngleRads, out=self.dihedraAngle)

        if hasattr(self, "gcb"):
            self._spec_glyCB()

    def _spec_glyCB(self) -> None:
        """Populate values for Gly C-beta."""
        Ca_Cb_Len = 1.53363
        if hasattr(self, "scale"):  # used for openscad output
            Ca_Cb_Len *= self.scale  # type: ignore

        for gcbd in self.gcb.values():  # gcb dict created by _create_edra
            cbak = gcbd[3]
            self.atomArrayValid[self.atomArrayIndex[cbak]] = False
            ric = cbak.ric
            rN, rCA, rC, rO = (
                ric.rak("N"),
                ric.rak("CA"),
                ric.rak("C"),
                ric.rak("O"),
            )
            gCBd = self.dihedra[gcbd]
            dndx = gCBd.ndx
            # generated dihedron is O-Ca-C-Cb
            # hedron2 is reversed: Cb-Ca-C (also h1 reversed: C-Ca-O)
            h2ndx = gCBd.hedron2.ndx
            self.hedraL12[h2ndx] = Ca_Cb_Len
            self.hedraAngle[h2ndx] = 110.17513
            self.hedraL23[h2ndx] = self.hedraL12[self.hedraNdx[(rCA, rC, rO)]]

            self.hAtoms_needs_update[gCBd.hedron2.ndx] = True
            for ak in gCBd.hedron2.atomkeys:
                self.atomArrayValid[self.atomArrayIndex[ak]] = False

            refval = self.dihedra.get((rN, rCA, rC, rO), None)
            if refval:
                angl = 122.68219 + self.dihedraAngle[refval.ndx]
                self.dihedraAngle[dndx] = angl if (angl <= 180.0) else angl - 360.0
            else:
                self.dihedraAngle[dndx] = 120

    @staticmethod
    def _write_mtx(fp: TextIO, mtx: np.array) -> None:
        fp.write("[ ")
        rowsStarted = False
        for row in mtx:
            if rowsStarted:
                fp.write(", [ ")
            else:
                fp.write("[ ")
                rowsStarted = True
            colsStarted = False
            for col in row:
                if colsStarted:
                    fp.write(", " + str(col))
                else:
                    fp.write(str(col))
                    colsStarted = True
            fp.write(" ]")  # close row
        fp.write(" ]")

    @staticmethod
    def _writeSCAD_dihed(
        fp: TextIO, d: "Dihedron", hedraNdx: Dict, hedraSet: Set[EKT]
    ) -> None:
        fp.write(
            "[ {:9.5f}, {}, {}, {}, ".format(
                d.angle,
                hedraNdx[d.h1key],
                hedraNdx[d.h2key],
                (1 if d.reverse else 0),
            )
        )
        fp.write(
            f"{0 if d.h1key in hedraSet else 1}, "
            f"{0 if d.h2key in hedraSet else 1}, "
        )
        fp.write(
            "    // {} [ {} -- {} ] {}\n".format(
                d.id,
                d.hedron1.id,
                d.hedron2.id,
                ("reversed" if d.reverse else ""),
            )
        )
        fp.write("        ")
        IC_Chain._write_mtx(fp, d.rcst)
        fp.write(" ]")  # close residue array of dihedra entry

    def _write_SCAD(self, fp: TextIO, backboneOnly: bool, start=None, fin=None) -> None:
        """Write self to file fp as OpenSCAD data matrices.

        See `OpenSCAD <https://www.openscad.org>`_.
        Works with :func:`.write_SCAD` and embedded OpenSCAD routines therein.
        """
        fp.write(f'   "{self.chain.id}", // chain id\n')

        # generate dict for all hedra to eliminate redundant references
        hedra = {}
        for ric in self.ordered_aa_ic_list:
            respos = ric.residue.id[1]
            if start is not None and respos < start - 1:
                # start-1 because rprev has some hedra for residue r
                continue
            if fin is not None and respos > fin:
                continue
            for k, h in ric.hedra.items():
                hedra[k] = h

        atomSet: Set[AtomKey] = set()
        bondDict: Dict = {}  # set()
        hedraSet: Set[EKT] = set()
        ndx = 0
        hedraNdx = {}

        for hk in sorted(hedra):
            hedraNdx[hk] = ndx
            ndx += 1

        # write residue dihedra table

        fp.write("   [  // residue array of dihedra")
        resNdx = {}
        dihedraNdx = {}
        ndx = 0
        chnStarted = False

        for ric in self.ordered_aa_ic_list:
            respos = ric.residue.id[1]
            if start is not None and respos < start:
                continue
            if fin is not None and respos > fin:
                continue

            if "O" not in ric.akc:
                if ric.lc != "G" and ric.lc != "A":
                    print(
                        "Unable to generate complete sidechain for "
                        f"{ric} {ric.lc} missing O atom"
                    )
            resNdx[ric] = ndx
            if chnStarted:
                fp.write("\n     ],")
            else:
                chnStarted = True
            fp.write(
                "\n     [ // "
                + str(ndx)
                + " : "
                + str(ric.residue.id)
                + " "
                + ric.lc
                + " backbone\n"
            )
            ndx += 1

            # assemble with no start position, return transform matrices
            ric.clear_transforms()

            # compute residue atom coords for no start position
            # dump results because only want rcst
            # IC_Chain.adbg = True
            ric.assemble(resetLocation=True)
            # IC_Chain.adbg = False

            ndx2 = 0
            started = False
            for i in range(1 if backboneOnly else 2):
                if i == 1:
                    cma = "," if started else ""
                    fp.write(
                        f"{cma}\n       // {ric.residue.id!s} {ric.lc} sidechain\n"
                    )
                started = False
                for dk, d in sorted(ric.dihedra.items()):
                    if d.h2key in hedraNdx and (
                        (i == 0 and d.is_backbone()) or (i == 1 and not d.is_backbone())
                    ):
                        if d.cic.dcsValid[d.ndx]:
                            if started:
                                fp.write(",\n")
                            else:
                                started = True
                            fp.write("      ")
                            IC_Chain._writeSCAD_dihed(fp, d, hedraNdx, hedraSet)
                            dihedraNdx[dk] = ndx2
                            hedraSet.add(d.h1key)
                            hedraSet.add(d.h2key)
                            ndx2 += 1
                        else:
                            print(
                                f"Atom missing for {d.id3}-{d.id32}, OpenSCAD"
                                f" chain may be discontiguous"
                            )
        fp.write("   ],")  # end of residue entry dihedra table
        fp.write("\n  ],\n")  # end of all dihedra table

        # write hedra table

        fp.write("   [  //hedra\n")
        for hk in sorted(hedra):
            hed = hedra[hk]
            fp.write("     [ ")
            fp.write(
                "{:9.5f}, {:9.5f}, {:9.5f}".format(
                    set_accuracy_95(hed.len12),
                    set_accuracy_95(hed.angle),
                    set_accuracy_95(hed.len23),
                )
            )
            atom_str = ""  # atom and bond state
            atom_done_str = ""  # create each only once
            akndx = 0
            for ak in hed.atomkeys:
                atm = ak.akl[AtomKey.fields.atm]
                res = ak.akl[AtomKey.fields.resname]
                # try first for generic backbone/Cbeta atoms
                ab_state_res = residue_atom_bond_state["X"]
                ab_state = ab_state_res.get(atm, None)
                if "H" == atm[0]:
                    ab_state = "Hsb"
                if ab_state is None:
                    # not found above, must be sidechain atom
                    ab_state_res = residue_atom_bond_state.get(res, None)
                    if ab_state_res is not None:
                        ab_state = ab_state_res.get(atm, "")
                    else:
                        ab_state = ""
                atom_str += ', "' + ab_state + '"'

                if ak in atomSet:
                    atom_done_str += ", 0"
                elif hk in hedraSet:
                    if (
                        hasattr(hed, "flex_female_1") or hasattr(hed, "flex_male_1")
                    ) and akndx != 2:
                        if akndx == 0:
                            atom_done_str += ", 0"
                        elif akndx == 1:
                            atom_done_str += ", 1"
                            atomSet.add(ak)
                    elif (
                        hasattr(hed, "flex_female_2") or hasattr(hed, "flex_male_2")
                    ) and akndx != 0:
                        if akndx == 2:
                            atom_done_str += ", 0"
                        elif akndx == 1:
                            atom_done_str += ", 1"
                            atomSet.add(ak)
                    else:
                        atom_done_str += ", 1"
                        atomSet.add(ak)
                else:
                    atom_done_str += ", 0"
                akndx += 1
            fp.write(atom_str)
            fp.write(atom_done_str)

            # specify bond options

            bond = []
            bond.append(hed.atomkeys[0].id + "-" + hed.atomkeys[1].id)
            bond.append(hed.atomkeys[1].id + "-" + hed.atomkeys[2].id)
            b0 = True
            for b in bond:
                wstr = ""
                if b in bondDict and bondDict[b] == "StdBond":
                    wstr = ", 0"
                elif hk in hedraSet:
                    bondType = "StdBond"
                    if b0:
                        if hasattr(hed, "flex_female_1"):
                            bondType = "FemaleJoinBond"
                        elif hasattr(hed, "flex_male_1"):
                            bondType = "MaleJoinBond"
                        elif hasattr(hed, "skinny_1"):
                            bondType = "SkinnyBond"
                        elif hasattr(hed, "hbond_1"):
                            bondType = "HBond"
                    else:
                        if hasattr(hed, "flex_female_2"):
                            bondType = "FemaleJoinBond"
                        elif hasattr(hed, "flex_male_2"):
                            bondType = "MaleJoinBond"
                        # elif hasattr(hed, 'skinny_2'):  # unused
                        #     bondType = 'SkinnyBond'
                        elif hasattr(hed, "hbond_2"):
                            bondType = "HBond"
                    if b in bondDict:
                        bondDict[b] = "StdBond"
                    else:
                        bondDict[b] = bondType
                    wstr = ", " + str(bondType)
                else:
                    wstr = ", 0"
                fp.write(wstr)
                b0 = False
            akl = hed.atomkeys[0].akl
            fp.write(
                ', "'
                + akl[AtomKey.fields.resname]
                + '", '
                + akl[AtomKey.fields.respos]
                + ', "'
                + hed.e_class
                + '"'
            )
            fp.write(" ], // " + str(hk) + "\n")
        fp.write("   ],\n")  # end of hedra table

        # write chain table

        self.atomArrayValid[:] = False
        self.internal_to_atom_coordinates()

        fp.write("\n[  // chain - world transform for each residue\n")
        chnStarted = False
        for ric in self.ordered_aa_ic_list:
            # rtm handle start / end
            respos = ric.residue.id[1]
            if start is not None and respos < start:
                continue
            if fin is not None and respos > fin:
                continue
            for k, h in ric.hedra.items():
                hedra[k] = h

            for NCaCKey in sorted(ric.NCaCKey):  # type: ignore
                mtr = None
                if 0 < len(ric.rprev):
                    acl = [self.atomArray[self.atomArrayIndex[ak]] for ak in NCaCKey]
                    mt, mtr = coord_space(acl[0], acl[1], acl[2], True)
                else:
                    mtr = np.identity(4, dtype=np.float64)
                if chnStarted:
                    fp.write(",\n")
                else:
                    chnStarted = True
                fp.write("     [ " + str(resNdx[ric]) + ', "' + str(ric.residue.id[1]))
                fp.write(ric.lc + '", //' + str(NCaCKey) + "\n")
                IC_Chain._write_mtx(fp, mtr)
                fp.write(" ]")
        fp.write("\n   ]\n")

    def distance_plot(
        self, filter: Optional[Union[np.ndarray, None]] = None
    ) -> np.ndarray:
        """Generate 2D distance plot from atomArray.

        Default is to calculate distances for all atoms.  To generate the
        classic C-alpha distance plot, pass a boolean mask array like::

            atmNameNdx = internal_coords.AtomKey.fields.atm
            CaSelect = [
                atomArrayIndex.get(k)
                for k in atomArrayIndex.keys()
                if k.akl[atmNameNdx] == "CA"
            ]
            plot = cic.distance_plot(CaSelect)

        Alternatively, this will select all backbone atoms::

            backboneSelect = [
                atomArrayIndex.get(k)
                for k in atomArrayIndex.keys()
                if k.is_backbone()
            ]

        :param [bool] filter: restrict atoms for calculation

        .. seealso::
            :meth:`.distance_to_internal_coordinates`, which requires the
            default all atom distance plot.

        """
        if filter is None:
            atomSet = self.atomArray
        else:
            atomSet = self.atomArray[filter]

        # create distance matrix without scipy
        # see https://jbencook.com/pairwise-distance-in-numpy/
        return np.linalg.norm(atomSet[:, None, :] - atomSet[None, :, :], axis=-1)

    def dihedral_signs(self) -> np.ndarray:
        """Get sign array (+1/-1) for each element of chain dihedraAngle array.

        Required for :meth:`.distance_to_internal_coordinates`
        """
        return np.sign(self.dihedraAngle)

    def distplot_to_dh_arrays(
        self, distplot: np.ndarray, dihedra_signs: np.ndarray
    ) -> None:
        """Load di/hedra distance arays from distplot.

        Fill :class:`IC_Chain` arrays hedraL12, L23, L13 and dihedraL14
        distance value arrays from input distplot, dihedra_signs array from
        input dihedra_signs.  Distplot and di/hedra distance arrays must index
        according to AtomKey mappings in :class:`IC_Chain` .hedraNdx and .dihedraNdx
        (created in :meth:`IC_Chain.init_edra`)

        Call :meth:`atom_to_internal_coordinates` (or at least :meth:`init_edra`)
        to generate a2ha_map and d2a_map before running this.

        Explcitly removed from :meth:`.distance_to_internal_coordinates` so
        user may populate these chain di/hedra arrays by other
        methods.
        """
        ha = self.a2ha_map.reshape(-1, 3)
        self.hedraL12 = distplot[ha[:, 0], ha[:, 1]]
        self.hedraL23 = distplot[ha[:, 1], ha[:, 2]]
        self.hedraL13 = distplot[ha[:, 0], ha[:, 2]]
        da = self.d2a_map
        self.dihedraL14 = distplot[da[:, 0], da[:, 3]]
        self.dihedra_signs = dihedra_signs

    def distance_to_internal_coordinates(
        self, resetAtoms: Optional[Union[bool, None]] = True
    ) -> None:
        """Compute chain di/hedra from from distance and chirality data.

        Distance properties on hedra L12, L23, L13 and dihedra L14 configured
        by :meth:`.distplot_to_dh_arrays` or alternative loader.

        dihedraAngles result is multiplied by dihedra_signs at final step
        recover chirality information lost in distance plot (mirror image of
        structure has same distances but opposite sign dihedral angles).

        Note that chain breaks will cause errors in rebuilt structure, use
        :meth:`.copy_initNCaCs` to avoid this

        Based on Blue, the Hedronometer's answer to `The dihedral angles of a tetrahedron
        in terms of its edge lengths <https://math.stackexchange.com/a/49340/972353>`_
        on `math.stackexchange.com <https://math.stackexchange.com/>`_.  See also:
        `"Heron-like Hedronometric Results for Tetrahedral Volume"
        <http://daylateanddollarshort.com/mathdocs/Heron-like-Results-for-Tetrahedral-Volume.pdf>`_.

        Other values from that analysis included here as comments for
        completeness:

        * oa = hedron1 L12 if reverse else hedron1 L23
        * ob = hedron1 L23 if reverse else hedron1 L12
        * ac = hedron2 L12 if reverse else hedron2 L23
        * ab = hedron1 L13 = law of cosines on OA, OB (hedron1 L12, L23)
        * oc = hedron2 L13 = law of cosines on OA, AC (hedron2 L12, L23)
        * bc = dihedron L14

        target is OA, the dihedral angle along edge oa.

        :param bool resetAtoms: default True.
            Mark all atoms in di/hedra and atomArray for updating by
            :meth:`.internal_to_atom_coordinates`.  Alternatvely set this to
            False and manipulate `atomArrayValid`, `dAtoms_needs_update` and
            `hAtoms_needs_update` directly to reduce computation.
        """  # noqa
        oa = self.hedraL12[self.dH1ndx]
        oa[self.dFwd] = self.hedraL23[self.dH1ndx][self.dFwd]
        ob = self.hedraL23[self.dH1ndx]
        ob[self.dFwd] = self.hedraL12[self.dH1ndx][self.dFwd]
        ac = self.hedraL12[self.dH2ndx]
        ac[self.dFwd] = self.hedraL23[self.dH2ndx][self.dFwd]
        ab = self.hedraL13[self.dH1ndx]
        oc = self.hedraL13[self.dH2ndx]
        bc = self.dihedraL14

        # Ws = (ab + ac + bc) / 2
        # Xs = (ob + bc + oc) / 2
        Ys = (oa + ac + oc) / 2
        Zs = (oa + ob + ab) / 2
        # Wsqr = Ws * (Ws - ab) * (Ws - ac) * (Ws - bc)
        # Xsqr = Xs * (Xs - ob) * (Xs - bc) * (Xs - oc)
        Ysqr = Ys * (Ys - oa) * (Ys - ac) * (Ys - oc)
        Zsqr = Zs * (Zs - oa) * (Zs - ob) * (Zs - ab)
        Hsqr = (
            4 * oa * oa * bc * bc - np.square((ob * ob + ac * ac) - (oc * oc + ab * ab))
        ) / 16
        """
        Jsqr = (
            4 * ob * ob * ac * ac
            - np.square((oc * oc + ab * ab) - (oa * oa + bc * bc))
        ) / 16
        Ksqr = (
            4 * oc * oc * ab * ab
            - np.square((oa * oa + bc * bc) - (ob * ob + ac * ac))
        ) / 16
        """

        Y = np.sqrt(Ysqr)
        Z = np.sqrt(Zsqr)
        # X = np.sqrt(Xsqr)
        # W = np.sqrt(Wsqr)

        cosOA = (Ysqr + Zsqr - Hsqr) / (2 * Y * Z)
        # cosOB = (Zsqr + Xsqr - Jsqr) / (2 * Z * X)
        # cosOC = (Xsqr + Ysqr - Ksqr) / (2 * X * Y)
        # cosBC = (Wsqr + Xsqr - Hsqr) / (2 * W * X)
        # cosCA = (Wsqr + Ysqr - Jsqr) / (2 * W * Y)
        # cosAB = (Wsqr + Zsqr - Ksqr) / (2 * W * Z)

        # OA =
        # compute dihedral angles
        # ensure cosOA is in range [-1,1] for arccos
        cosOA[cosOA < -1.0] = -1.0
        cosOA[cosOA > 1.0] = 1.0
        # without np.longdouble here a few OCCACB angles lose last digit match
        np.arccos(cosOA, out=self.dihedraAngleRads, dtype=np.longdouble)
        self.dihedraAngleRads *= self.dihedra_signs
        np.rad2deg(self.dihedraAngleRads, out=self.dihedraAngle)
        # OB = np.rad2deg(np.arccos(cosOB))
        # OC = np.rad2deg(np.arccos(cosOC))
        # BC = np.rad2deg(np.arccos(cosBC))
        # CA = np.rad2deg(np.arccos(cosCA))
        # AB = np.rad2deg(np.arccos(cosAB))

        # law of cosines for hedra angles
        np.rad2deg(
            np.arccos(
                (
                    np.square(self.hedraL12)
                    + np.square(self.hedraL23)
                    - np.square(self.hedraL13)
                )
                / (2 * self.hedraL12 * self.hedraL23)
            ),
            out=self.hedraAngle,
        )

        if resetAtoms:
            self.atomArrayValid[:] = False
            self.dAtoms_needs_update[:] = True
            self.hAtoms_needs_update[:] = True

    def copy_initNCaCs(self, other: "IC_Chain") -> None:
        """Copy atom coordinates for initNCaC atoms from other IC_Chain.

        Copies the coordinates and sets atomArrayValid flags True for initial
        NCaC and after any chain breaks.

        Needed for :meth:`.distance_to_internal_coordinates` if target has
        chain breaks (otherwise each fragment will start at origin).

        Also useful if copying internal coordinates from another chain.

        N.B. :meth:`IC_Residue.set_angle()` and :meth:`IC_Residue.set_length()`
        invalidate their relevant atoms, so apply them before calling this
        function.
        """
        ndx = [self.atomArrayIndex[ak] for iNCaC in other.initNCaCs for ak in iNCaC]
        self.atomArray[ndx] = other.atomArray[ndx]
        self.atomArrayValid[ndx] = True

    def make_extended(self):
        """Set all psi and phi angles to extended conformation (123, -104)."""
        for ric in self.ordered_aa_ic_list:
            ric.set_angle("psi", 123)
            ric.set_angle("phi", -104)


class IC_Residue:
    """Class to extend Biopython Residue with internal coordinate data.

    Parameters
    ----------
    parent: biopython Residue object this class extends

    Attributes
    ----------
    no_altloc: bool default False
        **Class** variable, disable processing of ALTLOC atoms if True, use
        only selected atoms.

    accept_atoms: tuple
        **Class** variable :data:`accept_atoms`, list of PDB atom names to use
        when generating internal coordinates.
        Default is::

            accept_atoms = accept_mainchain + accept_hydrogens

        to exclude hydrogens in internal coordinates and generated PDB files,
        override as::

            IC_Residue.accept_atoms = IC_Residue.accept_mainchain

        to get only mainchain atoms plus amide proton, use::

            IC_Residue.accept_atoms = IC_Residue.accept_mainchain + ('H',)

        to convert D atoms to H, set :data:`AtomKey.d2h` = True and use::

            IC_Residue.accept_atoms = (
                accept_mainchain + accept_hydrogens + accept_deuteriums
            )

        Note that `accept_mainchain = accept_backbone + accept_sidechain`.
        Thus to generate sequence-agnostic conformational data for e.g.
        structure alignment in dihedral angle space, use::

            IC_Residue.accept_atoms = accept_backbone

        or set gly_Cbeta = True and use::

            IC_Residue.accept_atoms = accept_backbone + ('CB',)

        Changing accept_atoms will cause the default `structure_rebuild_test` in
        :mod:`.ic_rebuild` to fail if some atoms are filtered (obviously).  Use
        the `quick=True` option to test only the coordinates of filtered atoms
        to avoid this.

        There is currently no option to output internal coordinates with D
        instead of H.

    accept_resnames: tuple
        **Class** variable :data:`accept_resnames`, list of 3-letter residue
        names for HETATMs to accept when generating internal coordinates from
        atoms.  HETATM sidechain will be ignored, but normal backbone atoms (N,
        CA, C, O, CB) will be included.  Currently only CYG, YCM and UNK;
        override at your own risk.  To generate sidechain, add appropriate
        entries to `ic_data_sidechains` in :mod:`.ic_data` and support in
        :meth:`IC_Chain.atom_to_internal_coordinates`.

    gly_Cbeta: bool default False
        **Class** variable :data:`gly_Cbeta`, override to True to generate
        internal coordinates for glycine CB atoms in
        :meth:`IC_Chain.atom_to_internal_coordinates` ::

            IC_Residue.gly_Cbeta = True

    pic_accuracy: str default "17.13f"
        **Class** variable :data:`pic_accuracy` sets accuracy for numeric values
        (angles, lengths) in .pic files.  Default set high to support mmCIF file
        accuracy in rebuild tests.  If you find rebuild tests fail with
        'ERROR -COORDINATES-' and verbose=True shows only small discrepancies,
        try raising this value (or lower it to 9.5 if only working with PDB
        format files).  ::

            IC_Residue.pic_accuracy = "9.5f"

    residue: Biopython Residue object reference
        The :class:`.Residue` object this extends
    hedra: dict indexed by 3-tuples of AtomKeys
        Hedra forming this residue
    dihedra: dict indexed by 4-tuples of AtomKeys
        Dihedra forming (overlapping) this residue
    rprev, rnext: lists of IC_Residue objects
        References to adjacent (bonded, not missing, possibly disordered)
        residues in chain
    atom_coords: AtomKey indexed dict of numpy [4] arrays
        **removed**
        Use AtomKeys and atomArrayIndex to build if needed
    ak_set: set of AtomKeys in dihedra
        AtomKeys in all dihedra overlapping this residue (see __contains__())
    alt_ids: list of char
        AltLoc IDs from PDB file
    bfactors: dict
        AtomKey indexed B-factors as read from PDB file
    NCaCKey: List of tuples of AtomKeys
        List of tuples of N, Ca, C backbone atom AtomKeys; usually only 1
        but more if backbone altlocs.
    is20AA: bool
        True if residue is one of 20 standard amino acids, based on
        Residue resname
    isAccept: bool
        True if is20AA or in accept_resnames below
    rbase: tuple
        residue position, insert code or none, resname (1 letter if standard
        amino acid)
    cic: IC_Chain default None
        parent chain :class:`IC_Chain` object
    scale: optional float
        used for OpenSCAD output to generate gly_Cbeta bond length

    Methods
    -------
    assemble(atomCoordsIn, resetLocation, verbose)
        Compute atom coordinates for this residue from internal coordinates
    get_angle()
        Return angle for passed key
    get_length()
        Return bond length for specified pair
    pick_angle()
        Find Hedron or Dihedron for passed key
    pick_length()
        Find hedra for passed AtomKey pair
    set_angle()
        Set angle for passed key (no position updates)
    set_length()
        Set bond length in all relevant hedra for specified pair
    bond_rotate(delta)
        adjusts related dihedra angles by delta, e.g. rotating psi (N-Ca-C-N)
        will adjust the adjacent N-Ca-C-O by the same amount to avoid clashes
    bond_set(angle)
        uses bond_rotate to set specified dihedral to angle and adjust related
        dihedra accordingly
    rak(atom info)
        cached AtomKeys for this residue
    """

    accept_resnames = ("CYG", "YCM", "UNK")
    """Add 3-letter residue name here for non-standard residues with
    normal backbone.  CYG included for test case 4LGY (1305 residue
    contiguous chain).  Safe to add more names for N-CA-C-O backbones, any
    more complexity will need additions to :data:`accept_atoms`,
    `ic_data_sidechains` in :mod:`.ic_data` and support in
    :meth:`IC_Chain.atom_to_internal_coordinates`"""

    _AllBonds: bool = False
    """For OpenSCAD output, generate explicit hedra covering all bonds.
    **Class** variable, whereas a PDB file just specifies atoms, OpenSCAD
    output for 3D printing needs all bonds specified explicitly - otherwise
    e.g. PHE rings will not be closed.  This variable is managed by the
    :func:`.SCADIO.write_SCAD` code."""

    no_altloc: bool = False
    """Set True to filter altloc atoms on input and only work with Biopython
    default Atoms"""

    gly_Cbeta: bool = False
    """Create beta carbons on all Gly residues.

    Setting this to True will generate internal coordinates for Gly C-beta
    carbons in :meth:`atom_to_internal_coordinates`.

    Data averaged from Sep 2019 Dunbrack cullpdb_pc20_res2.2_R1.0
    restricted to structures with amide protons.
    Please see

    `PISCES: A Protein Sequence Culling Server <https://dunbrack.fccc.edu/pisces/>`_

    'G. Wang and R. L. Dunbrack, Jr. PISCES: a protein sequence culling
    server. Bioinformatics, 19:1589-1591, 2003.'

    Ala avg rotation of OCCACB from NCACO query::

        select avg(g.rslt) as avg_rslt, stddev(g.rslt) as sd_rslt, count(*)
        from
        (select f.d1d, f.d2d,
        (case when f.rslt > 0 then f.rslt-360.0 else f.rslt end) as rslt
        from (select d1.angle as d1d, d2.angle as d2d,
        (d2.angle - d1.angle) as rslt from dihedron d1,
        dihedron d2 where d1.re_class='AOACACAACB' and
        d2.re_class='ANACAACAO' and d1.pdb=d2.pdb and d1.chn=d2.chn
        and d1.res=d2.res) as f) as g

    results::

        | avg_rslt          | sd_rslt          | count   |
        | -122.682194862932 | 5.04403040513919 | 14098   |
"""
    pic_accuracy: str = (
        "17.13f"  # output accuracy for angle and len values in .pic files
    )

    accept_backbone = (
        "N",
        "CA",
        "C",
        "O",
        "OXT",
    )
    accept_sidechain = (
        "CB",
        "CG",
        "CG1",
        "OG1",
        "OG",
        "SG",
        "CG2",
        "CD",
        "CD1",
        "SD",
        "OD1",
        "ND1",
        "CD2",
        "ND2",
        "CE",
        "CE1",
        "NE",
        "OE1",
        "NE1",
        "CE2",
        "OE2",
        "NE2",
        "CE3",
        "CZ",
        "NZ",
        "CZ2",
        "CZ3",
        "OD2",
        "OH",
        "CH2",
        "NH1",
        "NH2",
    )

    accept_mainchain = accept_backbone + accept_sidechain

    accept_hydrogens = (
        "H",
        "H1",
        "H2",
        "H3",
        "HA",
        "HA2",
        "HA3",
        "HB",
        "HB1",
        "HB2",
        "HB3",
        "HG2",
        "HG3",
        "HD2",
        "HD3",
        "HE2",
        "HE3",
        "HZ1",
        "HZ2",
        "HZ3",
        "HG11",
        "HG12",
        "HG13",
        "HG21",
        "HG22",
        "HG23",
        "HZ",
        "HD1",
        "HE1",
        "HD11",
        "HD12",
        "HD13",
        "HG",
        "HG1",
        "HD21",
        "HD22",
        "HD23",
        "NH1",
        "NH2",
        "HE",
        "HH11",
        "HH12",
        "HH21",
        "HH22",
        "HE21",
        "HE22",
        "HE2",
        "HH",
        "HH2",
    )
    accept_deuteriums = (
        "D",
        "D1",
        "D2",
        "D3",
        "DA",
        "DA2",
        "DA3",
        "DB",
        "DB1",
        "DB2",
        "DB3",
        "DG2",
        "DG3",
        "DD2",
        "DD3",
        "DE2",
        "DE3",
        "DZ1",
        "DZ2",
        "DZ3",
        "DG11",
        "DG12",
        "DG13",
        "DG21",
        "DG22",
        "DG23",
        "DZ",
        "DD1",
        "DE1",
        "DD11",
        "DD12",
        "DD13",
        "DG",
        "DG1",
        "DD21",
        "DD22",
        "DD23",
        "ND1",
        "ND2",
        "DE",
        "DH11",
        "DH12",
        "DH21",
        "DH22",
        "DE21",
        "DE22",
        "DE2",
        "DH",
        "DH2",
    )
    accept_atoms = accept_mainchain + accept_hydrogens
    """Change accept_atoms to restrict atoms processed. See :class:`IC_Residue`
    for usage."""

    def __init__(self, parent: "Residue") -> None:
        """Initialize IC_Residue with parent Biopython Residue.

        :param Residue parent: Biopython Residue object.
            The Biopython Residue this object extends
        """
        self.residue = parent
        self.cic: IC_Chain
        # dict of hedron objects indexed by hedron keys
        self.hedra: Dict[HKT, Hedron] = {}
        # dict of dihedron objects indexed by dihedron keys
        self.dihedra: Dict[DKT, Dihedron] = {}
        # cache of AtomKey results for rak()
        self.akc: Dict[Union[str, Atom], AtomKey] = {}
        # set of AtomKeys involved in dihedra, used by split_akl,
        # build_rak_cache.  Built by __init__ for XYZ (PDB coord) input,
        # _link_dihedra for PIC input
        self.ak_set: Set[AtomKey] = set()
        # reference to adjacent residues in chain
        self.rprev: List[IC_Residue] = []
        self.rnext: List[IC_Residue] = []
        # bfactors copied from PDB file
        self.bfactors: Dict[str, float] = {}
        self.alt_ids: Union[List[str], None] = None if IC_Residue.no_altloc else []
        self.is20AA = True
        self.isAccept = True
        # self.NCaCKey Set by _link_dihedra()
        # rbase = position, insert code or none, resname (1 letter if in 20)
        rid = parent.id
        rbase = [rid[1], rid[2] if " " != rid[2] else None, parent.resname]
        try:
            rbase[2] = protein_letters_3to1[rbase[2]]
        except KeyError:
            self.is20AA = False
            if rbase[2] not in self.accept_resnames:
                self.isAccept = False

        self.rbase = tuple(rbase)
        self.lc = rbase[2]

        if self.isAccept:
            for atom in parent.get_atoms():
                if hasattr(atom, "child_dict"):
                    if IC_Residue.no_altloc:
                        self._add_atom(atom.selected_child)
                    else:
                        for atm in atom.child_dict.values():
                            self._add_atom(atm)
                else:
                    self._add_atom(atom)
            if self.ak_set:
                # only for coordinate (pdb) input, _add_atom loads
                # init cache ready for atom_to_internal_coords
                self._build_rak_cache()

    def __deepcopy__(self, memo):
        """Deep copy implementation for IC_Residue."""
        existing = memo.get(id(self), False)
        if existing:
            return existing
        dup = type(self).__new__(self.__class__)
        memo[id(self)] = dup
        dup.__dict__.update(self.__dict__)  # later replace what is not static
        dup.cic = memo[id(self.cic)]
        dup.residue = memo[id(self.residue)]
        # still need to update: rnext, rprev, akc, ak_set, di/hedra
        return dup

    def __contains__(self, ak: "AtomKey") -> bool:
        """Return True if atomkey is in this residue."""
        if ak in self.ak_set:
            akl = ak.akl
            if (
                int(akl[0]) == self.rbase[0]
                and akl[1] == self.rbase[1]
                and akl[2] == self.rbase[2]
            ):
                return True
        return False

    def rak(self, atm: Union[str, Atom]) -> "AtomKey":
        """Cache calls to AtomKey for this residue."""
        try:
            ak = self.akc[atm]
        except KeyError:
            ak = self.akc[atm] = AtomKey(self, atm)
            if isinstance(atm, str):
                ak.missing = True
        return ak

    def _build_rak_cache(self) -> None:
        """Create explicit entries for for atoms so don't miss altlocs.

        This ensures that self.akc (atom key cache) has an entry for selected
        atom name (e.g. "CA") amongst any that have altlocs.  Without this,
        rak() on the other altloc atom first may result in the main atom being
        missed.
        """
        for ak in sorted(self.ak_set):
            atmName = ak.akl[3]
            if self.akc.get(atmName) is None:
                self.akc[atmName] = ak

    def _add_atom(self, atm: Atom) -> None:
        """Filter Biopython Atom with accept_atoms; set ak_set.

        Arbitrarily renames O' and O'' to O and OXT
        """
        if "O" == atm.name[0]:
            if "O'" == atm.name:
                atm.name = "O"
            elif "O''" == atm.name:
                atm.name = "OXT"

        if atm.name not in self.accept_atoms:
            # print('skip:', atm.name)
            return
        ak = self.rak(atm)  # passing Atom here not string
        self.ak_set.add(ak)

    def __repr__(self) -> str:
        """Print string is parent Residue ID."""
        return str(self.residue.full_id)

    def pretty_str(self) -> str:
        """Nice string for residue ID."""
        id = self.residue.id
        return f"{self.residue.resname} {id[0]}{id[1]!s}{id[2]}"

    def _link_dihedra(self, verbose: bool = False) -> None:
        """Housekeeping after loading all residues and dihedra.

        - Link dihedra to this residue
        - form id3_dh_index
        - form ak_set
        - set NCaCKey to be available AtomKeys

        called for loading PDB / atom coords
        """
        for dh in self.dihedra.values():
            dh.ric = self  # each dihedron can find its IC_Residue
            dh.cic = self.cic  # each dihedron can update chain dihedral angles
            self.ak_set.update(dh.atomkeys)
        for h in self.hedra.values():  # collect any atoms in orphan hedra
            self.ak_set.update(h.atomkeys)  # e.g. alternate CB path with no O
            h.cic = self.cic  # each hedron can update chain hedra

        # if loaded PIC data, akc not initialised yet
        if not self.akc:
            self._build_rak_cache()

        # initialise NCaCKey here:
        self.NCaCKey = []
        self.NCaCKey.extend(
            self.split_akl(
                (AtomKey(self, "N"), AtomKey(self, "CA"), AtomKey(self, "C"))
            )
        )

    def set_flexible(self) -> None:
        """For OpenSCAD, mark N-CA and CA-C bonds to be flexible joints.

        See :func:`.SCADIO.write_SCAD`
        """
        for h in self.hedra.values():
            if h.e_class == "NCAC":
                h.flex_female_1 = True
                h.flex_female_2 = True
            elif h.e_class.endswith("NCA"):
                h.flex_male_2 = True
            elif h.e_class.startswith("CAC") and h.atomkeys[1].akl[3] == "C":
                h.flex_male_1 = True
            elif h.e_class == "CBCAC":
                h.skinny_1 = True  # CA-CB bond interferes with flex join

    def set_hbond(self) -> None:
        """For OpenSCAD, mark H-N and C-O bonds to be hbonds (magnets).

        See :func:`.SCADIO.write_SCAD`
        """
        for h in self.hedra.values():
            if h.e_class == "HNCA":
                h.hbond_1 = True
            elif h.e_class == "CACO":
                h.hbond_2 = True

    def _default_startpos(self) -> Dict["AtomKey", np.array]:
        """Generate default N-Ca-C coordinates to build this residue from."""
        atomCoords = {}
        cic = self.cic
        dlist0 = [cic.id3_dh_index.get(akl, None) for akl in sorted(self.NCaCKey)]
        dlist1 = [d for d in dlist0 if d is not None]
        # https://stackoverflow.com/questions/11264684/flatten-list-of-lists
        dlist = [cic.dihedra[val] for sublist in dlist1 for val in sublist]
        # dlist = self.id3_dh_index[NCaCKey]
        for d in dlist:
            for i, a in enumerate(d.atomkeys):
                # atomCoords[a] = d.initial_coords[i]
                atomCoords[a] = cic.dAtoms[d.ndx][i]
                # cic.atomArray[cic.atomArrayIndex[a]] = atomCoords[a]
                # cic.atomArrayValid[cic.atomArrayIndex[a]] = True
        return atomCoords

    def _get_startpos(self) -> Dict["AtomKey", np.array]:
        """Find N-Ca-C coordinates to build this residue from."""
        # only used by assemble()
        startPos = {}
        cic = self.cic
        for ncac in self.NCaCKey:
            if np.all(cic.atomArrayValid[[cic.atomArrayIndex[ak] for ak in ncac]]):
                for ak in ncac:
                    startPos[ak] = cic.atomArray[cic.atomArrayIndex[ak]]
        if startPos == {}:
            startPos = self._default_startpos()
        return startPos

    def clear_transforms(self):
        """Invalidate dihedra coordinate space attributes before assemble().

        Coordinate space attributes are Dihedron.cst and .rcst, and
        :data:`IC_Chain.dCoordSpace`
        """
        for d in self.dihedra.values():
            self.cic.dcsValid[d.ndx] = False

    def assemble(
        self,
        resetLocation: bool = False,
        verbose: bool = False,
    ) -> Union[Dict["AtomKey", np.array], Dict[HKT, np.array], None]:
        """Compute atom coordinates for this residue from internal coordinates.

        This is the IC_Residue part of the :meth:`.assemble_residues_ser` serial
        version, see :meth:`.assemble_residues` for numpy vectorized approach
        which works at the :class:`IC_Chain` level.

        Join prepared dihedra starting from N-CA-C and N-CA-CB hedrons,
        computing protein space coordinates for backbone and sidechain atoms

        Sets forward and reverse transforms on each Dihedron to convert from
        protein coordinates to dihedron space coordinates for first three
        atoms (see :data:`IC_Chain.dCoordSpace`)

        Call :meth:`.init_atom_coords` to update any modified di/hedra before
        coming here, this only assembles dihedra into protein coordinate space.

        **Algorithm**

        Form double-ended queue, start with c-ca-n, o-c-ca, n-ca-cb, n-ca-c.

        if resetLocation=True, use initial coords from generating dihedron
        for n-ca-c initial positions (result in dihedron coordinate space)

        while queue not empty
            get 3-atom hedron key

            for each dihedron starting with hedron key (1st hedron of dihedron)

                if have coordinates for all 4 atoms already
                    add 2nd hedron key to back of queue
                else if have coordinates for 1st 3 atoms
                    compute forward and reverse transforms to take 1st 3 atoms
                    to/from dihedron initial coordinate space

                    use reverse transform to get position of 4th atom in
                    current coordinates from dihedron initial coordinates

                    add 2nd hedron key to back of queue
                else
                    ordering failed, put hedron key at back of queue and hope
                    next time we have 1st 3 atom positions (should not happen)

        loop terminates (queue drains) as hedron keys which do not start any
        dihedra are removed without action

        :param bool resetLocation: default False.
            - Option to ignore start location and orient so initial N-Ca-C
            hedron at origin.

        :returns:
            Dict of AtomKey -> homogeneous atom coords for residue in protein
            space relative to previous residue

            **Also** directly updates :data:`IC_Chain.atomArray` as
            :meth:`.assemble_residues` does.

        """
        # debug statements below still useful, commented for performance
        # dbg = False
        # if hasattr(IC_Chain, "adbg"):
        #    dbg = IC_Chain.adbg

        cic = self.cic
        dcsValid = cic.dcsValid
        aaValid = cic.atomArrayValid
        aaNdx = cic.atomArrayIndex
        aa = cic.atomArray

        if not self.ak_set:
            return None  # give up now if no atoms to work with

        NCaCKey = sorted(self.NCaCKey)
        rseqpos = self.rbase[0]

        # order of these startLst entries matters
        startLst = self.split_akl((self.rak("C"), self.rak("CA"), self.rak("N")))
        if "CB" in self.akc:
            startLst.extend(
                self.split_akl((self.rak("N"), self.rak("CA"), self.rak("CB")))
            )
        if "O" in self.akc:
            startLst.extend(
                self.split_akl((self.rak("O"), self.rak("C"), self.rak("CA")))
            )

        startLst.extend(NCaCKey)

        q = deque(startLst)
        # resnum = self.rbase[0]

        # get initial coords from previous residue or IC_Chain info
        # or default coords
        if resetLocation:
            # use N-CA-C initial coords from creating dihedral
            atomCoords = self._default_startpos()
        else:
            atomCoords = self._get_startpos()

        while q:  # deque is not empty
            """
            if dbg:
                print("assemble loop start q=", q)
            """
            h1k = cast(HKT, q.pop())
            dihedraKeys = cic.id3_dh_index.get(h1k, None)
            """
            if dbg:
                print(
                    "  h1k:",
                    h1k,
                    "len dihedra: ",
                    len(dihedraKeys) if dihedraKeys is not None else "None",
                )
            """
            if dihedraKeys is not None:
                for dk in dihedraKeys:
                    d = cic.dihedra[dk]
                    dseqpos = int(d.atomkeys[0].akl[AtomKey.fields.respos])
                    d.initial_coords = cic.dAtoms[d.ndx]
                    if 4 == len(d.initial_coords) and d.initial_coords[3] is not None:
                        # skip incomplete dihedron if don't have 4th atom due
                        # to missing input data
                        d_h2key = d.hedron2.atomkeys
                        ak = d.atomkeys[3]
                        """
                        if dbg:
                            print("    process", d, d_h2key, d.atomkeys)
                        """

                        acount = len([a for a in d.atomkeys if a in atomCoords])

                        if 4 == acount:
                            # dihedron already done, queue 2nd hedron key
                            if dseqpos == rseqpos:  # only this residue
                                q.appendleft(d_h2key)
                            """
                            if dbg:
                                print("    4- already done, append left")
                            """
                            if not dcsValid[d.ndx]:  # missing transform
                                # can happen for altloc atoms
                                # only needed for write_SCAD output
                                acs = [atomCoords[a] for a in h1k]
                                d.cst, d.rcst = coord_space(
                                    acs[0], acs[1], acs[2], True
                                )
                                dcsValid[d.ndx] = True
                        elif 3 == acount:
                            """
                            if dbg:
                                print("    3- call coord_space")
                            """

                            acs = np.asarray([atomCoords[a] for a in h1k])
                            d.cst, d.rcst = coord_space(acs[0], acs[1], acs[2], True)
                            dcsValid[d.ndx] = True
                            """
                            if dbg:
                                print("     acs:", acs.transpose())
                                print("cst", d.cst)
                                print("rcst", d.rcst)
                                print(
                                    "        initial_coords[3]=",
                                    d.initial_coords[3].transpose(),
                                )
                            """
                            acak3 = d.rcst.dot(d.initial_coords[3])
                            """
                            if dbg:
                                print("        acak3=", acak3.transpose())
                            """
                            atomCoords[ak] = acak3
                            aa[aaNdx[ak]] = acak3
                            aaValid[aaNdx[ak]] = True
                            """
                            if dbg:
                                print(
                                    "        3- finished, ak:",
                                    ak,
                                    "coords:",
                                    atomCoords[ak].transpose(),
                                )
                            """
                            if dseqpos == rseqpos:  # only this residue
                                q.appendleft(d_h2key)
                        else:
                            if verbose:
                                print("no coords to start", d)
                                print(
                                    [
                                        a
                                        for a in d.atomkeys
                                        if atomCoords.get(a, None) is not None
                                    ]
                                )
                    else:
                        if verbose:
                            print("no initial coords for", d)

        return atomCoords

    def split_akl(
        self,
        lst: Union[Tuple["AtomKey", ...], List["AtomKey"]],
        missingOK: bool = False,
    ) -> List[Tuple["AtomKey", ...]]:
        """Get AtomKeys for this residue (ak_set) for generic list of AtomKeys.

        Changes and/or expands a list of 'generic' AtomKeys (e.g. 'N, C, C') to
        be specific to this Residue's altlocs etc., e.g.
        '(N-Ca_A_0.3-C, N-Ca_B_0.7-C)'

        Given a list of AtomKeys for a Hedron or Dihedron,
          return:
                list of matching atomkeys that have id3_dh in this residue
                (ak may change if occupancy != 1.00)

            or
                multiple lists of matching atomkeys expanded for all atom altlocs

            or
                empty list if any of atom_coord(ak) missing and not missingOK

        :param list lst: list[3] or [4] of AtomKeys.
            Non-altloc AtomKeys to match to specific AtomKeys for this residue
        :param bool missingOK: default False, see above.
        """
        altloc_ndx = AtomKey.fields.altloc
        occ_ndx = AtomKey.fields.occ

        # step 1
        # given a list of AtomKeys
        #  form a new list of same atomkeys with coords or diheds in this residue
        #      plus lists of matching altloc atomkeys in coords or diheds
        edraLst: List[Tuple[AtomKey, ...]] = []
        altlocs = set()
        posnAltlocs: Dict["AtomKey", Set[str]] = {}
        akMap = {}
        for ak in lst:
            posnAltlocs[ak] = set()
            if (
                ak in self.ak_set
                and ak.akl[altloc_ndx] is None
                and ak.akl[occ_ndx] is None
            ):
                # simple case no altloc and exact match in set
                edraLst.append((ak,))  # tuple of ak
            else:
                ak2_lst = []
                for ak2 in self.ak_set:
                    if ak.altloc_match(ak2):
                        # print(key)
                        ak2_lst.append(ak2)
                        akMap[ak2] = ak
                        altloc = ak2.akl[altloc_ndx]
                        if altloc is not None:
                            altlocs.add(altloc)
                            posnAltlocs[ak].add(altloc)
                edraLst.append(tuple(ak2_lst))

        # step 2
        # check and finish for
        #   missing atoms
        #   simple case no altlocs
        # else form new AtomKey lists covering all altloc permutations
        maxc = 0
        for akl in edraLst:
            lenAKL = len(akl)
            if 0 == lenAKL and not missingOK:
                return []  # atom missing in atom_coords, cannot form object
            elif maxc < lenAKL:
                maxc = lenAKL
        if 1 == maxc:  # simple case no altlocs for any ak in list
            newAKL = []
            for akl in edraLst:
                if akl:  # may have empty lists if missingOK, do not append
                    newAKL.append(akl[0])
            return [tuple(newAKL)]
        else:
            new_edraLst = []
            for al in altlocs:
                # form complete new list for each altloc
                alhl = []
                for akl in edraLst:
                    lenAKL = len(akl)
                    if 0 == lenAKL:
                        continue  # ignore empty list from missingOK
                    if 1 == lenAKL:
                        alhl.append(akl[0])  # not all atoms will have altloc
                    # elif (lenAKL < maxc
                    #      and al not in posnAltlocs[akMap[akl[0]]]):
                    elif al not in posnAltlocs[akMap[akl[0]]]:
                        # this position has fewer altlocs than other positions
                        # and this position does not have this al,
                        # so just grab first to form angle as could be any
                        alhl.append(sorted(akl)[0])
                    else:
                        for ak in akl:
                            if ak.akl[altloc_ndx] == al:
                                alhl.append(ak)
                new_edraLst.append(tuple(alhl))

            # print(new_edraLst)
            return new_edraLst

    def _gen_edra(self, lst: Union[Tuple["AtomKey", ...], List["AtomKey"]]) -> None:
        """Populate hedra/dihedra given edron ID tuple.

        Given list of AtomKeys defining hedron or dihedron
          convert to AtomKeys with coordinates in this residue
          add appropriately to self.di/hedra, expand as needed atom altlocs

        :param list lst: tuple of AtomKeys.
            Specifies Hedron or Dihedron
        """
        for ak in lst:
            if ak.missing:
                return  # give up if atoms actually missing

        lenLst = len(lst)
        if 4 > lenLst:
            cdct, dct, obj = self.cic.hedra, self.hedra, Hedron
        else:
            cdct, dct, obj = self.cic.dihedra, self.dihedra, Dihedron  # type: ignore # noqa

        if isinstance(lst, List):
            tlst = tuple(lst)
        else:
            tlst = lst

        hl = self.split_akl(tlst)  # expand tlst with any altlocs
        # returns list of tuples

        for tnlst in hl:
            # do not add edron if split_akl() made something shorter
            if len(tnlst) == lenLst:
                # if edron already exists, then update not replace with new
                if tnlst not in cdct:
                    cdct[tnlst] = obj(tnlst)  # type: ignore
                if tnlst not in dct:
                    dct[tnlst] = cdct[tnlst]  # type: ignore

                dct[tnlst].needs_update = True  # type: ignore

    # @profile
    def _create_edra(self, verbose: bool = False) -> None:
        """Create IC_Chain and IC_Residue di/hedra for atom coordinates.

        AllBonds handled here.

        :param bool verbose: default False.
            Warn about missing N, Ca, C backbone atoms.
        """
        # on entry we have all Biopython Atoms loaded
        if not self.ak_set:
            return  # so give up if no atoms loaded for this residue

        sN, sCA, sC = self.rak("N"), self.rak("CA"), self.rak("C")
        if self.lc != "G":
            sCB = self.rak("CB")

        # first init di/hedra, AtomKey objects and atom_coords for di/hedra
        # which extend into next residue.

        if 0 < len(self.rnext) and self.rnext[0].ak_set:
            # atom_coords, hedra and dihedra for backbone dihedra
            # which reach into next residue
            for rn in self.rnext:
                nN, nCA, nC = rn.rak("N"), rn.rak("CA"), rn.rak("C")

                nextNCaC = rn.split_akl((nN, nCA, nC), missingOK=True)

                for tpl in nextNCaC:
                    for ak in tpl:
                        if ak in rn.ak_set:
                            self.ak_set.add(ak)
                        else:
                            for rn_ak in rn.ak_set:
                                if rn_ak.altloc_match(ak):
                                    self.ak_set.add(rn_ak)

                self._gen_edra((sN, sCA, sC, nN))  # psi
                self._gen_edra((sCA, sC, nN, nCA))  # omega i+1
                self._gen_edra((sC, nN, nCA, nC))  # phi i+1
                self._gen_edra((sCA, sC, nN))
                self._gen_edra((sC, nN, nCA))
                self._gen_edra((nN, nCA, nC))  # tau i+1

                # redundant next residue C-beta locator (alternate CB path)
                # otherwise missing O will cause no sidechain
                try:
                    nO = rn.akc["O"]  # noqa: F841
                except KeyError:
                    # not rn.rak here so don't trigger missing CB for Gly
                    nCB = rn.akc.get("CB", None)
                    if nCB is not None and nCB in rn.ak_set:
                        self.ak_set.add(nCB)
                        self._gen_edra((nN, nCA, nCB))
                        self._gen_edra((sC, nN, nCA, nCB))

        # if start of chain then need to init NCaC hedron as not in previous
        # residue
        if 0 == len(self.rprev):
            self._gen_edra((sN, sCA, sC))

        # now init di/hedra for standard backbone atoms independent of
        # neighbours
        backbone = ic_data_backbone
        for edra in backbone:
            # only need to build if this residue has all the atoms in the edra
            if all(atm in self.akc for atm in edra):
                r_edra = [self.rak(atom) for atom in edra]
                self._gen_edra(r_edra)  # [4] is label on some table entries

        # next init sidechain di/hedra
        if self.lc is not None:
            sidechain = ic_data_sidechains.get(self.lc, [])
            for edraLong in sidechain:
                edra = edraLong[0:4]  # [4] is label on some sidechain table entries
                # lots of H di/hedra can be avoided if don't have those atoms
                if all(atm in self.akc for atm in edra):
                    r_edra = [self.rak(atom) for atom in edra]
                    self._gen_edra(r_edra)
            if (
                IC_Residue._AllBonds
            ):  # openscad output needs all bond cylinders explicit
                sidechain = ic_data_sidechain_extras.get(self.lc, [])
                for edra in sidechain:
                    # test less useful here but avoids populating rak cache if
                    # possible
                    if all(atm in self.akc for atm in edra):
                        r_edra = [self.rak(atom) for atom in edra]
                        self._gen_edra(r_edra)

        # create di/hedra for gly Cbeta if needed, populate values later
        if self.gly_Cbeta and "G" == self.lc:
            # add C-beta for Gly
            self.ak_set.add(AtomKey(self, "CB"))
            sCB = self.rak("CB")
            sCB.missing = False  # was True because akc cache did not have entry
            self.cic.akset.add(sCB)

            # main orientation comes from O-C-Ca-Cb so make Cb-Ca-C hedron
            sO = self.rak("O")
            htpl = (sCB, sCA, sC)
            self._gen_edra(htpl)

            # generate dihedral based on N-Ca-C-O offset from db query above
            dtpl = (sO, sC, sCA, sCB)
            self._gen_edra(dtpl)
            d = self.dihedra[dtpl]
            d.ric = self
            d._set_hedra()

            # prepare to add new Gly CB atom(s)
            # in IC_Chain.atom_to_internal_coordinates()
            if not hasattr(self.cic, "gcb"):
                self.cic.gcb = {}
            self.cic.gcb[sCB] = dtpl

        # final processing of all dihedra just generated
        self._link_dihedra(verbose)  # re-run for new dihedra

        if verbose:
            # oAtom =
            self.rak("O")  # trigger missing flag if needed
            missing = []
            for akk, akv in self.akc.items():
                if isinstance(akk, str) and akv.missing:
                    missing.append(akv)
            if missing:
                chn = self.residue.parent
                chn_id = chn.id
                chn_len = len(chn.internal_coord.ordered_aa_ic_list)
                print(f"chain {chn_id} len {chn_len} missing atom(s): {missing}")

    # rtm
    atom_sernum = None
    atom_chain = None

    @staticmethod
    def _pdb_atom_string(atm: Atom, cif_extend: bool = False) -> str:
        """Generate PDB ATOM record.

        :param Atom atm: Biopython Atom object reference
        :param IC_Residue.atom_sernum: Class variable default None.
            override atom serial number if not None
        :param IC_Residue.atom_chain: Class variable default None.
            override atom chain id if not None
        """
        if 2 == atm.is_disordered():
            if IC_Residue.no_altloc:
                return IC_Residue._pdb_atom_string(atm.selected_child, cif_extend)
            s = ""
            for a in atm.child_dict.values():
                s += IC_Residue._pdb_atom_string(a, cif_extend)
            return s
        else:
            res = atm.parent
            chn = res.parent
            fmt = "{:6}{:5d} {:4}{:1}{:3} {:1}{:4}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}        {:>4}\n"
            if cif_extend:
                fmt = "{:6}{:5d} {:4}{:1}{:3} {:1}{:4}{:1}   {:10.5f}{:10.5f}{:10.5f}{:7.3f}{:6.2f}        {:>4}\n"
            s = (fmt).format(
                "ATOM",
                IC_Residue.atom_sernum
                if IC_Residue.atom_sernum is not None
                else atm.serial_number,
                atm.fullname,
                atm.altloc,
                res.resname,
                IC_Residue.atom_chain if IC_Residue.atom_chain is not None else chn.id,
                res.id[1],
                res.id[2],
                atm.coord[0],
                atm.coord[1],
                atm.coord[2],
                atm.occupancy,
                atm.bfactor,
                atm.element,
            )
            # print(s)
        return s

    # rtm
    def pdb_residue_string(self) -> str:
        """Generate PDB ATOM records for this residue as string.

        Convenience method for functionality not exposed in PDBIO.py.
        Increments :data:`IC_Residue.atom_sernum` if not None

        :param IC_Residue.atom_sernum: Class variable default None.
            Override and increment atom serial number if not None
        :param IC_Residue.atom_chain: Class variable.
            Override atom chain id if not None

        .. todo::
            move to PDBIO
        """
        str = ""
        atomArrayIndex = self.cic.atomArrayIndex
        bpAtomArray = self.cic.bpAtomArray
        respos = self.rbase[0]
        resposNdx = AtomKey.fields.respos

        for ak in sorted(self.ak_set):
            if int(ak.akl[resposNdx]) == respos:  # skip rnext atoms
                str += IC_Residue._pdb_atom_string(bpAtomArray[atomArrayIndex[ak]])
                if IC_Residue.atom_sernum is not None:
                    IC_Residue.atom_sernum += 1
        return str

    @staticmethod
    def _residue_string(res: "Residue") -> str:
        """Generate PIC Residue string.

        Enough to create Biopython Residue object without actual Atoms.

        :param Residue res: Biopython Residue object reference
        """
        segid = res.get_segid()
        if segid.isspace() or "" == segid:
            segid = ""
        else:
            segid = " [" + segid + "]"
        return str(res.get_full_id()) + " " + res.resname + segid + "\n"

    _pfDef = namedtuple(
        # general supersedes specific, so pomg + omg = omg, tau + hedra = hedra
        "_pfDef",
        [
            "psi",  # _b[0]
            "omg",
            "phi",
            "tau",  # tau hedron (N-Ca-C)
            "chi1",
            "chi2",
            "chi3",
            "chi4",
            "chi5",
            "pomg",  # _b[9] : proline omega
            "chi",  # chi1 | ... | chi5
            "classic_b",  # psi | phi | tau | pomg
            "classic",  # classic_b | chi
            "hedra",  # _b[10] : all hedra
            "primary",  # _b[11] : all primary dihedra
            "secondary",  # _b[12] : all secondary dihedra
            "all",  # hedra | primary | secondary
            "initAtoms",  # _b[13] : XYZ coordinates of initial Tau (N-Ca-C)
            "bFactors",  # _b[14]
        ],
    )

    _b = [1 << i for i in range(16)]
    _bChi = _b[4] | _b[5] | _b[6] | _b[7] | _b[8]
    _bClassB = _b[0] | _b[2] | _b[3] | _b[9]
    _bClass = _bClassB | _bChi
    _bAll = _b[10] | _b[11] | _b[12]

    pic_flags = _pfDef(
        _b[0],
        _b[1],
        _b[2],
        _b[3],
        _b[4],
        _b[5],
        _b[6],
        _b[7],
        _b[8],
        _b[9],
        _bChi,
        _bClassB,
        _bClass,
        _b[10],
        _b[11],
        _b[12],
        _bAll,
        _b[13],
        _b[14],
    )
    """Used by :func:`.PICIO.write_PIC` to control classes of values to be defaulted."""

    picFlagsDefault = pic_flags.all | pic_flags.initAtoms | pic_flags.bFactors
    """Default is all dihedra + initial tau atoms + bFactors."""

    picFlagsDict = pic_flags._asdict()
    """Dictionary of pic_flags values to use as needed."""

    def _write_pic_bfac(self, atm: Atom, s: str, col: int) -> Tuple[str, int]:
        ak = self.rak(atm)
        if 0 == col % 5:
            s += "BFAC:"
        s += " " + ak.id + " " + f"{atm.get_bfactor():6.2f}"
        col += 1
        if 0 == col % 5:
            s += "\n"
        return s, col

    def _write_PIC(
        self,
        pdbid: str = "0PDB",
        chainid: str = "A",
        picFlags: int = picFlagsDefault,
        hCut: Optional[Union[float, None]] = None,
        pCut: Optional[Union[float, None]] = None,
    ) -> str:
        """Write PIC format lines for this residue.

        See :func:`.PICIO.write_PIC`.

        :param str pdbid: PDB idcode string; default 0PDB
        :param str chainid: PDB Chain ID character; default A
        :param int picFlags: control details written to PIC file; see
            :meth:`.PICIO.write_PIC`
        :param float hCut: only write hedra with ref db angle std dev > this
            value; default None
        :param float pCut: only write primary dihedra with ref db angle
            std dev > this value; default None
        """
        pAcc = IC_Residue.pic_accuracy
        if pdbid is None:
            pdbid = "0PDB"
        if chainid is None:
            chainid = "A"
        icr = IC_Residue
        s = icr._residue_string(self.residue)

        if (
            picFlags & icr.pic_flags.initAtoms
            and 0 == len(self.rprev)  # no prev residue
            and hasattr(self, "NCaCKey")
            and self.NCaCKey is not None  # have valid NCacKey
            # N coords valid (e.g. not all 0.00)
            and not (np.all(self.residue["N"].coord == self.residue["N"].coord[0]))
        ):
            NCaChedron = self.pick_angle(self.NCaCKey[0])  # first tau
            if NCaChedron is not None:
                try:
                    ts = IC_Residue._pdb_atom_string(self.residue["N"], cif_extend=True)
                    ts += IC_Residue._pdb_atom_string(
                        self.residue["CA"], cif_extend=True
                    )
                    ts += IC_Residue._pdb_atom_string(
                        self.residue["C"], cif_extend=True
                    )
                    s += ts  # only if no exception: have all 3 atoms
                except KeyError:
                    pass

        base = pdbid + " " + chainid + " "

        cic = self.cic
        if picFlags & icr.pic_flags.hedra or picFlags & icr.pic_flags.tau:
            for h in sorted(self.hedra.values()):
                if (
                    not picFlags & icr.pic_flags.hedra  # not all hedra
                    and picFlags & icr.pic_flags.tau  # but yes tau hedron
                    and h.e_class != "NCAC"  # and is not tau
                ):
                    continue
                if hCut is not None:
                    hc = h.xrh_class if hasattr(h, "xrh_class") else h.e_class
                    if hc in hedra_defaults and hedra_defaults[hc][1] <= hCut:
                        continue
                hndx = h.ndx
                try:
                    s += (
                        base
                        + h.id
                        + " "
                        + f"{cic.hedraL12[hndx]:{pAcc}} {cic.hedraAngle[hndx]:{pAcc}} {cic.hedraL23[hndx]:{pAcc}}"
                        + "\n"
                    )
                except KeyError:
                    pass
        for d in sorted(self.dihedra.values()):
            if d.primary:
                if not picFlags & icr.pic_flags.primary:
                    # primary and not primary flag so keep checking filters
                    # db = d.bits()
                    if not picFlags & d.bits():
                        continue
            elif not picFlags & icr.pic_flags.secondary:
                continue  # secondary and flag not set -> skip

            if pCut is not None:
                if (
                    d.primary
                    and d.pclass in dihedra_primary_defaults
                    and dihedra_primary_defaults[d.pclass][1] <= pCut
                ):
                    continue
            try:
                s += base + d.id + " " + f"{cic.dihedraAngle[d.ndx]:{pAcc}}" + "\n"
            except KeyError:
                pass

        if picFlags & icr.pic_flags.bFactors:
            col = 0
            for a in sorted(self.residue.get_atoms()):
                if 2 == a.is_disordered():
                    if IC_Residue.no_altloc or self.alt_ids is None:
                        s, col = self._write_pic_bfac(a.selected_child, s, col)
                    else:
                        for atm in a.child_dict.values():
                            s, col = self._write_pic_bfac(atm, s, col)
                else:
                    s, col = self._write_pic_bfac(a, s, col)
            if 0 != col % 5:
                s += "\n"

        return s

    def _get_ak_tuple(self, ak_str: str) -> Optional[Tuple["AtomKey", ...]]:
        """Convert atom pair string to AtomKey tuple.

        :param str ak_str:
            Two atom names separated by ':', e.g. 'N:CA'
            Optional position specifier relative to self,
            e.g. '-1C:N' for preceding peptide bond.
        """
        AK = AtomKey
        S = self
        angle_key2 = []
        akstr_list = ak_str.split(":")
        lenInput = len(akstr_list)
        for a in akstr_list:
            m = self._relative_atom_re.match(a)
            if m:
                if m.group(1) == "-1":
                    if 0 < len(S.rprev):
                        angle_key2.append(AK(S.rprev[0], m.group(2)))
                elif m.group(1) == "1":
                    if 0 < len(S.rnext):
                        angle_key2.append(AK(S.rnext[0], m.group(2)))
                elif m.group(1) == "0":
                    angle_key2.append(self.rak(m.group(2)))
            else:
                angle_key2.append(self.rak(a))
        if len(angle_key2) != lenInput:
            return None
        return tuple(angle_key2)

    _relative_atom_re = re.compile(r"^(-?[10])([A-Z]+)$")

    def _get_angle_for_tuple(
        self, angle_key: EKT
    ) -> Optional[Union["Hedron", "Dihedron"]]:
        len_mkey = len(angle_key)
        rval: Optional[Union["Hedron", "Dihedron"]]
        if 4 == len_mkey:
            rval = self.dihedra.get(cast(DKT, angle_key), None)
        elif 3 == len_mkey:
            rval = self.hedra.get(cast(HKT, angle_key), None)
        else:
            return None
        return rval

    # @profile
    def pick_angle(
        self, angle_key: Union[EKT, str]
    ) -> Optional[Union["Hedron", "Dihedron"]]:
        """Get Hedron or Dihedron for angle_key.

        :param angle_key:
            - tuple of 3 or 4 AtomKeys
            - string of atom names ('CA') separated by :'s
            - string of [-1, 0, 1]<atom name> separated by ':'s. -1 is
              previous residue, 0 is this residue, 1 is next residue
            - psi, phi, omg, omega, chi1, chi2, chi3, chi4, chi5
            - tau (N-CA-C angle) see Richardson1981
            - tuples of AtomKeys is only access for alternate disordered atoms

        Observe that a residue's phi and omega dihedrals, as well as the hedra
        comprising them (including the N:Ca:C `tau` hedron), are stored in the
        n-1 di/hedra sets; this overlap is handled here, but may be an issue if
        accessing directly.

        The following print commands are equivalent (except for sidechains with
        non-carbon atoms for chi2)::

            ric = r.internal_coord
            print(
                r,
                ric.get_angle("psi"),
                ric.get_angle("phi"),
                ric.get_angle("omg"),
                ric.get_angle("tau"),
                ric.get_angle("chi2"),
            )
            print(
                r,
                ric.get_angle("N:CA:C:1N"),
                ric.get_angle("-1C:N:CA:C"),
                ric.get_angle("-1CA:-1C:N:CA"),
                ric.get_angle("N:CA:C"),
                ric.get_angle("CA:CB:CG:CD"),
            )

        See ic_data.py for detail of atoms in the enumerated sidechain angles
        and the backbone angles which do not span the peptide bond. Using 's'
        for current residue ('self') and 'n' for next residue, the spanning
        (overlapping) angles are::

                (sN, sCA, sC, nN)   # psi
                (sCA, sC, nN, nCA)  # omega i+1
                (sC, nN, nCA, nC)   # phi i+1
                (sCA, sC, nN)
                (sC, nN, nCA)
                (nN, nCA, nC)       # tau i+1

        :return: Matching Hedron, Dihedron, or None.
        """
        rval: Optional[Union["Hedron", "Dihedron"]] = None
        if isinstance(angle_key, tuple):
            rval = self._get_angle_for_tuple(angle_key)
            if rval is None and self.rprev:
                rval = self.rprev[0]._get_angle_for_tuple(angle_key)
        elif ":" in angle_key:
            angle_key = cast(EKT, self._get_ak_tuple(cast(str, angle_key)))
            if angle_key is None:
                return None
            rval = self._get_angle_for_tuple(angle_key)
            if rval is None and self.rprev:
                rval = self.rprev[0]._get_angle_for_tuple(angle_key)
        elif "psi" == angle_key:
            if 0 == len(self.rnext):
                return None
            rn = self.rnext[0]
            sN, sCA, sC = self.rak("N"), self.rak("CA"), self.rak("C")
            nN = rn.rak("N")
            rval = self.dihedra.get((sN, sCA, sC, nN), None)
        elif "phi" == angle_key:
            if 0 == len(self.rprev):
                return None
            rp = self.rprev[0]
            pC, sN, sCA = rp.rak("C"), self.rak("N"), self.rak("CA")
            sC = self.rak("C")
            rval = rp.dihedra.get((pC, sN, sCA, sC), None)
        elif "omg" == angle_key or "omega" == angle_key:
            if 0 == len(self.rprev):
                return None
            rp = self.rprev[0]
            pCA, pC, sN = rp.rak("CA"), rp.rak("C"), self.rak("N")
            sCA = self.rak("CA")
            rval = rp.dihedra.get((pCA, pC, sN, sCA), None)
        elif "tau" == angle_key:
            sN, sCA, sC = self.rak("N"), self.rak("CA"), self.rak("C")
            rval = self.hedra.get((sN, sCA, sC), None)
            if rval is None and 0 != len(self.rprev):
                rp = self.rprev[0]  # tau in prev residue for all but first
                rval = rp.hedra.get((sN, sCA, sC), None)
        elif angle_key.startswith("chi"):
            sclist = ic_data_sidechains.get(self.lc, None)
            if sclist is None:
                return None
            ndx = (2 * int(angle_key[-1])) - 1
            try:
                akl = sclist[ndx]
                if akl[4] == angle_key:
                    klst = [self.rak(a) for a in akl[0:4]]
                    tklst = cast(DKT, tuple(klst))
                    rval = self.dihedra.get(tklst, None)
                else:
                    return None
            except IndexError:
                return None

        return rval

    def get_angle(self, angle_key: Union[EKT, str]) -> Optional[float]:
        """Get dihedron or hedron angle for specified key.

        See :meth:`.pick_angle` for key specifications.
        """
        edron = self.pick_angle(angle_key)
        if edron:
            return edron.angle
        return None

    def set_angle(self, angle_key: Union[EKT, str], v: float, overlap=True):
        """Set dihedron or hedron angle for specified key.

        If angle is a `Dihedron` and `overlap` is True (default), overlapping
        dihedra are also changed as appropriate.  The overlap is a result of
        protein chain definitions in :mod:`.ic_data` and :meth:`_create_edra`
        (e.g. psi overlaps N-CA-C-O).

        Te default overlap=True is probably what you want for:
        `set_angle("chi1", val)`

        The default is probably NOT what you want when processing all dihedrals
        in a chain or residue (such as copying from another structure), as the
        overlaping dihedra will likely be in the set as well.

        N.B. setting e.g. PRO chi2 is permitted without error or warning!

        See :meth:`.pick_angle` for angle_key specifications.
        See :meth:`.bond_rotate` to change a dihedral by a number of degrees

        :param angle_key: angle identifier.
        :param float v: new angle in degrees (result adjusted to +/-180).
        :param bool overlap: default True.
            Modify overlapping dihedra as needed
        """
        edron = self.pick_angle(angle_key)
        if edron is None:
            return
        elif isinstance(edron, Hedron) or not overlap:
            edron.angle = v
        else:  # Dihedron, do overlap angles
            delta = Dihedron.angle_dif(edron.angle, v)
            self._do_bond_rotate(edron, delta)

    def _do_bond_rotate(self, base: "Dihedron", delta: float):
        """Find and modify related dihedra through id3_dh_index."""
        try:
            for dk in self.cic.id3_dh_index[base.id3]:
                # change all diheds with same first hedron
                dihed = self.cic.dihedra[dk]
                dihed.angle += delta  # +/- 180 handled in setter
                # for changed dihed, change any with reverse key 2nd hedron
                # so change N-Ca-C-N will change O-Ca-C-Cb
                try:
                    for d2rk in self.cic.id3_dh_index[dihed.id32[::-1]]:
                        self.cic.dihedra[d2rk].angle += delta
                except KeyError:
                    pass
        except AttributeError:
            raise RuntimeError("bond_rotate, bond_set only for dihedral angles")

    def bond_rotate(self, angle_key: Union[EKT, str], delta: float):
        """Rotate set of overlapping dihedrals by delta degrees.

        Changes a dihedral angle by a given delta, i.e.
        new_angle = current_angle + delta
        Values are adjusted so new_angle iwll be within +/-180.

        Changes overlapping dihedra as in :meth:`.set_angle`

        See :meth:`.pick_angle` for key specifications.
        """
        base = self.pick_angle(angle_key)
        if base is not None:
            self._do_bond_rotate(base, delta)

    def bond_set(self, angle_key: Union[EKT, str], val: float):
        """Set dihedron to val, update overlapping dihedra by same amount.

        Redundant to :meth:`.set_angle`, retained for compatibility.  Unlike
        :meth:`.set_angle` this is for dihedra only and no option to not update
        overlapping dihedra.

        See :meth:`.pick_angle` for key specifications.
        """
        base = self.pick_angle(angle_key)
        if base is not None:
            delta = Dihedron.angle_dif(base.angle, val)
            self._do_bond_rotate(base, delta)

    def pick_length(
        self, ak_spec: Union[str, BKT]
    ) -> Tuple[Optional[List["Hedron"]], Optional[BKT]]:
        """Get list of hedra containing specified atom pair.

        :param ak_spec:
            - tuple of two AtomKeys
            - string: two atom names separated by ':', e.g. 'N:CA' with
              optional position specifier relative to self, e.g. '-1C:N' for
              preceding peptide bond.  Position specifiers are -1, 0, 1.

        The following are equivalent::

            ric = r.internal_coord
            print(
                r,
                ric.get_length("0C:1N"),
            )
            print(
                r,
                None
                if not ric.rnext
                else ric.get_length((ric.rak("C"), ric.rnext[0].rak("N"))),
            )

        If atom not found on current residue then will look on rprev[0] to
        handle cases like Gly N:CA.  For finer control please access
        `IC_Chain.hedra` directly.

        :return: list of hedra containing specified atom pair as tuples of
                AtomKeys
        """
        rlst: List[Hedron] = []
        # if ":" in ak_spec:
        if isinstance(ak_spec, str):
            ak_spec = cast(BKT, self._get_ak_tuple(ak_spec))
        if ak_spec is None:
            return None, None
        for hed_key, hed_val in self.hedra.items():
            if all(ak in hed_key for ak in ak_spec):
                rlst.append(hed_val)
        # handle bonds stored on rprev, e.g. set backbone, read gly N:CA
        for rp in self.rprev:
            for hed_key, hed_val in rp.hedra.items():
                if all(ak in hed_key for ak in ak_spec):
                    rlst.append(hed_val)
        return rlst, ak_spec

    def get_length(self, ak_spec: Union[str, BKT]) -> Optional[float]:
        """Get bond length for specified atom pair.

        See :meth:`.pick_length` for ak_spec and details.
        """
        hed_lst, ak_spec2 = self.pick_length(ak_spec)
        if hed_lst is None or ak_spec2 is None:
            return None

        for hed in hed_lst:
            val = hed.get_length(ak_spec2)
            if val is not None:
                return val
        return None

    def set_length(self, ak_spec: Union[str, BKT], val: float) -> None:
        """Set bond length for specified atom pair.

        See :meth:`.pick_length` for ak_spec.
        """
        hed_lst, ak_spec2 = self.pick_length(ak_spec)
        if hed_lst is not None and ak_spec2 is not None:
            for hed in hed_lst:
                hed.set_length(ak_spec2, val)

    def applyMtx(self, mtx: np.array) -> None:
        """Apply matrix to atom_coords for this IC_Residue."""
        aa = self.cic.atomArray
        aai = self.cic.atomArrayIndex
        rpndx = AtomKey.fields.respos
        rp = str(self.rbase[0])
        aselect = [aai.get(k) for k in aai.keys() if k.akl[rpndx] == rp]
        aas = aa[aselect]
        # numpy will broadcast the transform matrix over all points if dot()
        # applied in this order
        aa[aselect] = aas.dot(mtx.transpose())
        """
        # slower way, one at a time
        for ak in sorted(self.ak_set):
            ndx = self.cic.atomArrayIndex[ak]
            self.cic.atomArray[ndx] = mtx.dot(self.cic.atomArray[ndx])
        """


class Edron:
    """Base class for Hedron and Dihedron classes.

    Supports rich comparison based on lists of AtomKeys.

    Attributes
    ----------
    atomkeys: tuple
        3 (hedron) or 4 (dihedron) :class:`.AtomKey` s defining this di/hedron
    id: str
        ':'-joined string of AtomKeys for this di/hedron
    needs_update: bool
        indicates di/hedron local atom_coords do NOT reflect current di/hedron
        angle and length values in hedron local coordinate space
    e_class: str
        sequence of atoms (no position or residue) comprising di/hedron
        for statistics
    re_class: str
        sequence of residue, atoms comprising di/hedron for statistics
    cre_class: str
        sequence of covalent radii classses comprising di/hedron for statistics
    edron_re: compiled regex (Class Attribute)
        A compiled regular expression matching string IDs for Hedron
        and Dihedron objects
    cic: IC_Chain reference
        Chain internal coords object containing this hedron
    ndx: int
        index into IC_Chain level numpy data arrays for di/hedra.
        Set in :meth:`IC_Chain.init_edra`
    rc: int
        number of residues involved in this edron

    Methods
    -------
    gen_key([AtomKey, ...] or AtomKey, ...) (Static Method)
        generate a ':'-joined string of AtomKey Ids
    is_backbone()
        Return True if all atomkeys atoms are N, Ca, C or O

    """

    # regular expression to capture hedron and dihedron specifications, as in
    #  .pic files
    edron_re = re.compile(
        # pdbid and chain id
        r"^(?P<pdbid>\w+)?\s(?P<chn>[\w|\s])?\s"
        # 3 atom specifiers for hedron
        r"(?P<a1>[\w\-\.]+):(?P<a2>[\w\-\.]+):(?P<a3>[\w\-\.]+)"
        # 4th atom specifier for dihedron
        r"(:(?P<a4>[\w\-\.]+))?"
        r"\s+"
        # len-angle-len for hedron
        r"(((?P<len12>\S+)\s+(?P<angle>\S+)\s+(?P<len23>\S+)\s*$)|"
        # dihedral angle for dihedron
        r"((?P<dihedral>\S+)\s*$))"
    )
    """ A compiled regular expression matching string IDs for Hedron and
    Dihedron objects"""

    @staticmethod
    def gen_key(lst: List["AtomKey"]) -> str:
        """Generate string of ':'-joined AtomKey strings from input.

        Generate '2_A_C:3_P_N:3_P_CA' from (2_A_C, 3_P_N, 3_P_CA)
        :param list lst: list of AtomKey objects
        """
        if 4 == len(lst):
            return f"{lst[0].id}:{lst[1].id}:{lst[2].id}:{lst[3].id}"
        else:
            return f"{lst[0].id}:{lst[1].id}:{lst[2].id}"

    @staticmethod
    def gen_tuple(akstr: str) -> Tuple:
        """Generate AtomKey tuple for ':'-joined AtomKey string.

        Generate (2_A_C, 3_P_N, 3_P_CA) from '2_A_C:3_P_N:3_P_CA'
        :param str akstr: string of ':'-separated AtomKey strings
        """
        return tuple([AtomKey(i) for i in akstr.split(":")])

    # @profile
    def __init__(self, *args: Union[List["AtomKey"], EKT], **kwargs: str) -> None:
        """Initialize Edron with sequence of AtomKeys.

        Acceptable input:

            [ AtomKey, ... ]  : list of AtomKeys
            AtomKey, ...      : sequence of AtomKeys as args
            {'a1': str, 'a2': str, ... }  : dict of AtomKeys as 'a1', 'a2' ...
        """
        atomkeys: List[AtomKey] = []
        for arg in args:
            if isinstance(arg, list):
                atomkeys = arg
            elif isinstance(arg, tuple):
                atomkeys = list(arg)
            else:
                if arg is not None:
                    atomkeys.append(arg)
        if [] == atomkeys and all(k in kwargs for k in ("a1", "a2", "a3")):
            atomkeys = [
                AtomKey(kwargs["a1"]),
                AtomKey(kwargs["a2"]),
                AtomKey(kwargs["a3"]),
            ]
            if "a4" in kwargs and kwargs["a4"] is not None:
                atomkeys.append(AtomKey(kwargs["a4"]))

        self.atomkeys = tuple(atomkeys)
        self.id = Edron.gen_key(atomkeys)
        self._hash = hash(self.atomkeys)

        # flag indicating that atom coordinates are up to date
        # (do not need to be recalculated from angle and or length values)
        self.needs_update = True

        # IC_Chain which contains this di/hedron
        self.cic: IC_Chain  # set in :meth:`IC_Residue._link_dihedra`

        # no residue or position, just atoms
        self.e_class = ""
        # same but residue specific
        self.re_class = ""
        self.cre_class = ""
        rset = set()  # what residues this involves

        atmNdx = AtomKey.fields.atm
        resNdx = AtomKey.fields.resname
        resPos = AtomKey.fields.respos
        icode = AtomKey.fields.icode

        for ak in atomkeys:
            akl = ak.akl
            self.e_class += akl[atmNdx]
            self.re_class += akl[resNdx] + akl[atmNdx]
            rset.add(akl[resPos] + (akl[icode] or ""))
            self.cre_class += ak.cr_class()

        self.rc = len(rset)

    def __deepcopy__(self, memo):
        """Deep copy implementation for Edron."""
        existing = memo.get(id(self), False)
        if existing:
            return existing
        dup = type(self).__new__(self.__class__)
        memo[id(self)] = dup
        dup.__dict__.update(self.__dict__)  # mostly static attribs
        dup.cic = memo[id(self.cic)]
        dup.atomkeys = copy.deepcopy(self.atomkeys, memo)
        return dup

    def __contains__(self, ak: "AtomKey") -> bool:
        """Return True if atomkey is in this edron."""
        return ak in self.atomkeys

    def is_backbone(self) -> bool:
        """Report True for contains only N, C, CA, O, H atoms."""
        return all(ak.is_backbone() for ak in self.atomkeys)

    def __repr__(self) -> str:
        """Tuple of AtomKeys is default repr string."""
        return str(self.atomkeys)

    def __hash__(self) -> int:
        """Hash calculated at init from atomkeys tuple."""
        return self._hash

    def _cmp(self, other: "Edron") -> Union[Tuple["AtomKey", "AtomKey"], bool]:
        """Comparison function ranking self vs. other; False on equal.

        Priority is lowest value for sort: psi < chi1.
        """
        for ak_s, ak_o in zip(self.atomkeys, other.atomkeys):
            if ak_s != ak_o:
                return ak_s, ak_o
        return False

    def __eq__(self, other: object) -> bool:
        """Test for equality."""
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.id == other.id

    def __ne__(self, other: object) -> bool:
        """Test for inequality."""
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.id != other.id

    def __gt__(self, other: object) -> bool:
        """Test greater than."""
        if not isinstance(other, type(self)):
            return NotImplemented
        rslt = self._cmp(other)
        if rslt:
            rslt = cast(Tuple[AtomKey, AtomKey], rslt)
            return rslt[0] > rslt[1]
        return False

    def __ge__(self, other: object) -> bool:
        """Test greater or equal."""
        if not isinstance(other, type(self)):
            return NotImplemented
        rslt = self._cmp(other)
        if rslt:
            rslt = cast(Tuple[AtomKey, AtomKey], rslt)
            return rslt[0] >= rslt[1]
        return True

    def __lt__(self, other: object) -> bool:
        """Test less than."""
        if not isinstance(other, type(self)):
            return NotImplemented
        rslt = self._cmp(other)
        if rslt:
            rslt = cast(Tuple[AtomKey, AtomKey], rslt)
            return rslt[0] < rslt[1]
        return False

    def __le__(self, other: object) -> bool:
        """Test less or equal."""
        if not isinstance(other, type(self)):
            return NotImplemented
        rslt = self._cmp(other)
        if rslt:
            rslt = cast(Tuple[AtomKey, AtomKey], rslt)
            return rslt[0] <= rslt[1]
        return True


class Hedron(Edron):
    """Class to represent three joined atoms forming a plane.

    Contains atom coordinates in local coordinate space: central atom
    at origin, one terminal atom on XZ plane, and the other on the +Z axis.
    Stored in two orientations, with the 3rd (forward) or first (reversed)
    atom on the +Z axis.  See :class:`Dihedron` for use of forward and
    reverse orientations.

    Attributes
    ----------
    len12: float
        distance between first and second atoms
    len23: float
        distance between second and third atoms
    angle: float
        angle (degrees) formed by three atoms in hedron
    xrh_class: string
        only for hedron spanning 2 residues, will have 'X' for residue
        contributing only one atom

    Methods
    -------
    get_length()
        get bond length for specified atom pair
    set_length()
        set bond length for specified atom pair
    angle(), len12(), len23()
        setters for relevant attributes (angle in degrees)
    """

    def __init__(self, *args: Union[List["AtomKey"], HKT], **kwargs: str) -> None:
        """Initialize Hedron with sequence of AtomKeys, kwargs.

        Acceptable input:
            As for Edron, plus optional 'len12', 'angle', 'len23'
            keyworded values.
        """
        super().__init__(*args, **kwargs)
        if self.rc == 2:  # hedron crosses residue boundary
            resPos = AtomKey.fields.respos
            icode = AtomKey.fields.icode
            resNdx = AtomKey.fields.resname
            atmNdx = AtomKey.fields.atm
            akl0, akl1 = self.atomkeys[0].akl, self.atomkeys[1].akl
            if akl0[resPos] != akl1[resPos] or akl0[icode] != akl1[icode]:
                self.xrh_class = "X" + self.re_class[1:]
            else:
                xrhc = ""
                for i in range(2):
                    xrhc += self.atomkeys[i].akl[resNdx] + self.atomkeys[i].akl[atmNdx]
                self.xrh_class = xrhc + "X" + self.atomkeys[2].akl[atmNdx]

    # __deepcopy__ covered by Edron superclass

    def __repr__(self) -> str:
        """Print string for Hedron object."""
        return (
            f"3-{self.id} {self.re_class} {self.len12!s} "
            f"{self.angle!s} {self.len23!s}"
        )

    @property
    def angle(self) -> float:
        """Get this hedron angle."""
        try:
            return self.cic.hedraAngle[self.ndx]
        except AttributeError:
            return 0.0

    def _invalidate_atoms(self):
        self.cic.hAtoms_needs_update[self.ndx] = True
        for ak in self.atomkeys:
            self.cic.atomArrayValid[self.cic.atomArrayIndex[ak]] = False

    @angle.setter
    def angle(self, angle_deg) -> None:
        """Set this hedron angle; sets needs_update."""
        self.cic.hedraAngle[self.ndx] = angle_deg
        self.cic.hAtoms_needs_update[self.ndx] = True
        self.cic.atomArrayValid[self.cic.atomArrayIndex[self.atomkeys[2]]] = False

    @property
    def len12(self):
        """Get first length for Hedron."""
        try:
            return self.cic.hedraL12[self.ndx]
        except AttributeError:
            return 0.0

    @len12.setter
    def len12(self, len):
        """Set first length for Hedron; sets needs_update."""
        self.cic.hedraL12[self.ndx] = len
        self.cic.hAtoms_needs_update[self.ndx] = True
        self.cic.atomArrayValid[self.cic.atomArrayIndex[self.atomkeys[1]]] = False
        self.cic.atomArrayValid[self.cic.atomArrayIndex[self.atomkeys[2]]] = False

    @property
    def len23(self) -> float:
        """Get second length for Hedron."""
        try:
            return self.cic.hedraL23[self.ndx]
        except AttributeError:
            return 0.0

    @len23.setter
    def len23(self, len):
        """Set second length for Hedron; sets needs_update."""
        self.cic.hedraL23[self.ndx] = len
        self.cic.hAtoms_needs_update[self.ndx] = True
        self.cic.atomArrayValid[self.cic.atomArrayIndex[self.atomkeys[2]]] = False

    def get_length(self, ak_tpl: BKT) -> Optional[float]:
        """Get bond length for specified atom pair.

        :param tuple ak_tpl: tuple of AtomKeys.
            Pair of atoms in this Hedron
        """
        if 2 > len(ak_tpl):
            return None
        if all(ak in self.atomkeys[:2] for ak in ak_tpl):
            return self.cic.hedraL12[self.ndx]
        if all(ak in self.atomkeys[1:] for ak in ak_tpl):
            return self.cic.hedraL23[self.ndx]
        return None

    def set_length(self, ak_tpl: BKT, newLength: float):
        """Set bond length for specified atom pair; sets needs_update.

        :param tuple .ak_tpl: tuple of AtomKeys
            Pair of atoms in this Hedron
        """
        if 2 > len(ak_tpl):
            raise TypeError(f"Require exactly 2 AtomKeys: {ak_tpl!s}")
        elif all(ak in self.atomkeys[:2] for ak in ak_tpl):
            self.cic.hedraL12[self.ndx] = newLength
        elif all(ak in self.atomkeys[1:] for ak in ak_tpl):
            self.cic.hedraL23[self.ndx] = newLength
        else:
            raise TypeError("%s not found in %s" % (str(ak_tpl), self))
        self._invalidate_atoms()


class Dihedron(Edron):
    """Class to represent four joined atoms forming a dihedral angle.

    Attributes
    ----------
    angle: float
        Measurement or specification of dihedral angle in degrees; prefer
        :meth:`IC_Residue.bond_set` to set
    hedron1, hedron2: Hedron object references
        The two hedra which form the dihedral angle
    h1key, h2key: tuples of AtomKeys
        Hash keys for hedron1 and hedron2
    id3,id32: tuples of AtomKeys
        First 3 and second 3 atoms comprising dihedron; hxkey orders may differ
    ric: IC_Residue object reference
        :class:`.IC_Residue` object containing this dihedral
    reverse: bool
        Indicates order of atoms in dihedron is reversed from order of atoms
        in hedra
    primary: bool
        True if this is psi, phi, omega or a sidechain chi angle
    pclass: string (primary angle class)
        re_class with X for adjacent residue according to nomenclature
        (psi, omega, phi)
    cst, rcst: numpy [4][4] arrays
        transformations to (cst) and from (rcst) Dihedron coordinate space
        defined with atom 2 (Hedron 1 center atom) at the origin.  Views on
        :data:`IC_Chain.dCoordSpace`.

    Methods
    -------
    angle()
        getter/setter for dihdral angle in degrees; prefer
        :meth:`IC_Residue.bond_set`
    bits()
        return :data:`IC_Residue.pic_flags` bitmask for dihedron psi, omega, etc
    """

    def __init__(self, *args: Union[List["AtomKey"], DKT], **kwargs: str) -> None:
        """Init Dihedron with sequence of AtomKeys and optional dihedral angle.

        Acceptable input:
            As for Edron, plus optional 'dihedral' keyworded angle value.
        """
        super().__init__(*args, **kwargs)

        # hedra making up this dihedron; set by self:_set_hedra()
        self.hedron1: Hedron  # = None
        self.hedron2: Hedron  # = None

        self.h1key: HKT  # = None
        self.h2key: HKT  # = None

        # h1, h2key above may be reversed; id3,2 will not be

        self.id3: HKT = cast(HKT, tuple(self.atomkeys[0:3]))
        self.id32: HKT = cast(HKT, tuple(self.atomkeys[1:4]))

        self._setPrimary()

        # IC_Residue object which includes this dihedron;
        # set by Residue:linkDihedra()
        self.ric: IC_Residue
        # order of atoms in dihedron is reversed from order of atoms in hedra
        self.reverse = False  # configured by :meth:`._set_hedra`

    def __repr__(self) -> str:
        """Print string for Dihedron object."""
        return f"4-{self.id!s} {self.re_class} {self.angle!s} {self.ric!s}"

    @staticmethod
    def _get_hedron(ic_res: IC_Residue, id3: HKT) -> Optional[Hedron]:
        """Find specified hedron on this residue or its adjacent neighbors."""
        hedron = ic_res.hedra.get(id3, None)
        if not hedron and 0 < len(ic_res.rprev):
            for rp in ic_res.rprev:
                hedron = rp.hedra.get(id3, None)
                if hedron is not None:
                    break
        if not hedron and 0 < len(ic_res.rnext):
            for rn in ic_res.rnext:
                hedron = rn.hedra.get(id3, None)
                if hedron is not None:
                    break
        return hedron

    def _setPrimary(self) -> bool:
        """Mark dihedra required for psi, phi, omega, chi and other angles."""
        # http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
        dhc = self.e_class
        if dhc == "NCACN":  # psi
            self.pclass = self.re_class[0:7] + "XN"
            self.primary = True
        elif dhc == "CACNCA":  # omg
            self.pclass = "XCAXC" + self.re_class[5:]
            self.primary = True
        elif dhc == "CNCAC":  # phi
            self.pclass = "XC" + self.re_class[2:]
            self.primary = True
        elif dhc == "CNCACB":  # alternate Cbeta locator
            self.altCB_class = "XC" + self.re_class[2:]
            self.primary = False
        elif dhc in primary_angles:
            self.primary = True
            self.pclass = self.re_class
        else:
            self.primary = False

    def _set_hedra(self) -> Tuple[bool, Hedron, Hedron]:
        """Work out hedra keys and set rev flag."""
        try:
            return self.rev, self.hedron1, self.hedron2
        except AttributeError:
            pass

        rev = False
        res = self.ric
        h1key = self.id3
        hedron1 = Dihedron._get_hedron(res, h1key)
        if not hedron1:
            rev = True
            h1key = cast(HKT, tuple(self.atomkeys[2::-1]))
            hedron1 = Dihedron._get_hedron(res, h1key)
            h2key = cast(HKT, tuple(self.atomkeys[3:0:-1]))
        else:
            h2key = self.id32

        if not hedron1:
            raise HedronMatchError(
                f"can't find 1st hedron for key {h1key} dihedron {self}"
            )

        hedron2 = Dihedron._get_hedron(res, h2key)

        if not hedron2:
            raise HedronMatchError(
                f"can't find 2nd hedron for key {h2key} dihedron {self}"
            )

        self.hedron1 = hedron1
        self.h1key = h1key
        self.hedron2 = hedron2
        self.h2key = h2key

        self.reverse = rev

        return rev, hedron1, hedron2

    @property
    def angle(self) -> float:
        """Get dihedral angle."""
        try:
            return self.cic.dihedraAngle[self.ndx]
        except AttributeError:
            try:
                return self._dihedral
            except AttributeError:
                return 360.0  # error value without type hint hassles

    @angle.setter
    def angle(self, dangle_deg_in: float) -> None:
        """Save new dihedral angle; sets needs_update.

        Faster to modify IC_Chain level arrays directly.

        This is probably not the routine you are looking for.  See
        :meth:`IC_Residue.set_angle` or :meth:`IC_Residue.bond_rotate` to change
        a dihedral angle along with its overlapping dihedra, i.e. without
        clashing atoms.

        N.B. dihedron (i-1)C-N-CA-CB is ignored if O exists.
        C-beta is by default placed using O-C-CA-CB, but O is missing
        in some PDB file residues, which means the sidechain cannot be
        placed.  The alternate CB path (i-1)C-N-CA-CB is provided to
        circumvent this, but if this is needed then it must be adjusted in
        conjunction with PHI ((i-1)C-N-CA-C) as they overlap.  This is handled
        by the `IC_Residue` routines above.

        :param float dangle_deg: new dihedral angle in degrees
        """
        if dangle_deg_in > 180.0:
            dangle_deg = dangle_deg_in - 360.0
        elif dangle_deg_in < -180.0:
            dangle_deg = dangle_deg_in + 360.0
        else:
            dangle_deg = dangle_deg_in

        self._dihedral = dangle_deg
        self.needs_update = True

        cic = self.cic
        dndx = self.ndx
        cic.dihedraAngle[dndx] = dangle_deg
        cic.dihedraAngleRads[dndx] = np.deg2rad(dangle_deg)
        cic.dAtoms_needs_update[dndx] = True
        cic.atomArrayValid[cic.atomArrayIndex[self.atomkeys[3]]] = False

    @staticmethod
    def angle_dif(a1: Union[float, np.ndarray], a2: Union[float, np.ndarray]):
        """Get angle difference between two +/- 180 angles.

        https://stackoverflow.com/a/36001014/2783487
        """
        return 180.0 - ((180.0 - a2) + a1) % 360.0

    @staticmethod
    def angle_avg(alst: List, in_rads: bool = False, out_rads: bool = False):
        """Get average of list of +/-180 angles.

        :param List alst: list of angles to average
        :param bool in_rads: input values are in radians
        :param bool out_rads: report result in radians
        """
        walst = alst if in_rads else np.deg2rad(alst)
        ravg = np.arctan2(np.sum(np.sin(walst)), np.sum(np.cos(walst)))
        return ravg if out_rads else np.rad2deg(ravg)

    @staticmethod
    def angle_pop_sd(alst: List, avg: float):
        """Get population standard deviation for list of +/-180 angles.

        should be sample std dev but avoid len(alst)=1 -> div by 0
        """
        return np.sqrt(np.sum(np.square(Dihedron.angle_dif(alst, avg))) / len(alst))

    def difference(self, other: "Dihedron") -> float:
        """Get angle difference between this and other +/- 180 angles."""
        return Dihedron.angle_dif(self.angle, other.angle)

    def bits(self) -> int:
        """Get :data:`IC_Residue.pic_flags` bitmasks for self is psi, omg, phi, pomg, chiX."""
        icr = IC_Residue
        if self.e_class == "NCACN":
            # i psi
            return icr.pic_flags.psi
        elif hasattr(self, "pclass") and self.pclass == "XCAXCPNPCA":
            # i+1 is pro so i+1 omg
            return icr.pic_flags.omg | icr.pic_flags.pomg
        elif self.e_class == "CACNCA":
            # i+1 omg
            return icr.pic_flags.omg
        elif self.e_class == "CNCAC":
            # i+1 phi
            return icr.pic_flags.phi
        else:
            # i chiX
            atmNdx = AtomKey.fields.atm
            scList = ic_data_sidechains.get(self.ric.lc)
            aLst = tuple(ak.akl[atmNdx] for ak in self.atomkeys)
            for e in scList:
                if len(e) != 5:  # only chi entries have label at [4]
                    continue
                if aLst == e[0:4]:
                    return icr.pic_flags.chi1 << (int(e[4][-1]) - 1)
        return 0


class AtomKey:
    """Class for dict keys to reference atom coordinates.

    AtomKeys capture residue and disorder information together, and
    provide a no-whitespace string key for .pic files.

    Supports rich comparison and multiple ways to instantiate.

    AtomKeys contain:
     residue position (respos), insertion code (icode), 1 or 3 char residue
     name (resname), atom name (atm), altloc (altloc), and occupancy (occ)

    Use :data:`AtomKey.fields` to get the index to the component of interest by
    name:

    Get C-alpha atoms from IC_Chain atomArray and atomArrayIndex with
    AtomKeys::

        atmNameNdx = internal_coords.AtomKey.fields.atm
        CaSelection = [
            atomArrayIndex.get(k)
            for k in atomArrayIndex.keys()
            if k.akl[atmNameNdx] == "CA"
        ]
        AtomArrayCa = atomArray[CaSelection]

    Get all phenylalanine atoms in a chain::

        resNameNdx = internal_coords.AtomKey.fields.resname
        PheSelection = [
            atomArrayIndex.get(k)
            for k in atomArrayIndex.keys()
            if k.akl[resNameNdx] == "F"
        ]
        AtomArrayPhe = atomArray[PheSelection]

    'resname' will be the uppercase 1-letter amino acid code if one of the 20
    standard residues, otherwise the supplied 3-letter code.  Supplied as input
    or read from .rbase attribute of :class:`IC_Residue`.

    Attributes
    ----------
    akl: tuple
        All six fields of AtomKey
    fieldNames: tuple (Class Attribute)
        Mapping of key index positions to names
    fields: namedtuple (Class Attribute)
        Mapping of field names to index positions.
    id: str
        '_'-joined AtomKey fields, excluding 'None' fields
    atom_re: compiled regex (Class Attribute)
        A compiled regular expression matching the string form of the key
    d2h: bool (Class Attribute) default False
        Convert D atoms to H on input if True; must also modify
        :data:`IC_Residue.accept_atoms`
    missing: bool default False
        AtomKey __init__'d from string is probably missing, set this flag to
        note the issue.  Set by :meth:`.IC_Residue.rak`
    ric: IC_Residue default None
        *If* initialised with IC_Residue, this references the IC_residue

    Methods
    -------
    altloc_match(other)
        Returns True if this AtomKey matches other AtomKey excluding altloc
        and occupancy fields
    is_backbone()
        Returns True if atom is N, CA, C, O or H
    atm()
        Returns atom name, e.g. N, CA, CB, etc.
    cr_class()
        Returns covalent radii class e.g. Csb

    """

    atom_re = re.compile(
        r"^(?P<respos>-?\d+)(?P<icode>[A-Za-z])?"
        r"_(?P<resname>[a-zA-Z]+)_(?P<atm>[A-Za-z0-9]+)"
        r"(?:_(?P<altloc>\w))?(?:_(?P<occ>-?\d\.\d+?))?$"
    )
    """Pre-compiled regular expression to match an AtomKey string."""

    _endnum_re = re.compile(r"\D+(\d+)$")

    # PDB altLoc = Character = [\w ] (any non-ctrl ASCII incl space)
    # PDB iCode = AChar = [A-Za-z]

    fieldNames = ("respos", "icode", "resname", "atm", "altloc", "occ")
    _fieldsDef = namedtuple(
        "_fieldsDef", ["respos", "icode", "resname", "atm", "altloc", "occ"]
    )
    fields = _fieldsDef(0, 1, 2, 3, 4, 5)
    """Use this namedtuple to access AtomKey fields.  See :class:`AtomKey`"""

    d2h = False
    """Set True to convert D Deuterium to H Hydrogen on input."""

    def __init__(
        self, *args: Union[IC_Residue, Atom, List, Dict, str], **kwargs: str
    ) -> None:
        """Initialize AtomKey with residue and atom data.

        Examples of acceptable input::

            (<IC_Residue>, 'CA', ...)    : IC_Residue with atom info
            (<IC_Residue>, <Atom>)       : IC_Residue with Biopython Atom
            ([52, None, 'G', 'CA', ...])  : list of ordered data fields
            (52, None, 'G', 'CA', ...)    : multiple ordered arguments
            ({respos: 52, icode: None, atm: 'CA', ...}) : dict with fieldNames
            (respos: 52, icode: None, atm: 'CA', ...) : kwargs with fieldNames
            52_G_CA, 52B_G_CA, 52_G_CA_0.33, 52_G_CA_B_0.33  : id strings
        """
        akl: List[Optional[str]] = []
        self.ric = None

        for arg in args:
            if isinstance(arg, str):
                if "_" in arg:
                    # AtomKey.icd["_"] += 1
                    # got atom key string, recurse with regex parse
                    m = self.atom_re.match(arg)
                    if m is not None:
                        if akl != []:  # [] != akl:
                            raise Exception(
                                "Atom Key init full key not first argument: " + arg
                            )
                        akl = list(map(m.group, AtomKey.fieldNames))
                else:
                    akl.append(arg)

            elif isinstance(arg, IC_Residue):
                if akl != []:
                    raise Exception("Atom Key init Residue not first argument")
                akl = list(arg.rbase)
                self.ric = arg
            elif isinstance(arg, Atom):
                if 3 != len(akl):
                    raise Exception("Atom Key init Atom before Residue info")
                akl.append(arg.name)
                if not IC_Residue.no_altloc:
                    altloc = arg.altloc
                    akl.append(altloc if altloc != " " else None)
                    occ = float(arg.occupancy)
                    akl.append(str(occ) if occ != 1.00 else None)
                else:
                    akl += [None, None]
            elif isinstance(arg, (list, tuple)):
                akl += arg
            elif isinstance(arg, dict):
                for k in AtomKey.fieldNames:
                    akl.append(arg.get(k, None))
            else:
                raise Exception("Atom Key init not recognised")

        # process kwargs, initialize occ and altloc to None
        for i in range(len(akl), 6):
            if len(akl) <= i:
                fld = kwargs.get(AtomKey.fieldNames[i])
                if fld is not None:
                    akl.append(fld)

        # tweak local akl to generate id string
        if isinstance(akl[0], Integral):
            akl[0] = str(akl[0])  # numeric residue position to string

        if self.d2h:
            atmNdx = AtomKey.fields.atm
            if akl[atmNdx][0] == "D":
                akl[atmNdx] = re.sub("D", "H", akl[atmNdx], count=1)

        self.id = "_".join(
            [
                "".join(filter(None, akl[:2])),
                str(akl[2]),  # exclude None
                "_".join(filter(None, akl[3:])),
            ]
        )

        akl += [None] * (6 - len(akl))

        self.akl = tuple(akl)
        self._hash = hash(self.akl)
        self.missing = False

    def __deepcopy__(self, memo):
        """Deep copy implementation for AtomKey."""
        # will fail if .ric not in memo
        existing = memo.get(id(self), False)
        if existing:
            return existing
        dup = type(self).__new__(self.__class__)
        memo[id(self)] = dup
        dup.__dict__.update(self.__dict__)  # all static attribs except .ric
        if self.ric is not None:
            dup.ric = memo[id(self.ric)]
        # deepcopy complete
        return dup

    def __repr__(self) -> str:
        """Repr string from id."""
        return self.id

    def __hash__(self) -> int:
        """Hash calculated at init from akl tuple."""
        return self._hash

    _backbone_sort_keys = {"N": 0, "CA": 1, "C": 2, "O": 3}

    _sidechain_sort_keys = {
        "CB": 1,
        "CG": 2,
        "CG1": 2,
        "OG": 2,
        "OG1": 2,
        "SG": 2,
        "CG2": 3,
        "CD": 4,
        "CD1": 4,
        "SD": 4,
        "OD1": 4,
        "ND1": 4,
        "CD2": 5,
        "ND2": 5,
        "OD2": 5,
        "CE": 6,
        "NE": 6,
        "CE1": 6,
        "OE1": 6,
        "NE1": 6,
        "CE2": 7,
        "OE2": 7,
        "NE2": 7,
        "CE3": 8,
        "CZ": 9,
        "CZ2": 9,
        "NZ": 9,
        "NH1": 10,
        "OH": 10,
        "CZ3": 10,
        "CH2": 11,
        "NH2": 11,
        "OXT": 12,
        "H": 13,
    }

    _greek_sort_keys = {"A": 0, "B": 1, "G": 2, "D": 3, "E": 4, "Z": 5, "H": 6}

    def altloc_match(self, other: "AtomKey") -> bool:
        """Test AtomKey match to other discounting occupancy and altloc."""
        if isinstance(other, type(self)):
            return self.akl[:4] == other.akl[:4]
        else:
            return NotImplemented

    def is_backbone(self) -> bool:
        """Return True if is N, C, CA, O, or H."""
        return self.akl[self.fields.atm] in ("N", "C", "CA", "O", "H")

    def atm(self) -> str:
        """Return atom name : N, CA, CB, O etc."""
        return self.akl[self.fields.atm]

    def cr_class(self) -> Union[str, None]:
        """Return covalent radii class for atom or None."""
        akl = self.akl
        atmNdx = self.fields.atm
        try:
            return residue_atom_bond_state["X"][akl[atmNdx]]
        except KeyError:
            try:
                resNdx = self.fields.resname
                return residue_atom_bond_state[akl[resNdx]][akl[atmNdx]]
            except KeyError:
                return "Hsb" if akl[atmNdx][0] == "H" else None

    # @profile
    def _cmp(self, other: "AtomKey") -> Tuple[int, int]:
        """Comparison function ranking self vs. other.

        Priority is lower value, i.e. (CA, CB) gives (0, 1) for sorting.
        """
        for i in range(6):
            s, o = self.akl[i], other.akl[i]
            if s != o:
                # insert_code, altloc can be None, deal with first
                if s is None and o is not None:
                    # no insert code before named insert code
                    return 0, 1
                elif o is None and s is not None:
                    return 1, 0
                # now we know s, o not None
                # s, o = cast(str, s), cast(str, o) # performance critical code

                if AtomKey.fields.atm != i:
                    # only sorting complications at atom level, occ.
                    # otherwise respos, insertion code will trigger
                    # before residue name
                    if AtomKey.fields.occ == i:
                        oi = int(float(s) * 100)
                        si = int(float(o) * 100)
                        return si, oi  # swap so higher occupancy comes first
                    elif AtomKey.fields.respos == i:
                        return int(s), int(o)
                    elif AtomKey.fields.resname == i:
                        sac, oac = (
                            self.akl[AtomKey.fields.altloc],
                            other.akl[AtomKey.fields.altloc],
                        )
                        if sac is not None:
                            if oac is not None:
                                return ord(sac), ord(oac)  # altloc over resname
                            else:  # sac has val and oac is None
                                return 1, 0
                        elif oac is not None:  # oac has val and sac is None
                            return 0, 1
                    # else:  # altloc
                    # fall through for altloc, resname with both altloc = None
                    return ord(s), ord(o)

                # atom names from here
                # backbone atoms before sidechain atoms

                sb = self._backbone_sort_keys.get(s, None)
                ob = self._backbone_sort_keys.get(o, None)
                if sb is not None and ob is not None:
                    return sb, ob
                elif sb is not None and ob is None:
                    return 0, 1
                elif sb is None and ob is not None:
                    return 1, 0
                # finished backbone and backbone vs. sidechain atoms

                # sidechain vs sidechain, sidechain vs H
                ss = self._sidechain_sort_keys.get(s, None)
                os = self._sidechain_sort_keys.get(o, None)
                if ss is not None and os is not None:
                    return ss, os
                elif ss is not None and os is None:
                    return 0, 1
                elif ss is None and os is not None:
                    return 1, 0

                # amide single 'H' captured above in sidechain sort
                # now 'complex'' hydrogens after sidechain
                s0, s1, o0, o1 = s[0], s[1], o[0], o[1]
                s1d, o1d = s1.isdigit(), o1.isdigit()
                # if "H" == s0 == o0: # breaks cython
                if ("H" == s0) and ("H" == o0):
                    if (s1 == o1) or (s1d and o1d):
                        enmS = self._endnum_re.findall(s)
                        enmO = self._endnum_re.findall(o)
                        if (enmS != []) and (enmO != []):
                            return int(enmS[0]), int(enmO[0])
                        elif enmS == []:
                            return 0, 1
                        else:
                            return 1, 0
                    elif s1d:
                        return 0, 1
                    elif o1d:
                        return 1, 0
                    else:
                        return (
                            self._greek_sort_keys[s1],
                            self._greek_sort_keys[o1],
                        )
                return int(s), int(o)  # raise exception?
        return 1, 1

    def __ne__(self, other: object) -> bool:
        """Test for inequality."""
        if isinstance(other, type(self)):
            return self.akl != other.akl
        else:
            return NotImplemented

    def __eq__(self, other: object) -> bool:  # type: ignore
        """Test for equality."""
        if isinstance(other, type(self)):
            return self.akl == other.akl
        else:
            return NotImplemented

    def __gt__(self, other: object) -> bool:
        """Test greater than."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] > rslt[1]
        else:
            return NotImplemented

    def __ge__(self, other: object) -> bool:
        """Test greater or equal."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] >= rslt[1]
        else:
            return NotImplemented

    def __lt__(self, other: object) -> bool:
        """Test less than."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] < rslt[1]
        else:
            return NotImplemented

    def __le__(self, other: object) -> bool:
        """Test less or equal."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] <= rslt[1]
        else:
            return NotImplemented


def set_accuracy_95(num: float) -> float:
    """Reduce floating point accuracy to 9.5 (xxxx.xxxxx).

    Used by :class:`IC_Residue` class writing PIC and SCAD
    files.

    :param float num: input number
    :returns: float with specified accuracy
    """
    # return round(num, 5)  # much slower
    return float(f"{num:9.5f}")


# internal coordinates construction Exceptions
class HedronMatchError(Exception):
    """Cannot find hedron in residue for given key."""


class MissingAtomError(Exception):
    """Missing atom coordinates for hedron or dihedron."""
