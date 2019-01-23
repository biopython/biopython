# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for chopping up (dicing) a structure.

This module is used internally by the Bio.PDB.extract() function.
It extracts a section of a chain from start to end and saves
the file in PDB format

Example
--------
Dice 2BEG.pdb, extract chain "B" from residues 18 to 20.
>>> from Bio.PDB import Dice
>>> from Bio.PDB import PDBParser
>>> parser = PDBParser()
>>> structure = parser.get_structure("scr", "PDB/2BEG.pdb")
>>> Dice.extract(structure, "B", 18, 20, "/tmp/2BEG_diced.pdb")

The file  2BEG_diced.pdb will contain a section of chain "B" from residue
18 to 20, three residues; VAL, PHE, PHE and their respective atoms,
total 29 atoms.
"""

import re
import warnings

from Bio.PDB.PDBIO import PDBIO
from Bio import BiopythonWarning

_hydrogen = re.compile("[123 ]*H.*")


class ChainSelector(object):
    """Only accepts residues with right chainid, between start and end.

    Remove hydrogens, waters and ligands. Only use model 0 by default.
    """

    def __init__(self, chain_id, start, end, model_id=0):
        """Initialize the class."""
        self.chain_id = chain_id
        self.start = start
        self.end = end
        self.model_id = model_id

    def accept_model(self, model):
        """Verify if model match the model identifier."""
        # model - only keep model 0
        if model.get_id() == self.model_id:
            return 1
        return 0

    def accept_chain(self, chain):
        """Verify if chain match chain identifier."""
        if chain.get_id() == self.chain_id:
            return 1
        return 0

    def accept_residue(self, residue):
        """Verify if a residue sequence is between the start and end sequence."""
        # residue - between start and end
        hetatm_flag, resseq, icode = residue.get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode != " ":
            warnings.warn("WARNING: Icode %s at position %s"
                          % (icode, resseq), BiopythonWarning)
        if self.start <= resseq <= self.end:
            return 1
        return 0

    def accept_atom(self, atom):
        """Verify if atoms are not Hydrogen."""
        # atoms - get rid of hydrogens
        name = atom.get_id()
        if _hydrogen.match(name):
            return 0
        else:
            return 1


def extract(structure, chain_id, start, end, filename):
    """Write out selected portion to filename."""
    sel = ChainSelector(chain_id, start, end)
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)
