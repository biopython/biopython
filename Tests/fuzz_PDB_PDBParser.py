#!/usr/bin/python3
"""Fuzzes the PDB parser."""

import warnings
import atheris
import sys

with atheris.instrument_imports():
    from Bio import BiopythonWarning
    from Bio.PDB import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning


def TestOneInput(data):
    warnings.filterwarnings("ignore", category=PDBConstructionWarning)
    fdp = atheris.FuzzedDataProvider(data)
    number_of_lines = fdp.ConsumeIntInRange(1, 1000)
    max_line_length = 1000
    lines = [
        fdp.ConsumeString(fdp.ConsumeIntInRange(0, max_line_length))
        for _ in range(1, number_of_lines)
    ]
    parser = PDBParser()
    parser.header = None
    parser.trailer = None
    # Make a StructureBuilder instance (pass id of structure as parameter)
    parser.structure_builder.init_structure(fdp.ConsumeString(16))
    parser._parse(lines)
    parser.structure_builder.set_header(parser.header)
    # Return the Structure instance
    structure = parser.structure_builder.get_structure()
    atoms = structure.get_atoms()
    residues = structure.get_residues()


atheris.Setup(sys.argv, TestOneInput)
atheris.Fuzz()
