#!/bin/python
"""Fuzzes the PDB.Superimposer api."""
import atheris
import sys
import numpy

from numpy.linalg import LinAlgError

with atheris.instrument_imports():
    from Bio.PDB import Superimposer, Selection
    from Bio.PDB import Atom
    from Bio.Data import IUPACData


def _fuzzed_atom(fdp):
    max_str = 1024
    atom_name = fdp.PickValueInList(list(IUPACData.atom_weights.keys()))
    return Atom.Atom(
        name=atom_name,
        coord=numpy.array(fdp.ConsumeFloatList(3)),
        bfactor=fdp.ConsumeProbability(),
        occupancy=fdp.ConsumeProbability(),
        altloc=fdp.ConsumeString(fdp.ConsumeIntInRange(0, max_str)),
        fullname=atom_name,
        element=atom_name.upper(),
        pqr_charge=fdp.ConsumeFloat(),
        radius=fdp.ConsumeFloat(),
        serial_number=fdp.ConsumeString(fdp.ConsumeIntInRange(0, max_str)),
    )


def _fuzzed_atom_list(fdp, len):
    return [_fuzzed_atom(fdp) for _ in range(0, len)]


def TestOneInput(data):
    fdp = atheris.FuzzedDataProvider(data)
    # It is possible that len(fixed) != len(moving), this is intentional with
    # the goal of catching vulnerabilities in the underlying library.
    max_list = 1024
    len = fdp.ConsumeIntInRange(0, max_list)
    fixed = _fuzzed_atom_list(fdp, len)
    moving = _fuzzed_atom_list(fdp, len)
    try:
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        sup.apply(moving)
    except ZeroDivisionError:
        # Probably just an invalid atom?
        pass
    except LinAlgError:
        pass


atheris.Setup(sys.argv, TestOneInput)
atheris.Fuzz()
