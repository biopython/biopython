# Copyright (C) 2006, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Wrappers for PSEA, a program for secondary structure assignment.

See this citation for P-SEA, PMID: 9183534

Labesse G, Colloc'h N, Pothier J, Mornon J-P:  P-SEA: a new efficient
assignment of secondary structure from C_alpha.
Comput Appl Biosci 1997 , 13:291-295

ftp://ftp.lmcp.jussieu.fr/pub/sincris/software/protein/p-sea/
"""

import tempfile  # <-- Add this at the very top of the file!
import shutil  # <-- Add this at the very top too!
import os
import subprocess

from Bio.PDB.Polypeptide import is_aa


from pathlib import Path


def run_psea(fname, verbose=False):
    """Run PSEA and return the output."""

    # We create a secure temporary directory
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Determine the expected output name
        # If input is "protein.pdb", PSEA will name output "protein.sea"
        input_file = Path(fname)
        output_name = input_file.stem + ".sea"

        # Run PSEA, but set the 'cwd' (current working directory) to our temp folder
        p = subprocess.run(["psea", fname], capture_output=True, text=True, cwd=tmp_dir)

        if verbose:
            print(p.stdout)

        temp_output = tmp_path / output_name

        if p.returncode == 0 and temp_output.exists():
            # IMPORTANT: Since the temp folder will be deleted,
            # we should probably read the data or move the file
            # to a permanent location before the 'with' block ends.
            return temp_output.read_text()
        else:
            raise RuntimeError(f"Error running p-sea: {p.stderr}")


def psea(pname):
    """Parse PSEA output file."""
    fname = run_psea(pname)
    start = 0
    ss = ""
    with open(fname) as fp:
        for line in fp:
            if line[0:6] == ">p-sea":
                start = 1
                continue
            if not start:
                continue
            if line[0] == "\n":
                break
            ss = ss + line[0:-1]
    return ss


def psea2HEC(pseq):
    """Translate PSEA secondary structure string into HEC."""
    seq = []
    for ss in pseq:
        if ss == "a":
            n = "H"
        elif ss == "b":
            n = "E"
        elif ss == "c":
            n = "C"
        seq.append(n)
    return seq


def annotate(m, ss_seq):
    """Apply secondary structure information to residues in model."""
    c = m.get_list()[0]
    all = c.get_list()
    residues = []
    # Now remove HOH etc.
    for res in all:
        if is_aa(res):
            residues.append(res)
    L = len(residues)
    if not L == len(ss_seq):
        raise ValueError("Length mismatch %i %i" % (L, len(ss_seq)))
    for i in range(L):
        residues[i].xtra["SS_PSEA"] = ss_seq[i]
    # subprocess.call(["rm", fname])


class PSEA:
    """Define PSEA class.

    PSEA object is a wrapper to PSEA program for secondary structure assignment.
    """

    def __init__(self, model, filename):
        """Initialize the class."""
        ss_seq = psea(filename)
        ss_seq = psea2HEC(ss_seq)
        annotate(model, ss_seq)
        self.ss_seq = ss_seq

    def get_seq(self):
        """Return secondary structure string."""
        return self.ss_seq
