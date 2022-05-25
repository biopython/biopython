# Copyright 2017 by Valentin Vareskic (valentin.vareskic@gmail.com).
# All rights reserved. This code is part of the Biopython distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.

"""Tests for PDB PSEA."""

import unittest
from ftplib import FTP
from Bio.PDB.PSEA import run_psea


# ftp://ftp.lmcp.jussieu.fr/pub/sincris/software/protein/p-sea/
def setup():
    with FTP("ftp.lmcp.jussieu.fr") as ftp:
        ftp.login()
        ftp.dir()


# TODO Create unittest for run_psea(fname)
def test_run_psea():
    pass


if __name__ == "__main__":
    setup()
# TODO Create unittest psea(pname)
# TODO Create unittest psea2HEC(pseq)
# TODO Create unittest annotate(m, ss_seq)
# TODO Create unittest class PSEA
