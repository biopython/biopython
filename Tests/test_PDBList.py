# Copyright 2016 by Jacek Smietanski.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Testing access to the PDB over the internet."""

import contextlib
import os
import shutil
import tempfile
import unittest

import requires_internet
requires_internet.check()

# We want to test this module:
from Bio.PDB.PDBList import PDBList


class TestPBDListGetList(unittest.TestCase):
    """Test methods responsible for getting lists of entries."""

    def test_get_recent_changes(self):
        """Tests the Bio.PDB.PDBList.get_recent_changes method."""
        pdblist = PDBList(obsolete_pdb="unimportant")  # obsolete_pdb declared to prevent from creating the "obsolete" directory
        url = pdblist.pdb_server + '/pub/pdb/data/status/latest/added.pdb'
        entries = pdblist.get_status_list(url)
        self.assertIsNotNone(entries)

    def test_get_all_entries(self):
        """Tests the Bio.PDB.PDBList.get_all_entries method."""
        pdblist = PDBList(obsolete_pdb="unimportant")  # obsolete_pdb declared to prevent from creating the "obsolete" directory
        entries = pdblist.get_all_entries()
        # As number of entries constantly grow, test checks if a certain number was exceeded
        self.assertTrue(len(entries) > 100000)

    def test_get_all_obsolete(self):
        """Tests the Bio.PDB.PDBList.get_all_obsolete method."""
        pdblist = PDBList(obsolete_pdb="unimportant")  # obsolete_pdb declared to prevent from creating the "obsolete" directory
        entries = pdblist.get_all_obsolete()
        # As number of obsolete entries constantly grow, test checks if a certain number was exceeded
        self.assertTrue(len(entries) > 3000)


class TestPDBListGetStructure(unittest.TestCase):
    """Test methods responsible for getting structures."""

    @contextlib.contextmanager
    def make_temp_directory(self, directory):
        temp_dir = tempfile.mkdtemp(dir=directory)
        try:
            yield temp_dir
        finally:
            shutil.rmtree(temp_dir)

    def check(self, structure, filename, file_format, obsolete=False, pdir=None):
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, filename)
            if pdir:
                pdir = os.path.join(tmp, pdir)
            pdblist.retrieve_pdb_file(structure, obsolete=obsolete, pdir=pdir, file_format=file_format)
            self.assertTrue(os.path.isfile(path))
            os.remove(path)

    def test_retrieve_pdb_file_small_pdb(self):
        """Tests retrieving the small molecule in pdb format."""
        structure = "127d"
        self.check(structure, os.path.join(structure[1:3], "pdb%s.ent" % structure), "pdb")

    def test_retrieve_pdb_file_large_pdb(self):
        """Tests retrieving the bundle for large molecule in pdb-like format."""
        structure = "3k1q"
        self.check(structure, os.path.join(structure[1:3], "%s-pdb-bundle.tar" % structure), "bundle")

    def test_retrieve_pdb_file_obsolete_pdb(self):
        """Tests retrieving the obsolete molecule in pdb format"""
        structure = "347d"
        self.check(structure, os.path.join("obsolete", structure[1:3], "pdb%s.ent" % structure), "pdb", obsolete=True)

    def test_retrieve_pdb_file_obsolete_mmcif(self):
        """Tests retrieving the obsolete molecule in mmcif format."""
        structure = "347d"
        self.check(structure, os.path.join("obsolete", structure[1:3], "%s.cif" % structure), "mmCif", obsolete=True)

    def test_retrieve_pdb_file_mmcif(self):
        """Tests retrieving the (non-obsolete) molecule in mmcif format."""
        structure = "127d"
        self.check(structure, os.path.join(structure[1:3], "%s.cif" % structure), "mmCif")

    def test_retrieve_pdb_file_obsolete_xml(self):
        """Tests retrieving the obsolete molecule in mmcif format."""
        structure = "347d"
        self.check(structure, os.path.join("obsolete", structure[1:3], "%s.xml" % structure), "xml", obsolete=True)

    def test_retrieve_pdb_file_xml(self):
        """Tests retrieving the (non obsolete) molecule in xml format."""
        structure = "127d"
        self.check(structure, os.path.join(structure[1:3], "%s.xml" % structure), "xml")

    def test_retrieve_pdb_file_mmtf(self):
        """Tests retrieving the molecule in mmtf format."""
        structure = "127d"
        self.check(structure, os.path.join(structure[1:3], "%s.mmtf" % structure), "mmtf")

    def test_double_retrieve(self):
        """Tests retrieving the same file to different directories."""
        structure = "127d"
        self.check(structure, os.path.join("a", "%s.cif" % structure), "mmCif", pdir="a")
        self.check(structure, os.path.join("b", "%s.cif" % structure), "mmCif", pdir="b")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
