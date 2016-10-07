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
    """
    Test methods responsible for getting lists of entries.
    """

    def test_get_recent_changes(self):
        """
        Tests the Bio.PDB.PDBList.get_recent_changes method.
        """
        pdblist = PDBList(obsolete_pdb="unimpotrant")  # obsolete_pdb declared to prevent from creating the "obsolete" directory
        url = pdblist.pdb_server + '/pub/pdb/data/status/latest/added.pdb'
        entries = pdblist.get_status_list(url)
        self.assertIsNotNone(entries)

    def test_get_all_entries(self):
        """
        Tests the Bio.PDB.PDBList.get_all_entries method.
        """
        pdblist = PDBList(obsolete_pdb="unimpotrant")  # obsolete_pdb declared to prevent from creating the "obsolete" directory
        entries = pdblist.get_all_entries()
        # As number of entries constantly grow, test checks if a certain number was exceeded
        self.assertTrue(len(entries) > 100000)

    def test_get_all_obsolete(self):
        """
        Tests the Bio.PDB.PDBList.get_all_obsolete method.
        """
        pdblist = PDBList(obsolete_pdb="unimpotrant")  # obsolete_pdb declared to prevent from creating the "obsolete" directory
        entries = pdblist.get_all_obsolete()
        # As number of obsolete entries constantly grow, test checks if a certain number was exceeded
        self.assertTrue(len(entries) > 3000)


class TestPDBListGetStructure(unittest.TestCase):
    """
    Test methods responsible for getting structures.
    """

    @contextlib.contextmanager
    def make_temp_directory(self, dir):
        temp_dir = tempfile.mkdtemp(dir=dir)
        try:
            yield temp_dir
        finally:
            shutil.rmtree(temp_dir)

    def test_retrieve_pdb_file_small_pdb(self):
        """
        Tests retrieving the small molecule in pdb format
        """
        structure = "127d"
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, structure[1:3], "pdb%s.ent" % structure)
            pdblist.retrieve_pdb_file(structure, file_format="pdb")
            self.assertTrue(os.path.isfile(path))
            os.remove(path)

    def test_retrieve_pdb_file_large_pdb(self):
        """
        Tests retrieving the bundle for large molecule in pdb-like format
        """
        structure = "3k1q"
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, structure[1:3], "%s-pdb-bundle.tar" % structure)
            pdblist.retrieve_pdb_file(structure, file_format="bundle")
            self.assertTrue(os.path.isfile(path))
            os.remove(path)

    def test_retrieve_pdb_file_obsolete_pdb(self):
        """
        Tests retrieving the obsolete molecule in pdb format
        """
        structure = "347d"
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, "obsolete", structure[1:3], "pdb%s.ent" % structure)
            pdblist.retrieve_pdb_file(structure, obsolete=True, file_format="pdb")
            self.assertTrue(os.path.isfile(path))
            os.remove(path)

    def test_retrieve_pdb_file_obsolete_mmcif(self):
        """
        Tests retrieving the obsolete molecule in mmcif format
        """
        structure = "347d"
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, "obsolete", structure[1:3], "%s.cif" % structure)
            pdblist.retrieve_pdb_file(structure, obsolete=True, file_format="mmCif")
            self.assertTrue(os.path.isfile(path))
            os.remove(path)

    def test_retrieve_pdb_file_mmcif(self):
        """
        Tests retrieving the (non-obsolete) molecule in mmcif format
        """
        structure = "127d"
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, structure[1:3], "%s.cif" % structure)
            pdblist.retrieve_pdb_file(structure, file_format="mmCif")
            self.assertTrue(os.path.isfile(path))
            os.remove(path)

    def test_retrieve_pdb_file_obsolete_xml(self):
        """
        Tests retrieving the obsolete molecule in mmcif format
        """
        structure = "347d"
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, "obsolete", structure[1:3], "%s.xml" % structure)
            pdblist.retrieve_pdb_file(structure, obsolete=True, file_format="xml")
            self.assertTrue(os.path.isfile(path))
            os.remove(path)

    def test_retrieve_pdb_file_xml(self):
        """
        Tests retrieving the (non obsolete) molecule in xml format
        """
        structure = "127d"
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, structure[1:3], "%s.xml" % structure)
            pdblist.retrieve_pdb_file(structure, file_format="xml")
            self.assertTrue(os.path.isfile(path))
            os.remove(path)

    def test_retrieve_pdb_file_mmtf(self):
        """
        Tests retrieving the molecule in mmtf format
        """
        structure = "127d"
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp, obsolete_pdb=os.path.join(tmp, "obsolete"))
            path = os.path.join(tmp, structure[1:3], "%s.mmtf" % structure)
            pdblist.retrieve_pdb_file(structure, file_format="mmtf")
            self.assertTrue(os.path.isfile(path))
            os.remove(path)
