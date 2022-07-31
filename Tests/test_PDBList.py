# Copyright 2016 by Jacek Smietanski.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Testing access to the PDB over the internet."""

import os
import pathlib
import unittest
import unittest.mock

import requires_internet

from Bio.PDB.PDBList import PDBFileFormat
from Bio.PDB.PDBList import PDBList
import helpers

requires_internet.check()


class TestPBDListGetList(unittest.TestCase):
    """Test methods responsible for getting lists of entries."""

    def test_get_recent_changes(self):
        """Tests the Bio.PDB.PDBList.get_recent_changes method."""
        # obsolete_pdb declared to prevent from creating the "obsolete" directory
        pdblist = PDBList(obsolete_pdb="unimportant")
        url = pdblist.pdb_server + "/pub/pdb/data/status/latest/added.pdb"
        entries = pdblist.get_status_list(url)
        self.assertIsNotNone(entries)

    def test_get_all_entries(self):
        """Tests the Bio.PDB.PDBList.get_all_entries method."""
        # obsolete_pdb declared to prevent from creating the "obsolete" directory
        pdblist = PDBList(obsolete_pdb="unimportant")
        entries = pdblist.get_all_entries()
        # As number of entries constantly grow, test checks if a certain number was
        # exceeded
        self.assertGreater(len(entries), 100000)

    def test_get_all_obsolete(self):
        """Tests the Bio.PDB.PDBList.get_all_obsolete method."""
        # obsolete_pdb declared to prevent from creating the "obsolete" directory
        pdblist = PDBList(obsolete_pdb="unimportant")
        entries = pdblist.get_all_obsolete()
        # As number of obsolete entries constantly grow, test checks if a certain number
        # was exceeded
        self.assertGreater(len(entries), 3000)


class TestPDBListGetStructure(unittest.TestCase):
    """Test methods responsible for getting structures."""

    @helpers.temporary_directory()
    def check(
        self,
        structure,
        filename,
        file_format,
        obsolete=False,
        pdir=None,
        temporary_directory=None,
    ):
        pdb_list = PDBList(pdb=temporary_directory)
        pdir = str(pathlib.Path(pdb_list.local_pdb, pdir)) if pdir else None
        pdb_list.retrieve_pdb_file(
            structure, obsolete=obsolete, pdir=pdir, file_format=file_format
        )
        pdb_filepath = pathlib.Path(pdb_list.local_pdb, filename)
        self.assertTrue(pdb_filepath.is_file())
        pdb_filepath.unlink()

    def test_retrieve_pdb_file_small_pdb(self):
        """Tests retrieving the small molecule in pdb format."""
        structure = "127d"
        self.check(
            structure, os.path.join(structure[1:3], f"pdb{structure}.ent"), "pdb"
        )

    def test_retrieve_pdb_file_large_pdb(self):
        """Tests retrieving the bundle for large molecule in pdb-like format."""
        structure = "3k1q"
        self.check(
            structure,
            os.path.join(structure[1:3], f"{structure}-pdb-bundle.tar"),
            "bundle",
        )

    def test_retrieve_pdb_file_obsolete_pdb(self):
        """Tests retrieving the obsolete molecule in pdb format."""
        structure = "347d"
        self.check(
            structure,
            os.path.join("obsolete", structure[1:3], f"pdb{structure}.ent"),
            "pdb",
            obsolete=True,
        )

    def test_retrieve_pdb_file_obsolete_mmcif(self):
        """Tests retrieving the obsolete molecule in mmcif format."""
        structure = "347d"
        self.check(
            structure,
            os.path.join("obsolete", structure[1:3], f"{structure}.cif"),
            "mmCif",
            obsolete=True,
        )

    def test_retrieve_pdb_file_mmcif(self):
        """Tests retrieving the (non-obsolete) molecule in mmcif format."""
        structure = "127d"
        self.check(structure, os.path.join(structure[1:3], f"{structure}.cif"), "mmCif")

    def test_retrieve_pdb_file_obsolete_xml(self):
        """Tests retrieving the obsolete molecule in mmcif format."""
        structure = "347d"
        self.check(
            structure,
            os.path.join("obsolete", structure[1:3], f"{structure}.xml"),
            "xml",
            obsolete=True,
        )

    def test_retrieve_pdb_file_xml(self):
        """Tests retrieving the (non obsolete) molecule in xml format."""
        structure = "127d"
        self.check(structure, os.path.join(structure[1:3], f"{structure}.xml"), "xml")

    def test_retrieve_pdb_file_mmtf(self):
        """Tests retrieving the molecule in mmtf format."""
        structure = "127d"
        self.check(structure, os.path.join(structure[1:3], f"{structure}.mmtf"), "mmtf")

    def test_double_retrieve(self):
        """Tests retrieving the same file to different directories."""
        structure = "127d"
        self.check(structure, os.path.join("a", f"{structure}.cif"), "mmCif", pdir="a")
        self.check(structure, os.path.join("b", f"{structure}.cif"), "mmCif", pdir="b")


class TestPDBListUpdatePDB(unittest.TestCase):
    """Test PDBList.update_pdb."""

    @helpers.temporary_directory()
    def test_update_pdb_new_file(self, temporary_directory=None):
        """Test that the new files are retrieved locally."""
        pdb_list = PDBList(pdb=temporary_directory)
        with unittest.mock.patch(
            "Bio.PDB.PDBList.get_recent_changes", return_value=(["127d"], [], [])
        ):
            pdb_list.update_pdb(PDBFileFormat.PDB)
            self.assertTrue(pathlib.Path(pdb_list.local_pdb, "27/pdb127d.ent").exists())

    @helpers.temporary_directory()
    def test_update_pdb_modified_file(self, temporary_directory=None):
        """Test that the modified files are overridden locally."""
        pdb_list = PDBList(pdb=temporary_directory)
        directory = "27"
        filename = "127d.cif"
        outdated_pdb_file_content = "FAKE PDB FILE"

        pathlib.Path(pdb_list.local_pdb, directory).mkdir()
        with open(f"{pdb_list.local_pdb}/{directory}/{filename}", "w+") as file_stream:
            file_stream.write(outdated_pdb_file_content)

        with unittest.mock.patch(
            "Bio.PDB.PDBList.get_recent_changes",
            return_value=([], [filename.split(".")[0]], []),
        ):
            pdb_list.update_pdb(PDBFileFormat.MMCIF)

        with open(f"{pdb_list.local_pdb}/{directory}/{filename}", "r") as file_stream:
            self.assertNotEqual(outdated_pdb_file_content, file_stream.read())

    @helpers.temporary_directory()
    def test_update_pdb_obsolete_file(self, temporary_directory=None):
        """Test that the obsolete files are moved to the local obsolete directory."""
        outdated_pdb_file_content = "FAKE PDB FILE"
        directory = "47"
        filename = "347d.xml"
        pdb_list = PDBList(pdb=temporary_directory)

        pathlib.Path(pdb_list.local_pdb, directory).mkdir()
        with open(f"{pdb_list.local_pdb}/{directory}/{filename}", "w+") as file_stream:
            file_stream.write(outdated_pdb_file_content)

        with unittest.mock.patch(
            "Bio.PDB.PDBList.get_recent_changes",
            return_value=([], [], [filename.split(".")[0]]),
        ):
            pdb_list.update_pdb(PDBFileFormat.XML)

        self.assertFalse(pathlib.Path(pdb_list.local_pdb, directory, filename).exists())
        self.assertTrue(
            pathlib.Path(pdb_list.obsolete_pdb, directory, filename).exists()
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
