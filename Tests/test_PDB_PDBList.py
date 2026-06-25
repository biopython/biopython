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
from unittest import mock

import requires_internet

# We want to test this module:
from Bio.PDB.PDBList import PDBList

requires_internet.check()


class TestPBDListGetList(unittest.TestCase):
    """Test methods responsible for getting lists of entries."""

    def test_default_server_includes_pdb_path_prefix(self):
        """The default server should point at the wwPDB data root."""
        pdblist = PDBList(obsolete_pdb="unimportant")
        self.assertEqual(pdblist.pdb_server, PDBList.DEFAULT_SERVER)

    def test_mirror_server_urls(self):
        """Custom servers should only need the mirror-specific data root."""
        ebi = PDBList(server="https://ftp.ebi.ac.uk/pub/databases/pdb")
        status_url = ebi.pdb_server + "/data/status/latest/added.pdb"
        self.assertEqual(
            status_url,
            "https://ftp.ebi.ac.uk/pub/databases/pdb/data/status/latest/added.pdb",
        )

        pdbj = PDBList(server="https://ftp.pdbj.org/pub/pdb")
        structure_url = (
            pdbj.pdb_server + "/data/structures/divided/pdb/27/pdb127d.ent.gz"
        )
        self.assertEqual(
            structure_url,
            "https://ftp.pdbj.org/pub/pdb/data/structures/divided/pdb/27/pdb127d.ent.gz",
        )

    def test_server_trailing_slash_is_stripped(self):
        pdblist = PDBList(
            server="https://ftp.ebi.ac.uk/pub/databases/pdb/", obsolete_pdb="unimportant"
        )
        self.assertEqual(pdblist.pdb_server, "https://ftp.ebi.ac.uk/pub/databases/pdb")

    def test_get_recent_changes(self):
        """Tests the Bio.PDB.PDBList.get_recent_changes method."""
        # obsolete_pdb declared to prevent from creating the "obsolete" directory
        pdblist = PDBList(obsolete_pdb="unimportant")
        url = pdblist.pdb_server + "/data/status/latest/added.pdb"
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

    def test_get_all_assemblies(self):
        """Tests the Bio.PDB.PDBList.get_all_assemblies method."""
        # obsolete_pdb declared to prevent from creating the "obsolete" directory
        pdblist = PDBList(obsolete_pdb="unimportant")
        entries = pdblist.get_all_assemblies()
        # As number of obsolete entries constantly grow, test checks if a certain number
        # was exceeded
        self.assertGreater(len(entries), 100000)


class TestPDBListURLConstruction(unittest.TestCase):
    """Check download URLs without hitting the network."""

    @staticmethod
    def _fake_urlretrieve(_url, filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, "wb") as handle:
            handle.write(b"")
        return filename, {}

    def test_retrieve_pdb_file_builds_rcsb_structure_url(self):
        pdblist = PDBList(pdb=".", obsolete_pdb="unimportant", verbose=False)
        with tempfile.TemporaryDirectory() as tmp:
            pdblist.local_pdb = tmp
            with mock.patch(
                "Bio.PDB.PDBList.urlretrieve", side_effect=self._fake_urlretrieve
            ) as urlretrieve:
                with mock.patch("Bio.PDB.PDBList.gzip.open") as gzip_open:
                    gzip_open.return_value.__enter__.return_value = iter([b""])
                    pdblist.retrieve_pdb_file("127d", file_format="pdb")
            url = urlretrieve.call_args[0][0]
        self.assertEqual(
            url,
            "https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/27/pdb127d.ent.gz",
        )

    def test_retrieve_pdb_file_builds_ebi_structure_url(self):
        pdblist = PDBList(
            server="https://ftp.ebi.ac.uk/pub/databases/pdb",
            pdb=".",
            obsolete_pdb="unimportant",
            verbose=False,
        )
        with tempfile.TemporaryDirectory() as tmp:
            pdblist.local_pdb = tmp
            with mock.patch(
                "Bio.PDB.PDBList.urlretrieve", side_effect=self._fake_urlretrieve
            ) as urlretrieve:
                with mock.patch("Bio.PDB.PDBList.gzip.open") as gzip_open:
                    gzip_open.return_value.__enter__.return_value = iter([b""])
                    pdblist.retrieve_pdb_file("127d", file_format="mmCif")
            url = urlretrieve.call_args[0][0]
        self.assertEqual(
            url,
            "https://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/mmCIF/27/127d.cif.gz",
        )

    def test_get_all_entries_uses_data_root(self):
        pdblist = PDBList(
            server="https://ftp.pdbj.org/pub/pdb", obsolete_pdb="unimportant", verbose=False
        )
        with mock.patch("Bio.PDB.PDBList.urlopen") as urlopen:
            urlopen.return_value.__enter__.return_value.readlines.return_value = [
                b"HEADER\n",
                b"----\n",
            ]
            pdblist.get_all_entries()
        url = urlopen.call_args[0][0]
        self.assertEqual(
            url, "https://ftp.pdbj.org/pub/pdb/derived_data/index/entries.idx"
        )

    def test_get_recent_changes_uses_data_root(self):
        pdblist = PDBList(
            server="https://ftp.ebi.ac.uk/pub/databases/pdb",
            obsolete_pdb="unimportant",
            verbose=False,
        )
        with mock.patch.object(PDBList, "get_status_list", return_value=["1abc"]) as get_list:
            pdblist.get_recent_changes()
        urls = [call.args[0] for call in get_list.call_args_list]
        self.assertEqual(
            urls[0],
            "https://ftp.ebi.ac.uk/pub/databases/pdb/data/status/latest/added.pdb",
        )
        self.assertEqual(
            urls[1],
            "https://ftp.ebi.ac.uk/pub/databases/pdb/data/status/latest/modified.pdb",
        )
        self.assertEqual(
            urls[2],
            "https://ftp.ebi.ac.uk/pub/databases/pdb/data/status/latest/obsolete.pdb",
        )


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
            pdblist = PDBList(pdb=tmp)
            path = os.path.join(tmp, filename)
            if pdir:
                pdir = os.path.join(tmp, pdir)
            pdblist.retrieve_pdb_file(
                structure, obsolete=obsolete, pdir=pdir, file_format=file_format
            )
            self.assertTrue(os.path.isfile(path))

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

    def test_double_retrieve_structure(self):
        """Tests retrieving the same file to different directories."""
        structure = "127d"
        self.check(structure, os.path.join("a", f"{structure}.cif"), "mmCif", pdir="a")
        self.check(structure, os.path.join("b", f"{structure}.cif"), "mmCif", pdir="b")

    def test_invalid_file_format(self):
        pdb_list = PDBList()
        with self.assertRaises(ValueError) as context:
            pdb_list.retrieve_pdb_file("127d", file_format="invalid")

        self.assertEqual(
            "Specified file_format invalid does not exist or is not supported. Please use one of the "
            "following: pdb, mmCif, xml, mmtf, bundle.",
            str(context.exception),
        )

    def test_retrieve_pdb_file_not_existing(self):
        """Tests retrieving a non-existent molecule - returns None and prints error message."""
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp)
            with unittest.patch("sys.stderr") as mock_stderr:
                result = pdblist.retrieve_pdb_file("zzzz", file_format="pdb")
                self.assertIsNone(result)
                self.assertTrue(
                    any(
                        "Desired structure not found or download failed." in str(call)
                        for call in mock_stderr.write.call_args_list
                    ),
                    "Error message not printed for non-existent molecule",
                )

    def test_retrieve_pdb_file_bad_url(self):
        """Tests retrieving with bad server URL - returns None and prints error message."""
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(server=" http://something.wrong ", pdb=tmp)
            with unittest.patch("sys.stderr") as mock_stderr:
                result = pdblist.retrieve_pdb_file("127d", file_format="pdb")
                self.assertIsNone(result)
                self.assertTrue(
                    any(
                        "Desired structure not found or download failed." in str(call)
                        for call in mock_stderr.write.call_args_list
                    ),
                    "Error message not printed for bad server URL",
                )


class TestPDBListGetAssembly(unittest.TestCase):
    """Test methods responsible for getting assemblies."""

    @contextlib.contextmanager
    def make_temp_directory(self, directory):
        temp_dir = tempfile.mkdtemp(dir=directory)
        try:
            yield temp_dir
        finally:
            shutil.rmtree(temp_dir)

    def check(self, structure, assembly_num, filename, file_format, pdir=None):
        with self.make_temp_directory(os.getcwd()) as tmp:
            pdblist = PDBList(pdb=tmp)
            path = os.path.join(tmp, filename)
            if pdir:
                pdir = os.path.join(tmp, pdir)
            pdblist.retrieve_assembly_file(
                structure, assembly_num, pdir=pdir, file_format=file_format
            )
            self.assertTrue(os.path.isfile(path))

    def test_retrieve_assembly_file_mmcif(self):
        """Tests retrieving a small assembly in mmCif format."""
        structure = "127d"
        assembly_num = "1"
        self.check(
            structure,
            assembly_num,
            os.path.join(structure[1:3], f"{structure}-assembly{assembly_num}.cif"),
            "mmCif",
        )

    def test_retrieve_assembly_file_pdb(self):
        """Tests retrieving a small assembly in pdb format."""
        structure = "127d"
        assembly_num = "1"
        self.check(
            structure,
            assembly_num,
            os.path.join(structure[1:3], f"{structure}.pdb{assembly_num}"),
            "pdb",
        )

    def test_double_retrieve_assembly(self):
        """Tests retrieving the same file to different directories."""
        structure = "127d"
        assembly_num = "1"
        self.check(
            structure,
            assembly_num,
            os.path.join("a", f"{structure}-assembly{assembly_num}.cif"),
            "mmCif",
            pdir="a",
        )
        self.check(
            structure,
            assembly_num,
            os.path.join("b", f"{structure}-assembly{assembly_num}.cif"),
            "mmCif",
            pdir="b",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
