# Copyright 2016 by Jacek Smietanski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Testing access to the PDB over the internet."""

import unittest
import os
import shutil
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
        pdblist = PDBList()
        url = pdblist.pdb_server + '/pub/pdb/data/status/latest/added.pdb'
        entries = pdblist.get_status_list(url)
        self.assertIsNotNone(entries)

    def test_get_all_entries(self):
        """
        Tests the Bio.PDB.PDBList.get_all_entries method.
        """
        pdblist = PDBList()
        entries = pdblist.get_all_entries()
        # As number of entries constantly grow, test checks if a certain number was exceeded
        self.assertTrue(len(entries) > 100000)

    def test_get_all_obsolete(self):
        """
        Tests the Bio.PDB.PDBList.get_all_obsolete method.
        """
        pdblist = PDBList()
        entries = pdblist.get_all_obsolete()
        # As number of obsolete entries constantly grow, test checks if a certain number was exceeded
        self.assertTrue(len(entries) > 3000)


class TestPDBListGetStructure(unittest.TestCase):
    """
    Test methods responsible for getting structures.
    """

    def test_retrieve_pdb_file(self):
        """
        Tests the process of retrieving files from PDB database.
        """
        # This test downloads several large files from external database so it may run a while.

        # Create temporary directory for files downloaded while testing.
        basedir = os.path.join(os.getcwd(), "pdb-test")
        try:
            os.mkdir(basedir)
        except FileExistsError:
            pass

        pdblist = PDBList(pdb = basedir, obsolete_pdb=os.path.join(basedir, "obsolete"))
        small = "1esy"
        obslte = "347d"
        large = "3k1q"
        nonexistent = "0000"
        formats = {"pdb": "pdb%s.ent", "mmCif": "%s.cif", "xml": "%s.xml", "mmtf": "%s.mmtf",
                   "bundle": "%s-pdb-bundle.tar"}
        for file_format in formats.keys():
            for structure in (small, large, obslte, nonexistent):
                for overwrite in (False, True):
                    for pdir in (None, os.path.join(basedir, "test-pdir")):
                        print(file_format, structure, "overwrite:", overwrite, "pdir:", pdir)
                        if pdir is None:
                            if structure == obslte and format != 'mmtf':
                                obsolete = True
                                path = os.path.join(basedir, "obsolete", structure[1:3], formats[file_format] % structure)
                            else:
                                obsolete = False
                                path = os.path.join(basedir, structure[1:3], formats[file_format] % structure)
                        else:  # pdir specified
                            path = os.path.join(pdir, formats[file_format] % structure)
                        pdblist.retrieve_pdb_file(structure, file_format, overwrite, obsolete, pdir)
                        exists = os.path.isfile(path)
                        error_msg = "error with " + structure + " " + file_format + " overwrite=" + str(overwrite) + \
                                    " obsolete=" + str(obsolete) + " pdir:" + str(pdir)
                        if structure == nonexistent or \
                           (structure == large and file_format == "pdb") or \
                           (file_format == "bundle" and structure in (small, obslte)) or \
                           (file_format) == "mmtf" and structure==obslte:
                            self.assertFalse(exists, msg=error_msg)
                        else:
                            self.assertTrue(exists, msg=error_msg)
        # Delete temporary directory
        shutil.rmtree(basedir)
