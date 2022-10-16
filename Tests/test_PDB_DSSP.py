# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
# Revisions copyright 2020 Joao Rodrigues. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# Merged related test files into one, by Joao Rodrigues (2020)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB.DSSP submodule."""

import re
import subprocess
import unittest
import warnings

try:
    import numpy  # noqa: F401
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None


from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import DSSP, make_dssp_dict


VERSION_2_2_0 = (2, 2, 0)


def parse_dssp_version(version_string):
    """Parse the DSSP version into a tuple from the tool output."""
    match = re.search(r"\s*([\d.]+)", version_string)
    if match:
        version = match.group(1)
    return tuple(map(int, version.split(".")))


def will_it_float(s):  # well played, whoever this was :)
    """Convert the input into a float if it is a number.

    If the input is a string, the output does not change.
    """
    try:
        return float(s)
    except ValueError:
        return s


class DSSP_tool_test(unittest.TestCase):
    """Test calling DSSP from Bio.PDB."""

    @classmethod
    def setUpClass(cls):

        cls.dssp_version = (0, 0, 0)
        is_dssp_available = False
        # Check if DSSP is installed
        quiet_kwargs = {"stdout": subprocess.PIPE, "stderr": subprocess.STDOUT}
        try:
            try:
                # Newer versions of DSSP
                version_string = subprocess.check_output(
                    ["dssp", "--version"], text=True
                )
                cls.dssp_version = parse_dssp_version(version_string)
                is_dssp_available = True
            except subprocess.CalledProcessError:
                # Older versions of DSSP
                subprocess.check_call(["dssp", "-h"], **quiet_kwargs)
                is_dssp_available = True
        except OSError:
            try:
                version_string = subprocess.check_output(
                    ["mkdssp", "--version"], text=True
                )
                cls.dssp_version = parse_dssp_version(version_string)
                is_dssp_available = True
            except OSError:
                pass

        if not is_dssp_available:
            raise unittest.SkipTest(
                "Install dssp if you want to use it from Biopython."
            )

        cls.pdbparser = PDBParser()
        cls.cifparser = MMCIFParser()

    def test_dssp(self):
        """Test DSSP generation from PDB."""
        pdbfile = "PDB/2BEG.pdb"
        model = self.pdbparser.get_structure("2BEG", pdbfile)[0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # silence DSSP warnings
            dssp = DSSP(model, pdbfile)
        self.assertEqual(len(dssp), 130)

    # Only run mmCIF tests if DSSP version installed supports mmcif
    def test_dssp_with_mmcif_file(self):
        """Test DSSP generation from MMCIF."""
        if self.dssp_version < VERSION_2_2_0:
            self.skipTest("Test requires DSSP version 2.2.0 or greater")

        pdbfile = "PDB/4ZHL.cif"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # silence all warnings
            model = self.cifparser.get_structure("4ZHL", pdbfile)[0]
            dssp = DSSP(model, pdbfile)
        self.assertEqual(len(dssp), 257)

    def test_dssp_with_mmcif_file_and_nonstandard_residues(self):
        """Test DSSP generation from MMCIF with non-standard residues."""
        if self.dssp_version < VERSION_2_2_0:
            self.skipTest("Test requires DSSP version 2.2.0 or greater")

        pdbfile = "PDB/1AS5.cif"
        model = self.cifparser.get_structure("1AS5", pdbfile)[0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # silence DSSP warnings
            dssp = DSSP(model, pdbfile)
        self.assertEqual(len(dssp), 24)

    def test_dssp_with_mmcif_file_and_different_chain_ids(self):
        """Test DSSP generation from MMCIF which has different label and author chain IDs."""
        if self.dssp_version < VERSION_2_2_0:
            self.skipTest("Test requires DSSP version 2.2.0 or greater")

        pdbfile = "PDB/1A7G.cif"
        model = self.cifparser.get_structure("1A7G", pdbfile)[0]
        dssp = DSSP(model, pdbfile)
        self.assertEqual(len(dssp), 82)
        self.assertEqual(dssp.keys()[0][0], "E")


class DSSP_test(unittest.TestCase):
    """Tests for DSSP parsing etc which don't need the binary tool."""

    def test_DSSP_file(self):
        """Test parsing of pregenerated DSSP."""
        dssp, keys = make_dssp_dict("PDB/2BEG.dssp")
        self.assertEqual(len(dssp), 130)

    def test_DSSP_noheader_file(self):
        """Test parsing of pregenerated DSSP missing header information."""
        # New DSSP prints a line containing only whitespace and "."
        dssp, keys = make_dssp_dict("PDB/2BEG_noheader.dssp")
        self.assertEqual(len(dssp), 130)

    def test_DSSP_hbonds(self):
        """Test parsing of DSSP hydrogen bond information."""
        dssp, keys = make_dssp_dict("PDB/2BEG.dssp")

        dssp_indices = {v[5] for v in dssp.values()}
        hb_indices = set()

        # The integers preceding each hydrogen bond energy (kcal/mol) in the
        # "N-H-->O    O-->H-N    N-H-->O    O-->H-N" dssp output columns are
        # relative dssp indices. Therefore, "hb_indices" contains the absolute
        # dssp indices of residues participating in (provisional) h-bonds. Note
        # that actual h-bonds are typically determined by an energetic
        # threshold.
        for val in dssp.values():
            hb_indices |= {val[5] + x for x in (val[6], val[8], val[10], val[12])}

        # Check if all h-bond partner indices were successfully parsed.
        self.assertEqual((dssp_indices & hb_indices), hb_indices)

    def test_DSSP_in_model_obj(self):
        """All elements correctly added to xtra attribute of input model object."""
        p = PDBParser()
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        # Read the DSSP data into the pdb object:
        _ = DSSP(m, "PDB/2BEG.dssp", "dssp", "Sander", "DSSP")
        # Now compare the xtra attribute of the pdb object
        # residue by residue with the pre-computed values:
        i = 0
        with open("PDB/dssp_xtra_Sander.txt") as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    # Split the pre-computed values into a list:
                    xtra_list_ref = ref_lines[i].rstrip().split("\t")
                    # Then convert each element to float where possible:
                    xtra_list_ref = list(map(will_it_float, xtra_list_ref))
                    # The xtra attribute is a dict.
                    # To compare with the pre-computed values first sort according to keys:
                    xtra_itemts = sorted(
                        res.xtra.items(), key=lambda s: s[0]
                    )  # noqa: E731
                    # Then extract the list of xtra values for the residue
                    # and convert to floats where possible:
                    xtra_list = [t[1] for t in xtra_itemts]
                    xtra_list = list(map(will_it_float, xtra_list))
                    # The reason for converting to float is, that casting a float to a string in python2.6
                    # will include fewer decimals than python3 and an assertion error will be thrown.
                    self.assertEqual(xtra_list, xtra_list_ref)
                    i += 1

    def test_DSSP_RSA(self):
        """Tests the usage of different ASA tables."""
        # Tests include Sander/default, Wilke and Miller
        p = PDBParser()
        # Sander/default:
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        # Read the DSSP data into the pdb object:
        _ = DSSP(m, "PDB/2BEG.dssp", "dssp", "Sander", "DSSP")
        # Then compare the RASA values for each residue with the pre-computed values:
        i = 0
        with open("PDB/Sander_RASA.txt") as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    rasa_ref = float(ref_lines[i].rstrip())
                    rasa = float(res.xtra["EXP_DSSP_RASA"])
                    self.assertAlmostEqual(rasa, rasa_ref)
                    i += 1

        # Wilke (procedure similar as for the Sander values above):
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        _ = DSSP(m, "PDB/2BEG.dssp", "dssp", "Wilke", "DSSP")
        i = 0
        with open("PDB/Wilke_RASA.txt") as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    rasa_ref = float(ref_lines[i].rstrip())
                    rasa = float(res.xtra["EXP_DSSP_RASA"])
                    self.assertAlmostEqual(rasa, rasa_ref)
                    i += 1

        # Miller (procedure similar as for the Sander values above):
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        _ = DSSP(m, "PDB/2BEG.dssp", "dssp", "Miller", "DSSP")
        i = 0
        with open("PDB/Miller_RASA.txt") as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    rasa_ref = float(ref_lines[i].rstrip())
                    rasa = float(res.xtra["EXP_DSSP_RASA"])
                    self.assertAlmostEqual(rasa, rasa_ref)
                    i += 1

        # Ahmad (procedure similar as for the Sander values above):
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        _ = DSSP(m, "PDB/2BEG.dssp", "dssp", "Ahmad", "DSSP")
        i = 0
        with open("PDB/Ahmad_RASA.txt") as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    rasa_ref = float(ref_lines[i].rstrip())
                    rasa = float(res.xtra["EXP_DSSP_RASA"])
                    self.assertAlmostEqual(rasa, rasa_ref)
                    i += 1


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
