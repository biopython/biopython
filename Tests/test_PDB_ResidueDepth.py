# Copyright 2018 by Francesco Gastaldello. All rights reserved.
#
# Converted by Francesco Gastaldello from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# Merged related test files into one, by Joao Rodrigues (2020)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.ResidueDepth module."""

import subprocess
import unittest
import warnings

from Bio.PDB import MMCIFParser, PDBParser, ResidueDepth
from Bio import MissingExternalDependencyError
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.ResidueDepth import _get_atom_radius


class MSMS_tests(unittest.TestCase):
    """Test calling MSMS via Bio.PDB.ResidueDepth."""

    @classmethod
    def setUpClass(cls):
        # Check if MSMS is installed
        try:
            v = subprocess.check_output(
                ["msms", "-h"], text=True, stderr=subprocess.STDOUT
            )
        except OSError:
            raise unittest.SkipTest(
                "Install MSMS if you want to use it from Biopython."
            )

        cls.pdbparser = PDBParser()
        cls.cifparser = MMCIFParser()

    def check_msms(self, prot_file, first_100_residues):
        """Wrap calls to MSMS and the respective tests."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            s = self.pdbparser.get_structure("X", prot_file)

        model = s[0]
        rd = ResidueDepth(model)
        residues = []
        for item in rd.property_list[:100]:
            residues.append(item[0].get_resname())
        self.assertEqual("".join(residues), first_100_residues)

    # def test_ResidueDepth_2XHE(self):
    #     self.check_msms(
    #         "PDB/2XHE.pdb",
    #         "HISMETSERLEULYSSERALAVALLYSTHRVALLEUTHRASNSERLEUARGSERVALALAASPGLYGLYASPTR"
    #         "PLYSVALLEUVALVALASPLYSPROALALEUARGMETILESERGLUCYSALAARGMETSERGLUILELEUASPL"
    #         "EUGLYVALTHRVALVALGLUASPVALSERLYSGLNARGLYSVALLEUPROGLNPHEHISGLYVALTYRPHEILE"
    #         "GLUPROTHRGLUGLUASNLEUASPTYRVALILEARGASPPHEALAASPARGTHRPROTHRTYRGLUALAALAHI"
    #         "SLEU",
    #     )

    def test_ResidueDepth_2BEG(self):
        self.check_msms(
            "PDB/2BEG.pdb",
            "LEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLYVALVALIL"
            "EALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLYVALV"
            "ALILEALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLY"
            "VALVALILEALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGL"
            "YGLY",
        )

    def test_ResidueDepth_1LCD(self):
        self.check_msms(
            "PDB/1LCD.pdb",
            "METLYSPROVALTHRLEUTYRASPVALALAGLUTYRALAGLYVALSERTYRGLNTHRVALSERARGVALVALAS"
            "NGLNALASERHISVALSERALALYSTHRARGGLULYSVALGLUALAALAMETALAGLULEUASNTYRILEPROA"
            "SNARG",
        )

    def test_ResidueDepth_1A8O(self):
        self.check_msms(
            "PDB/1A8O.pdb",
            "MSEASPILEARGGLNGLYPROLYSGLUPROPHEARGASPTYRVALASPARGPHETYRLYSTHRLEUARGALAGL"
            "UGLNALASERGLNGLUVALLYSASNTRPMSETHRGLUTHRLEULEUVALGLNASNALAASNPROASPCYSLYST"
            "HRILELEULYSALALEUGLYPROGLYALATHRLEUGLUGLUMSEMSETHRALACYSGLNGLY",
        )


class ResidueDepth_tests(unittest.TestCase):
    """Tests for Bio.PDB.ResidueDepth, except for running MSMS itself."""

    def test_pdb_to_xyzr(self):
        """Test generation of xyzr (atomic radii) file."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            p = PDBParser(PERMISSIVE=1)
            structure = p.get_structure("example", "PDB/1A8O.pdb")

        # Read radii produced with original shell script
        with open("PDB/1A8O.xyzr") as handle:
            msms_radii = []
            for line in handle:
                fields = line.split()
                radius = float(fields[3])
                msms_radii.append(radius)

        model = structure[0]
        biopy_radii = []
        for atom in model.get_atoms():
            biopy_radii.append(_get_atom_radius(atom, rtype="united"))
        self.assertEqual(msms_radii, biopy_radii)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
