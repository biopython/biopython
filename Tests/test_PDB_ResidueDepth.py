# Copyright 2018 by Francesco Gastaldello. All rights reserved.
#
# Converted by Francesco Gastaldello from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.ResidueDepth module."""

import unittest
import warnings
from Bio.PDB import PDBParser, ResidueDepth
from Bio import MissingExternalDependencyError
from Bio.PDB.PDBExceptions import PDBConstructionWarning

msms_exe = None
from Bio._py3k import getoutput
try:
    output = getoutput("msms -h")
    if output.startswith("Usage : msms parameters"):
        msms_exe = "msms"
except OSError:
    pass

if not msms_exe:
    raise MissingExternalDependencyError(
        "Install MSMS if you want to use it in Biopython.")


class ResidueDepthTests(unittest.TestCase):
    """Test ResidueDepth module."""

    def check_msms(self, prot_file, first_100_residues):
        p = PDBParser()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            s = p.get_structure("X", prot_file)
        model = s[0]
        rd = ResidueDepth(model)
        res_chain = ''
        for item in rd.property_list[:100]:
            res_chain = res_chain + item[0].get_resname()
        self.assertEqual(res_chain, first_100_residues)

    def test_ResidueDepth_2XHE(self):
        self.check_msms('PDB/2XHE.pdb', 'HISMETSERLEULYSSERALAVALLYSTHRVALLEUTH'
                                        'RASNSERLEUARGSERVALALAASPGLYGLYASPTRPL'
                                        'YSVALLEUVALVALASPLYSPROALALEUARGMETILE'
                                        'SERGLUCYSALAARGMETSERGLUILELEUASPLEUGL'
                                        'YVALTHRVALVALGLUASPVALSERLYSGLNARGLYSV'
                                        'ALLEUPROGLNPHEHISGLYVALTYRPHEILEGLUPRO'
                                        'THRGLUGLUASNLEUASPTYRVALILEARGASPPHEAL'
                                        'AASPARGTHRPROTHRTYRGLUALAALAHISLEU')

    def test_ResidueDepth_2BEG(self):
        self.check_msms('PDB/2BEG.pdb', 'LEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGL'
                                        'YALAILEILEGLYLEUMETVALGLYGLYVALVALILEA'
                                        'LALEUVALPHEPHEALAGLUASPVALGLYSERASNLYS'
                                        'GLYALAILEILEGLYLEUMETVALGLYGLYVALVALIL'
                                        'EALALEUVALPHEPHEALAGLUASPVALGLYSERASNL'
                                        'YSGLYALAILEILEGLYLEUMETVALGLYGLYVALVAL'
                                        'ILEALALEUVALPHEPHEALAGLUASPVALGLYSERAS'
                                        'NLYSGLYALAILEILEGLYLEUMETVALGLYGLY')

    def test_ResidueDepth_1LCD(self):
        self.check_msms('PDB/1LCD.pdb', 'METLYSPROVALTHRLEUTYRASPVALALAGLUTYRAL'
                                        'AGLYVALSERTYRGLNTHRVALSERARGVALVALASNG'
                                        'LNALASERHISVALSERALALYSTHRARGGLULYSVAL'
                                        'GLUALAALAMETALAGLULEUASNTYRILEPROASNAR'
                                        'G')

    def test_ResidueDepth_1A8O(self):
        self.check_msms('PDB/1A8O.pdb', 'MSEASPILEARGGLNGLYPROLYSGLUPROPHEARGAS'
                                        'PTYRVALASPARGPHETYRLYSTHRLEUARGALAGLUG'
                                        'LNALASERGLNGLUVALLYSASNTRPMSETHRGLUTHR'
                                        'LEULEUVALGLNASNALAASNPROASPCYSLYSTHRIL'
                                        'ELEULYSALALEUGLYPROGLYALATHRLEUGLUGLUM'
                                        'SEMSETHRALACYSGLNGLY')


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
