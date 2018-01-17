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
from Bio.PDB import PDBParser, ResidueDepth
from Bio import MissingExternalDependencyError


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

    def test_ResidueDepth_2XHE(self):
        """Test on module that calculate residue depth via MSMS on protein structures."""
        prot_file = 'PDB/2XHE.pdb'
        p = PDBParser()
        s = p.get_structure("X", prot_file)
        model = s[0]
        rd = ResidueDepth(model, prot_file)
        res_chain = ''
        for item in rd.property_list[:100]:
            res_chain = res_chain + item[0].get_resname()
        self.assertEqual(res_chain, """HISMETSERLEULYSSERALAVALLYSTHRVALLEUTHRASNSERLEUARGSERVALALAASPGLYGLYASPTRPLYSVALLEUVALVALASPLYSPROALALEUARGMETILESERGLUCYSALAARGMETSERGLUILELEUASPLEUGLYVALTHRVALVALGLUASPVALSERLYSGLNARGLYSVALLEUPROGLNPHEHISGLYVALTYRPHEILEGLUPROTHRGLUGLUASNLEUASPTYRVALILEARGASPPHEALAASPARGTHRPROTHRTYRGLUALAALAHISLEU""")

    def test_ResidueDepth_2BEG(self):
        prot_file = 'PDB/2BEG.pdb'
        p = PDBParser()
        s = p.get_structure("X", prot_file)
        model = s[0]
        rd = ResidueDepth(model, prot_file)
        res_chain = ''
        for item in rd.property_list[:100]:
            res_chain = res_chain + item[0].get_resname()
        self.assertEqual(res_chain, """LEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLYVALVALILEALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLYVALVALILEALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLYVALVALILEALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLY""")

    def test_ResidueDepth_1LCD(self):
        prot_file = 'PDB/1LCD.pdb'
        p = PDBParser()
        s = p.get_structure("X", prot_file)
        model = s[0]
        rd = ResidueDepth(model, prot_file)
        res_chain = ''
        for item in rd.property_list[:100]:
            res_chain = res_chain + item[0].get_resname()
        self.assertEqual(res_chain, """METLYSPROVALTHRLEUTYRASPVALALAGLUTYRALAGLYVALSERTYRGLNTHRVALSERARGVALVALASNGLNALASERHISVALSERALALYSTHRARGGLULYSVALGLUALAALAMETALAGLULEUASNTYRILEPROASNARG""")

    def test_ResidueDepth_1A8O(self):
        prot_file = 'PDB/1A8O.pdb'
        p = PDBParser()
        s = p.get_structure("X", prot_file)
        model = s[0]
        rd = ResidueDepth(model, prot_file)
        res_chain = ''
        for item in rd.property_list[:100]:
            res_chain = res_chain + item[0].get_resname()
        self.assertEqual(res_chain, """MSEASPILEARGGLNGLYPROLYSGLUPROPHEARGASPTYRVALASPARGPHETYRLYSTHRLEUARGALAGLUGLNALASERGLNGLUVALLYSASNTRPMSETHRGLUTHRLEULEUVALGLNASNALAASNPROASPCYSLYSTHRILELEULYSALALEUGLYPROGLYALATHRLEUGLUGLUMSEMSETHRALACYSGLNGLY""")


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
