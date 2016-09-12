import unittest
import warnings
from Bio.PDB.mmtf import MMTFParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

class ParseMicrohet(unittest.TestCase):
    """Just parse some Microheterogenous files"""

    def test_cif(self):
        """Parse MMTF file """
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            structure = MMTFParser.get_structure("PDB/1EJG.mmtf")

    def test_mmtf(self):
        """Parse mmCIF file """
        with warnings.catch_warnings():
	    mmcif_parser = MMCIFParser()
            warnings.simplefilter('ignore', PDBConstructionWarning)
            structure = mmcif_parser.get_structure("MICR","PDB/1EJG.cif")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

