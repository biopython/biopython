
import unittest
from Bio import RNA_Structure
from Bio.RNA_Structure.API import API_RNA_STRAND
from Bio.RNA_Structure.API import API_NDB
from Bio.RNA_Structure.RBP_score import RBP_score
from Bio.RNA_Structure import RNAfold
from Bio import MissingExternalDependencyError
import os
import sys

if sys.platform == "win32":
    # TODO
    raise MissingExternalDependencyError("Testing this on Windows is not implemented yet")
else:
    from Bio._py3k import getoutput
    output = getoutput("RNAfold --version")
    if "not found" not in output:
        rna_fold = "RNAfold"

if not rna_fold:
    raise MissingExternalDependencyError(
        "Install RNAfold if you want to use RNAfold from Biopython.")



class API_PROJEKTfunctions(unittest.TestCase):
    """ The internet connection is necessary  """
    def test_getsequence(self):
        example = API_NDB.via_sequence(pdb_id = "5SWE")
        assert example.get_seq_record() is not None

    def test_download_fasta_sequence(self):
        example = API_NDB.via_sequence(pdb_id = "5SWE")
        result = example.download_fasta_sequence()
        expected = "Fasta file is ready"
        self.assertEqual(result, expected)
        os.remove("5SWE_sequence.fasta")

    def test_metadata_to_file(self):
        example = API_NDB.via_sequence(pdb_id = "5SWE")
        result = example.metadata_to_file()
        expected = "File with metadata is ready"
        self.assertEqual(result, expected)
        os.remove("report_5SWE")

    def test_structure_download(self):
        example = API_NDB.via_sequence(pdb_id = "5SWE")
        result = example.download_pdb_structure()
        expected = "PDB file 5SWE is ready"
        self.assertEqual(result, expected)
        os.remove("5swe.pdb")

    def test_download_database(self):
        example = API_NDB.via_sequence(pdb_id = "5SWE")
        result = example.download_database()
        expected = "NDB Database was updated and converted to csv file"
        self.assertEqual(result, expected)

class ApiRNASTRANDfunctions(unittest.TestCase):

    def test_searchbysequence(self):
        example = API_RNA_STRAND.RNA_STRAND('UAAGCCCUA')
        assert example.search_by_sequence() is not None

    def test_download_database(self):
        example = API_RNA_STRAND.RNA_STRAND('UAAGCCCUA')
        result = example.download_database()
        expected = "RNA STRAND database downloaded to given localizaton."
        self.assertEqual(result, expected)
        os.remove("RNA_STRAND_data.tar.gz")

    def test_class_instance(self):
        example = API_RNA_STRAND.RNA_STRAND('UAAGCCCUA')
        self.assertIsInstance(example, API_RNA_STRAND.RNA_STRAND)

class Calculations(unittest.TestCase):

    def test_rmsd_calculation(self):
        example = RNA_Structure.rmsd_calculation("RNA_Structure/5swe_TEST.pdb", "RNA_Structure/5swe_TEST.pdb")
        expected = 'Normal RMSD: 0.0'
        self.assertEqual(example, expected)

    def test_rbp_calculation(self):
        example = RBP_score.struct_comparison('RNA_Structure/TMR_00273_structure_test','RNA_Structure/TMR_00200_structure_test', 1)
        result = example.rbp_score()
        expected = 1.0
        self.assertEqual(result, expected)
class RNA_Fold_wrapper(unittest.TestCase):

    def  test_rna_fold(self):
        cline = RNAfold.RNAfoldCommandLine(energyModel=0, dangles=2,maxBPspan=1)
        self.assertIsInstance(cline, RNAfold.RNAfoldCommandLine)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
