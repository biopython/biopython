# Copyright 2012 by Eric Talevich.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import warnings

try:
    import numpy
    from numpy import dot  # Missing on PyPy's micronumpy
    del dot
    # We don't need this (?) but Bio.PDB imports it automatically :(
    from numpy.linalg import svd, det # Missing in PyPy 2.0 numpypy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use PDB formats with SeqIO.")

from Bio import SeqIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)


class TestPdbSeqres(unittest.TestCase):
    def test_seqres_parse(self):
        """Parse a multi-chain PDB by SEQRES entries.

        Reference:
        http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=2BEG
        """
        chains = list(SeqIO.parse('PDB/2BEG.pdb', 'pdb-seqres'))
        self.assertEqual(len(chains), 5)
        actual_seq = 'DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA'
        for chain, chn_id in zip(chains, 'ABCDE'):
            self.assertEqual(chain.id, '2BEG:' + chn_id)
            self.assertEqual(chain.annotations['chain'], chn_id)
            self.assertEqual(str(chain.seq), actual_seq)

    def test_seqres_read(self):
        """Read a single-chain PDB by SEQRES entries.

        Reference:
        http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=1A8O
        """
        chain = SeqIO.read('PDB/1A8O.pdb', 'pdb-seqres')
        self.assertEqual(chain.id, '1A8O:A')
        self.assertEqual(chain.annotations['chain'], 'A')
        self.assertEqual(str(chain.seq),
                         'MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTIL'
                         'KALGPGATLEEMMTACQG')

    def test_seqres_missing(self):
        """Parse a PDB with no SEQRES entries."""
        chains = list(SeqIO.parse('PDB/1MOT.pdb', 'pdb-seqres'))
        self.assertEqual(len(chains), 0)


class TestPdbAtom(unittest.TestCase):
    def test_atom_parse(self):
        """Parse a multi-chain PDB by ATOM entries.

        Reference:
        http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=2BEG
        """
        chains = list(SeqIO.parse('PDB/2BEG.pdb', 'pdb-atom'))
        self.assertEqual(len(chains), 5)
        actual_seq = 'LVFFAEDVGSNKGAIIGLMVGGVVIA'
        for chain, chn_id in zip(chains, 'ABCDE'):
            self.assertEqual(chain.id, '2BEG:' + chn_id)
            self.assertEqual(chain.annotations['chain'], chn_id)
            self.assertEqual(str(chain.seq), actual_seq)

        chains = list(SeqIO.parse('PDB/2XHE.pdb', 'pdb-atom'))
        actual_seq = 'DRLSRLRQMAAENQXXXXXXXXXXXXXXXXXXXXXXXPEPFMADFFNRVK'\
                     'RIRDNIEDIEQAIEQVAQLHTESLVAVSKEDRDRLNEKLQDTMARISALG'\
                     'NKIRADLKQIEKENKRAQQEGTFEDGTVSTDLRIRQSQHSSLSRKFVKVM'\
                     'TRYNDVQAENKRRYGENVARQCRVVEPSLSDDAIQKVIEHGXXXXXXXXX'\
                     'XXXXXXXXNEIRDRHKDIQQLERSLLELHEMFTDMSTLVASQGEMIDRIE'\
                     'FSVEQSHNYV'
        self.assertEqual(str(chains[1].seq), actual_seq)


    def test_atom_read(self):
        """Read a single-chain PDB by ATOM entries.

        Reference:
        http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=1A8O
        """
        chain = SeqIO.read('PDB/1A8O.pdb', 'pdb-atom')
        self.assertEqual(chain.id, '1A8O:A')
        self.assertEqual(chain.annotations['chain'], 'A')
        self.assertEqual(str(chain.seq),
                         'MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTIL'
                         'KALGPGATLEEMMTACQG')
        chain = SeqIO.read('PDB/a_structure.pdb', 'pdb-atom')
        self.assertEqual(chain.id, '????:A')
        self.assertEqual(chain.annotations['chain'], 'A')
        self.assertEqual(str(chain.seq), 'E')

    def test_atom_noheader(self):
        """Parse a PDB with no HEADER line."""
        warnings.simplefilter('ignore', UserWarning)
        chains = list(SeqIO.parse('PDB/1MOT.pdb', 'pdb-atom'))
        self.assertEqual(len(chains), 1)
        self.assertEqual(str(chains[0].seq), 'APARVGLGITTVLTMTTQSSGSRASLPK')


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
