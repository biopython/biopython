# Copyright 2012 by Eric Talevich.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

from Bio import SeqIO


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

        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

