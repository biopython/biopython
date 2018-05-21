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
    from numpy.linalg import svd, det  # Missing in PyPy 2.0 numpypy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use PDB formats with SeqIO.")

from Bio import SeqIO
from Bio import BiopythonParserWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning


def SeqresTestGenerator(extension, parser):
    """Test factory for tests reading SEQRES (or similar) records.

    This is a factory returning a parameterised superclass for tests reading
    sequences from the sequence records of structure files.

    Arguments:
        extension:
            The extension of the files to read from the ``PDB`` directory (e.g.
            ``pdb`` or ``cif``).
        parser:
            The name of the SeqIO parser to use (e.g. ``pdb-atom``).
    """

    class SeqresTests(unittest.TestCase):
        """Use "parser" to parse sequence records from a structure file.

        Args:
            parser (str): Name of the parser used by SeqIO.
            extension (str): Extension of the files to parse.
        """

        def test_seqres_parse(self):
            """Parse a multi-chain PDB by SEQRES entries.

            Reference:
            http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=2BEG
            """
            chains = list(SeqIO.parse('PDB/2BEG.' + extension, parser))
            self.assertEqual(len(chains), 5)
            actual_seq = 'DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA'
            for chain, chn_id in zip(chains, 'ABCDE'):
                self.assertEqual(chain.id, '2BEG:' + chn_id)
                self.assertEqual(chain.annotations['chain'], chn_id)
                self.assertEqual(str(chain.seq), actual_seq)

        def test_seqres_read(self):
            """Read a single-chain structure by sequence entries.

            Reference:
            http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=1A8O
            """
            chain = SeqIO.read('PDB/1A8O.' + extension, parser)
            self.assertEqual(chain.id, '1A8O:A')
            self.assertEqual(chain.annotations['chain'], 'A')
            self.assertEqual(str(chain.seq),
                             'MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPD'
                             'CKTILKALGPGATLEEMMTACQG')

        def test_seqres_missing(self):
            """Parse a PDB with no SEQRES entries."""
            chains = list(SeqIO.parse('PDB/a_structure.' + extension, parser))
            self.assertEqual(len(chains), 0)
    return SeqresTests


class TestPdbSeqres(SeqresTestGenerator("pdb", "pdb-seqres")):
    """Test pdb-seqres SeqIO driver."""
    pass


class TestCifSeqres(SeqresTestGenerator("cif", "cif-seqres")):
    """Test cif-seqres SeqIO driver."""
    pass


def AtomTestGenerator(extension, parser):
    """Test factory for tests reading ATOM (or similar) records.

    See SeqresTestGenerator for more information.
    """

    class AtomTests(unittest.TestCase):
        def test_atom_parse(self):
            """Parse a multi-chain structure by ATOM entries.

            Reference:
            http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=2BEG
            """
            chains = list(SeqIO.parse('PDB/2BEG.' + extension, parser))
            self.assertEqual(len(chains), 5)
            actual_seq = 'LVFFAEDVGSNKGAIIGLMVGGVVIA'
            for chain, chn_id in zip(chains, 'ABCDE'):
                self.assertEqual(chain.id, '2BEG:' + chn_id)
                self.assertEqual(chain.annotations['chain'], chn_id)
                self.assertEqual(chain.annotations['model'], 0)
                self.assertEqual(str(chain.seq), actual_seq)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", PDBConstructionWarning)
                chains = list(SeqIO.parse('PDB/2XHE.' + extension, parser))
            actual_seq = 'DRLSRLRQMAAENQXXXXXXXXXXXXXXXXXXXXXXXPEPFMADFFNRVK'\
                         'RIRDNIEDIEQAIEQVAQLHTESLVAVSKEDRDRLNEKLQDTMARISALG'\
                         'NKIRADLKQIEKENKRAQQEGTFEDGTVSTDLRIRQSQHSSLSRKFVKVM'\
                         'TRYNDVQAENKRRYGENVARQCRVVEPSLSDDAIQKVIEHGXXXXXXXXX'\
                         'XXXXXXXXNEIRDRHKDIQQLERSLLELHEMFTDMSTLVASQGEMIDRIE'\
                         'FSVEQSHNYV'
            self.assertEqual(str(chains[1].seq), actual_seq)

        def test_atom_read(self):
            """Read a single-chain structure by ATOM entries.

            Reference:
            http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=1A8O
            """
            chain = SeqIO.read('PDB/1A8O.' + extension, parser)
            self.assertEqual(chain.id, '1A8O:A')
            self.assertEqual(chain.annotations['chain'], 'A')
            self.assertEqual(chain.annotations['model'], 0)
            self.assertEqual(str(chain.seq),
                             'MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTIL'
                             'KALGPGATLEEMMTACQG')

    return AtomTests


class TestPdbAtom(AtomTestGenerator('pdb', 'pdb-atom')):
    """Test pdb-atom SeqIO driver."""

    def test_atom_noheader(self):
        """Parse a PDB with no HEADER line."""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            warnings.simplefilter("ignore", BiopythonParserWarning)
            chains = list(SeqIO.parse('PDB/1LCD.pdb', 'pdb-atom'))

        self.assertEqual(len(chains), 1)
        self.assertEqual(str(chains[0].seq), 'MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNR')

    def test_atom_read_noheader(self):
        """Read a single-chain PDB without a header by ATOM entries."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            warnings.simplefilter("ignore", BiopythonParserWarning)
            chain = SeqIO.read('PDB/a_structure.pdb', 'pdb-atom')
        self.assertEqual(chain.id, '????:A')
        self.assertEqual(chain.annotations['chain'], 'A')
        self.assertEqual(str(chain.seq), 'E')


class TestCifAtom(AtomTestGenerator('cif', 'cif-atom')):
    """Test cif-atom SeqIO driver."""

    def test_atom_read_noheader(self):
        """Read a single-chain CIF without a header by ATOM entries."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            warnings.simplefilter("ignore", BiopythonParserWarning)
            chain = SeqIO.read('PDB/a_structure.cif', 'cif-atom')
        self.assertEqual(chain.id, '????:A')
        self.assertEqual(chain.annotations['chain'], 'A')
        self.assertEqual(str(chain.seq),
                         'MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTIL'
                         'KALGPGATLEEMMTACQG')


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
