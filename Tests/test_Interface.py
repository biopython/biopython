# Copyright 2011 by Mikael Trellet.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PDB module."""
import os
import unittest
import warnings
from StringIO import StringIO

from Bio.PDB import PDBParser
from Bio.PDB import InterfaceBuilder

class Interface(unittest.TestCase):
    " Test Interface Module "

    def setUp(self):
        "Setup"
    
        P = PDBParser(QUIET=1)
        structure = P.get_structure('test', 'PDB/2WFU.pdb')[0]
        self.structure = structure

    def test_InterfaceBuilder(self):
        """InterfaceBuilder Test"""

        interfaceBuilder = InterfaceBuilder.InterfaceBuilder(self.structure)
        interface = interfaceBuilder.get_interface()

        self.assertEqual(interface.level, 'I')
        self.assertEqual(interface.id, 'Interface_AB')
        self.assertEqual(len(interface), 9)

        interface_ids = [(' ', 19, ' '), (' ', 19, ' '), (' ', 20, ' '), (' ', 20, ' '), (' ', 21, ' '), (' ', 21, ' '), (' ', 22, ' '), (' ', 22, ' '), (' ', 23, ' ')]
        # Indirect way of testing because of formatting problems otherwise
        self.assertEqual(sorted([res.id for res in interface.get_list()]), interface_ids)

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
