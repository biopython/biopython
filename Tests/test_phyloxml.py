#!/usr/bin/env python

"""Unit tests for the Bio.PhyloXML module."""

import unittest
import warnings

from Bio import PhyloXML


class ParseTests(unittest.TestCase):
    pass


# -------------------------------------------------------------

if __name__ == '__main__':
    # Hide warnings from the user
    warnings.simplefilter('ignore')
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

