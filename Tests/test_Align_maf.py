# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.maf module."""
import unittest
from io import StringIO


from Bio.Align import maf


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.stockholm."
    ) from None


class TestAlign_reading(unittest.TestCase):
    def test_reading1(self):
        """Test parsing ucsc_mm9_chr10."""
        path = "MAF/ucsc_mm9_chr10.maf"
        alignments = maf.AlignmentIterator(path)
        alignment = next(alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
