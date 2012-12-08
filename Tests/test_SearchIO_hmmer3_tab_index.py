# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO hmmer3-tab indexing."""

import os
import unittest

from search_tests_common import CheckRaw, CheckIndex


class Hmmer3TabRawCases(CheckRaw):

    fmt = 'hmmer3-tab'

    def test_hmmer3tab_30_multiple_first(self):
        """Test hmmer3-tab raw string retrieval, HMMER 3.0, multiple queries, first (tab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_001.out')
        raw = """Globin               PF00042.17 gi|4885477|ref|NP_005359.1| -              6e-21   74.6   0.3   9.2e-21   74.0   0.2   1.3   1   0   0   1   1   1   1 Globin
"""
        self.check_raw(filename, "gi|4885477|ref|NP_005359.1|", raw)

    def test_hmmer3tab_30_multiple_middle(self):
        """Test hmmer3-tab raw string retrieval, HMMER 3.0, multiple queries, middle (tab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_001.out')
        raw = """Ig_3                 PF13927.1  gi|126362951:116-221 -            1.4e-09   38.2   0.4   2.1e-09   37.6   0.3   1.3   1   0   0   1   1   1   1 Immunoglobulin domain
Ig_2                 PF13895.1  gi|126362951:116-221 -            3.5e-05   23.7   0.1   4.3e-05   23.4   0.1   1.1   1   0   0   1   1   1   1 Immunoglobulin domain
"""
        self.check_raw(filename, "gi|126362951:116-221", raw)

    def test_hmmer3tab_30_multiple_last(self):
        """Test hmmer3-tab raw string retrieval, HMMER 3.0, multiple queries, last (tab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_001.out')
        raw = """Pou                  PF00157.12 gi|125490392|ref|NP_038661.2| -              7e-37  124.8   0.5   1.4e-36  123.9   0.3   1.5   1   0   0   1   1   1   1 Pou domain - N-terminal to homeobox domain
Homeobox             PF00046.24 gi|125490392|ref|NP_038661.2| -            2.1e-18   65.5   1.1   4.1e-18   64.6   0.7   1.5   1   0   0   1   1   1   1 Homeobox domain
HTH_31               PF13560.1  gi|125490392|ref|NP_038661.2| -              0.012   15.6   0.0      0.16   12.0   0.0   2.2   2   0   0   2   2   2   0 Helix-turn-helix domain
Homeobox_KN          PF05920.6  gi|125490392|ref|NP_038661.2| -              0.039   13.5   0.0     0.095   12.3   0.0   1.6   1   0   0   1   1   1   0 Homeobox KN domain
DUF521               PF04412.8  gi|125490392|ref|NP_038661.2| -               0.14   10.5   0.1      0.26    9.6   0.1   1.4   1   0   0   1   1   1   0 Protein of unknown function (DUF521)
"""
        self.check_raw(filename, "gi|125490392|ref|NP_038661.2|", raw)

    def test_hmmer3tab_30_single(self):
        """Test hmmer3-tab raw string retrieval, HMMER 3.0, single query (tab_30_hmmscan_004.out)"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_004.out')
        raw = """Ig_3                 PF13927.1  gi|126362951:116-221 -            1.4e-09   38.2   0.4   2.1e-09   37.6   0.3   1.3   1   0   0   1   1   1   1 Immunoglobulin domain
Ig_2                 PF13895.1  gi|126362951:116-221 -            3.5e-05   23.7   0.1   4.3e-05   23.4   0.1   1.1   1   0   0   1   1   1   1 Immunoglobulin domain
"""
        self.check_raw(filename, "gi|126362951:116-221", raw)


class Hmmer3TabIndexCases(CheckIndex):

    fmt = 'hmmer3-tab'

    def test_hmmer3tab_30_hmmscan_001(self):
        """Test hmmer3-tab indexing, HMMER 3.0, multiple queries"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_001.out')
        self.check_index(filename, self.fmt)

    def test_hmmer3tab_30_hmmscan_002(self):
        """Test hmmer3-tab indexing, HMMER 3.0, single query, no hits"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_002.out')
        self.check_index(filename, self.fmt)

    def test_hmmer3tab_30_hmmscan_003(self):
        """Test hmmer3-tab indexing, HMMER 3.0, single query, multiple hits"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_003.out')
        self.check_index(filename, self.fmt)

    def test_hmmer3tab_30_hmmscan_004(self):
        """Test hmmer3-tab indexing, HMMER 3.0, single query, no alignments"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_004.out')
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
