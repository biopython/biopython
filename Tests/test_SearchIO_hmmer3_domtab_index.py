# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO hmmer3-domtab indexing."""

import os
import unittest

from search_tests_common import CheckRaw, CheckIndex


class HmmerDomtabRawCases(CheckRaw):

    fmt = 'hmmscan3-domtab'

    def test_hmmerdomtab_30_multiple_first(self):
        """Test hmmscan-domtab raw string retrieval, HMMER 3.0, multiple queries, first (domtab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_001.out')
        raw = """Globin               PF00042.17   108 gi|4885477|ref|NP_005359.1| -            154     6e-21   74.6   0.3   1   1   6.7e-25   9.2e-21   74.0   0.2     1   107     7   112     7   113 0.97 Globin
"""
        self.check_raw(filename, "gi|4885477|ref|NP_005359.1|", raw)

    def test_hmmerdomtab_30_multiple_middle(self):
        """Test hmmscan-domtab raw string retrieval, HMMER 3.0, multiple queries, middle (domtab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_001.out')
        raw = """Ig_3                 PF13927.1     75 gi|126362951:116-221 -            106   1.4e-09   38.2   0.4   1   1     3e-13   2.1e-09   37.6   0.3     1    73     9    84     9    88 0.94 Immunoglobulin domain
Ig_2                 PF13895.1     80 gi|126362951:116-221 -            106   3.5e-05   23.7   0.1   1   1   6.2e-09   4.3e-05   23.4   0.1     1    80     9   104     9   104 0.71 Immunoglobulin domain
"""
        self.check_raw(filename, "gi|126362951:116-221", raw)

    def test_hmmerdomtab_30_multiple_last(self):
        """Test hmmscan-domtab raw string retrieval, HMMER 3.0, multiple queries, last (domtab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_001.out')
        raw = """Pou                  PF00157.12    75 gi|125490392|ref|NP_038661.2| -            352     7e-37  124.8   0.5   1   1     5e-40   1.4e-36  123.9   0.3     3    75   133   205   131   205 0.97 Pou domain - N-terminal to homeobox domain
Homeobox             PF00046.24    57 gi|125490392|ref|NP_038661.2| -            352   2.1e-18   65.5   1.1   1   1   1.5e-21   4.1e-18   64.6   0.7     1    57   224   280   224   280 0.98 Homeobox domain
HTH_31               PF13560.1     64 gi|125490392|ref|NP_038661.2| -            352     0.012   15.6   0.0   1   2   5.7e-05      0.16   12.0   0.0     1    35   141   181   141   184 0.96 Helix-turn-helix domain
HTH_31               PF13560.1     64 gi|125490392|ref|NP_038661.2| -            352     0.012   15.6   0.0   2   2      0.19   5.2e+02    0.8   0.0    39    62   245   268   243   270 0.86 Helix-turn-helix domain
Homeobox_KN          PF05920.6     40 gi|125490392|ref|NP_038661.2| -            352     0.039   13.5   0.0   1   1   3.5e-05     0.095   12.3   0.0     7    39   244   276   241   277 0.91 Homeobox KN domain
DUF521               PF04412.8    400 gi|125490392|ref|NP_038661.2| -            352      0.14   10.5   0.1   1   1   9.4e-05      0.26    9.6   0.1   273   334   221   280   197   294 0.77 Protein of unknown function (DUF521)
"""
        self.check_raw(filename, "gi|125490392|ref|NP_038661.2|", raw)

    def test_hmmerdomtab_30_single(self):
        """Test hmmscan-domtab raw string retrieval, HMMER 3.0, single query (domtab_30_hmmscan_004.out)"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_004.out')
        raw = """Ig_3                 PF13927.1     75 gi|126362951:116-221 -            106   1.4e-09   38.2   0.4   1   1     3e-13   2.1e-09   37.6   0.3     1    73     9    84     9    88 0.94 Immunoglobulin domain
Ig_2                 PF13895.1     80 gi|126362951:116-221 -            106   3.5e-05   23.7   0.1   1   1   6.2e-09   4.3e-05   23.4   0.1     1    80     9   104     9   104 0.71 Immunoglobulin domain
"""
        self.check_raw(filename, "gi|126362951:116-221", raw)


class HmmerDomtabIndexCases(CheckIndex):

    def test_hmmerdomtab_30_hmmscan_001(self):
        """Test hmmscan-domtab indexing, HMMER 3.0, multiple queries"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_001.out')
        self.check_index(filename, 'hmmscan3-domtab')

    def test_hmmerdomtab_30_hmmscan_002(self):
        """Test hmmscan-domtab indexing, HMMER 3.0, single query, no hits"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_002.out')
        self.check_index(filename, 'hmmscan3-domtab')

    def test_hmmerdomtab_30_hmmscan_003(self):
        """Test hmmscan-domtab indexing, HMMER 3.0, single query, multiple hits"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_003.out')
        self.check_index(filename, 'hmmscan3-domtab')

    def test_hmmerdomtab_30_hmmscan_004(self):
        """Test hmmscan-domtab indexing, HMMER 3.0, single query, no alignments"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_004.out')
        self.check_index(filename, 'hmmscan3-domtab')

    def test_hmmerdomtab_30_hmmsearch_001(self):
        """Test hmmsearch-domtab indexing, HMMER 3.0, single query, no alignments"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmsearch_001.out')
        self.check_index(filename, 'hmmsearch3-domtab')


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
