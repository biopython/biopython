# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Tests for SearchIO InfernalIO infernal-tab indexing"""

import os
import unittest

from search_tests_common import CheckIndex
from search_tests_common import CheckRaw

FMT = "infernal-tab"


class InfernalTabRawCases(CheckRaw):
    fmt = FMT

    def test_infernal_tab_single(self):
        """Test infernal-tab raw string retrieval, cmsearch, single query (U2_Yeast)."""
        filename = os.path.join("Infernal", "cmsearch_114_U2_Yeast.tbl")
        raw = """ENA|BK006936|BK006936.2 -         U2                   RF00004    cm        1      193   681858   681747      -    no    1 0.33   0.1   98.7   5.9e-20 !   TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.
"""
        self.check_raw(filename, "U2", raw)

    def test_infernal_tab_multiple_first(self):
        """Test infernal-tab raw string retrieval, cmsearch, multiple queries, first."""
        filename = os.path.join("Infernal", "cmscan_115_IRES_5S_U2_Yeast.tbl")
        raw = """U2                   RF00004   ENA|BK006935|BK006935.2 -          cm        1      193    52929    53083      +    no    1 0.44   0.0   13.5      0.91 ?   U2 spliceosomal RNA
U2                   RF00004   ENA|BK006935|BK006935.2 -          cm        1      193   196571   196389      -    no    1 0.33   5.3   12.8       1.3 ?   U2 spliceosomal RNA
"""
        self.check_raw(filename, "ENA|BK006935|BK006935.2", raw)

    def test_infernal_tab_multiple_middle(self):
        """Test infernal-tab raw string retrieval, cmsearch, multiple queries, middle."""
        filename = os.path.join("Infernal", "cmscan_115_IRES_5S_U2_Yeast.tbl")
        raw = """U2                   RF00004   ENA|BK006936|BK006936.2 -          cm        1      193   681858   681747      -    no    1 0.33   0.1   98.7   1.2e-20 !   U2 spliceosomal RNA
"""
        self.check_raw(filename, "ENA|BK006936|BK006936.2", raw)

    def test_infernal_tab_multiple_last(self):
        """Test infernal-tab raw string retrieval, cmsearch, multiple queries, last."""
        filename = os.path.join("Infernal", "cmscan_115_IRES_5S_U2_Yeast.tbl")
        raw = """5S_rRNA              RF00001   ENA|BK006937|BK006937.2 -          cm        1      119      761      644      -    no    1 0.41   0.3   14.1       2.4 ?   5S ribosomal RNA
U2                   RF00004   ENA|BK006937|BK006937.2 -          cm        1      193   229986   229885      -    no    1 0.32   0.1   11.1       4.7 ?   U2 spliceosomal RNA
"""
        self.check_raw(filename, "ENA|BK006937|BK006937.2", raw)

    def test_infernal_tab_multiple_fmt_2(self):
        """Test infernal-tab raw string retrieval, cmsearch, multiple queries, fmt 2."""
        filename = os.path.join("Infernal", "cmscan_115_IRES_5S_U2_Yeast_fmt_2.tbl")
        raw = """1    U2                   RF00004   ENA|BK006936|BK006936.2 -         -          cm        1      193   681858   681747      -    no    1 0.33   0.1   98.7   1.2e-20  !   *       -      -      -      -      -      -     193  813184 U2 spliceosomal RNA
"""
        self.check_raw(filename, "ENA|BK006936|BK006936.2", raw)


class InfernalTabIndexCases(CheckIndex):
    def test_infernal_tab_1q(self):
        """Test infernal-tab indexing, cmsearch, one query, one hit."""
        filename = os.path.join("Infernal", "cmsearch_114_U2_Yeast.tbl")
        self.check_index(filename, FMT)

    def test_infernal_tab_1q_0m(self):
        """Test infernal-tab indexing, cmsearch, single query, no hits."""
        filename = os.path.join("Infernal", "cmsearch_114_IRES_Yeast.tbl")
        self.check_index(filename, FMT)

    def test_infernal_tab_1q_mm(self):
        """Test infernal-tab indexing, cmsearch, single query, multiple hits."""
        filename = os.path.join("Infernal", "cmsearch_114_5S_Yeast.tbl")
        self.check_index(filename, FMT)

    def test_infernal_tab_mq_mm(self):
        """Test infernal-tab indexing, cmscan, multiple query, multiple matches."""
        filename = os.path.join("Infernal", "cmscan_115_IRES_5S_U2_Yeast.tbl")
        self.check_index(filename, FMT)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
