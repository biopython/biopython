# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Test for Blast records."""

import unittest

from Bio.Blast.Record import HSP


class TestHsp(unittest.TestCase):
    def test_str(self):
        # Test empty instance
        self.assertEqual(
            str(HSP()),
            "Score <unknown> (<unknown> bits), expectation <unknown>, alignment length <unknown>",
        )

        # Test instance with non-default attributes
        hsp = HSP()
        hsp.score = 1.0
        hsp.bits = 2.0
        hsp.expect = 3.0
        hsp.align_length = 4
        # Ignore trailing whitespace in output
        self.assertEqual(
            "\n".join(line.strip() for line in str(hsp).split("\n")),
            """Score 1 (2 bits), expectation 3.0e+00, alignment length 4
Query:    None  None

Sbjct:    None  None""",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
