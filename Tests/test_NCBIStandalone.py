# Copyright 1999-2009 by Jeffrey Chang and Michiel de Hoon.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from Bio.Blast import NCBIStandalone


class TestNCBIStandalone(unittest.TestCase):

    def test_command_injection(self):
        #Check the simple detection of command injection,
        for func in [NCBIStandalone.blastall,
                     NCBIStandalone.blastpgp,
                     NCBIStandalone.rpsblast]:
            try:
                handle = func("/somewhere/blast", "blastz", "nr",
                              "/tmp/example.fasta",
                              expectation=10**-4,
                              matrix="IDENTITY -F 0; cat /etc/passwd'")
            except ValueError, e:
                self.assertEqual(str(e), "Rejecting suspicious argument for matrix")
                #Good
            else:
                self.fail("Attempted command injection not caught!")

    def test_pipe_redirection(self):
        #Now check something similar using pipe redirection
        for func in [NCBIStandalone.blastall,
                     NCBIStandalone.blastpgp,
                     NCBIStandalone.rpsblast]:
            try:
                handle = func("/somewhere/blast", "blastz", "nr",
                              "/tmp/example.fasta",
                              nprocessors=4,
                              expectation="0.001",
                              filter= "F > /etc/passwd'")
            except ValueError, e:
                self.assertEqual(str(e), "Rejecting suspicious argument for filter")
                #Good
            else:
                self.fail("Attempted output redirection not caught!")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
