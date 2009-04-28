# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import unittest

from Bio.Application import generic_run
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
from Bio import AlignIO

#################################################################

muscle_exe = None
if sys.platform=="win32" :
    raise MissingExternalDependencyError("Testing with MUSCLE not implemented on Windows yet")
else :
    import commands
    output = commands.getoutput("muscle")
    if "not found" not in output and "MUSCLE" in output.upper() :
        muscle_exe = "muscle"

if not muscle_exe :
    raise MissingExternalDependencyError(\
        "Install MUSCLE if you want to use the Bio.Align.Applications wrapper.")

#################################################################

class SimpleAlignTest(unittest.TestCase) :
    """Simple MUSCLE tests."""

    def test_simple(self) :
        """Simple muscle call"""
        input_file = "Fasta/f002"
        self.assert_(os.path.isfile(input_file))
        records = list(SeqIO.parse(open(input_file),"fasta"))
        #Prepare the command...
        cline = MuscleCommandline(muscle_exe)
        cline.set_parameter("input", input_file)
        #Preserve input record order (makes checking output easier)
        cline.set_parameter("stable")
        #Use clustal output
        cline.set_parameter("clwstrict")
        #TODO - Fix the trailing space!
        self.assertEqual(str(cline).rstrip(), "muscle -in Fasta/f002 -clwstrict -stable")
        result, out_handle, err_handle = generic_run(cline)
        align = AlignIO.read(out_handle, "clustal")
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align) :
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq))

    def test_long(self) :
        """Simple muscle call using long file."""
        #Create a large input file by converting some of another example file
        temp_large_fasta_file = "temp_cw_prot.fasta"
        handle = open(temp_large_fasta_file, "w")
        records = list(SeqIO.parse(open("NBRF/Cw_prot.pir", "rU"), "pir"))[:40]
        SeqIO.write(records, handle, "fasta")
        handle.close()

        #Prepare the command...
        cline = MuscleCommandline(muscle_exe)
        cline.set_parameter("input", temp_large_fasta_file)
        #Preserve input record order
        cline.set_parameter("stable")
        #Use fast options
        cline.set_parameter("maxiters", 1)
        #cline.set_parameter("diags")
        #Use clustal output
        cline.set_parameter("clwstrict")
        #TODO - Fix the trailing space!
        self.assertEqual(str(cline).rstrip(), "muscle -in temp_cw_prot.fasta -maxiters 1 -clwstrict -stable")
        result, out_handle, err_handle = generic_run(cline)
        align = AlignIO.read(out_handle, "clustal")
        self.assertEqual(len(records), len(align))
        for old, new in zip(records, align) :
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq))
        os.remove(temp_large_fasta_file)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
