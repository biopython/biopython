# Copyright 2008-2011 by Peter Cock.  All rights reserved.
# Revisions copyright 2012 by Christian Brueffer.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import MissingExternalDependencyError

import sys
import os
import shutil
import unittest
from Bio import SeqIO
from Bio import AlignIO
from Bio._py3k import getoutput
from Bio.Align.Applications import Guidance2Commandline
from Bio.Application import ApplicationError

#################################################################

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

guidance_exe = None
if sys.platform == "win32":
    raise MissingExternalDependencyError("Testing with Guidance2 not implemented on Windows yet")

try:
    output = getoutput("guidance")
    if output.startswith("USAGE: --seqFile <seqFile>"):
        guidance_exe = "guidance"
except FileNotFoundError:
    pass

if not guidance_exe:
    raise MissingExternalDependencyError(
        "Install guidance and put it in your PATH if you want to use Guidance2 from Biopython.")

##################################################################################################


class Guidance2TestCase(unittest.TestCase):
    """
    --seqFile ADGRB1.ffn --msaProgram CLUSTALW --seqType codon --outDir ~/Guidance2/data
        --bootstraps 20 --seqCutoff 0.63 --colCutoff 0.9 --outOrder as_input --dataset ADGRB1
    """

    def setUp(self):
        self.input = {"NA_seqFile": "Guidance2/Seqs/ADGRB1.ffn",
                      "AA_seqFile": "Guidance2/Seqs/ADGRB1.faa"}
        self.output = {"outDir": "Guidance2/Output"}

    def tearDown(self):
        outDir = self.output["outDir"]
        for item in os.listdir(outDir):
            if "Results" in item:
                continue
            else:
                item_path = os.path.join(outDir, item)
                if os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                if os.path.isfile(item_path):
                    os.remove(item_path)

    def standard_test_procedure(self, cline=Guidance2Commandline):
        """Standard testing procedure used by all tests."""

        output_path = cline.outDir

        # input_records = SeqIO.to_dict(SeqIO.parse(cline.seqFile, "fasta"))
        # self.assertEqual(str(eval(repr(cline))), str(cline))
        # output, error = cline()
        # self.assertTrue(not output or output.strip().startswith("CLUSTAL"))
        #
        # # Test if ClustalOmega executed successfully.
        # self.assertTrue(error.strip() == "" or
        #                 error.startswith("WARNING: Sequence type is DNA.") or
        #                 error.startswith("WARNING: DNA alignment is still experimental."))
        #
        # # Check the output...
        # align = AlignIO.read(cline.outfile, "clustal")
        # output_records = SeqIO.to_dict(SeqIO.parse(cline.outfile, "clustal"))
        # self.assertEqual(len(set(input_records.keys())), len(set(output_records.keys())))
        # for record in align:
        #     self.assertEqual(str(record.seq), str(output_records[record.id].seq))
        #
        # # TODO - Try and parse this with Bio.Nexus?
        # if cline.guidetree_out:
        #     self.assertTrue(os.path.isfile(cline.guidetree_out))

