#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Run a set of PyUnit-based regression tests.

This will find all modules whose name is "test_*.py" in the test
directory, and run them.  Various command line options provide
additional facilities.

Command line options::

    --help        -- show usage info
    --offline     -- skip tests which require internet access
    -g;--generate -- write the output file for a test instead of comparing it.
                     The name of the test to write the output for must be
                     specified.
    -v;--verbose  -- run tests with higher verbosity (does not affect our
                     print-and-compare style unit tests).
    <test_name>   -- supply the name of one (or more) tests to be run.
                     The .py file extension is optional.
    doctest       -- run the docstring tests.

By default, all tests are run.
"""

from __future__ import print_function

# standard modules
import sys
import os
import re
import getopt
import time
import traceback
import unittest
import doctest
import distutils.util
import gc
from io import BytesIO


# Note, we want to be able to call run_tests.py BEFORE
# Biopython is installed, so we can't use this:
# from Bio._py3k import StringIO
try:
    from StringIO import StringIO  # Python 2 (byte strings)
except ImportError:
    from io import StringIO  # Python 3 (unicode strings)


try:
    import numpy
    try:
        # NumPy 1.14 changed repr output breaking our doctests,
        # request the legacy 1.13 style
        numpy.set_printoptions(legacy="1.13")
    except TypeError:
        # Old Numpy, output should be fine as it is :)
        # TypeError: set_printoptions() got an unexpected keyword argument 'legacy'
        pass
except ImportError:
    numpy = None


def is_pypy():
    import platform
    try:
        if platform.python_implementation() == 'PyPy':
            return True
    except AttributeError:
        # New in Python 2.6, not in Jython yet either
        pass
    return False


# The default verbosity (not verbose)
VERBOSITY = 0

# This is the list of modules containing docstring tests.
# If you develop docstring tests for other modules, please add
# those modules here. Please sort names alphabetically.
DOCTEST_MODULES = [
    "Bio",
    "Bio._utils",
    "Bio.Affy",
    "Bio.Align",
    "Bio.Align._aligners",
    "Bio.Align.AlignInfo",
    "Bio.Align.Applications._Clustalw",
    "Bio.Align.Applications._ClustalOmega",
    "Bio.Align.Applications._Dialign",
    "Bio.Align.Applications._MSAProbs",
    "Bio.Align.Applications._Mafft",
    "Bio.Align.Applications._Muscle",
    "Bio.Align.Applications._Probcons",
    "Bio.Align.Applications._Prank",
    "Bio.Align.Applications._TCoffee",
    "Bio.AlignIO",
    "Bio.AlignIO.ClustalIO",
    "Bio.AlignIO.EmbossIO",
    "Bio.AlignIO.FastaIO",
    "Bio.AlignIO.Interfaces",
    "Bio.AlignIO.MafIO",
    "Bio.AlignIO.NexusIO",
    "Bio.AlignIO.PhylipIO",
    "Bio.AlignIO.StockholmIO",
    "Bio.Alphabet",
    "Bio.Alphabet.IUPAC",
    "Bio.Alphabet.Reduced",
    "Bio.Application",
    "Bio.bgzf",
    "Bio.Blast",
    "Bio.Blast.Applications",
    "Bio.Blast.NCBIWWW",
    "Bio.Blast.NCBIXML",
    "Bio.Blast.ParseBlastTable",
    "Bio.Blast.Record",
    "Bio.CAPS",
    "Bio.codonalign",
    "Bio.codonalign.chisq",
    "Bio.codonalign.codonalignment",
    "Bio.codonalign.codonalphabet",
    "Bio.codonalign.codonseq",
    "Bio.Compass",
    "Bio.Crystal",
    "Bio.Data",
    "Bio.Data.IUPACData",
    "Bio.Data.SCOPData",
    "Bio.Emboss",
    "Bio.Emboss.Applications",
    "Bio.Emboss.Primer3",
    "Bio.Emboss.PrimerSearch",
    "Bio.Entrez.Parser",
    "Bio.ExPASy.Enzyme",
    "Bio.ExPASy.Prodoc",
    "Bio.ExPASy.Prosite",
    "Bio.ExPASy.ScanProsite",
    "Bio.FSSP",
    "Bio.FSSP.fssp_rec",
    "Bio.FSSP.FSSPTools",
    "Bio.GenBank",
    "Bio.GenBank.Record",
    "Bio.GenBank.Scanner",
    "Bio.GenBank.utils",
    "Bio.Geo",
    "Bio.Geo.Record",
    "Bio.Graphics",
    "Bio.Graphics.BasicChromosome",
    "Bio.Graphics.ColorSpiral",
    "Bio.Graphics.Comparative",
    "Bio.Graphics.DisplayRepresentation",
    "Bio.Graphics.Distribution",
    "Bio.Graphics.GenomeDiagram._AbstractDrawer",
    "Bio.Graphics.GenomeDiagram._CircularDrawer",
    "Bio.Graphics.GenomeDiagram._Colors",
    "Bio.Graphics.GenomeDiagram._CrossLink",
    "Bio.Graphics.GenomeDiagram._Diagram",
    "Bio.Graphics.GenomeDiagram._Feature",
    "Bio.Graphics.GenomeDiagram._FeatureSet",
    "Bio.Graphics.GenomeDiagram._Graph",
    "Bio.Graphics.GenomeDiagram._GraphSet",
    "Bio.Graphics.GenomeDiagram._LinearDrawer",
    "Bio.Graphics.GenomeDiagram._Track",
    "Bio.Graphics.KGML_vis",
    "Bio.HMM",
    "Bio.HMM.DynamicProgramming",
    "Bio.HMM.MarkovModel",
    "Bio.HMM.Trainer",
    "Bio.HMM.Utilities",
    "Bio.Index",
    "Bio.KEGG",
    "Bio.KEGG.Compound",
    "Bio.KEGG.Enzyme",
    "Bio.KEGG.Gene",
    "Bio.KEGG.KGML",
    "Bio.KEGG.KGML.KGML_parser",
    "Bio.KEGG.KGML.KGML_pathway",
    "Bio.KEGG.Map",
    "Bio.KEGG.REST",
    "Bio.motifs",
    "Bio.motifs.alignace",
    "Bio.motifs.applications",
    "Bio.motifs.applications._xxmotif",
    "Bio.motifs.jaspar",
    "Bio.motifs.matrix",
    "Bio.motifs.thresholds",
    "Bio.motifs.transfac",
    "Bio.Nexus",
    "Bio.Nexus.Nexus",
    "Bio.Nexus.Nodes",
    "Bio.Nexus.StandardData",
    "Bio.Nexus.Trees",
    "Bio.NMR",
    "Bio.NMR.NOEtools",
    "Bio.NMR.xpktools",
    "Bio.pairwise2",
    "Bio.Pathway",
    "Bio.Pathway.Rep.Graph",
    "Bio.Pathway.Rep",
    "Bio.Pathway.Rep.MultiGraph",
    "Bio.pairwise2",
    "Bio.Phylo",
    "Bio.Phylo.Applications",
    "Bio.Phylo.Applications._Phyml",
    "Bio.Phylo.Applications._Raxml",
    "Bio.Phylo.BaseTree",
    "Bio.Phylo.CDAOIO",
    "Bio.Phylo._cdao_owl",
    "Bio.Phylo.CDAO",
    "Bio.Phylo.Consensus",
    "Bio.Phylo.NewickIO",
    "Bio.Phylo.Newick",
    "Bio.Phylo.NeXMLIO",
    "Bio.Phylo.NeXML",
    "Bio.Phylo.NexusIO",
    "Bio.Phylo.PAML.baseml",
    "Bio.Phylo.PAML.chi2",
    "Bio.Phylo.PAML.codeml",
    "Bio.Phylo.PAML",
    "Bio.Phylo.PAML._paml",
    "Bio.Phylo.PAML._parse_baseml",
    "Bio.Phylo.PAML._parse_codeml",
    "Bio.Phylo.PAML._parse_yn00",
    "Bio.Phylo.PAML.yn00",
    "Bio.Phylo.PhyloXMLIO",
    "Bio.Phylo.PhyloXML",
    "Bio.PopGen.GenePop.Controller",
    "Bio.PopGen.GenePop.EasyController",
    "Bio.PopGen.GenePop.FileParser",
    "Bio.PopGen.GenePop",
    "Bio.PopGen.GenePop.LargeFileParser",
    "Bio.PopGen",
    "Bio._py3k",
    "Bio.Restriction.RanaConfig",
    "Bio.SCOP.Cla",
    "Bio.SCOP.Des",
    "Bio.SCOP.Dom",
    "Bio.SCOP.Hie",
    "Bio.SCOP",
    "Bio.SCOP.Raf",
    "Bio.SCOP.Residues",
    "Bio.SearchIO",
    "Bio.SearchIO.BlastIO.blast_tab",
    "Bio.SearchIO.BlastIO.blast_text",
    "Bio.SearchIO.BlastIO.blast_xml",
    "Bio.SearchIO.BlastIO",
    "Bio.SearchIO.BlatIO",
    "Bio.SearchIO.ExonerateIO.exonerate_cigar",
    "Bio.SearchIO.ExonerateIO.exonerate_text",
    "Bio.SearchIO.ExonerateIO.exonerate_vulgar",
    "Bio.SearchIO.ExonerateIO",
    "Bio.SearchIO.FastaIO",
    "Bio.SearchIO.HmmerIO._base",
    "Bio.SearchIO.HmmerIO.hmmer2_text",
    "Bio.SearchIO.HmmerIO.hmmer3_domtab",
    "Bio.SearchIO.HmmerIO.hmmer3_tab",
    "Bio.SearchIO.HmmerIO.hmmer3_text",
    "Bio.SearchIO.HmmerIO",
    "Bio.SearchIO._index",
    "Bio.SearchIO.InterproscanIO",
    "Bio.SearchIO.InterproscanIO.interproscan_xml",
    "Bio.SearchIO._legacy",
    "Bio.SearchIO._legacy.NCBIStandalone",
    "Bio.SearchIO._legacy.ParserSupport",
    "Bio.SearchIO._model._base",
    "Bio.SearchIO._model.hit",
    "Bio.SearchIO._model.hsp",
    "Bio.SearchIO._model",
    "Bio.SearchIO._model.query",
    "Bio.SearchIO._utils",
    "Bio.Seq",
    "Bio.SeqFeature",
    "Bio.SeqIO",
    "Bio.SeqIO.AbiIO",
    "Bio.SeqIO.AceIO",
    "Bio.SeqIO._convert",
    "Bio.SeqIO.FastaIO",
    "Bio.SeqIO.IgIO",
    "Bio.SeqIO._index",
    "Bio.SeqIO.InsdcIO",
    "Bio.SeqIO.Interfaces",
    "Bio.SeqIO.PhdIO",
    "Bio.SeqIO.PirIO",
    "Bio.SeqIO.QualityIO",
    "Bio.SeqIO.SeqXmlIO",
    "Bio.SeqIO.SffIO",
    "Bio.SeqIO.SwissIO",
    "Bio.SeqIO.TabIO",
    "Bio.SeqIO.UniprotIO",
    "Bio.SeqRecord",
    "Bio.Sequencing.Ace",
    "Bio.Sequencing.Applications._bwa",
    "Bio.Sequencing.Applications",
    "Bio.Sequencing.Applications._Novoalign",
    "Bio.Sequencing.Applications._samtools",
    "Bio.Sequencing",
    "Bio.Sequencing.Phd",
    "Bio.SeqUtils",
    "Bio.SeqUtils.CheckSum",
    "Bio.SeqUtils.CodonUsageIndices",
    "Bio.SeqUtils.CodonUsage",
    "Bio.SeqUtils.IsoelectricPoint",
    "Bio.SeqUtils.lcc",
    "Bio.SeqUtils.MeltingTemp",
    "Bio.SeqUtils.ProtParamData",
    "Bio.SeqUtils.ProtParam",
    "BioSQL",
    "BioSQL.BioSeq",
    "BioSQL.DBUtils",
    "BioSQL.Loader",
    "Bio.Statistics",
    "Bio.SubsMat.MatrixInfo",
    "Bio.SwissProt",
    "Bio.triefind",
    "Bio.UniGene",
    "Bio.UniProt",
    "Bio.UniProt.GOA",
    "Bio._utils",
    "Bio.Wise",
    "Bio.Wise.dnal",
    "Bio.Wise.psw",
]
# Silently ignore any doctests for modules requiring numpy!
if numpy:
    DOCTEST_MODULES.extend([
        "Bio.Affy.CelFile",
        "Bio.KDTree",
        "Bio.KDTree.KDTree",
        "Bio.kNN",
        "Bio.LogisticRegression",
        "Bio.MarkovModel",
        "Bio.MaxEntropy",
        "Bio.NaiveBayes",
        "Bio.PDB.Chain",
        "Bio.PDB.Dice",
        "Bio.PDB.HSExposure",
        "Bio.PDB.MMCIF2Dict",
        "Bio.PDB.MMCIFParser",
        "Bio.PDB.mmtf.DefaultParser",
        "Bio.PDB.mmtf",
        "Bio.PDB.Model",
        "Bio.PDB.NACCESS",
        "Bio.PDB.NeighborSearch",
        "Bio.PDB.parse_pdb_header",
        "Bio.PDB.PDBExceptions",
        "Bio.PDB.PDBList",
        "Bio.PDB.PDBParser",
        "Bio.PDB.Polypeptide",
        "Bio.PDB.PSEA",
        "Bio.PDB.QCPSuperimposer",
        "Bio.PDB.Residue",
        "Bio.PDB.Selection",
        "Bio.PDB.StructureAlignment",
        "Bio.PDB.StructureBuilder",
        "Bio.PDB.Structure",
        "Bio.PDB.Superimposer",
        "Bio.PDB.Vector",
        "Bio.phenotype.pm_fitting",
        "Bio.SeqIO.PdbIO",
        "Bio.Statistics.lowess",
        "Bio.SVDSuperimposer",
    ])


try:
    import sqlite3
    del sqlite3
except ImportError:
    # Missing on Jython or Python 2.4
    DOCTEST_MODULES.remove("Bio.SeqIO")
    DOCTEST_MODULES.remove("Bio.SearchIO")

# Skip Bio.Seq doctest under Python 2, see http://bugs.python.org/issue7490
# Can't easily write exceptions with consistent class name in python 2 and 3
if sys.version_info[0] == 2:
    DOCTEST_MODULES.remove("Bio.Seq")


# Skip Bio.bgzf doctest for broken gzip, see http://bugs.python.org/issue17666
def _have_bug17666():
    """Debug function to check if Python's gzip is broken (PRIVATE).

    Checks for http://bugs.python.org/issue17666 expected in Python 2.7.4,
    3.2.4 and 3.3.1 only.
    """
    if os.name == 'java':
        # Jython not affected
        return False
    import gzip
    # Would like to use byte literal here:
    bgzf_eof = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC" + \
               "\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    if sys.version_info[0] >= 3:
        import codecs
        bgzf_eof = codecs.latin_1_encode(bgzf_eof)[0]
    handle = gzip.GzipFile(fileobj=BytesIO(bgzf_eof))
    try:
        data = handle.read()
        handle.close()
        assert not data, "Should be zero length, not %i" % len(data)
        return False
    except TypeError as err:
        # TypeError: integer argument expected, got 'tuple'
        handle.close()
        return True


if _have_bug17666():
    DOCTEST_MODULES.remove("Bio.bgzf")

SYSTEM_LANG = os.environ.get('LANG', 'C')  # Cache this


def main(argv):
    """Run tests, return number of failures (integer)."""
    # insert our paths in sys.path:
    # ../build/lib.*
    # ..
    # Q. Why this order?
    # A. To find the C modules (which are in ../build/lib.*/Bio)
    # Q. Then, why ".."?
    # A. Because Martel may not be in ../build/lib.*
    test_path = sys.path[0] or "."
    source_path = os.path.abspath("%s/.." % test_path)
    sys.path.insert(1, source_path)
    build_path = os.path.abspath("%s/../build/lib.%s-%s" % (
        test_path, distutils.util.get_platform(), sys.version[:3]))
    if os.access(build_path, os.F_OK):
        sys.path.insert(1, build_path)

    # Using "export LANG=C" (which should work on Linux and similar) can
    # avoid problems detecting optional command line tools on
    # non-English OS (we may want 'command not found' in English).
    # HOWEVER, we do not want to change the default encoding which is
    # rather important on Python 3 with unicode.
    # lang = os.environ['LANG']

    # get the command line options
    try:
        opts, args = getopt.getopt(argv, 'gv', ["generate", "verbose",
                                                "doctest", "help", "offline"])
    except getopt.error as msg:
        print(msg)
        print(__doc__)
        return 2

    verbosity = VERBOSITY

    # deal with the options
    for opt, _ in opts:
        if opt == "--help":
            print(__doc__)
            return 0
        if opt == "--offline":
            print("Skipping any tests requiring internet access")
            # This is a bit of a hack...
            import requires_internet
            requires_internet.check.available = False
            # Monkey patch for urlopen()
            import Bio._py3k

            def dummy_urlopen(url):
                raise RuntimeError("Internal test suite error, attempting to use internet despite --offline setting")

            Bio._py3k.urlopen = dummy_urlopen
        if opt == "-g" or opt == "--generate":
            if len(args) > 1:
                print("Only one argument (the test name) needed for generate")
                print(__doc__)
                return 2
            elif len(args) == 0:
                print("No test name specified to generate output for.")
                print(__doc__)
                return 2
            # strip off .py if it was included
            if args[0][-3:] == ".py":
                args[0] = args[0][:-3]

            test = ComparisonTestCase(args[0])
            test.generate_output()
            return 0

        if opt == "-v" or opt == "--verbose":
            verbosity = 2

    # deal with the arguments, which should be names of tests to run
    for arg_num in range(len(args)):
        # strip off the .py if it was included
        if args[arg_num][-3:] == ".py":
            args[arg_num] = args[arg_num][:-3]

    print("Python version: %s" % sys.version)
    print("Operating system: %s %s" % (os.name, sys.platform))

    # run the tests
    runner = TestRunner(args, verbosity)
    return runner.run()


class ComparisonTestCase(unittest.TestCase):
    """Run a print-and-compare test and check it against expected output."""

    def __init__(self, name, output=None):
        """Initialize with the test to run.

        Arguments:
            - name - The name of the test. The expected output should be
              stored in the file output/name.
            - output - The output that was generated when this test was run.

        """
        unittest.TestCase.__init__(self)
        self.name = name
        self.output = output

    def shortDescription(self):
        return self.name

    def runTest(self):
        # check the expected output to be consistent with what
        # we generated
        outputdir = os.path.join(TestRunner.testdir, "output")
        outputfile = os.path.join(outputdir, self.name)
        try:
            if sys.version_info[0] >= 3:
                # Python 3 problem: Can't use utf8 on output/test_geo
                # due to micro (\xb5) and degrees (\xb0) symbols
                # Also universal new lines mode deprecated on Python 3
                expected = open(outputfile, encoding="latin")
            else:
                expected = open(outputfile, "rU")
        except IOError:
            self.fail("Warning: Can't open %s for test %s" %
                      (outputfile, self.name))

        self.output.seek(0)
        # first check that we are dealing with the right output
        # the first line of the output file is the test name
        expected_test = expected.readline().strip()

        if expected_test != self.name:
            expected.close()
            raise ValueError("\nOutput:   %s\nExpected: %s"
                             % (self.name, expected_test))

        # Track the line number. Starts at 1 to account for the output file
        # header line.
        line_number = 1

        # now loop through the output and compare it to the expected file
        while True:
            expected_line = expected.readline()
            output_line = self.output.readline()
            line_number += 1

            # stop looping if either of the info handles reach the end
            if (not expected_line) or (not output_line):
                # make sure both have no information left
                assert expected_line == '', "Unread: %s" % expected_line
                assert output_line == '', "Extra output: %s" % output_line
                break

            # normalize the newlines in the two lines
            expected_line = expected_line.strip("\r\n")
            output_line = output_line.strip("\r\n")

            # if the line is a doctest or PyUnit time output like:
            # Ran 2 tests in 0.285s
            # ignore it, so we don't have problems with different running times
            if re.compile("^Ran [0-9]+ tests? in ").match(expected_line):
                pass
            # otherwise make sure the two lines are the same
            elif expected_line != output_line:
                expected.close()
                raise ValueError("\nOutput  : %s\nExpected: %s\n%s line %s"
                                 % (repr(output_line), repr(expected_line),
                                    outputfile, line_number))

        expected.close()

    def generate_output(self):
        """Generate the golden output for the specified test."""
        outputdir = os.path.join(TestRunner.testdir, "output")
        outputfile = os.path.join(outputdir, self.name)

        output_handle = open(outputfile, 'w')

        # write the test name as the first line of the output
        output_handle.write(self.name + "\n")

        # remember standard out so we can reset it after we are done
        save_stdout = sys.stdout
        try:
            # write the output from the test into a string
            sys.stdout = output_handle
            __import__(self.name)
        finally:
            output_handle.close()
            # return standard out to its normal setting
            sys.stdout = save_stdout


class TestRunner(unittest.TextTestRunner):

    if __name__ == '__main__':
        file = sys.argv[0]
    else:
        file = __file__
    testdir = os.path.abspath(os.path.dirname(file) or os.curdir)

    def __init__(self, tests=None, verbosity=0):
        """Initialise test runner.

        If not tests are specified, we run them all,
        including the doctests.

        Defaults to running without any verbose logging.
        """
        if tests is None:
            self.tests = []
        else:
            self.tests = tests
        if not self.tests:
            # Make a list of all applicable test modules.
            names = os.listdir(TestRunner.testdir)
            for name in names:
                if name[:5] == "test_" and name[-3:] == ".py":
                    self.tests.append(name[:-3])
            self.tests.sort()
            self.tests.append("doctest")
        if "doctest" in self.tests:
            self.tests.remove("doctest")
            self.tests.extend(DOCTEST_MODULES)
        stream = StringIO()
        unittest.TextTestRunner.__init__(self, stream,
                                         verbosity=verbosity)

    def runTest(self, name):
        from Bio import MissingExternalDependencyError
        result = self._makeResult()
        output = StringIO()
        # Restore the language and thus default encoding (in case a prior
        # test changed this, e.g. to help with detecting command line tools)
        global SYSTEM_LANG
        os.environ['LANG'] = SYSTEM_LANG
        # Always run tests from the Tests/ folder where run_tests.py
        # should be located (as we assume this with relative paths etc)
        os.chdir(self.testdir)
        try:
            stdout = sys.stdout
            sys.stdout = output
            if name.startswith("test_"):
                sys.stderr.write("%s ... " % name)
                # It's either a unittest or a print-and-compare test
                loader = unittest.TestLoader()
                suite = loader.loadTestsFromName(name)
                if hasattr(loader, "errors") and loader.errors:
                    # New in Python 3.5, don't always get an exception anymore
                    # Instead this is a list of error messages as strings
                    for msg in loader.errors:
                        if "Bio.MissingExternalDependencyError: " in msg or \
                                "Bio.MissingPythonDependencyError: " in msg:
                            # Remove the traceback etc
                            msg = msg[msg.find("Bio.Missing"):]
                            msg = msg[msg.find("Error: "):]
                            sys.stderr.write("skipping. %s\n" % msg)
                            return True
                    # Looks like a real failure
                    sys.stderr.write("loading tests failed:\n")
                    for msg in loader.errors:
                        sys.stderr.write("%s\n" % msg)
                    return False
                if suite.countTestCases() == 0:
                    # This is a print-and-compare test instead of a
                    # unittest-type test.
                    test = ComparisonTestCase(name, output)
                    suite = unittest.TestSuite([test])
            else:
                # It's a doc test
                sys.stderr.write("%s docstring test ... " % name)
                module = __import__(name, fromlist=name.split("."))
                suite = doctest.DocTestSuite(module,
                                             optionflags=doctest.ELLIPSIS)
                del module
            suite.run(result)
            if self.testdir != os.path.abspath("."):
                sys.stderr.write("FAIL\n")
                result.stream.write(result.separator1 + "\n")
                result.stream.write("ERROR: %s\n" % name)
                result.stream.write(result.separator2 + "\n")
                result.stream.write("Current directory changed\n")
                result.stream.write("Was: %s\n" % self.testdir)
                result.stream.write("Now: %s\n" % os.path.abspath("."))
                os.chdir(self.testdir)
                if not result.wasSuccessful():
                    result.printErrors()
                return False
            elif result.wasSuccessful():
                sys.stderr.write("ok\n")
                return True
            else:
                sys.stderr.write("FAIL\n")
                result.printErrors()
            return False
        except MissingExternalDependencyError as msg:
            # Seems this isn't always triggered on Python 3.5,
            # exception messages can be in loader.errors instead.
            sys.stderr.write("skipping. %s\n" % msg)
            return True
        except Exception as msg:
            # This happened during the import
            sys.stderr.write("ERROR\n")
            result.stream.write(result.separator1 + "\n")
            result.stream.write("ERROR: %s\n" % name)
            result.stream.write(result.separator2 + "\n")
            result.stream.write(traceback.format_exc())
            return False
        except KeyboardInterrupt as err:
            # Want to allow this, and abort the test
            # (see below for special case)
            raise err
        except:  # noqa: B901
            # This happens in Jython with java.lang.ClassFormatError:
            # Invalid method Code length ...
            sys.stderr.write("ERROR\n")
            result.stream.write(result.separator1 + "\n")
            result.stream.write("ERROR: %s\n" % name)
            result.stream.write(result.separator2 + "\n")
            result.stream.write(traceback.format_exc())
            return False
        finally:
            sys.stdout = stdout
            # Running under PyPy we were leaking file handles...
            gc.collect()

    def run(self):
        """Run tests, return number of failures (integer)."""
        failures = 0
        start_time = time.time()
        for test in self.tests:
            ok = self.runTest(test)
            if not ok:
                failures += 1
        total = len(self.tests)
        stop_time = time.time()
        time_taken = stop_time - start_time
        sys.stderr.write(self.stream.getvalue())
        sys.stderr.write('-' * 70 + "\n")
        sys.stderr.write("Ran %d test%s in %.3f seconds\n" %
                         (total, total != 1 and "s" or "", time_taken))
        sys.stderr.write("\n")
        if failures:
            sys.stderr.write("FAILED (failures = %d)\n" % failures)
        return failures


if __name__ == "__main__":
    errors = main(sys.argv[1:])
    if errors:
        # Doing a sys.exit(...) isn't nice if run from IDLE...
        sys.exit(1)
