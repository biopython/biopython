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
    -v;--verbose  -- run tests with higher verbosity
    <test_name>   -- supply the name of one (or more) tests to be run.
                     The .py file extension is optional.
    doctest       -- run the docstring tests.

<test_name> and doctest, if supplied, must come after the other options.
By default, all tests are run.
"""


# standard modules
import sys
import os
import getopt
import time
import traceback
import unittest
import doctest
import distutils.util
import gc
from io import BytesIO
from pkgutil import iter_modules
from setuptools import find_packages
from io import StringIO

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


# The default verbosity (not verbose)
VERBOSITY = 0

# Bio.Alphabet has been removed from Biopython. Importing it will raise an
# informative ImportError.
EXCLUDE_DOCTEST_MODULES = ["Bio.Alphabet"]

# Exclude modules with online activity
# They are not excluded by default, use --offline to exclude them
ONLINE_DOCTEST_MODULES = [
    "Bio.Entrez",
    "Bio.ExPASy",
    "Bio.TogoWS",
]

# Silently ignore any doctests for modules requiring numpy!
if numpy is None:
    EXCLUDE_DOCTEST_MODULES.extend(
        [
            "Bio.Affy.CelFile",
            "Bio.Align",
            "Bio.Align.substitution_matrices",
            "Bio.Cluster",
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
            "Bio.phenotype",
            "Bio.phenotype.parse",
            "Bio.phenotype.phen_micro",
            "Bio.phenotype.pm_fitting",
            "Bio.SeqIO.PdbIO",
            "Bio.SVDSuperimposer",
        ]
    )


try:
    import sqlite3

    del sqlite3
except ImportError:
    # May be missing on self-compiled Python
    EXCLUDE_DOCTEST_MODULES.append("Bio.SeqIO")
    EXCLUDE_DOCTEST_MODULES.append("Bio.SearchIO")


def find_modules(path):
    modules = set()
    for pkg in find_packages(path):
        modules.add(pkg)
        pkgpath = path + "/" + pkg.replace(".", "/")
        for info in iter_modules([pkgpath]):
            if not info.ispkg:
                modules.add(pkg + "." + info.name)
    return modules


SYSTEM_LANG = os.environ.get("LANG", "C")  # Cache this


def main(argv):
    """Run tests, return number of failures (integer)."""
    # Using "export LANG=C" (which should work on Linux and similar) can
    # avoid problems detecting optional command line tools on
    # non-English OS (we may want 'command not found' in English).
    # HOWEVER, we do not want to change the default encoding which is
    # rather important on Python 3 with unicode.
    # lang = os.environ['LANG']

    # get the command line options
    try:
        opts, args = getopt.getopt(
            argv, "gv", ["generate", "verbose", "doctest", "help", "offline"]
        )
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
            EXCLUDE_DOCTEST_MODULES.extend(ONLINE_DOCTEST_MODULES)
            # This is a bit of a hack...
            import requires_internet

            requires_internet.check.available = False
            # Monkey patch for urlopen()
            import urllib.request

            def dummy_urlopen(url):
                raise RuntimeError(
                    "Internal test suite error, attempting to use internet despite --offline setting"
                )

            urllib.request.urlopen = dummy_urlopen

        if opt == "-v" or opt == "--verbose":
            verbosity = 2

    # deal with the arguments, which should be names of tests to run
    for arg_num in range(len(args)):
        # strip off the .py if it was included
        if args[arg_num][-3:] == ".py":
            args[arg_num] = args[arg_num][:-3]

    print(f"Python version: {sys.version}")
    print(f"Operating system: {os.name} {sys.platform}")

    # run the tests
    runner = TestRunner(args, verbosity)
    return runner.run()


class TestRunner(unittest.TextTestRunner):

    if __name__ == "__main__":
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
            modules = find_modules(self.testdir + "/..")
            modules.difference_update(set(EXCLUDE_DOCTEST_MODULES))
            self.tests.extend(sorted(modules))
        stream = StringIO()
        unittest.TextTestRunner.__init__(self, stream, verbosity=verbosity)

    def runTest(self, name):
        from Bio import MissingPythonDependencyError
        from Bio import MissingExternalDependencyError

        result = self._makeResult()
        output = StringIO()
        # Restore the language and thus default encoding (in case a prior
        # test changed this, e.g. to help with detecting command line tools)
        global SYSTEM_LANG
        os.environ["LANG"] = SYSTEM_LANG
        # Always run tests from the Tests/ folder where run_tests.py
        # should be located (as we assume this with relative paths etc)
        os.chdir(self.testdir)
        try:
            stdout = sys.stdout
            sys.stdout = output
            if name.startswith("test_"):
                # It's a unittest
                sys.stderr.write(f"{name} ... ")
                loader = unittest.TestLoader()
                suite = loader.loadTestsFromName(name)
                if hasattr(loader, "errors") and loader.errors:
                    # New in Python 3.5, don't always get an exception anymore
                    # Instead this is a list of error messages as strings
                    for msg in loader.errors:
                        if (
                            "Bio.MissingExternalDependencyError: " in msg
                            or "Bio.MissingPythonDependencyError: " in msg
                        ):
                            # Remove the traceback etc
                            msg = msg[msg.find("Bio.Missing") :]
                            msg = msg[msg.find("Error: ") :]
                            sys.stderr.write(f"skipping. {msg}\n")
                            return True
                    # Looks like a real failure
                    sys.stderr.write("loading tests failed:\n")
                    for msg in loader.errors:
                        sys.stderr.write(f"{msg}\n")
                    return False
                if suite.countTestCases() == 0:
                    raise RuntimeError(f"No tests found in {name}")
            else:
                # It's a doc test
                sys.stderr.write(f"{name} docstring test ... ")
                try:
                    module = __import__(name, fromlist=name.split("."))
                except MissingPythonDependencyError:
                    sys.stderr.write("skipped, missing Python dependency\n")
                    return True
                except ImportError as e:
                    sys.stderr.write("FAIL, ImportError\n")
                    result.stream.write(f"ERROR while importing {name}: {e}\n")
                    result.printErrors()
                    return False
                suite = doctest.DocTestSuite(module, optionflags=doctest.ELLIPSIS)
                del module
            suite.run(result)
            if self.testdir != os.path.abspath("."):
                sys.stderr.write("FAIL\n")
                result.stream.write(result.separator1 + "\n")
                result.stream.write(f"ERROR: {name}\n")
                result.stream.write(result.separator2 + "\n")
                result.stream.write("Current directory changed\n")
                result.stream.write(f"Was: {self.testdir}\n")
                result.stream.write(f"Now: {os.path.abspath('.')}\n")
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
            sys.stderr.write(f"skipping. {msg}\n")
            return True
        except Exception as msg:
            # This happened during the import
            sys.stderr.write("ERROR\n")
            result.stream.write(result.separator1 + "\n")
            result.stream.write(f"ERROR: {name}\n")
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
        sys.stderr.write("-" * 70 + "\n")
        sys.stderr.write(
            "Ran %d test%s in %.3f seconds\n"
            % (total, total != 1 and "s" or "", time_taken)
        )
        sys.stderr.write("\n")
        if failures:
            sys.stderr.write("FAILED (failures = %d)\n" % failures)
        return failures


if __name__ == "__main__":
    errors = main(sys.argv[1:])
    if errors:
        # Doing a sys.exit(...) isn't nice if run from IDLE...
        sys.exit(1)
