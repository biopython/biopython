#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Run a set of PyUnit-based regression tests.

This will find all modules whose name is "test_*.py" in the test
directory, and run them.  Various command line options provide
additional facilities.

Command line options:

--help        -- show usage info
-g;--generate -- write the output file for a test instead of comparing it.
                 The name of the  test to write the output for must be
                 specified.
<test_name>   -- supply the name of one (or more) tests to be run.
                 The .py file extension is optional.
doctest       -- run the docstring tests.
By default, all tests are run.
"""

# This is the list of modules containing docstring tests.
# If you develop docstring tests for other modules, please add
# those modules here.
DOCTEST_MODULES = ["Bio.Seq",
                   "Bio.SeqRecord",
                   "Bio.SeqIO",
                   "Bio.SeqIO.PhdIO",
                   "Bio.SeqIO.QualityIO",
                   "Bio.SeqIO.AceIO",
                   "Bio.SeqUtils",
                   "Bio.Align.Generic",
                   "Bio.AlignIO",
                   "Bio.AlignIO.StockholmIO",
                   "Bio.Application",
                   "Bio.Clustalw",
                   "Bio.KEGG.Compound",
                   "Bio.KEGG.Enzyme",
                   "Bio.Wise",
                   "Bio.Wise.psw",
                   "Bio.Motif",
                  ]
#Silently ignore any doctests for modules requiring numpy!
try:
    import numpy
    DOCTEST_MODULES.extend(["Bio.Statistics.lowess"])
except ImportError:
    pass

# standard modules
import sys
import cStringIO
import os
import re
import getopt
import time
import traceback
import unittest
import doctest
import distutils.util



def main(argv):
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
    
    # get the command line options
    try:
        opts, args = getopt.getopt(argv, 'g', ["generate", "doctest", "help"])
    except getopt.error, msg:
        print msg
        print __doc__
        return 2

    # deal with the options
    for o, a in opts:
        if o == "--help":
            print __doc__
            return 0
        if o == "-g" or o == "--generate":
            if len(args) > 1:
                print "Only one argument (the test name) needed for generate"
                print __doc__
                return 2
            elif len(args) == 0:
                print "No test name specified to generate output for."
                print __doc__
                return 2
            # strip off .py if it was included
            if args[0][-3:] == ".py":
                args[0] = args[0][:-3]

            test = ComparisonTestCase(args[0])
            test.generate_output()
            return 0

    # deal with the arguments, which should be names of tests to run
    for arg_num in range(len(args)):
        # strip off the .py if it was included
        if args[arg_num][-3:] == ".py":
            args[arg_num] = args[arg_num][:-3]

    # run the tests
    runner = TestRunner(args)
    runner.run()


class ComparisonTestCase(unittest.TestCase):
    """Run a print-and-compare test and compare its output against expected output.
    """

    def __init__(self, name, output=None):
        """Initialize with the test to run.

        Arguments:
        o name - The name of the test. The expected output should be
          stored in the file output/name.
        o output - The output that was generated when this test was run.
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
            expected = open(outputfile, 'r')
        except IOError:
            self.fail("Warning: Can't open %s for test %s" % (outputfile, self.name))

        self.output.seek(0)
        # first check that we are dealing with the right output
        # the first line of the output file is the test name
        expected_test = expected.readline().strip()

        assert expected_test == self.name, "\nOutput:   %s\nExpected: %s" % \
               (self.name, expected_test)

        # now loop through the output and compare it to the expected file
        while True:
            expected_line = expected.readline()
            output_line = self.output.readline()

            # stop looping if either of the info handles reach the end
            if not(expected_line) or not(output_line):
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
            else:
                assert expected_line == output_line, \
                      "\nOutput  : %s\nExpected: %s" \
                      % (repr(output_line), repr(expected_line))

    def generate_output(self):
        """Generate the golden output for the specified test.
        """
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
    testdir = os.path.dirname(file) or os.curdir

    def __init__(self, tests=[]):
        # if no tests were specified to run, we run them all
        # including the doctests
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
        stream = cStringIO.StringIO()
        unittest.TextTestRunner.__init__(self, stream, verbosity=0)

    def runTest(self, name):
        from Bio import MissingExternalDependencyError
        result = self._makeResult()
        output = cStringIO.StringIO()
        # Run the actual test inside a try/except to catch import errors.
        # Have to do a nested try because try/except/except/finally requires
        # python 2.5+
        try:
            try :
                stdout = sys.stdout
                sys.stdout = output
                if name.startswith("test_"):
                    sys.stderr.write("%s ... " % name)
                    #It's either a unittest or a print-and-compare test
                    suite = unittest.TestLoader().loadTestsFromName(name)
                    if suite.countTestCases()==0:
                        # This is a print-and-compare test instead of a
                        # unittest-type test.
                        test = ComparisonTestCase(name, output)
                        suite = unittest.TestSuite([test])
                else :
                    #It's a doc test
                    sys.stderr.write("%s docstring test ... " % name)
                    #Can't use fromlist=name.split(".") until python 2.5+
                    module = __import__(name, None, None, name.split("."))
                    suite = doctest.DocTestSuite(module)
                    del module
                suite.run(result)
                if result.wasSuccessful():
                    sys.stderr.write("ok\n")
                    return True
                else:
                    sys.stderr.write("FAIL\n")
                    result.printErrors()
                return False
            except MissingExternalDependencyError, msg:
                sys.stderr.write("skipping. %s\n" % msg)
                return True
            except Exception, msg:
                # This happened during the import
                sys.stderr.write("ERROR\n")
                result.stream.write(result.separator1+"\n")
                result.stream.write("ERROR: %s\n" % name)
                result.stream.write(result.separator2+"\n")
                result.stream.write(traceback.format_exc())
                return False
            except KeyboardInterrupt, err :
                # Want to allow this, and abort the test
                # (see below for special case)
                raise err
            except :
                # This happens in Jython with java.lang.ClassFormatError:
                # Invalid method Code length ...
                sys.stderr.write("ERROR\n")
                result.stream.write(result.separator1+"\n")
                result.stream.write("ERROR: %s\n" % name)
                result.stream.write(result.separator2+"\n")
                result.stream.write(traceback.format_exc())
                return False
        finally:
            sys.stdout = stdout

    def run(self):
        failures = 0
        startTime = time.time()
        for test in self.tests:
            ok = self.runTest(test)
            if not ok:
                failures += 1
        total = len(self.tests)
        stopTime = time.time()
        timeTaken = stopTime - startTime
        sys.stderr.write(self.stream.getvalue())
        sys.stderr.write('-' * 70 + "\n")
        sys.stderr.write("Ran %d test%s in %.3f seconds\n" %
                            (total, total != 1 and "s" or "", timeTaken))
        sys.stderr.write("\n")
        if failures:
            sys.stderr.write("FAILED (failures = %d)\n" % failures)


if __name__ == "__main__":
    #Don't do a sys.exit(...) as it isn't nice if run from IDLE.
    main(sys.argv[1:])
