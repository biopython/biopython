#!/usr/bin/env python
"""Run a set of PyUnit-based regression tests.

This will find all modules whose name is "test_*.py" in the test
directory, and run them.  Various command line options provide
additional facilities.

Command line options:

-g;--generate -- write the output file for a test instead of comparing it.
                 A test to write the output for must be specified.
--no-gui      -- Obsolete option (there used to be a GUI for the tests).
--help        -- show usage info
<test_name>   -- supply the name of one (or more) tests to be run.
                 The .py file extension is optional.
"""
# standard modules
import sys
import cStringIO
import os
import re
import sys
import getopt

# PyUnit
import unittest

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
        opts, args = getopt.getopt(argv, 'g',
				   ["generate", "no-gui", "help"])
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

            generate_output(args[0])
            return 0

    # deal with the arguments, which should be names of tests to run
    for arg_num in range(len(args)):
        # strip off the .py if it was included
        if args[arg_num][-3:] == ".py":
            args[arg_num] = args[arg_num][:-3]

    # run the tests
    test_suite = testing_suite(args)
    runner = unittest.TextTestRunner(verbosity = 2)
    runner.run(test_suite)

def testing_suite(tests_to_run = []):
    all_tests = tests_to_run
    # if no tests were specified to run, we run them all
    if len(all_tests) == 0:
        all_tests = findtests()

    test_suite = unittest.TestSuite()

    for test in all_tests:
        # add the test to the test suite
        test_suite.addTest(RegressionTest(test))

    return test_suite

class RegressionTest(unittest.TestCase):
    """Run a unit test and compare its output against expected output.
    """
    def __init__(self, test_name):
        """Initialize with the test to run.

        Arguments:
        o test_name - The name of the test. This should be able to be
        imported so that '__import__(test_name)' will run the test.
        """
        unittest.TestCase.__init__(self)
        self.test_name = test_name

    def __str__(self):
        return self.shortDescription()

    def shortDescription(self):
        return self.test_name

    def runTest(self):
        """Run the actual test inside a try/except to catch import errors.
        """
        from Bio import MissingExternalDependencyError
        try:
            self.runSafeTest()
        except MissingExternalDependencyError, msg:
            print "skipping.", msg

    def runSafeTest(self):
        generated_output = ''
        output = cStringIO.StringIO()

        # remember standard out so we can reset it after we are done
        save_stdout = sys.stdout
        try:
            # write the output from the test into a string
            sys.stdout = output
            cur_test = __import__(self.test_name)

            # some tests may be runnable by run_tests()
            try:        
                cur_test.run_tests([])
            except AttributeError:
                pass

            generated_output = output.getvalue()
        finally:
            # return standard out to its normal setting
            sys.stdout = save_stdout

        # get the expected output
        testdir = findtestdir()
        outputdir = os.path.join(testdir, "output")
        outputfile = os.path.join(outputdir, self.test_name)

        try:
            expected_handle = open(outputfile, 'r')
            output_handle = cStringIO.StringIO(generated_output)
            # check the expected output to be consistent with what
            # we generated
            compare_output(self.test_name, output_handle,
                           expected_handle)
        except IOError:
            raise IOError, "Warning: Can't open %s for test %s" % \
                  (outputfile, self.test_name)

def findtests():
    """Return a list of all applicable test modules."""
    testdir = findtestdir()
    names = os.listdir(testdir)
    tests = []
    for name in names:
        if name[:5] == "test_" and name[-3:] == ".py":
	    tests.append(name[:-3])
    tests.sort()
    return tests

def findtestdir():
    if __name__ == '__main__':
        file = sys.argv[0]
    else:
        file = __file__
    testdir = os.path.dirname(file) or os.curdir
    return testdir

def generate_output(test_name):
    """Generate the golden output for the specified test.
    """
    testdir = findtestdir()
    outputdir = os.path.join(testdir, "output")
    outputfile = os.path.join(outputdir, test_name)

    output_handle = open(outputfile, 'w')

    # write the test name as the first line of the output
    output_handle.write(test_name + "\n")

    # remember standard out so we can reset it after we are done
    save_stdout = sys.stdout
    try:
        # write the output from the test into a string
        sys.stdout = output_handle
        cur_test = __import__(test_name)

        # run tests that use run_tests() to signify running
        try:
            cur_test.run_tests([])
        except AttributeError:
            pass
    finally:
        output_handle.close()
        # return standard out to its normal setting
        sys.stdout = save_stdout

def compare_output(test_name, output_handle, expected_handle):
    """Compare output from a test to the expected output.

    Arguments:

    o test_name - The name of the test we are running.

    o output_handle - A handle to all of the output generated by running
    a test.

    o expected_handle - A handle to the expected output for a test.
    """
    # first check that we are dealing with the right output
    # the first line of the output file is the test name
    expected_test = expected_handle.readline().strip()

    assert expected_test == test_name, "\nOutput:   %s\nExpected: %s" % \
           (test_name, expected_test)

    # now loop through the output and compare it to the expected file
    while 1:
        expected_line = expected_handle.readline()
        output_line = output_handle.readline()

        # stop looping if either of the info handles reach the end
        if not(expected_line) or not(output_line):
            # make sure both have no information left
            assert expected_line == '', "Unread: %s" % expected_line
            assert output_line == '', "Extra output: %s" % output_line
            break

        # normalize the newlines in the two lines
        expected_line = convert_newlines(expected_line)
        output_line = convert_newlines(output_line)

        # normalize the internal newlines inside strings
        expected_line = convert_string_newlines(expected_line)
        output_line = convert_string_newlines(output_line)

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

def convert_newlines(line):
    """Convert all end of line characters to '\n'

    Ths deals with differences in newline conventions between different
    platforms by converting all newlines to the simple \n.
    """
    newlines = ["\r\n", "\r"]

    for newline_to_replace in newlines:
        line = line.replace(newline_to_replace, "\n")

    return line
        
def convert_string_newlines(line):
    """Convert all newlines inside strings in the given line into '\\n'.

    This helps deal with the problem between python2.1 and older pythons
    in how line breaks are generated inside strings. Older python versions
    used '\012', and 2.1 uses '\n'.
    """
    # two slashes are used since we are dealing with newlines inside strings
    # where they are "escaped" with the extra '\' 
    newlines = ["\\012"]
        
    for newline_to_replace in newlines:
        line = line.replace(newline_to_replace, "\\n")

    return line
        
if __name__ == "__main__":
    #Don't do a sys.exit(...) as it isn't nice if run from IDLE.
    main(sys.argv[1:])

