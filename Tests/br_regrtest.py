#! /usr/bin/env python

# Regression test code based on the standard Python 1.5 regr_test,
# (add Python copyright statement here) and with additions made by
# Bioreason Inc. (add Bioreason copyright here, and under LGPL).  More
# additions by Andrew Dalke.

"""Regression test.

This will find all modules whose name is "test_*" in the test
directory, and run them.  Various command line options provide
additional facilities.

Command line options:

-v: verbose -- run tests in verbose mode with output to stdout
-q: quiet -- don't print anything except if a test fails
-g: generate -- write the output file for a test instead of comparing it
-x: exclude -- arguments are tests to *exclude*
-p <name>=<directory>: treat the directory as the source for a package of the
    given name

If non-option arguments are present, they are names for tests to run,
unless -x is given, in which case they are names for tests not to run.
If no test names are given, all tests are run.

-v is incompatible with -g and does not compare test output files.
"""

import sys, string, os
import getopt, traceback, imp
try:
    import test_support
except ImportError:
    # looks like Python 1.5.2 changed the location on me
    test_support = __import__("test/test_support")

SUCCEEDED = 1
FAILED = 0
ERROR = -1

class TestResult:
    """This class stores the results of a test, including any error msgs."""
    def __init__(self, status=1, info=None):
        """Initialize a new instance.
        `status' indicates the completion status of the test:
           SUCCEEDED == 1 => test succeeded
           FAILED == 0 => test failed
           ERROR == -1 => test could not be executed, usually due to import errors

        `info' should be an object containing extra information about the
        test result, e.g. an error message string."""

        self.status = status
        self.info = info


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'vgqxp:',
				   ["verbose", "quiet", "generate",
				    "exclude", "package=", "help",
				    "version"])
    except getopt.error, msg:
        print msg
        print __doc__
        return 2
    verbose = 0
    quiet = 0
    generate = 0
    exclude = 0
    for o, a in opts:
        if o == '-v' or o == "--verbose": verbose = verbose+1
        if o == '-q' or o == "--quiet": quiet = 1; verbose = 0
        if o == '-g' or o == "--generate": generate = 1
        if o == '-x' or o == "--exclude": exclude = 1
	if o == '-p' or o == "--package":
	    pos = string.find(a, "=")
	    if pos == -1:
		print "package option `-p' is missing `='"
		return 2
	    module_name = a[:pos]
	    dir_name = a[pos+1:]
	    imp.load_module(module_name, None, dir_name,
			    ("", "r", imp.PKG_DIRECTORY) )
	if o == "--help":
	    print __doc__
	    return 0
	if o == "--version":
	    print "br_regrtest 1.0"
	    return 0
		
		
    if generate and verbose:
        print "-g and -v don't go together!"
        return 2
    good = []
    bad = []
    skipped = []
    for i in range(len(args)):
        # Strip trailing ".py" from arguments
        if args[i][-3:] == '.py':
            args[i] = args[i][:-3]
    if exclude:
        nottests[:0] = args
        args = []
    tests = args or findtests()
    test_support.verbose = verbose      # Tell tests to be moderately quiet
    for test in tests:
        if not quiet:
            print test
        result = runtest(test, generate, verbose)
        if result.status == SUCCEEDED:
            good.append(test)
        elif result.status == FAILED:
            bad.append(test)
        else:
            if not quiet:
                print "test", test, "skipped:", result.info
            skipped.append(test)
    if good and not quiet:
        if not bad and not skipped and len(good) > 1:
            print "All",
        print count(len(good), "test"), "OK."
    if bad:
        print count(len(bad), "test"), "failed:",
        print string.join(bad)
    if skipped and not quiet:
        print count(len(skipped), "test"), "skipped:",
        print string.join(skipped)
    return len(bad) > 0

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

def runtest(test, generate, verbose):
    test_support.unload(test)
    testdir = findtestdir()
    outputdir = os.path.join(testdir, "output")
    outputfile = os.path.join(outputdir, test)
    try:
        if generate:
            cfp = open(outputfile, "w")
        elif verbose:
            cfp = sys.stdout
        else:
            cfp = Compare(outputfile)
    except IOError:
        cfp = None
        print "Warning: can't open", outputfile

    try:
        save_stdout = sys.stdout
        try:
            if cfp:
                sys.stdout = cfp
                print test              # Output file starts with test name
            __import__(test, globals(), locals(), [])
	    if cfp and not generate and not verbose:
		cfp.close()
        finally:
            sys.stdout = save_stdout
    except ImportError, msg:
        return TestResult(ERROR, msg)
    except KeyboardInterrupt, v:
        raise KeyboardInterrupt, v, sys.exc_info()[2]
    except test_support.TestFailed, msg:
        print "test", test, "failed --", msg
        return TestResult(FAILED, msg)
    except:
        type, value = sys.exc_info()[:2]
        print "test", test, "crashed --", type, ":", value
        if verbose:
            traceback.print_exc(file=sys.stdout)
        return TestResult(FAILED, (type, value))
    else:
        return TestResult(SUCCEEDED)

def findtestdir():
    if __name__ == '__main__':
        file = sys.argv[0]
    else:
        file = __file__
    testdir = os.path.dirname(file) or os.curdir
    return testdir

# fancy things up to say "5 tests failed" or "1 test passed"
def count(n, word):
    if n == 1:
        return "%d %s" % (n, word)
    else:
        return "%d %ss" % (n, word)

class Compare:

    def __init__(self, filename):
        self.fp = open(filename, 'r')

    def write(self, data):
        expected = self.fp.read(len(data))
        if data <> expected:
            raise test_support.TestFailed, \
                    'Writing: '+`data`+', expected: '+`expected`

    def writelines(self, listoflines):
        map(self.write, listoflines)

    def flush(self):
        pass

    def close(self):
        leftover = self.fp.read()
        if leftover:
            raise test_support.TestFailed, 'Unread: '+`leftover`
        self.fp.close()

    def isatty(self):
        return 0

if __name__ == '__main__':
    sys.exit(main())
