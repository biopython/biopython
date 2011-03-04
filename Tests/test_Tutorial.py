import unittest
import doctest
import os
import sys

if sys.version_info[0] >= 3:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "This test doesn't work on Python 3 yet (need to call 2to3).")

tutorial = os.path.join(os.path.dirname(sys.argv[0]), "../Doc/Tutorial.tex")

def _extract(handle):
    line = handle.readline()
    if line != "\\begin{verbatim}\n":
        raise ValueError("Any '%doctest' or '%cont-doctest' line should be followed by '\\begin{verbatim}'")
    lines = []
    while True:
        line = handle.readline()
        if not line:
            if lines:
                print "".join(lines[:30])
                raise ValueError("Didn't find end of test starting: %r", lines[0])
            else:
                raise ValueError("Didn't find end of test!")
        elif line.startswith("\end{verbatim}"):
            break
        else:
            lines.append(line)
    return lines
    
def extract_doctests(latex_filename):
    """Scans LaTeX file and pulls out marked doctests as strings."""
    handle = open(latex_filename, "rU")
    in_test = False
    count = 0
    lines = []
    while True:
        line = handle.readline()
        if not line:
            #End of file
            break
        elif line.startswith("%cont-doctest"):
            lines.extend(_extract(handle))
        elif line.startswith("%doctest"):
            if lines:
                if not lines[0].startswith(">>> "):
                    raise ValueError("Should start '>>> ' not %r" % lines[0])
                yield name, "".join(lines)
                count += 1
                lines = []
            try:
                name = line.split(None,1)[1].strip()
            except:
                name = "unnamed_test_%i" % count
            lines = _extract(handle)
    handle.close()
    if lines:
        if not lines[0].startswith(">>> "):
            raise ValueError("Should start '>>> ' not %r" % lines[0])
        yield name, "".join(lines)
        count += 1
    #yield "dummy", ">>> 2 + 2\n5\n"

class TutorialDocTestHolder(object):
    """Python doctests extracted from the Biopython Tutorial."""
    pass


#Create dummy methods on the object purely to hold doctests
for name, example in extract_doctests(tutorial):
    #print name, repr(doctest)
    def funct(n, d):
        method = lambda x : None
        method.__doc__ = "%s\n\n%s\n" % (n, d)
        return method
    setattr(TutorialDocTestHolder,
            "doctest_%s" % name.replace(" ","_"),
            funct(name, example))
    del funct


#This is a TestCase class so it is found by run_tests.py
class TutorialTestCase(unittest.TestCase):
    """Python doctests extracted from the Biopython Tutorial."""
    #Single method to be invoked by run_tests.py
    def test_doctests(self):
        """Run tutorial doctests."""
        runner = doctest.DocTestRunner()
        failures = []
        for test in doctest.DocTestFinder().find(TutorialDocTestHolder):
            failed, success = runner.run(test)
            if failed:
                name = test.name
                assert name.startswith("TutorialDocTestHolder.doctest_")
                failures.append(name[30:])
                #raise ValueError("Tutorial doctest %s failed" % test.name[30:])
        if failures:
            raise ValueError("%i Tutorial doctests failed: %s" % \
                             (len(failures), ", ".join(failures)))


#This is to run the doctests if the script is called directly:
if __name__ == "__main__":
    doctest.testmod()
