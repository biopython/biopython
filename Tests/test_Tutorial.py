import unittest
import doctest
import os
import sys

if sys.version_info[0] >= 3:
    from lib2to3 import refactor
    rt = refactor.RefactoringTool(refactor.get_fixers_from_package("lib2to3.fixes"))
    assert rt.refactor_docstring(">>> print 2+2\n4\n", "example") == \
           ">>> print(2+2)\n4\n"
    
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
    line_number = 0
    in_test = False
    lines = []
    while True:
        line = handle.readline()
        line_number += 1
        if not line:
            #End of file
            break
        elif line.startswith("%cont-doctest"):
            x = _extract(handle)
            lines.extend(x)
            line_number += len(x) + 2
        elif line.startswith("%doctest"):
            if lines:
                if not lines[0].startswith(">>> "):
                    raise ValueError("Should start '>>> ' not %r" % lines[0])
                yield name, "".join(lines)
                lines = []
            try:
                name = line.split(None,1)[1].strip()
                assert not name
            except:
                pass
            name = "test_from_line_%05i" % line_number
            x = _extract(handle)
            lines.extend(x)
            line_number += len(x) + 2
    handle.close()
    if lines:
        if not lines[0].startswith(">>> "):
            raise ValueError("Should start '>>> ' not %r" % lines[0])
        yield name, "".join(lines)
    #yield "dummy", ">>> 2 + 2\n5\n"

class TutorialDocTestHolder(object):
    """Python doctests extracted from the Biopython Tutorial."""
    pass


#Create dummy methods on the object purely to hold doctests
for name, example in extract_doctests(tutorial):
    if sys.version_info[0] >= 3:
        example = rt.refactor_docstring(example, name)
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
    print "Runing Tutorial doctests..."
    import doctest
    tests = doctest.testmod()
    if tests[0]:
        #Note on Python 2.5+ can use tests.failed rather than tests[0]
        raise RuntimeError("%i/%i tests failed" % tests)
    print "Tests done"
