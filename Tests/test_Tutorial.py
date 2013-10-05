# Copyright 2011-2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This will apply to all the doctests too:
from __future__ import print_function

import unittest
import doctest
import os
import sys
import warnings
from Bio import BiopythonExperimentalWarning

warnings.simplefilter('ignore', BiopythonExperimentalWarning)

if sys.version_info[0] >= 3:
    from lib2to3 import refactor
    fixers = refactor.get_fixers_from_package("lib2to3.fixes")
    fixers.remove("lib2to3.fixes.fix_print") # Already using print function
    rt = refactor.RefactoringTool(fixers)
    assert rt.refactor_docstring(">>> print(2+2)\n4\n", "example1") == \
                                 ">>> print(2+2)\n4\n"
    assert rt.refactor_docstring('>>> print("Two plus two is", 2+2)\n'
                                 'Two plus two is 4\n', "example2") == \
                                 '>>> print("Two plus two is", 2+2)\nTwo plus two is 4\n'


tutorial = os.path.join(os.path.dirname(sys.argv[0]), "../Doc/Tutorial.tex")
if not os.path.isfile(tutorial) and sys.version_info[0] >= 3:
    tutorial = os.path.join(os.path.dirname(sys.argv[0]), "../../../Doc/Tutorial.tex")
if not os.path.isfile(tutorial):
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError("Could not find ../Doc/Tutorial.tex file")

tutorial_base = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "../Doc/"))
original_path = os.path.abspath(".")


def _extract(handle):
    line = handle.readline()
    if line != "\\begin{verbatim}\n":
        raise ValueError("Any '%doctest' or '%cont-doctest' line should be followed by '\\begin{verbatim}'")
    lines = []
    while True:
        line = handle.readline()
        if not line:
            if lines:
                print("".join(lines[:30]))
                raise ValueError("Didn't find end of test starting: %r", lines[0])
            else:
                raise ValueError("Didn't find end of test!")
        elif line.startswith("\end{verbatim}"):
            break
        else:
            lines.append(line)
    return lines


def extract_doctests(latex_filename):
    """Scans LaTeX file and pulls out marked doctests as strings.

    This is a generator, yielding one tuple per doctest.
    """
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
                yield name, "".join(lines), folder, deps
                lines = []
            deps = [x.strip() for x in line.split()[1:]]
            if deps:
                folder = deps[0]
                deps = deps[1:]
            else:
                folder = ""
            name = "test_from_line_%05i" % line_number
            x = _extract(handle)
            lines.extend(x)
            line_number += len(x) + 2
    handle.close()
    if lines:
        if not lines[0].startswith(">>> "):
            raise ValueError("Should start '>>> ' not %r" % lines[0])
        yield name, "".join(lines), folder, deps
    #yield "dummy", ">>> 2 + 2\n5\n"


class TutorialDocTestHolder(object):
    """Python doctests extracted from the Biopython Tutorial."""
    pass

def check_deps(dependencies):
    missing = []
    for dep in dependencies:
        assert dep.startswith("lib:"), dep
        lib = dep[4:]
        try:
            tmp = __import__(lib)
            del tmp
        except ImportError:
            missing.append(lib)
    return missing

#Create dummy methods on the object purely to hold doctests
missing_deps = set()
for name, example, folder, deps in extract_doctests(tutorial):
    missing = check_deps(deps)
    if missing:
        missing_deps.update(missing)
        continue

    if sys.version_info[0] >= 3:
        example = ">>> from __future__ import print_function\n" + example
        example = rt.refactor_docstring(example, name)

    def funct(n, d, f):
        global tutorial_base
        method = lambda x : None
        if f:
            p = os.path.join(tutorial_base, f)
            method.__doc__ = "%s\n\n>>> import os\n>>> os.chdir(%r)\n%s\n" \
                           % (n, p, d)
        else:
            method.__doc__ = "%s\n\n%s\n" % (n, d)
        method._folder = f
        return method

    setattr(TutorialDocTestHolder,
            "doctest_%s" % name.replace(" ", "_"),
            funct(name, example, folder))
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
            raise ValueError("%i Tutorial doctests failed: %s" %
                             (len(failures), ", ".join(failures)))

    def tearDown(self):
        global original_path
        os.chdir(original_path)


#This is to run the doctests if the script is called directly:
if __name__ == "__main__":
    if missing_deps:
        print("Skipping tests needing the following:")
        for dep in sorted(missing_deps):
            print(" - %s" % dep)
    print("Running Tutorial doctests...")
    import doctest
    tests = doctest.testmod()
    if tests[0]:
        #Note on Python 2.5+ can use tests.failed rather than tests[0]
        raise RuntimeError("%i/%i tests failed" % tests)
    print("Tests done")
