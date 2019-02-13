# Copyright 2011-2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This script looks for entries in the LaTeX source for the
# Biopython Tutorial which can be turned into Python doctests,
# e.g.
#
# %doctest
# \begin{verbatim}
# >>> from Bio.Alphabet import generic_dna
# >>> from Bio.Seq import Seq
# >>> len("ACGT")
# 4
# \end{verbatim}
#
# Code snippets can be extended using a similar syntax, which
# will create a single combined doctest:
#
# %cont-doctest
# \begin{verbatim}
# >>> Seq("ACGT") == Seq("ACGT", generic_dna)
# True
# \end{verbatim}
#
# The %doctest line also supports a relative working directory,
# and listing multiple Python dependencies as lib:XXX which will
# ensure "import XXX" works before using the test. e.g.
#
# %doctest examples lib:numpy lib:scipy
#
# Additionally after the path, special keyword 'internet' is
# used to flag online tests.
#
# Note if using lib:XXX or special value 'internet' you must
# include a relative path to the working directory, just use '.'
# for the default path, e.g.
#
# %doctest . lib:reportlab
#
# %doctest . internet
#
# TODO: Adding bin:XXX for checking binary XXX is on $PATH?
#
# See also "Writing doctests in the Tutorial" in the Tutorial
# itself.


# This future import will apply to all the doctests too:
from __future__ import print_function
from Bio._py3k import _universal_read_mode

import unittest
import doctest
import os
import sys
import warnings
from Bio import BiopythonExperimentalWarning, MissingExternalDependencyError

# This is the same mechanism used for run_tests.py --offline
# to skip tests requiring the network.
import requires_internet
try:
    requires_internet.check()
    online = True
except MissingExternalDependencyError:
    online = False
if "--offline" in sys.argv:
    # Allow manual override via "python test_Tutorial.py --offline"
    online = False

warnings.simplefilter('ignore', BiopythonExperimentalWarning)

if sys.version_info[0] >= 3:
    from lib2to3 import refactor
    from lib2to3.pgen2.tokenize import TokenError
    fixers = refactor.get_fixers_from_package("lib2to3.fixes")
    fixers.remove("lib2to3.fixes.fix_print")  # Already using print function
    rt = refactor.RefactoringTool(fixers)
    assert rt.refactor_docstring(
        ">>> print(2+2)\n4\n", "example1") == \
        ">>> print(2+2)\n4\n"
    assert rt.refactor_docstring(
        '>>> print("Two plus two is", 2+2)\n'
        'Two plus two is 4\n', "example2") == \
        '>>> print("Two plus two is", 2+2)\nTwo plus two is 4\n'

# Cache this to restore the cwd at the end of the tests
original_path = os.path.abspath(".")

if os.path.basename(sys.argv[0]) == "test_Tutorial.py":
    # sys.argv[0] will be (relative) path to test_Turorial.py - use this to allow, e.g.
    # [base]$ python Tests/test_Tutorial.py
    # [Tests/]$ python test_Tutorial.py
    tutorial_base = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "../Doc/"))
    tutorial = os.path.join(tutorial_base, "Tutorial.tex")
else:
    # Probably called via run_tests.py so current directory should (now) be Tests/
    # but may have been changed by run_tests.py so can't infer from sys.argv[0] with e.g.
    # [base]$ python Tests/run_tests.py test_Tutorial
    tutorial_base = os.path.abspath("../Doc/")
    tutorial = os.path.join(tutorial_base, "Tutorial.tex")
if not os.path.isfile(tutorial):
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError("Could not find ../Doc/Tutorial.tex file")

# Build a list of all the Tutorial LaTeX files:
files = [tutorial]
for latex in os.listdir(os.path.join(tutorial_base, "Tutorial/")):
    if latex.startswith("chapter_") and latex.endswith(".tex"):
        files.append(os.path.join(tutorial_base, "Tutorial", latex))


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
        elif line.startswith("\\end{verbatim}"):
            break
        else:
            lines.append(line)
    return lines


def extract_doctests(latex_filename):
    """Scan LaTeX file and pull out marked doctests as strings.

    This is a generator, yielding one tuple per doctest.
    """
    base_name = os.path.splitext(os.path.basename(latex_filename))[0]
    deps = ""
    folder = ""
    with open(latex_filename, _universal_read_mode) as handle:
        line_number = 0
        lines = []
        name = None
        while True:
            line = handle.readline()
            line_number += 1
            if not line:
                # End of file
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
                name = "test_%s_line_%05i" % (base_name, line_number)
                x = _extract(handle)
                lines.extend(x)
                line_number += len(x) + 2
    if lines:
        if not lines[0].startswith(">>> "):
            raise ValueError("Should start '>>> ' not %r" % lines[0])
        yield name, "".join(lines), folder, deps
    # yield "dummy", ">>> 2 + 2\n5\n"


class TutorialDocTestHolder(object):
    """Python doctests extracted from the Biopython Tutorial."""

    pass


def check_deps(dependencies):
    """Check 'lib:XXX' and 'internet' dependencies are met."""
    missing = []
    for dep in dependencies:
        if dep == "internet":
            if not online:
                missing.append("internet")
        else:
            assert dep.startswith("lib:"), dep
            lib = dep[4:]
            try:
                tmp = __import__(lib)
                del tmp
            except ImportError:
                missing.append(lib)
    return missing


# Create dummy methods on the object purely to hold doctests
missing_deps = set()
for latex in files:
    # print("Extracting doctests from %s" % latex)
    for name, example, folder, deps in extract_doctests(latex):
        missing = check_deps(deps)
        if missing:
            missing_deps.update(missing)
            continue

        if sys.version_info[0] >= 3:
            example = ">>> from __future__ import print_function\n" + example
            try:
                example = rt.refactor_docstring(example, name)
            except TokenError:
                raise ValueError("Problem with %s:\n%s" % (name, example))

        def funct(n, d, f):
            global tutorial_base
            method = lambda x: None
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


# This is a TestCase class so it is found by run_tests.py
class TutorialTestCase(unittest.TestCase):
    """Python doctests extracted from the Biopython Tutorial."""

    # Single method to be invoked by run_tests.py
    def test_doctests(self):
        """Run tutorial doctests."""
        # TODO: Remove IGNORE_EXCEPTION_DETAIL once drop Python 2 support
        runner = doctest.DocTestRunner(optionflags=doctest.IGNORE_EXCEPTION_DETAIL)
        failures = []
        for test in doctest.DocTestFinder().find(TutorialDocTestHolder):
            failed, success = runner.run(test)
            if failed:
                name = test.name
                assert name.startswith("TutorialDocTestHolder.doctest_")
                failures.append(name[30:])
                # raise ValueError("Tutorial doctest %s failed" % test.name[30:])
        if failures:
            raise ValueError("%i Tutorial doctests failed: %s" %
                             (len(failures), ", ".join(failures)))

    def tearDown(self):
        global original_path
        os.chdir(original_path)
        # files currently don't get created during test with python3.5 and pypy
        # remove files created from chapter_phylo.tex
        delete_phylo_tutorial = [
            "examples/tree1.nwk", "examples/other_trees.nwk"
            ]
        for file in delete_phylo_tutorial:
            if os.path.exists(os.path.join(tutorial_base, file)):
                os.remove(os.path.join(tutorial_base, file))
        # remove files created from chapter_cluster.tex
        tutorial_cluster_base = os.path.abspath("../Tests/")
        delete_cluster_tutorial = [
            "Cluster/cyano_result.atr", "Cluster/cyano_result.cdt",
            "Cluster/cyano_result.gtr", "Cluster/cyano_result_K_A2.kag",
            "Cluster/cyano_result_K_G5.kgg", "Cluster/cyano_result_K_G5_A2.cdt"
            ]
        for file in delete_cluster_tutorial:
            if os.path.exists(os.path.join(tutorial_cluster_base, file)):
                os.remove(os.path.join(tutorial_cluster_base, file))


# This is to run the doctests if the script is called directly:
if __name__ == "__main__":
    if missing_deps:
        print("Skipping tests needing the following:")
        for dep in sorted(missing_deps):
            print(" - %s" % dep)
    print("Running Tutorial doctests...")
    tests = doctest.testmod()
    if tests.failed:
        raise RuntimeError("%i/%i tests failed" % tests)
    print("Tests done")
