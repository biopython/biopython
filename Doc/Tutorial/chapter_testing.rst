.. _`chapter:testing`:

The Biopython testing framework
===============================

Biopython has a regression testing framework (the file ``run_tests.py``)
based on `unittest <https://docs.python.org/3/library/unittest.html>`__,
the standard unit testing framework for Python. Providing comprehensive
tests for modules is one of the most important aspects of making sure
that the Biopython code is as bug-free as possible before going out. It
also tends to be one of the most undervalued aspects of contributing.
This chapter is designed to make running the Biopython tests and writing
good test code as easy as possible. Ideally, every module that goes into
Biopython should have a test (and should also have documentation!). All
our developers, and anyone installing Biopython from source, are
strongly encouraged to run the unit tests.

Running the tests
-----------------

When you download the Biopython source code, or check it out from our
source code repository, you should find a subdirectory call ``Tests``.
This contains the key script ``run_tests.py``, lots of individual
scripts named ``test_XXX.py``, and lots of other subdirectories which
contain input files for the test suite.

As part of building and installing Biopython you will typically run the
full test suite at the command line from the Biopython source top level
directory using the following:

.. code:: console

   $ cd Tests
   $ python run_tests.py

This is actually equivalent to going to the ``Tests`` subdirectory and
running:

.. code:: console

   $ python run_tests.py

You’ll often want to run just some of the tests, and this is done like
this:

.. code:: console

   $ python run_tests.py test_SeqIO.py test_AlignIO.py

When giving the list of tests, the ``.py`` extension is optional, so you
can also just type:

.. code:: console

   $ python run_tests.py test_SeqIO test_AlignIO

To run the docstring tests (see section :ref:`sec:doctest`), you can
use

.. code:: console

   $ python run_tests.py doctest

You can also skip any tests which have been setup with an explicit
online component by adding ``--offline``, e.g.

.. code:: console

   $ python run_tests.py --offline

By default, ``run_tests.py`` runs all tests, including the docstring
tests.

If an individual test is failing, you can also try running it directly,
which may give you more information.

Tests based on Python’s standard ``unittest`` framework will
``import unittest`` and then define ``unittest.TestCase`` classes, each
with one or more sub-tests as methods starting with ``test_`` which
check some specific aspect of the code.

Running the tests using Tox
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like most Python projects, you can also use
`Tox <https://tox.readthedocs.org/en/latest/>`__ to run the tests on
multiple Python versions, provided they are already installed in your
system.

We do not provide the configuration ``tox.ini`` file in our code base
because of difficulties pinning down user-specific settings (e.g.
executable names of the Python versions). You may also only be
interested in testing Biopython only against a subset of the Python
versions that we support.

If you are interested in using Tox, you could start with the example
``tox.ini`` shown below:

.. code:: text

   [tox]
   envlist = pypy,py38,py39

   [testenv]
   changedir = Tests
   commands = {envpython} run_tests.py --offline
   deps =
       numpy
       reportlab

Using the template above, executing ``tox`` will test your Biopython
code against PyPy, Python 3.8 and 3.9. It assumes that those Pythons’
executables are named “python3.8“ for Python 3.8, and so on.

Writing tests
-------------

Let’s say you want to write some tests for a module called ``Biospam``.
This can be a module you wrote, or an existing module that doesn’t have
any tests yet. In the examples below, we assume that ``Biospam`` is a
module that does simple math.

Each Biopython test consists of a script containing the test itself, and
optionally a directory with input files used by the test:

#. ``test_Biospam.py`` – The actual test code for your module.

#. ``Biospam`` [optional]– A directory where any necessary input files
   will be located. If you have any output files that should be manually
   reviewed, output them here (but this is discouraged) to prevent
   clogging up the main Tests directory. In general, use a temporary
   file/folder.

Any script with a ``test_`` prefix in the ``Tests`` directory will be
found and run by ``run_tests.py``. Below, we show an example test script
``test_Biospam.py``. If you put this script in the Biopython ``Tests``
directory, then ``run_tests.py`` will find it and execute the tests
contained in it:

.. code:: console

   $ python run_tests.py
   test_Ace ... ok
   test_AlignIO ... ok
   test_BioSQL ... ok
   test_BioSQL_SeqIO ... ok
   test_Biospam ... ok
   test_CAPS ... ok
   test_Clustalw ... ok
   ...
   ----------------------------------------------------------------------
   Ran 107 tests in 86.127 seconds

Writing a test using ``unittest``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``unittest``-framework has been included with Python since version
2.1, and is documented in the Python Library Reference (which I know you
are keeping under your pillow, as recommended). There is also `online
documentation for
unittest <https://docs.python.org/3/library/unittest.html>`__. If you
are familiar with the ``unittest`` system (or something similar like the
nose test framework), you shouldn’t have any trouble. You may find
looking at the existing examples within Biopython helpful too.

Here’s a minimal ``unittest``-style test script for ``Biospam``, which
you can copy and paste to get started:

.. code:: python

   import unittest
   from Bio import Biospam


   class BiospamTestAddition(unittest.TestCase):
       def test_addition1(self):
           result = Biospam.addition(2, 3)
           self.assertEqual(result, 5)

       def test_addition2(self):
           result = Biospam.addition(9, -1)
           self.assertEqual(result, 8)


   class BiospamTestDivision(unittest.TestCase):
       def test_division1(self):
           result = Biospam.division(3.0, 2.0)
           self.assertAlmostEqual(result, 1.5)

       def test_division2(self):
           result = Biospam.division(10.0, -2.0)
           self.assertAlmostEqual(result, -5.0)


   if __name__ == "__main__":
       runner = unittest.TextTestRunner(verbosity=2)
       unittest.main(testRunner=runner)

In the division tests, we use ``assertAlmostEqual`` instead of
``assertEqual`` to avoid tests failing due to roundoff errors; see the
``unittest`` chapter in the Python documentation for details and for
other functionality available in ``unittest`` (`online
reference <https://docs.python.org/3/library/unittest.html>`__).

These are the key points of ``unittest``-based tests:

-  Test cases are stored in classes that derive from
   ``unittest.TestCase`` and cover one basic aspect of your code

-  You can use methods ``setUp`` and ``tearDown`` for any repeated code
   which should be run before and after each test method. For example,
   the ``setUp`` method might be used to create an instance of the
   object you are testing, or open a file handle. The ``tearDown``
   should do any “tidying up”, for example closing the file handle.

-  The tests are prefixed with ``test_`` and each test should cover one
   specific part of what you are trying to test. You can have as many
   tests as you want in a class.

-  At the end of the test script, you can use

   .. code:: python

      if __name__ == "__main__":
          runner = unittest.TextTestRunner(verbosity=2)
          unittest.main(testRunner=runner)

   to execute the tests when the script is run by itself (rather than
   imported from ``run_tests.py``). If you run this script, then you’ll
   see something like the following:

   .. code:: console

      $ python test_BiospamMyModule.py
      test_addition1 (__main__.TestAddition) ... ok
      test_addition2 (__main__.TestAddition) ... ok
      test_division1 (__main__.TestDivision) ... ok
      test_division2 (__main__.TestDivision) ... ok

      ----------------------------------------------------------------------
      Ran 4 tests in 0.059s

      OK

-  To indicate more clearly what each test is doing, you can add
   docstrings to each test. These are shown when running the tests,
   which can be useful information if a test is failing.

   .. code:: python

      import unittest
      from Bio import Biospam


      class BiospamTestAddition(unittest.TestCase):
          def test_addition1(self):
              """An addition test"""
              result = Biospam.addition(2, 3)
              self.assertEqual(result, 5)

          def test_addition2(self):
              """A second addition test"""
              result = Biospam.addition(9, -1)
              self.assertEqual(result, 8)


      class BiospamTestDivision(unittest.TestCase):
          def test_division1(self):
              """Now let's check division"""
              result = Biospam.division(3.0, 2.0)
              self.assertAlmostEqual(result, 1.5)

          def test_division2(self):
              """A second division test"""
              result = Biospam.division(10.0, -2.0)
              self.assertAlmostEqual(result, -5.0)


      if __name__ == "__main__":
          runner = unittest.TextTestRunner(verbosity=2)
          unittest.main(testRunner=runner)

   Running the script will now show you:

   .. code:: console

      $ python test_BiospamMyModule.py
      An addition test ... ok
      A second addition test ... ok
      Now let's check division ... ok
      A second division test ... ok

      ----------------------------------------------------------------------
      Ran 4 tests in 0.001s

      OK

If your module contains docstring tests (see section
:ref:`sec:doctest`), you *may* want to include those in the tests to
be run. You can do so as follows by modifying the code under
``if __name__ == "__main__":`` to look like this:

.. code:: python

   if __name__ == "__main__":
       unittest_suite = unittest.TestLoader().loadTestsFromName("test_Biospam")
       doctest_suite = doctest.DocTestSuite(Biospam)
       suite = unittest.TestSuite((unittest_suite, doctest_suite))
       runner = unittest.TextTestRunner(sys.stdout, verbosity=2)
       runner.run(suite)

This is only relevant if you want to run the docstring tests when you
execute ``python test_Biospam.py`` if it has some complex run-time
dependency checking.

In general instead include the docstring tests by adding them to the
``run_tests.py`` as explained below.

.. _`sec:doctest`:

Writing doctests
----------------

Python modules, classes and functions support built-in documentation
using docstrings. The `doctest
framework <https://docs.python.org/3/library/doctest.html>`__ (included
with Python) allows the developer to embed working examples in the
docstrings, and have these examples automatically tested.

Currently only part of Biopython includes doctests. The ``run_tests.py``
script takes care of running the doctests. For this purpose, at the top
of the ``run_tests.py`` script is a manually compiled list of modules to
skip, important where optional external dependencies which may not be
installed (e.g. the Reportlab and NumPy libraries). So, if you’ve added
some doctests to the docstrings in a Biopython module, in order to have
them excluded in the Biopython test suite, you must update
``run_tests.py`` to include your module. Currently, the relevant part of
``run_tests.py`` looks as follows:

.. code:: python

   # Following modules have historic failures. If you fix one of these
   # please remove here!
   EXCLUDE_DOCTEST_MODULES = [
       "Bio.PDB",
       "Bio.PDB.AbstractPropertyMap",
       "Bio.Phylo.Applications._Fasttree",
       "Bio.Phylo._io",
       "Bio.Phylo.TreeConstruction",
       "Bio.Phylo._utils",
   ]

   # Exclude modules with online activity
   # They are not excluded by default, use --offline to exclude them
   ONLINE_DOCTEST_MODULES = ["Bio.Entrez", "Bio.ExPASy", "Bio.TogoWS"]

   # Silently ignore any doctests for modules requiring numpy!
   if numpy is None:
       EXCLUDE_DOCTEST_MODULES.extend(
           [
               "Bio.Affy.CelFile",
               "Bio.Cluster",
               # ...
           ]
       )

Note that we regard doctests primarily as documentation, so you should
stick to typical usage. Generally complicated examples dealing with
error conditions and the like would be best left to a dedicated unit
test.

Note that if you want to write doctests involving file parsing, defining
the file location complicates matters. Ideally use relative paths
assuming the code will be run from the ``Tests`` directory, see the
``Bio.SeqIO`` doctests for an example of this.

To run the docstring tests only, use

.. code:: console

   $ python run_tests.py doctest

Note that the doctest system is fragile and care is needed to ensure
your output will match on all the different versions of Python that
Biopython supports (e.g. differences in floating point numbers).

.. _`sec:doctest-tutorial`:

Writing doctests in the Tutorial
--------------------------------

This Tutorial you are reading has a lot of code snippets, which are
often formatted like a doctest. We have our own system in file
``test_Tutorial.py`` to allow tagging code snippets in the Tutorial
source to be run as Python doctests. This works by adding special
``.. doctest`` comment lines before each Python Console (pycon) block,
e.g.

.. code:: rst

   .. doctest

   .. code:: pycon

      >>> from Bio.Seq import Seq
      >>> s = Seq("ACGT")
      >>> len(s)
      4

Often code examples are not self-contained, but continue from the
previous Python block. Here we use the magic comment ``.. cont-doctest``
as shown here:

.. code:: rst

   .. cont-doctest

   .. code:: pycon

      >>> s == "ACGT"
      True

The special ``.. doctest`` comment line can take a working directory
(relative to the ``Doc/`` folder) to use if you have any example data
files, e.g. ``.. doctest examples`` will use the ``Doc/examples`` folder,
while ``.. doctest ../Tests/GenBank`` will use the ``Tests/GenBank``
folder.

After the directory argument, you can specify any Python dependencies
which must be present in order to run the test by adding ``lib:XXX`` to
indicate ``import XXX`` must work, e.g. ``.. doctest examples lib:numpy``

You can run the Tutorial doctests via:

.. code:: console

   $ python test_Tutorial.py

or:

.. code:: console

   $ python run_tests.py test_Tutorial.py
