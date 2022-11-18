.. image:: https://img.shields.io/pypi/v/biopython.svg?logo=pypi
   :alt: Biopython on the Python Package Index (PyPI)
   :target: https://pypi.python.org/pypi/biopython
.. image:: https://img.shields.io/conda/vn/conda-forge/biopython.svg?logo=conda-forge
   :alt: Biopython on the Conda package conda-forge channel
   :target: https://anaconda.org/conda-forge/biopython
.. image:: https://results.pre-commit.ci/badge/github/biopython/biopython/master.svg
   :target: https://results.pre-commit.ci/latest/github/biopython/biopython/master
   :alt: pre-commit.ci status
.. image:: https://img.shields.io/circleci/build/github/biopython/biopython.svg?logo=circleci
   :alt: Linux testing with CircleCI
   :target: https://app.circleci.com/pipelines/github/biopython/biopython
.. image:: https://img.shields.io/appveyor/ci/biopython/biopython/master.svg?logo=appveyor
   :alt: Windows testing with AppVeyor
   :target: https://ci.appveyor.com/project/biopython/biopython/history
.. image:: https://img.shields.io/github/workflow/status/biopython/biopython/Basic%20Checks?logo=github-actions
   :alt: GitHub workflow status
   :target: https://github.com/biopython/biopython/actions
.. image:: https://img.shields.io/codecov/c/github/biopython/biopython/master.svg?logo=codecov
   :alt: Test coverage on CodeCov
   :target: https://codecov.io/github/biopython/biopython/
.. image:: http://depsy.org/api/package/pypi/biopython/badge.svg
   :alt: Research software impact on Depsy
   :target: http://depsy.org/package/python/biopython

.. image:: https://github.com/biopython/biopython/raw/master/Doc/images/biopython_logo_m.png
   :alt: The Biopython Project
   :target: http://biopython.org

Biopython README file
=====================

The Biopython Project is an international association of developers of freely
available Python tools for computational molecular biology.

Our user-centric documentation is hosted on https://biopython.org including
our `API Documentation <https://biopython.org/docs/latest/api/>`_ and the main
`Biopython Tutorial and Cookbook
<http://biopython.org/DIST/docs/tutorial/Tutorial.html>`_
(`PDF <http://biopython.org/DIST/docs/tutorial/Tutorial.pdf>`_).

This README file is intended primarily for people interested in working
with the Biopython source code, either one of the releases from the
http://biopython.org website, or from our repository on GitHub
https://github.com/biopython/biopython

The `NEWS <https://github.com/biopython/biopython/blob/master/NEWS.rst>`_
file summarises the changes in each release of Biopython.

The Biopython package is open source software made available under generous
terms. Please see the `LICENSE
<https://github.com/biopython/biopython/blob/master/LICENSE.rst>`_ file for
further details.

If you use Biopython in work contributing to a scientific publication, we ask
that you cite our application note (below) or one of the module specific
publications (listed on our website):

Cock, P.J.A. et al. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3
https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878


For the impatient
=================

Python includes the package management system "pip" which should allow you to
install Biopython (and its dependency NumPy if needed), upgrade or uninstall
with just one terminal command::

    pip install biopython
    pip install --upgrade biopython
    pip uninstall biopython

Since Biopython 1.70 we have provided pre-compiled binary wheel packages on
PyPI for Linux, Mac OS X and Windows. This means pip install should be quick,
and not require a compiler.

As a developer or potential contributor, you may wish to download, build and
install Biopython yourself. This is described below.


Python Requirements
===================

We currently recommend using Python 3.10 from http://www.python.org

Biopython is currently supported and tested on the following Python
implementations:

- Python 3.7, 3.8, 3.9, 3.10 and 3.11 -- see http://www.python.org

- PyPy3.7 v7.3.5 -- or later, see http://www.pypy.org


Optional Dependencies
=====================

Biopython requires NumPy (see http://www.numpy.org) which will be installed
automatically if you install Biopython with pip (see below for compiling
Biopython yourself).

Depending on which parts of Biopython you plan to use, there are a number of
other optional Python dependencies, which can be installed later if needed:

- ReportLab, see http://www.reportlab.com/opensource/ (optional)
  This package is only used in ``Bio.Graphics``, so if you do not need this
  functionality, you will not need to install this package.

- matplotlib, see http://matplotlib.org/ (optional)
  ``Bio.Phylo`` uses this package to plot phylogenetic trees.

- networkx, see https://networkx.github.io/ (optional) and
  pygraphviz or pydot, see https://pygraphviz.github.io/ and
  http://code.google.com/p/pydot/ (optional)
  These packages are used for certain niche functions in ``Bio.Phylo``.

- rdflib, see https://github.com/RDFLib/rdflib (optional)
  This package is used in the CDAO parser under ``Bio.Phylo``.

- psycopg2, see http://initd.org/psycopg/ (optional) or
  PyGreSQL (pgdb), see http://www.pygresql.org/ (optional)
  These packages are used by ``BioSQL`` to access a PostgreSQL database.

- MySQL Connector/Python, see http://dev.mysql.com/downloads/connector/python/
  This package is used by ``BioSQL`` to access a MySQL database, and is
  supported on PyPy too.

- mysqlclient, see https://github.com/PyMySQL/mysqlclient-python (optional)
  This is a fork of the older MySQLdb and is used by ``BioSQL`` to access a
  MySQL database. It is supported by PyPy.

In addition there are a number of useful third party tools you may wish to
install such as standalone NCBI BLAST, EMBOSS or ClustalW.


Installation From Source
========================

We recommend using the pre-compiled binary wheels available on PyPI using::

    pip install biopython

However, if you need to compile Biopython yourself, the following are required
at compile time:

- Python including development header files like ``python.h``, which on Linux
  are often not installed by default (trying looking for and installing a
  package named ``python-dev`` or ``python-devel`` as well as the ``python``
  package).

- Appropriate C compiler for your version of Python, for example GCC on Linux,
  MSVC on Windows. For Mac OS X, or as it is now branded, macOS, use Apple's
  command line tools, which can be installed with the terminal command::

      xcode-select --install

  This will offer to install Apple's XCode development suite - you can, but it
  is not needed and takes a lot of disk space.

Then either download and decompress our source code, or fetch it using git.
Now change directory to the Biopython source code folder and run::

    pip install -e .
    python setup.py test
    sudo python setup.py install

Substitute ``python`` with your specific version if required, for example
``python3``, or ``pypy3``.

To exclude tests that require an internet connection (and which may take a
long time), use the ``--offline`` option::

    python setup.py test --offline

If you need to do additional configuration, e.g. changing the install
directory prefix, please type ``python setup.py``.


Testing
=======

Biopython includes a suite of regression tests to check if everything is
running correctly. To run the tests, go to the biopython source code
directory and type::

    pip install -e .
    python setup.py test

If you want to skip the online tests (which is recommended when doing repeated
testing), use::

    python setup.py test --offline

Do not panic if you see messages warning of skipped tests::

    test_DocSQL ... skipping. Install MySQLdb if you want to use Bio.DocSQL.

This most likely means that a package is not installed.  You can
ignore this if it occurs in the tests for a module that you were not
planning on using.  If you did want to use that module, please install
the required dependency and re-run the tests.

Some of the tests may fail due to network issues, this is often down to
chance or a service outage. If the problem does not go away on
re-running the tests, you can use the ``--offline`` option.

There is more testing information in the Biopython Tutorial & Cookbook.


Experimental code
=================

Biopython 1.61 introduced a new warning, ``Bio.BiopythonExperimentalWarning``,
which is used to mark any experimental code included in the otherwise
stable Biopython releases. Such 'beta' level code is ready for wider
testing, but still likely to change, and should only be tried by early
adopters in order to give feedback via the biopython-dev mailing list.

We'd expect such experimental code to reach stable status within one or two
releases, at which point our normal policies about trying to preserve
backwards compatibility would apply.


Bugs
====

While we try to ship a robust package, bugs inevitably pop up.  If you are
having problems that might be caused by a bug in Biopython, it is possible
that it has already been identified. Update to the latest release if you are
not using it already, and retry. If the problem persists, please search our
bug database and our mailing lists to see if it has already been reported
(and hopefully fixed), and if not please do report the bug. We can't fix
problems we don't know about ;)

Issue tracker: https://github.com/biopython/biopython/issues

If you suspect the problem lies within a parser, it is likely that the data
format has changed and broken the parsing code.  (The text BLAST and GenBank
formats seem to be particularly fragile.)  Thus, the parsing code in
Biopython is sometimes updated faster than we can build Biopython releases.
You can get the most recent parser by pulling the relevant files (e.g. the
ones in ``Bio.SeqIO`` or ``Bio.Blast``) from our git repository. However, be
careful when doing this, because the code in github is not as well-tested
as released code, and may contain new dependencies.

In any bug report, please let us know:

1. Which operating system and hardware (32 bit or 64 bit) you are using
2. Python version
3. Biopython version (or git commit/date)
4. Traceback that occurs (the full error message)

And also ideally:

5. Example code that breaks
6. A data file that causes the problem


Contributing, Bug Reports
=========================

Biopython is run by volunteers from all over the world, with many types of
backgrounds. We are always looking for people interested in helping with code
development, web-site management, documentation writing, technical
administration, and whatever else comes up.

If you wish to contribute, please first read `CONTRIBUTING.rst
<https://github.com/biopython/biopython/blob/master/CONTRIBUTING.rst>`_ here,
visit our web site http://biopython.org and join our mailing list:
http://biopython.org/wiki/Mailing_lists


Distribution Structure
======================

- ``README.rst``  -- This file.
- ``NEWS.rst``    -- Release notes and news.
- ``LICENSE.rst`` -- What you can do with the code.
- ``CONTRIB.rst`` -- An (incomplete) list of people who helped Biopython in
  one way or another.
- ``CONTRIBUTING.rst`` -- An overview about how to contribute to Biopython.
- ``DEPRECATED.rst`` -- Contains information about modules in Biopython that
  were removed or no longer recommended for use, and how to update code that
  uses those modules.
- ``MANIFEST.in`` -- Configures which files to include in releases.
- ``setup.py``    -- Installation file.
- ``Bio/``        -- The main code base code.
- ``BioSQL/``     -- Code for using Biopython with BioSQL databases.
- ``Doc/``        -- Documentation.
- ``Scripts/``    -- Miscellaneous, possibly useful, standalone scripts.
- ``Tests/``      -- Regression testing code including sample data files.
