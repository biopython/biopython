.. image:: https://img.shields.io/pypi/v/biopython.svg
   :alt: Biopython on the Python Package Index (PyPI)
   :target: https://pypi.python.org/pypi/biopython
.. image:: https://img.shields.io/travis/biopython/biopython/master.svg
   :alt: Testing with TravisCI
   :target: https://travis-ci.org/biopython/biopython/branches
.. image:: https://img.shields.io/codecov/c/github/biopython/biopython/master.svg
   :alt: TravisCI test coverage
   :target: https://codecov.io/github/biopython/biopython/

.. image:: https://github.com/biopython/biopython/raw/master/Doc/images/biopython.jpg
   :alt: The Biopython Project
   :target: http://biopython.org
   :width: 100%

Biopython README file
=====================

The Biopython Project is an international association of developers of freely
available Python tools for computational molecular biology.

Our user-centric documentation is hosted on http://biopython.org including
the main Biopython Tutorial and Cookbook:

* HTML - http://biopython.org/DIST/docs/tutorial/Tutorial.html
* PDF - http://biopython.org/DIST/docs/tutorial/Tutorial.pdf

This README file is intended primarily for people interested in working
with the Biopython source code, either one of the releases from the
http://biopython.org website, or from our repository on GitHub
https://github.com/biopython/biopython

This Biopython package is open source software made available under generous
terms. Please see the ``LICENSE.rst`` file for further details.

If you use Biopython in work contributing to a scientific publication, we ask
that you cite our application note (below) or one of the module specific
publications (listed on our website):

Cock, P.J.A. et al. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878


For the impatient
=================

Windows users are currently recommended to use the installation packages provided
on our website, http://biopython.org -- further instructions are given below.

Python 2.7.9 onwards, and Python 3.4 onwards, include the package management
system "pip" which should allow you to install Biopython with just::

    pip install numpy
    pip install biopython

Otherwise you may have to build and install Biopython. Download and unzip the
source code, go to this directory at the command line, and type::

    python setup.py build
    python setup.py test
    sudo python setup.py install

Here you can replace ``python`` with a specific version, e.g. ``python3.5``.


Python Requirements
===================

We currently recommend using Python 3.6 from http://www.python.org

Biopython is currently supported and tested on the following Python
implementations:

- Python 2.7, 3.3, 3.4, 3.5, 3.6 -- see http://www.python.org

  This is the primary development platform for Biopython.

- PyPy v5.7 and also PyPy3.5 v5.7 beta -- see http://www.pypy.org

  Aside from ``Bio.trie`` (which does not compile as ``marshal.h`` is
  currently missing under PyPy), everything should work. Older versions
  of PyPy mostly work too.

- Jython 2.7 -- see http://www.jython.org

  We provide limited support for Jython, but aside from ``Bio.Restriction``,
  modules with C code, or dependent on SQLite3 or NumPy, everything should
  work. There are some known issues with test failures which have not yet
  been resolved.

Please note that support for Python 3.3 is deprecated as of Biopython 1.67.
Biopython 1.68 was our final release to support Python 2.6.


Dependencies
============

Depending on which parts of Biopython you plan to use, there are a number of
other optional Python dependencies, which (other than NumPy) can be installed
after Biopython.

- NumPy, see http://www.numpy.org (optional, but strongly recommended)
  This package is only used in the computationally-oriented modules. It is
  required for ``Bio.Cluster``, ``Bio.PDB`` and a few other modules.  If you
  think you might need these modules, then please install NumPy first BEFORE
  installing Biopython. The older Numeric library is no longer supported in
  Biopython.

- ReportLab, see http://www.reportlab.com/opensource/ (optional)
  This package is only used in ``Bio.Graphics``, so if you do not need this
  functionality, you will not need to install this package.  You can install
  it later if needed.

- matplotlib, see http://matplotlib.org/ (optional)
  ``Bio.Phylo`` uses this package to plot phylogenetic trees. As with
  ReportLab, you can install this at any time to enable the plotting
  functionality.

- networkx, see http://networkx.lanl.gov/ (optional) and
  pygraphviz or pydot, see http://networkx.lanl.gov/pygraphviz/ and
  http://code.google.com/p/pydot/ (optional)
  These packages are used for certain niche functions in ``Bio.Phylo``.
  Again, they are only needed to enable these functions and can be installed
  later if needed.

- rdflib, see https://github.com/RDFLib/rdflib (optional)
  This package is used in the CDAO parser under ``Bio.Phylo``, and can be
  installed as needed.

- psycopg2, see http://initd.org/psycopg/ (optional) or
  PyGreSQL (pgdb), see http://www.pygresql.org/ (optional)
  These packages are used by ``BioSQL`` to access a PostgreSQL database.

- MySQL Connector/Python, see http://dev.mysql.com/downloads/connector/python/
  This package is used by ``BioSQL`` to access a MySQL database, and is
  supported on Python 2 and 3 and PyPy too.

- MySQLdb, see http://sourceforge.net/projects/mysql-python (optional)
  This is an older alternative package used by ``BioSQL`` to access a MySQL
  database, but it is not available for Python 3 or PyPy.

Note that some of these libraries are not available for Jython or PyPy,
and not all are available for Python 3 yet either.

In addition there are a number of useful third party tools you may wish to
install such as standalone NCBI BLAST, EMBOSS or ClustalW.


Installation
============

First, **make sure that Python is installed correctly**. Second, we
recommend you install NumPy (see above). Then install Biopython.

Windows users should use the appropriate provided installation package
from our website (each is specific to a different Python version).

Installation from source should be as simple as going to the Biopython
source code directory, and typing::

    python setup.py build
    python setup.py test
    sudo python setup.py install

Substitute `python` with your specific version, for example `python3`,
`jython` or `pypy`.

If you need to do additional configuration, e.g. changing the base
directory, please type ``python setup.py``, or see the documentation here:

* HTML - http://biopython.org/DIST/docs/install/Installation.html
* PDF - http://biopython.org/DIST/docs/install/Installation.pdf


Testing
=======

Biopython includes a suite of regression tests to check if everything is
running correctly. To run the tests, go to the biopython source code
directory and type::

    python setup.py build
    python setup.py test

Do not panic if you see messages warning of skipped tests::

    test_DocSQL ... skipping. Install MySQLdb if you want to use Bio.DocSQL.

This most likely means that a package is not installed.  You can
ignore this if it occurs in the tests for a module that you were not
planning on using.  If you did want to use that module, please install
the required dependency and re-run the tests.

Some of the tests may fail due to network issues, this is often down to
chance or a service outage. If the problem does not go away on
re-running the tests, it is possible to run only the offline tests.

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

* Old issue tracker: https://redmine.open-bio.org/projects/biopython
* Current issue tracker: https://github.com/biopython/biopython/issues

If you suspect the problem lies within a parser, it is likely that the data
format has changed and broken the parsing code.  (The text BLAST and GenBank
formats seem to be particularly fragile.)  Thus, the parsing code in
Biopython is sometimes updated faster than we can build Biopython releases.
You can get the most recent parser by pulling the relevant files (e.g. the
ones in ``Bio.SeqIO`` or ``Bio.Blast``) from our git repository. However, be
careful when doing this, because the code in github is not as well-tested
as released code, and may contain new dependencies.

Finally, you can send a bug report to the bug database or the mailing list at
biopython@biopython.org (subscription required).  In the bug report, please
let us know:

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

If you wish to contribute, please visit the web site http://biopython.org
and join our mailing list: http://biopython.org/wiki/Mailing_lists


Distribution Structure
======================

- ``README.rst``  -- This file.
- ``NEWS.rst``    -- Release notes and news.
- ``LICENSE.rst`` -- What you can do with the code.
- ``CONTRIB.rst`` -- An (incomplete) list of people who helped Biopython in
  one way or another.
- ``DEPRECATED.rst`` -- Contains information about modules in Biopython that are
  removed or no longer recommended for use, and how to update code that uses
  those modules.
- ``MANIFEST.in`` -- Tells distutils what files to distribute.
- ``setup.py``    -- Installation file.
- ``Bio/``        -- The main code base code.
- ``BioSQL/``     -- Code for using Biopython with BioSQL databases.
- ``Doc/``        -- Documentation.
- ``Scripts/``    -- Miscellaneous, possibly useful, standalone scripts.
- ``Tests/``      -- Regression testing code including sample data files.
