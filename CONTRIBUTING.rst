Contributing to Biopython
=========================

We welcome pull requests to fix bugs or add new features. Please read
http://biopython.org/wiki/Contributing for a general overview about
contributing - this file is primarily concerned with the practicalities.


Licensing
---------

Biopython is moving towards dual licensing the code. In your git commit and/or
pull request, please explicitly state that you agree to your contributions
being dual licensed under *both* our original "Biopython License Agreement"
and the more widely used "3-Clause BSD License" (see our LICENSE file for more
details).


Git Usage
---------

We have a git introduction online at http://biopython.org/wiki/GitUsage

If you are planning to make a pull request, start by creating a new branch
with a short but descriptive name (rather than using your master branch).


Coding Conventions
------------------

Biopython tries to follow the coding conventions laid out in PEP8 and PEP257,
with the notable exception of existing module names which are not lower case.

- http://www.python.org/dev/peps/pep-0008/
- http://www.python.org/dev/peps/pep-0257/

For docstrings (Python's in-code documentation), in addition to PEP257 we are
using reStructuredText (RST) markup language which allows basic formatting
like *italics* and **bold** once rendered into HTML webpages for our online
API documentation. We also use reStructuredText for files like ``README.rst``.

To facilitate style checking, we make use of ``pre-commit`` to automatically
run various ``flake8`` plugins as well as ``black`` on the Python code, and
``restructuredtext-lint`` (also known as ``rst-lint``) to check the RST files
too. We also use continuous integration on GitHub to run these checks, but we
strongly suggest you install and run ``pre-commit`` in your local machine (see
below).


Testing
-------

Any new feature or functionality will not be accepted without tests. Likewise
for any bug fix, we encourage including an additional test. See the testing
chapter in the Biopython Tutorial for more information on our test framework:
http://biopython.org/DIST/docs/tutorial/Tutorial.html


Local Development
-----------------

As mentioned above, to simplify your contributions, we provide a `pre-commit
<https://pre-commit.com/>`_ configuration. Pre-commit is a Python package which
automatically hooks into git. When you run "git commit" within the biopython
repository, it will automatically run various fast checks (including ``black``,
``flake8``, ``rst-lint`` and ``doc8``) before the commit happens. To install
this, run::

    $ pip install pre-commit

    # Activate pre-commit for biopython
    $ cd biopython-repository
    $ pre-commit install


Local Testing
-------------

Please always run the full test suite locally before submitting a pull
request, e.g.::

    $ python setup.py build
    $ python setup.py test
    $ git commit ...

Have a look at the `related chapter <http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter%3Atesting>`_ in the documentation for more details.

Continuous Integration
----------------------

Once you submit your pull request on GitHub, several automated tests should be
run via continuous integration services, and their results reported on the pull
request. These will run most of the Biopython tests (although not with all the
optional dependencies included), plus also style checks using ``pre-commit``
(also used for git pre-commit checks, see above).

**The continuous integration checks must pass before your pull request will be
merged.**

The continuous integration tests collect test coverage information via
CodeCov: https://codecov.io/github/biopython/biopython/

Ideally the CodeCov checks will also pass, but we currently do not insist on
this when reviewing pull requests.

Contributing to the Biopython Tutorial
--------------------------------------

For instructions, see `the Biopython Tutorial README <Doc/README.rst>`_
