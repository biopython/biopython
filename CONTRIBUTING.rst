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

To facilitate style checking, we make use of ``pre-commit`` to automatically
run various ``flake8`` plugins as well as ``black``. We use the continuous
integration service TravisCI to run these checks, but we strongly suggest you
install and run ``pre-commit`` in your local machine (see below).

For docstrings (Python's in-code documentation), in addition to PEP257 we are
using reStructuredText (RST) markup language which allows basic formatting
like *italics* and **bold** once rendered into HTML webpages for our online
API documentation.

We also use reStructuredText for files like ``README.rst``. You can run the
reStructuredText checks with the ``restructuredtext-lint`` tool (also known as
``rst-lint``)::

    $ pip install restructuredtext_lint
    $ rst-lint --level warning *.rst


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
repository, it will automatically run various fast checks (including flake8 and
black) before the commit happens. In order to install it, run::

    $ pip install pre-commit

    # Activate pre-commit for biopython
    $ cd biopython-repository
    $ pre-commit install


Local Testing
-------------

Please always run the style checks (see above) and the full test suite on
your local computer before submitting a pull request, e.g.::

    $ git commit Bio/XXX.py Tests/test_XXX.py  -m "Fixed bug 123"
    $ python setup.py build
    $ python setup.py test

If you have multiple versions of Python installed, ideally test them all
(the Python tool ``tox`` can be helpful here).


Continuous Integration
----------------------

Once you submit your pull request on GitHub, several automated tests should
be run, and their results reported on the pull request.

We use TravisCI to run most of the Biopython tests (although currently only
under Linux, and not with all the optional dependencies included), plus also
check Python coding style using the ``flake8`` tool with lots of plugins, and
reStructuredText files using ``rst-lint`` and ``doc8``.
https://travis-ci.org/biopython/biopython/branches

We use AppVeyor to run most of the tests under Windows (although currently
without any optional dependencies).
https://ci.appveyor.com/project/biopython/biopython/history

**The TravisCI and AppVeyor checks must pass before your pull request will
be merged.**

Some of the TravisCI and AppVeyor runs collect test coverage information via
CodeCov: https://codecov.io/github/biopython/biopython/

Ideally the CodeCov checks will also pass, but we currently do not insist
on this when reviewing pull requests.

Contributing to the Biopython Tutorial
--------------------------------------

For instructions, see `the Biopython Tutorial README <Doc/README.rst>`_
