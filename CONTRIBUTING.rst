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

For docstrings (Python's in-code documentation) in addition to PEP257 we are
using reStructuredText (RST) markup language which allows basic formatting
like *italics* and **bold** once rendered into HTML webpages for our online
API documentation.

We use the continuous integration service TravisCI to enforce some of these
checks, so if you are making a contribution it is best to check this locally.

You can run the tool ``flake8`` directly as follows. This includes the checks
from  ``pycodestyle`` (formerly ``pep8``) and uses the ``pydocstyle`` plugin
(``flake8-docstrings``) for checking the docstrings within the Python code::

    $ pip install flake8 flake8-docstrings flake8-blind-except
    $ flake8 Bio/
    $ flake8 BioSQL/
    $ flake8 Tests/
    $ flake8 Scripts/
    $ flake8 Doc/examples/

Each of these folders has its own flake8 settings file, in the long term we
would like all the code to follow the same strict settings.

You can run the reStructuredText checks with the ``restructuredtext-lint``
tool (also known as ``rst-lint``)::

    $ pip install restructuredtext_lint
    $ rst-lint --level warning *.rst

More simply we currently suggest you install the git pre-commit hook described
here which will check our basic coding conventions as you work:
https://github.com/biopython/biopython/issues/493



Testing
-------

Any new feature or functionality will not be accepted without tests. Likewise
for any bug fix, we encourage including an additional test. See the testing
chapter in the Biopython Tutorial for more information on our test framework:
http://biopython.org/DIST/docs/tutorial/Tutorial.html


Local Testing
-------------

Please always run the style checks (see above) and the full test suite on
your local computer before submitting a pull request, e.g.::

    $ git commit Bio/XXX.py Tests/test_XXX.py  -m "Fixed bug 123"
    $ python3.5 setup.py build
    $ python3.5 setup.py test

If you have multiple versions of Python installed, ideally test them all
(the Python tool ``tox`` can be helpful here).


Continous Integration
---------------------

Once you submit your pull request on GitHub, several automated tests should
be run, and their results reported on the pull request.

We use TravisCI to run most of the Biopython tests (although currently only
under Linux, and not with all the optional dependencies included), plus also
check Python coding style using the ``flake8`` tool with the ``pydocstyle``
pluging (``flake8-docstrings``), and reStructuredText using ``rst-lint``.
These checks must pass before your pull request will be merged.
https://travis-ci.org/biopython/biopython

Some of the TravisCI runs collect test coverage information via CodeCov:
https://codecov.io/github/biopython/biopython/

We also are using QuantifiedCode which provides an independent code analysis
which can catch important errors. Some of the issues it flags are subjective,
so in general we do not require this pass before merging a pull request.
https://www.quantifiedcode.com/app/project/gh:biopython:biopython
