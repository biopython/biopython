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

We have adopted the code formating tool ``black`` to which automates a lot of
coding style choices like line breaks and spacing. This applies to all our
code, except the ``Tests/`` folder which we are still gradually updating.

We use the tool ``flake8`` for automating our code style checks, together with
various plugins which can be installed as follows::

    $ pip install flake8 flake8-docstrings flake8-black flake8-blind-except flake8-rst-docstrings flake8-comprehensions flake8-quotes flake8-bugbear flake8-implicit-str-concat

We use the continuous integration service TravisCI to run these ``flake8``
checks, but it is better to run them locally before making a pull request.
We currently strongly suggest you then install the ``flake8`` git pre-commit
hook which will check our basic coding conventions as you work::

    $ flake8 --install-hook git
    $ git config --bool flake8.strict true

For testing, you can also run ``flake8`` directly from the command line as
follows::

    $ flake8

For docstrings (Python's in-code documentation) in addition to PEP257 we are
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
