Contributing to the Biopython Tutorial
======================================

Installing the requirements
---------------------------

In the Biopython directory, install the required packages for building the documentation with::

    pip install -r .circleci/requirements-sphinx.txt

ReStructered Text Markup Language
---------------------------------

The tutorials are written in `ReStructered Text <https://en.m.wikipedia.org/wiki/ReStructuredText>`_ format.
Read the `Sphix reST Primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#hyperlinks>`_ to learn how to use this format.

Building the documentation
--------------------------

Build the documentation in HTML by running ``make html`` in the Docs folder.
Make sure that you have installed necessary requirements (see instructions above).

To view the generated documentation, open ``Doc/_build/html/index.html`` in your browser.

Build a PDF of the documentation using ``make latexpdf``.
This command depends on `Latexmk <https://mg.readthedocs.io/latexmk.html>`_ and LaTeX tooling.

On Ubuntu, you can install these with the following commands::

    sudo apt install texlive-latex-base texlive-latex-extra python-pygments

On Fedora (tested for Workstation 33)::

    dnf install texlive-scheme-basic texlive-preprint \
                texlive-comment texlive-minted

Testing code examples
------------------------

The Biopython tutorial uses its own system for testing code examples, which is located in
``Tests/test_Tutorial.py``. You can read more about this system in the testing
chapter of the Biopython Tutorial.
