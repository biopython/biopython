Contributing to the Biopython Tutorial
======================================

Installing the tools you need
-----------------------------

On Ubuntu, you can install these with the following commands::

    sudo apt install hevea texlive-latex-base texlive-latex-extra python-pygments

Formatting code examples
------------------------

Code examples should be formatted using minted. Also, the Biopython tutorial
uses its own system for testing code examples, which is located in
``test_Tutorial.py``. You can read more about this system in the testing
chapter of the Biopython Tutorial.

Here is an example of Python code that uses the console to show output::

    %doctest path-to-folder
    \begin{minted}{pycon}
    >>> from Bio import SeqIO
    \end{minted}

These examples should use %doctest to verify that the output is correct.

Here is an example of Python code that does not show output::

    \begin{minted}{python}
    from Bio import SeqIO
    \end{minted}


Formatting text
---------------

Format headings as follows:

Chapter titles::

    \chapter{Title of Chapter}
    \label{chapter:chapterlabel}

Section titles::

    \section{Title of Section}
    \label{sec:sectionlabel}

Subsection titles::

    \subsection{Title of Subsection}

When referring to code in the middle of a paragraph, format it as follows::

    \verb|variable_name|

This will render ``variable_name`` in a monospace font. Within the pipe
characters, underscores will be interpreted literally, so they do not need
to be escaped.

Building the documentation
--------------------------

Build the documentation by running "make" in the Docs folder.

Once the documentation has been generated, you can inspect either Tutorial.pdf
or Tutorial.html in the Docs directory to see if the output is correct.
