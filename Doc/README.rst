Contributing to the Biopython Tutorial
======================================

Installing the tools you need
-----------------------------

On Ubuntu, you can install these with the following commands::

    sudo apt install hevea texlive-latex-base texlive-latex-extra python-pygments

Formatting code examples
------------------------

Code examples should be formatted using minted. Also, the Biopython tutorial uses its own system for testing code examples, which is located in ``test_Tutorial.py``. You can read more about this system in the testing chapter of the Biopython Tutorial.

Here is an example::

    %doctest path-to-folder
    \begin{minted}{python}
    >>> from Bio import SeqIO
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

    \verb|ObjectName|

This will render ObjectName in a monospace font.

Building the documentation
--------------------------

Build the documentation by running "make" in the Docs folder.

Ideally, the documentation should be generated without a hitch, but some errors may cause the program to stop and ask for user input. Press R to continue building the documentation.

If you would like to ignore errors at the command-line, then you can use::

    yes R | make

However, these errors may provide valuable information if there are mistakes in your own LaTeX code.

Once the documentation has been generated, you can inspect either Tutorial.pdf or Tutorial.html in the Docs directory to see if the output is correct.
