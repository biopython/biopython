.. _`chapter:appendix`:

Appendix: Useful stuff about Python
===================================

If you haven’t spent a lot of time programming in Python, many questions
and problems that come up in using Biopython are often related to Python
itself. This section tries to present some ideas and code that come up
often (at least for us!) while using the Biopython libraries. If you
have any suggestions for useful pointers that could go here, please
contribute!

.. _`sec:appendix-handles`:

What the heck is a handle?
--------------------------

Handles are mentioned quite frequently throughout this documentation,
and are also fairly confusing (at least to me!). Basically, you can
think of a handle as being a “wrapper” around text information.

Handles provide (at least) two benefits over plain text information:

#. They provide a standard way to deal with information stored in
   different ways. The text information can be in a file, or in a string
   stored in memory, or the output from a command line program, or at
   some remote website, but the handle provides a common way of dealing
   with information in all of these formats.

#. They allow text information to be read incrementally, instead of all
   at once. This is really important when you are dealing with huge text
   files which would use up all of your memory if you had to load them
   all.

Handles can deal with text information that is being read (e. g. reading
from a file) or written (e. g. writing information to a file). In the
case of a “read” handle, commonly used functions are ``read()``, which
reads the entire text information from the handle, and ``readline()``,
which reads information one line at a time. For “write” handles, the
function ``write()`` is regularly used.

The most common usage for handles is reading information from a file,
which is done using the built-in Python function ``open``. Here, we
handle to the file ``m_cold.fasta`` which you can download
`here <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/m_cold.fasta>`__
(or find included in the Biopython source code as
``Doc/examples/m_cold.fasta``).

.. code:: pycon

   >>> handle = open("m_cold.fasta", "r")
   >>> handle.readline()
   ">gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum ...\n"

Handles are regularly used in Biopython for passing information to
parsers. For example, since Biopython 1.54 the main functions in
``Bio.SeqIO`` and ``Bio.AlignIO`` have allowed you to use a filename
instead of a handle:

.. code:: python

   from Bio import SeqIO

   for record in SeqIO.parse("m_cold.fasta", "fasta"):
       print(record.id, len(record))

On older versions of Biopython you had to use a handle, e.g.

.. code:: python

   from Bio import SeqIO

   handle = open("m_cold.fasta", "r")
   for record in SeqIO.parse(handle, "fasta"):
       print(record.id, len(record))
   handle.close()

This pattern is still useful - for example suppose you have a gzip
compressed FASTA file you want to parse:

.. code:: python

   import gzip
   from Bio import SeqIO

   handle = gzip.open("m_cold.fasta.gz", "rt")
   for record in SeqIO.parse(handle, "fasta"):
       print(record.id, len(record))
   handle.close()

With our parsers for plain text files, it is essential to use gzip in
text mode (the default is binary mode).

See Section :ref:`sec:SeqIO_compressed` for more
examples like this, including reading bzip2 compressed files.

Creating a handle from a string
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One useful thing is to be able to turn information contained in a string
into a handle. The following example shows how to do this using
``StringIO`` from the Python standard library:

.. doctest

.. code:: pycon

   >>> my_info = "A string\n with multiple lines."
   >>> print(my_info)
   A string
    with multiple lines.
   >>> from io import StringIO
   >>> my_info_handle = StringIO(my_info)
   >>> first_line = my_info_handle.readline()
   >>> print(first_line)
   A string
   <BLANKLINE>
   >>> second_line = my_info_handle.readline()
   >>> print(second_line)
    with multiple lines.
