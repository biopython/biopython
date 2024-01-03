# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Biopython interface for sequence search program outputs.

The SearchIO submodule provides parsers, indexers, and writers for outputs from
various sequence search programs. It provides an API similar to SeqIO and
AlignIO, with the following main functions: ``parse``, ``read``, ``to_dict``,
``index``, ``index_db``, ``write``, and ``convert``.

SearchIO parses a search output file's contents into a hierarchy of four nested
objects: QueryResult, Hit, HSP, and HSPFragment. Each of them models a part of
the search output file:

    - QueryResult represents a search query. This is the main object returned
      by the input functions and it contains all other objects.
    - Hit represents a database hit,
    - HSP represents high-scoring alignment region(s) in the hit,
    - HSPFragment represents a contiguous alignment within the HSP

In addition to the four objects above, SearchIO is also tightly integrated with
the SeqRecord objects (see SeqIO) and MultipleSeqAlignment objects (see
AlignIO). SeqRecord objects are used to store the actual matching hit and query
sequences, while MultipleSeqAlignment objects store the alignment between them.

A detailed description of these objects' features and their example usages are
available in their respective documentations.


Input
=====
The main function for parsing search output files is Bio.SearchIO.parse(...).
This function parses a given search output file and returns a generator object
that yields one QueryResult object per iteration.

``parse`` takes two arguments: 1) a file handle or a filename of the input file
(the search output file) and 2) the format name.

    >>> from Bio import SearchIO
    >>> for qresult in SearchIO.parse('Blast/mirna.xml', 'blast-xml'):
    ...     print("%s %s" % (qresult.id, qresult.description))
    ...
    33211 mir_1
    33212 mir_2
    33213 mir_3

SearchIO also provides the Bio.SearchIO.read(...) function, which is intended
for use on search output files containing only one query. ``read`` returns one
QueryResult object and will raise an exception if the source file contains more
than one query:

    >>> qresult = SearchIO.read('Blast/xml_2226_blastp_004.xml', 'blast-xml')
    >>> print("%s %s" % (qresult.id, qresult.description))
    ...
    gi|11464971:4-101 pleckstrin [Mus musculus]

    >>> SearchIO.read('Blast/mirna.xml', 'blast-xml')
    Traceback (most recent call last):
    ...
    ValueError: ...

For accessing search results of large output files, you may use the indexing
functions Bio.SearchIO.index(...) or Bio.SearchIO.index_db(...). They have a
similar interface to their counterparts in SeqIO and AlignIO, with the addition
of optional, format-specific keyword arguments.


Output
======
SearchIO has write support for several formats, accessible from the
Bio.SearchIO.write(...) function. This function returns a tuple of four
numbers: the number of QueryResult, Hit, HSP, and HSPFragment written::

    qresults = SearchIO.parse('Blast/mirna.xml', 'blast-xml')
    SearchIO.write(qresults, 'results.tab', 'blast-tab')
    <stdout> (3, 239, 277, 277)

Note that different writers may require different attribute values of the
SearchIO objects. This limits the scope of writable search results to search
results possessing the required attributes.

For example, the writer for HMMER domain table output requires
the conditional e-value attribute from each HSP object, among others. If you
try to write to the HMMER domain table format and your HSPs do not have this
attribute, an exception will be raised.


Conversion
==========
SearchIO provides a shortcut function Bio.SearchIO.convert(...) to convert a
given file into another format. Under the hood, ``convert`` simply parses a given
output file and writes it to another using the ``parse`` and ``write`` functions.

Note that the same restrictions found in Bio.SearchIO.write(...) apply to the
convert function as well.


Conventions
===========
The main goal of creating SearchIO is to have a common, easy to use interface
across different search output files. As such, we have also created some
conventions / standards for SearchIO that extend beyond the common object model.
These conventions apply to all files parsed by SearchIO, regardless of their
individual formats.

Python-style sequence coordinates
---------------------------------

When storing sequence coordinates (start and end values), SearchIO uses
the Python-style slice convention: zero-based and half-open intervals. For
example, if in a BLAST XML output file the start and end coordinates of an
HSP are 10 and 28, they would become 9 and 28 in SearchIO. The start
coordinate becomes 9 because Python indices start from zero, while the end
coordinate remains 28 as Python slices omit the last item in an interval.

Beside giving you the benefits of standardization, this convention also
makes the coordinates usable for slicing sequences. For example, given a
full query sequence and the start and end coordinates of an HSP, one can
use the coordinates to extract part of the query sequence that results in
the database hit.

When these objects are written to an output file using
SearchIO.write(...), the coordinate values are restored to their
respective format's convention. Using the example above, if the HSP would
be written to an XML file, the start and end coordinates would become 10
and 28 again.

Sequence coordinate order
-------------------------

Some search output formats reverse the start and end coordinate sequences
according to the sequence's strand.

In SearchIO, start coordinates are always smaller than the end
coordinates, regardless of their originating strand. This ensures
consistency when using the coordinates to slice full sequences.

Note that this coordinate order convention is only enforced in the
HSPFragment level. If an HSP object has several HSPFragment objects, each
individual fragment will conform to this convention. But the order of the
fragments within the HSP object follows what the search output file uses.

Similar to the coordinate style convention, the start and end coordinates'
order are restored to their respective formats when the objects are
written using Bio.SearchIO.write(...).

Frames and strand values
------------------------

SearchIO only allows -1, 0, 1 and None as strand values. For frames, the
only allowed values are integers from -3 to 3 (inclusive) and None. Both
of these are standard Biopython conventions.


Supported Formats
=================
Below is a list of search program output formats supported by SearchIO.

Support for parsing, indexing, and writing:

 - blast-tab        - BLAST+ tabular output. Both variants without comments
                      (-m 6 flag) and with comments (-m 7 flag) are supported.
 - blast-xml        - BLAST+ XML output.
 - blat-psl         - The default output of BLAT (PSL format). Variants with or
                      without header are both supported. PSLX (PSL + sequences)
                      is also supported.
 - hmmer3-tab       - HMMER3 table output.
 - hmmer3-domtab    - HMMER3 domain table output. When using this format, the
                      program name has to be specified. For example, for parsing
                      hmmscan output, the name would be 'hmmscan-domtab'.

Support for parsing and indexing:

 - exonerate-text   - Exonerate plain text output.
 - exonerate-vulgar - Exonerate vulgar line.
 - exonerate-cigar  - Exonerate cigar line.
 - fasta-m10        - Bill Pearson's FASTA -m 10 output.
 - hmmer3-text      - HMMER3 regular text output format. Supported HMMER3
                      subprograms are hmmscan, hmmsearch, and phmmer.
 - hmmer2-text      - HMMER2 regular text output format. Supported HMMER2
                      subprograms are hmmpfam, hmmsearch.

Support for parsing:

 - hhsuite2-text    - HHSUITE plain text output.

Each of these formats have different keyword arguments available for use with
the main SearchIO functions. More details and examples are available in each
of the format's documentation.

"""

from Bio.File import as_handle
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.SearchIO._utils import get_processor


__all__ = ("read", "parse", "to_dict", "index", "index_db", "write", "convert")


# dictionary of supported formats for parse() and read()
_ITERATOR_MAP = {
    "blast-tab": ("BlastIO", "BlastTabParser"),
    "blast-xml": ("BlastIO", "BlastXmlParser"),
    "blat-psl": ("BlatIO", "BlatPslParser"),
    "exonerate-cigar": ("ExonerateIO", "ExonerateCigarParser"),
    "exonerate-text": ("ExonerateIO", "ExonerateTextParser"),
    "exonerate-vulgar": ("ExonerateIO", "ExonerateVulgarParser"),
    "fasta-m10": ("FastaIO", "FastaM10Parser"),
    "hhsuite2-text": ("HHsuiteIO", "Hhsuite2TextParser"),
    "hhsuite3-text": ("HHsuiteIO", "Hhsuite2TextParser"),
    "hmmer2-text": ("HmmerIO", "Hmmer2TextParser"),
    "hmmer3-text": ("HmmerIO", "Hmmer3TextParser"),
    "hmmer3-tab": ("HmmerIO", "Hmmer3TabParser"),
    # for hmmer3-domtab, the specific program is part of the format name
    # as we need it distinguish hit / target coordinates
    "hmmscan3-domtab": ("HmmerIO", "Hmmer3DomtabHmmhitParser"),
    "hmmsearch3-domtab": ("HmmerIO", "Hmmer3DomtabHmmqueryParser"),
    "interproscan-xml": ("InterproscanIO", "InterproscanXmlParser"),
    "phmmer3-domtab": ("HmmerIO", "Hmmer3DomtabHmmqueryParser"),
}

# dictionary of supported formats for index()
_INDEXER_MAP = {
    "blast-tab": ("BlastIO", "BlastTabIndexer"),
    "blast-xml": ("BlastIO", "BlastXmlIndexer"),
    "blat-psl": ("BlatIO", "BlatPslIndexer"),
    "exonerate-cigar": ("ExonerateIO", "ExonerateCigarIndexer"),
    "exonerate-text": ("ExonerateIO", "ExonerateTextIndexer"),
    "exonerate-vulgar": ("ExonerateIO", "ExonerateVulgarIndexer"),
    "fasta-m10": ("FastaIO", "FastaM10Indexer"),
    "hmmer2-text": ("HmmerIO", "Hmmer2TextIndexer"),
    "hmmer3-text": ("HmmerIO", "Hmmer3TextIndexer"),
    "hmmer3-tab": ("HmmerIO", "Hmmer3TabIndexer"),
    "hmmscan3-domtab": ("HmmerIO", "Hmmer3DomtabHmmhitIndexer"),
    "hmmsearch3-domtab": ("HmmerIO", "Hmmer3DomtabHmmqueryIndexer"),
    "phmmer3-domtab": ("HmmerIO", "Hmmer3DomtabHmmqueryIndexer"),
}

# dictionary of supported formats for write()
_WRITER_MAP = {
    "blast-tab": ("BlastIO", "BlastTabWriter"),
    "blast-xml": ("BlastIO", "BlastXmlWriter"),
    "blat-psl": ("BlatIO", "BlatPslWriter"),
    "hmmer3-tab": ("HmmerIO", "Hmmer3TabWriter"),
    "hmmscan3-domtab": ("HmmerIO", "Hmmer3DomtabHmmhitWriter"),
    "hmmsearch3-domtab": ("HmmerIO", "Hmmer3DomtabHmmqueryWriter"),
    "phmmer3-domtab": ("HmmerIO", "Hmmer3DomtabHmmqueryWriter"),
}


def parse(handle, format=None, **kwargs):
    """Iterate over search tool output file as QueryResult objects.

    Arguments:
     - handle - Handle to the file, or the filename as a string.
     - format - Lower case string denoting one of the supported formats.
     - kwargs - Format-specific keyword arguments.

    This function is used to iterate over each query in a given search output
    file:

    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('Blast/mirna.xml', 'blast-xml')
    >>> qresults
    <generator object ...>
    >>> for qresult in qresults:
    ...     print("Search %s has %i hits" % (qresult.id, len(qresult)))
    ...
    Search 33211 has 100 hits
    Search 33212 has 44 hits
    Search 33213 has 95 hits

    Depending on the file format, ``parse`` may also accept additional keyword
    argument(s) that modifies the behavior of the format parser. Here is a
    simple example, where the keyword argument enables parsing of a commented
    BLAST tabular output file:

    >>> from Bio import SearchIO
    >>> for qresult in SearchIO.parse('Blast/mirna.tab', 'blast-tab', comments=True):
    ...     print("Search %s has %i hits" % (qresult.id, len(qresult)))
    ...
    Search 33211 has 100 hits
    Search 33212 has 44 hits
    Search 33213 has 95 hits

    """
    # get the iterator object and do error checking
    iterator = get_processor(format, _ITERATOR_MAP)

    # HACK: force BLAST XML decoding to use utf-8
    handle_kwargs = {}
    if format == "blast-xml":
        handle_kwargs["encoding"] = "utf-8"

    # and start iterating
    with as_handle(handle, **handle_kwargs) as source_file:
        generator = iterator(source_file, **kwargs)
        yield from generator


def read(handle, format=None, **kwargs):
    """Turn a search output file containing one query into a single QueryResult.

     - handle - Handle to the file, or the filename as a string.
     - format - Lower case string denoting one of the supported formats.
     - kwargs - Format-specific keyword arguments.

    ``read`` is used for parsing search output files containing exactly one query:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('Blast/xml_2226_blastp_004.xml', 'blast-xml')
    >>> print("%s %s" % (qresult.id, qresult.description))
    ...
    gi|11464971:4-101 pleckstrin [Mus musculus]

    If the given handle has no results, an exception will be raised:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('Blast/tab_2226_tblastn_002.txt', 'blast-tab')
    Traceback (most recent call last):
    ...
    ValueError: No query results found in handle

    Similarly, if the given handle has more than one result, an exception will
    be raised:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('Blast/tab_2226_tblastn_001.txt', 'blast-tab')
    Traceback (most recent call last):
    ...
    ValueError: More than one query result found in handle

    Like ``parse``, ``read`` may also accept keyword argument(s) depending on the
    search output file format.

    """
    query_results = parse(handle, format, **kwargs)

    try:
        query_result = next(query_results)
    except StopIteration:
        raise ValueError("No query results found in handle") from None
    try:
        next(query_results)
        raise ValueError("More than one query result found in handle")
    except StopIteration:
        pass

    return query_result


def to_dict(qresults, key_function=None):
    """Turn a QueryResult iterator or list into a dictionary.

     - qresults     - Iterable returning QueryResult objects.
     - key_function - Optional callback function which when given a
                      QueryResult object should return a unique key for the
                      dictionary. Defaults to using .id of the result.

    This function enables access of QueryResult objects from a single search
    output file using its identifier.

    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('Blast/wnts.xml', 'blast-xml')
    >>> search_dict = SearchIO.to_dict(qresults)
    >>> list(search_dict)
    ['gi|195230749:301-1383', 'gi|325053704:108-1166', ..., 'gi|53729353:216-1313']
    >>> search_dict['gi|156630997:105-1160']
    QueryResult(id='gi|156630997:105-1160', 5 hits)

    By default, the dictionary key is the QueryResult's string ID. This may be
    changed by supplying a callback function that returns the desired identifier.
    Here is an example using a function that removes the 'gi|' part in the
    beginning of the QueryResult ID.

    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('Blast/wnts.xml', 'blast-xml')
    >>> key_func = lambda qresult: qresult.id.split('|')[1]
    >>> search_dict = SearchIO.to_dict(qresults, key_func)
    >>> list(search_dict)
    ['195230749:301-1383', '325053704:108-1166', ..., '53729353:216-1313']
    >>> search_dict['156630997:105-1160']
    QueryResult(id='gi|156630997:105-1160', 5 hits)

    Note that the callback function does not change the QueryResult's ID value.
    It only changes the key value used to retrieve the associated QueryResult.

    As this function loads all QueryResult objects into memory, it may be
    unsuitable for dealing with files containing many queries. In that case, it
    is recommended that you use either ``index`` or ``index_db``.

    Since Python 3.7, the default dict class maintains key order, meaning
    this dictionary will reflect the order of records given to it. For
    CPython and PyPy, this was already implemented for Python 3.6, so
    effectively you can always assume the record order is preserved.
    """

    def _default_key_function(rec):
        return rec.id

    if key_function is None:
        key_function = _default_key_function

    qdict = {}
    for qresult in qresults:
        key = key_function(qresult)
        if key in qdict:
            raise ValueError("Duplicate key %r" % key)
        qdict[key] = qresult
    return qdict


def index(filename, format=None, key_function=None, **kwargs):
    """Indexes a search output file and returns a dictionary-like object.

     - filename     - string giving name of file to be indexed
     - format       - Lower case string denoting one of the supported formats.
     - key_function - Optional callback function which when given a
                      QueryResult should return a unique key for the dictionary.
     - kwargs       - Format-specific keyword arguments.

    Index returns a pseudo-dictionary object with QueryResult objects as its
    values and a string identifier as its keys. The function is mainly useful
    for dealing with large search output files, as it enables access to any
    given QueryResult object much faster than using parse or read.

    Index works by storing in-memory the start locations of all queries in a
    file. When a user requests access to the query, this function will jump
    to its start position, parse the whole query, and return it as a
    QueryResult object:

    >>> from Bio import SearchIO
    >>> search_idx = SearchIO.index('Blast/wnts.xml', 'blast-xml')
    >>> search_idx
    SearchIO.index('Blast/wnts.xml', 'blast-xml', key_function=None)
    >>> sorted(search_idx)
    ['gi|156630997:105-1160', 'gi|195230749:301-1383', ..., 'gi|53729353:216-1313']
    >>> search_idx['gi|195230749:301-1383']
    QueryResult(id='gi|195230749:301-1383', 5 hits)
    >>> search_idx.close()

    If the file is BGZF compressed, this is detected automatically. Ordinary
    GZIP files are not supported:

    >>> from Bio import SearchIO
    >>> search_idx = SearchIO.index('Blast/wnts.xml.bgz', 'blast-xml')
    >>> search_idx
    SearchIO.index('Blast/wnts.xml.bgz', 'blast-xml', key_function=None)
    >>> search_idx['gi|195230749:301-1383']
    QueryResult(id='gi|195230749:301-1383', 5 hits)
    >>> search_idx.close()

    You can supply a custom callback function to alter the default identifier
    string. This function should accept as its input the QueryResult ID string
    and return a modified version of it.

    >>> from Bio import SearchIO
    >>> key_func = lambda id: id.split('|')[1]
    >>> search_idx = SearchIO.index('Blast/wnts.xml', 'blast-xml', key_func)
    >>> search_idx
    SearchIO.index('Blast/wnts.xml', 'blast-xml', key_function=<function <lambda> at ...>)
    >>> sorted(search_idx)
    ['156630997:105-1160', ..., '371502086:108-1205', '53729353:216-1313']
    >>> search_idx['156630997:105-1160']
    QueryResult(id='gi|156630997:105-1160', 5 hits)
    >>> search_idx.close()

    Note that the callback function does not change the QueryResult's ID value.
    It only changes the key value used to retrieve the associated QueryResult.

    """
    if not isinstance(filename, str):
        raise TypeError("Need a filename (not a handle)")

    from Bio.File import _IndexedSeqFileDict

    proxy_class = get_processor(format, _INDEXER_MAP)
    repr = f"SearchIO.index({filename!r}, {format!r}, key_function={key_function!r})"
    return _IndexedSeqFileDict(
        proxy_class(filename, **kwargs), key_function, repr, "QueryResult"
    )


def index_db(index_filename, filenames=None, format=None, key_function=None, **kwargs):
    """Indexes several search output files into an SQLite database.

     - index_filename - The SQLite filename.
     - filenames    - List of strings specifying file(s) to be indexed, or when
                      indexing a single file this can be given as a string.
                      (optional if reloading an existing index, but must match)
     - format       - Lower case string denoting one of the supported formats.
                      (optional if reloading an existing index, but must match)
     - key_function - Optional callback function which when given a
                      QueryResult identifier string should return a unique
                      key for the dictionary.
     - kwargs       - Format-specific keyword arguments.

    The ``index_db`` function is similar to ``index`` in that it indexes the start
    position of all queries from search output files. The main difference is
    instead of storing these indices in-memory, they are written to disk as an
    SQLite database file. This allows the indices to persist between Python
    sessions. This enables access to any queries in the file without any
    indexing overhead, provided it has been indexed at least once.

    >>> from Bio import SearchIO
    >>> idx_filename = ":memory:" # Use a real filename, this is in RAM only!
    >>> db_idx = SearchIO.index_db(idx_filename, 'Blast/mirna.xml', 'blast-xml')
    >>> sorted(db_idx)
    ['33211', '33212', '33213']
    >>> db_idx['33212']
    QueryResult(id='33212', 44 hits)
    >>> db_idx.close()

    ``index_db`` can also index multiple files and store them in the same
    database, making it easier to group multiple search files and access them
    from a single interface.

    >>> from Bio import SearchIO
    >>> idx_filename = ":memory:" # Use a real filename, this is in RAM only!
    >>> files = ['Blast/mirna.xml', 'Blast/wnts.xml']
    >>> db_idx = SearchIO.index_db(idx_filename, files, 'blast-xml')
    >>> sorted(db_idx)
    ['33211', '33212', '33213', 'gi|156630997:105-1160', ..., 'gi|53729353:216-1313']
    >>> db_idx['33212']
    QueryResult(id='33212', 44 hits)
    >>> db_idx.close()

    One common example where this is helpful is if you had a large set of
    query sequences (say ten thousand) which you split into ten query files
    of one thousand sequences each in order to run as ten separate BLAST jobs
    on a cluster. You could use ``index_db`` to index the ten BLAST output
    files together for seamless access to all the results as one dictionary.

    Note that ':memory:' rather than an index filename tells SQLite to hold
    the index database in memory. This is useful for quick tests, but using
    the Bio.SearchIO.index(...) function instead would use less memory.

    BGZF compressed files are supported, and detected automatically. Ordinary
    GZIP compressed files are not supported.

    See also Bio.SearchIO.index(), Bio.SearchIO.to_dict(), and the Python module
    glob which is useful for building lists of files.
    """
    # cast filenames to list if it's a string
    # (can we check if it's a string or a generator?)
    if isinstance(filenames, str):
        filenames = [filenames]

    from Bio.File import _SQLiteManySeqFilesDict

    repr = f"SearchIO.index_db({index_filename!r}, filenames={filenames!r}, {format!r}, key_function={key_function!r})"

    def proxy_factory(format, filename=None):
        """Given a filename returns proxy object, else boolean if format OK."""
        if filename:
            return get_processor(format, _INDEXER_MAP)(filename, **kwargs)
        else:
            return format in _INDEXER_MAP

    return _SQLiteManySeqFilesDict(
        index_filename, filenames, proxy_factory, format, key_function, repr
    )


def write(qresults, handle, format=None, **kwargs):
    """Write QueryResult objects to a file in the given format.

     - qresults - An iterator returning QueryResult objects or a single
                  QueryResult object.
     - handle   - Handle to the file, or the filename as a string.
     - format   - Lower case string denoting one of the supported formats.
     - kwargs   - Format-specific keyword arguments.

    The ``write`` function writes QueryResult object(s) into the given output
    handle / filename. You can supply it with a single QueryResult object or an
    iterable returning one or more QueryResult objects. In both cases, the
    function will return a tuple of four values: the number of QueryResult, Hit,
    HSP, and HSPFragment objects it writes to the output file::

        from Bio import SearchIO
        qresults = SearchIO.parse('Blast/mirna.xml', 'blast-xml')
        SearchIO.write(qresults, 'results.tab', 'blast-tab')
        <stdout> (3, 239, 277, 277)

    The output of different formats may be adjusted using the format-specific
    keyword arguments. Here is an example that writes BLAT PSL output file with
    a header::

        from Bio import SearchIO
        qresults = SearchIO.parse('Blat/psl_34_001.psl', 'blat-psl')
        SearchIO.write(qresults, 'results.tab', 'blat-psl', header=True)
        <stdout> (2, 13, 22, 26)

    """
    # turn qresults into an iterator if it's a single QueryResult object
    if isinstance(qresults, QueryResult):
        qresults = iter([qresults])
    else:
        qresults = iter(qresults)

    # get the writer object and do error checking
    writer_class = get_processor(format, _WRITER_MAP)

    # write to the handle
    with as_handle(handle, "w") as target_file:
        writer = writer_class(target_file, **kwargs)
        # count how many qresults, hits, and hsps
        qresult_count, hit_count, hsp_count, frag_count = writer.write_file(qresults)

    return qresult_count, hit_count, hsp_count, frag_count


def convert(in_file, in_format, out_file, out_format, in_kwargs=None, out_kwargs=None):
    """Convert between two search output formats, return number of records.

     - in_file    - Handle to the input file, or the filename as string.
     - in_format  - Lower case string denoting the format of the input file.
     - out_file   - Handle to the output file, or the filename as string.
     - out_format - Lower case string denoting the format of the output file.
     - in_kwargs  - Dictionary of keyword arguments for the input function.
     - out_kwargs - Dictionary of keyword arguments for the output function.

    The convert function is a shortcut function for ``parse`` and ``write``. It has
    the same return type as ``write``. Format-specific arguments may be passed to
    the convert function, but only as dictionaries.

    Here is an example of using ``convert`` to convert from a BLAST+ XML file
    into a tabular file with comments::

        from Bio import SearchIO
        in_file = 'Blast/mirna.xml'
        in_fmt = 'blast-xml'
        out_file = 'results.tab'
        out_fmt = 'blast-tab'
        out_kwarg = {'comments': True}
        SearchIO.convert(in_file, in_fmt, out_file, out_fmt, out_kwargs=out_kwarg)
        <stdout> (3, 239, 277, 277)

    Given that different search output file provide different statistics and
    different level of details, the convert function is limited only to
    converting formats that have the same statistics and for conversion to
    formats with the same level of detail, or less.

    For example, converting from a BLAST+ XML output to a HMMER table file
    is not possible, as these are two search programs with different kinds of
    statistics. In theory, you may provide the necessary values required by the
    HMMER table file (e.g. conditional e-values, envelope coordinates, etc).
    However, these values are likely to hold little meaning as they are not true
    HMMER-computed values.

    Another example is converting from BLAST+ XML to BLAST+ tabular file. This
    is possible, as BLAST+ XML provide all the values necessary to create a
    BLAST+ tabular file. However, the reverse conversion may not be possible.
    There are more details covered in the XML file that are not found in a
    tabular file (e.g. the lambda and kappa values)

    """
    if in_kwargs is None:
        in_kwargs = {}
    if out_kwargs is None:
        out_kwargs = {}

    qresults = parse(in_file, in_format, **in_kwargs)
    return write(qresults, out_file, out_format, **out_kwargs)


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
