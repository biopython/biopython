# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Biopython interface for sequence search program outputs.

The SearchIO submodule provides parsers, indexers, and writers for outputs from
various sequence search programs. It provides an API similar to SeqIO and
AlignIO, with the following main functions: `parse`, `read`, `to_dict`, `index`,
`index_db`, `write`, and `convert`.

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
sequences, while MultipleSeqAlignment objects stores the alignment between them.

A detailed description of these objects' features and their example usages are
available in their respective documentations.


Input
=====
The main function used to parse search output files is Bio.SearchIO.parse(...).
This function reads in a given search output file and returns a generator that
yields one QueryResult object per iteration.

`parse` takes two arguments: 1) a file handle or a filename of the input file
(the search output) and 2) a string of one of the supported formats.

    >>> from Bio import SearchIO
    >>> for qresult in SearchIO.parse('Blast/mirna.xml', 'blast-xml'):
    ...     print qresult.id, qresult.description
    ...
    33211 mir_1
    33212 mir_2
    33213 mir_3

SearchIO also provides the Bio.SearchIO.read(...) function, which is intended
for use on search output files containing only one query. `read` returns one
QueryResult object and will raise an exception if the source file contains more
than one queries:

    >>> qresult = SearchIO.read('Blast/xml_2226_blastp_004.xml', 'blast-xml')
    >>> print qresult.id, qresult.description
    ...
    gi|11464971:4-101 pleckstrin [Mus musculus]

For accessing search results of large output files, you can use the indexing
functions Bio.SearchIO.index(...) or Bio.SearchIO.index_db(...). They have a
similar interface to their counterparts in SeqIO and AlignIO, with the addition
of optional, format-specific keyword arguments.


Output
======
SearchIO has writing support for several formats, accessible from the
Bio.SearchIO.write(...) function. This function returns a tuple of four
numbers: the number of QueryResult, Hit, HSP, and HSPFragment written:

    qresults = SearchIO.parse('Blast/mirna.xml', 'blast-xml')
    SearchIO.write(qresults, 'results.tab', 'blast-tab')
    <stdout> (3, 239, 277, 277)

Note that different writers may require different attribute values of the
SearchIO objects. This limits the scope of writable search results to search
results that have the required attributes.

For example, the writer for HMMER domain table output requires
the conditional e-value attribute from each HSP object. If you try to write
to the HMMER domain table format and your HSPs do not have this attribute,
an exception will be raised.


Conversion
==========
SearchIO provides a shortcut function Bio.SearchIO.convert(...) to convert a
given file into another format. Under the hood, `convert` simply parses a given
output file and writes it to another using the `parse` and `write` functions.

Note that the same restrictions found in Bio.SearchIO.write(...) applies to the
convert function as well.


Conventions
===========
The main goal of creating SearchIO is to have a common, easy to use interface
across different search output files. As such, we have also created some
conventions / standards for SearchIO that extend beyond the common object model.
You can expect these to apply to all files parsed by SearchIO, regardless of
their individual formats.

    * Python-style sequence coordinates.

      When storing sequence coordinates (start and end values), SearchIO uses
      the Python-style slice convention: zero-based and half-open intervals. For
      example, if in a BLAST XML output file the start and end coordinates of an
      HSP are 10 and 28, they would become 9 and 28 in SearchIO. The start
      coordinate becomes 9 because Python indices start from zero, while the end
      coordinate remains 28 because Python slices omit the last item in an
      interval.

      Besides giving you the benefits of standardization, this convention also
      makes the coordinates usable for slicing sequences. For example, given a
      full query sequence and the start and end coordinates of an HSP, one can
      use the coordinates to extract part of the query sequence that results in
      the database hit.

      When these objects are written to an output file using
      SearchIO.write(...), the coordinate values are restored to their
      respective format's convention. Using the example above, if the HSP would
      be written to an XML file, the start and end coordinates would become 10
      and 28 again.

    * Sequence coordinates' order.

      Some search output format reverses the start and end coordinate sequences
      according to the sequence's strand. For example, in BLAST plain text
      format if the matching strand lies in the minus orientation, then the
      start coordinate will always be bigger than the end coordinate.

      In SearchIO, the start coordinate is always smaller than the end
      coordinate, regardless of strand. This ensures consistency when using the
      coordinates to slice full sequences.

      Note that this coordinate order convention is only enforced in the
      HSPFragment level. If an HSP object has several HSPFragment objects, each
      individual fragment will conform to this convention. But the order of the
      fragments within the HSP object follows what the search output file uses.

      Similar to the coordinate style convention, the start and end coordinates'
      order are restored to their respective formats when the objects are
      written using Bio.SearchIO.write(...).

    * Frames and strand values.

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
 - exonerate-text   - Exonerate cigar line.
 - fasta-m10        - Bill Pearson's FASTA -m 10 output.
 - hmmer3-text      - HMMER regular text output format. Supported HMMER
                      subprograms are hmmscan, hmmsearch, and phmmer.

Support for parsing:

 - blast-text       - BLAST+ plain text output.

Each of these formats have different keyword arguments available for use with
the main SearchIO functions. More details and examples are available in each
of the format's documentation.

"""

# For using with statement in Python 2.5 or Jython
from __future__ import with_statement

__docformat__ = 'epytext en'

from Bio.File import as_handle
from Bio.SearchIO._objects import QueryResult, Hit, HSP, HSPFragment


__all__ = ['read', 'parse', 'to_dict', 'index', 'index_db', 'write', 'convert']


# dictionary of supported formats for parse() and read()
_ITERATOR_MAP = {
        'blast-tab': ('BlastIO', 'BlastTabParser'),
        'blast-text': ('BlastIO', 'BlastTextParser'),
        'blast-xml': ('BlastIO', 'BlastXmlParser'),
        'blat-psl': ('BlatIO', 'BlatPslParser'),
        'exonerate-cigar': ('ExonerateIO', 'ExonerateCigarParser'),
        'exonerate-text': ('ExonerateIO', 'ExonerateTextParser'),
        'exonerate-vulgar': ('ExonerateIO', 'ExonerateVulgarParser'),
        'fasta-m10': ('FastaIO', 'FastaM10Parser'),
        'hmmer3-text': ('HmmerIO', 'Hmmer3TextParser'),
        'hmmer3-tab': ('HmmerIO', 'Hmmer3TabParser'),
        # for hmmer3-domtab, the specific program is part of the format name
        # as we need it distinguish hit / target coordinates
        'hmmscan3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmhitParser'),
        'hmmsearch3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryParser'),
        'phmmer3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryParser'),
}

# dictionary of supported formats for index()
_INDEXER_MAP = {
        'blast-tab': ('BlastIO', 'BlastTabIndexer'),
        'blast-xml': ('BlastIO', 'BlastXmlIndexer'),
        'blat-psl': ('BlatIO', 'BlatPslIndexer'),
        'exonerate-cigar': ('ExonerateIO', 'ExonerateCigarIndexer'),
        'exonerate-text': ('ExonerateIO', 'ExonerateTextIndexer'),
        'exonerate-vulgar': ('ExonerateIO', 'ExonerateVulgarIndexer'),
        'fasta-m10': ('FastaIO', 'FastaM10Indexer'),
        'hmmer3-text': ('HmmerIO', 'Hmmer3TextIndexer'),
        'hmmer3-tab': ('HmmerIO', 'Hmmer3TabIndexer'),
        'hmmscan3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmhitIndexer'),
        'hmmsearch3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryIndexer'),
        'phmmer3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryIndexer'),
}

# dictionary of supported formats for write()
_WRITER_MAP = {
        'blast-tab': ('BlastIO', 'BlastTabWriter'),
        'blast-xml': ('BlastIO', 'BlastXmlWriter'),
        'blat-psl': ('BlatIO', 'BlatPslWriter'),
        'hmmer3-tab': ('HmmerIO', 'Hmmer3TabWriter'),
        'hmmscan3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmhitWriter'),
        'hmmsearch3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryWriter'),
        'phmmer3-domtab': ('HmmerIO', 'Hmmer3DomtabHmmqueryWriter'),
}


def _get_handler(format, mapping):
    """Returns the object to handle the given format according to the mapping.

    Arguments:
    format -- Lower case string denoting one of the supported formats.
    mapping -- Dictionary of format and object name mapping.

    """
    # map file format to iterator name
    try:
        obj_info = mapping[format]
    except KeyError:
        # handle the errors with helpful messages
        if format is None:
            raise ValueError("Format required (lower case string)")
        elif not isinstance(format, basestring):
            raise TypeError("Need a string for the file format (lower case)")
        elif format != format.lower():
            raise ValueError("Format string '%s' should be lower case" %
                    format)
        elif mapping == _WRITER_MAP and format in _ITERATOR_MAP:
            raise ValueError("Reading format '%s' is supported, but not "
                    "writing" % format)
        else:
            raise ValueError("Unknown format '%s'. Supported formats are "
                    "'%s'" % (format, "', '".join(mapping.keys())))

    mod_name, obj_name = obj_info
    mod = __import__('Bio.SearchIO.%s' % mod_name, fromlist=[''])

    return getattr(mod, obj_name)


def parse(handle, format=None, **kwargs):
    """Turns a search output file into a generator that yields QueryResult
    objects.

    Arguments:
    handle -- Handle to the file, or the filename as a string.
    format -- Lower case string denoting one of the supported formats.
    kwargs -- Format-specific keyword arguments.

    This function is used to iterate over each query in a given search output
    file:

    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('Blast/mirna.xml', 'blast-xml')
    >>> qresults
    <generator object ...>
    >>> for qresult in qresults:
    ...     print "Search %s has %i hits" % (qresult.id, len(qresult))
    ...
    Search 33211 has 100 hits
    Search 33212 has 44 hits
    Search 33213 has 95 hits

    Depending on the file format, parse may also take additional keyword
    argument(s) that modifies the behavior of the format parser. Here is a
    simple example, where the keyword argument enables parsing of a commented
    BLAST tabular output file:

    >>> from Bio import SearchIO
    >>> for qresult in SearchIO.parse('Blast/wnts.tab', 'blast-tab', comments=True):
    ...     print "Search %s has %i hits" % (qresult.id, len(qresult))
    ...
    Search gi|195230749:301-1383 has 5 hits
    Search gi|325053704:108-1166 has 5 hits
    Search gi|156630997:105-1160 has 5 hits
    Search gi|371502086:108-1205 has 5 hits
    Search gi|53729353:216-1313 has 5 hits

    """
    # get the iterator object and do error checking
    iterator = _get_handler(format, _ITERATOR_MAP)

    # and start iterating
    with as_handle(handle, 'rU') as source_file:
        generator = iterator(source_file, **kwargs)

        for qresult in generator:
            yield qresult


def read(handle, format=None, **kwargs):
    """Turns a search output file into a single QueryResult.

    Arguments:
    handle -- Handle to the file, or the filename as a string.
    format -- Lower case string denoting one of the supported formats.
    kwargs -- Format-specific keyword arguments.

    The read function is used for parsing search output files containing exactly
    one query:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('Blast/xml_2226_blastp_004.xml', 'blast-xml')
    >>> print qresult.id, qresult.description
    ...
    gi|11464971:4-101 pleckstrin [Mus musculus]

    If the given handle has no results, an exception will be raised:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('Blast/tab_2226_tblastn_002.txt', 'blast-tab')
    Traceback (most recent call last):
    ...
    ValueError: No query results found in handle

    Similarly, if the given handle has more than one results, an exception will
    also be raised:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('Blast/tab_2226_tblastn_001.txt', 'blast-tab')
    Traceback (most recent call last):
    ...
    ValueError: More than one query results found in handle

    Like parse, read may also accept keyword argument(s) depending on the search
    file format.

    """
    generator = parse(handle, format, **kwargs)

    try:
        first = generator.next()
    except StopIteration:
        raise ValueError("No query results found in handle")
    else:
        try:
            second = generator.next()
        except StopIteration:
            second = None

    if second is not None:
        raise ValueError("More than one query results found in handle")

    return first


def to_dict(qresults, key_function=lambda rec: rec.id):
    """Turns a QueryResult iterator or list into a dictionary.

    Arguments:
    qresults -- Iterable returning QueryResult objects.
    key_function -- Optional callback function which when given a
                    QueryResult object should return a unique key for the
                    dictionary.

    This function enables access of QueryResult objects in a single search
    output file using its identifier.

    >>> from Bio import SearchIO
    >>> qresults = SearchIO.parse('Blast/wnts.xml', 'blast-xml')
    >>> search_dict = SearchIO.to_dict(qresults)
    >>> sorted(search_dict.keys())
    ['gi|156630997:105-1160', ..., 'gi|371502086:108-1205', 'gi|53729353:216-1313']
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
    >>> sorted(search_dict.keys())
    ['156630997:105-1160', ..., '371502086:108-1205', '53729353:216-1313']
    >>> search_dict['156630997:105-1160']
    QueryResult(id='gi|156630997:105-1160', 5 hits)

    Note that the callback function does not change the QueryResult's ID value.
    It only changes the key value used to retrieve the associated QueryResult.

    As this function loads all QueryResult objects into memory, it may be
    unsuitable for dealing with files containing many queries. In that case, it
    is recommended that you use either SearchIO.index(...) or
    SearchIO.index_db(...)

    """
    qdict = {}
    for qresult in qresults:
        key = key_function(qresult)
        if key in qdict:
            raise ValueError("Duplicate key %r" % key)
        qdict[key] = qresult
    return qdict


def index(handle, format=None, key_function=None, **kwargs):
    """Indexes a search output file and returns a dictionary-like object.

    Arguments:
    handle -- Handle to the file, or the filename as a string.
    format -- Lower case string denoting one of the supported formats.
    key_function -- Optional callback function which when given a
                    QueryResult should return a unique key for the dictionary.
    kwargs -- Format-specific keyword arguments.

    Index returns a pseudo-dictionary object with QueryResult objects as its
    values and a string identifier as its keys. The function is mainly useful
    for dealing with large search output files, as it enables access to any
    given QueryResult object much faster than using parse or read.

    Index works by creating storing in-memory the start locations of all
    QueryResult objects in a file. When a user requested access to that
    QueryResult, this function will jump into that position, parse the query
    directly, and returns a QueryResult object:

    >>> from Bio import SearchIO
    >>> search_idx = SearchIO.index('Blast/wnts.xml', 'blast-xml')
    >>> search_idx
    SearchIO.index('Blast/wnts.xml', 'blast-xml', key_function=None)
    >>> sorted(search_idx.keys())
    ['gi|156630997:105-1160', 'gi|195230749:301-1383', ..., 'gi|53729353:216-1313']
    >>> search_idx['gi|195230749:301-1383']
    QueryResult(id='gi|195230749:301-1383', 5 hits)

    You can supply a custom callback function to alter the default identifier
    string. This function should accept as its input the QueryResult ID string
    and returns a modified version of it.

    >>> from Bio import SearchIO
    >>> key_func = lambda id: id.split('|')[1]
    >>> search_idx = SearchIO.index('Blast/wnts.xml', 'blast-xml', key_func)
    >>> search_idx
    SearchIO.index('Blast/wnts.xml', 'blast-xml', key_function=<function <lambda> at ...>)
    >>> sorted(search_idx.keys())
    ['156630997:105-1160', ..., '371502086:108-1205', '53729353:216-1313']
    >>> search_idx['156630997:105-1160']
    QueryResult(id='gi|156630997:105-1160', 5 hits)

    Note that the callback function does not change the QueryResult's ID value.
    It only changes the key value used to retrieve the associated QueryResult.

    """
    # check if handle type is correct
    if not isinstance(handle, basestring):
        raise TypeError("Handle must be a string of filename")

    from Bio.SearchIO._index import _IndexedSearch
    return _IndexedSearch(handle, format, key_function, **kwargs)


def index_db(index_filename, filenames=None, format=None,
        key_function=None, **kwargs):
    """Indexes several search output files into an SQLite database.

    Arguments:
    index_filename -- The SQLite filename.
    filenames -- List of strings specifying file(s) to be indexed, or when
                 indexing a single file this can be given as a string.
                 (optional if reloading an existing index, but must match)
    format -- Lower case string denoting one of the supported formats.
              (optional if reloading an existing index, but must match)
    key_function -- Optional callback function which when given a
                    QueryResult identifier string should return a unique
                    key for the dictionary.
    kwargs -- Format-specific keyword arguments.

    The index_db function is similar to Bio.SearchIO.index(...) in that it
    indexes the start position of all queries from search output files. The main
    difference is instead of storing these indices in-memory, they are written
    into a flat SQLite database. This allows the indices to persist between
    Python sessions, so the QueryResult objects in the source file can be
    accessed without any indexing overhead.

    >>> from Bio import SearchIO
    >>> db_idx = SearchIO.index_db(':memory:', 'Blast/mirna.xml', 'blast-xml')
    >>> sorted(db_idx.keys())
    ['33211', '33212', '33213']
    >>> db_idx['33212']
    QueryResult(id='33212', 44 hits)

    index_db can also index multiple files and store them in the same database,
    making it easier to group multiple search files and access them from a
    single interface.

    >>> from Bio import SearchIO
    >>> files = ['Blast/mirna.xml', 'Blast/wnts.xml']
    >>> db_idx = SearchIO.index_db(':memory:', files, 'blast-xml')
    >>> sorted(db_idx.keys())
    ['33211', '33212', '33213', 'gi|156630997:105-1160', ..., 'gi|53729353:216-1313']
    >>> db_idx['33212']
    QueryResult(id='33212', 44 hits)

    """
    # cast filenames to list if it's a string
    # (can we check if it's a string or a generator?)
    if isinstance(filenames, basestring):
        filenames = [filenames]

    from Bio.SearchIO._index import _DbIndexedSearch
    return _DbIndexedSearch(index_filename, filenames, format, key_function,
            **kwargs)


def write(qresults, handle, format=None, **kwargs):
    """Writes QueryResult objects to a file in the given format.

    Arguments:
    qresults -- An iterator returning QueryResult objects or a single
                QueryResult object.
    handle -- Handle to the file, or the filename as a string.
    format -- Lower case string denoting one of the supported formats.
    kwargs -- Format-specific keyword arguments.

    The write function writes QueryResult object(s) into a the given output
    handle / filename. You can supply it with a single QueryResult object or an
    iterable returning one or more QueryResult objects. In both cases, the
    function will return a tuple of four values: the number of QueryResult, Hit,
    HSP, and HSPFragment objects it writes to the output file.

    from Bio import SearchIO
    qresults = SearchIO.parse('Blast/mirna.xml', 'blast-xml')
    SearchIO.write(qresults, 'results.tab', 'blast-tab')
    <stdout> (3, 239, 277, 277)

    The output of different formats may be adjusted using the format-specific
    keyword arguments. Here is an example that writes BLAT PSL output file with
    a header:

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
    writer_class = _get_handler(format, _WRITER_MAP)

    # write to the handle
    with as_handle(handle, 'w') as target_file:
        writer = writer_class(target_file, **kwargs)
        # count how many qresults, hits, and hsps
        qresult_count, hit_count, hsp_count, frag_count = \
                writer.write_file(qresults)

    return qresult_count, hit_count, hsp_count, frag_count


def convert(in_file, in_format, out_file, out_format, in_kwargs=None,
        out_kwargs=None):
    """Convert between two search output formats, return number of records.

    Arguments:
    in_file -- Handle to the input file, or the filename as string.
    in_format -- Lower case string denoting the format of the input file.
    out_file -- Handle to the output file, or the filename as string.
    out_format -- Lower case string denoting the format of the output file.
    in_kwargs --  Dictionary of keyword arguments for the input function.
    out_kwargs -- Dictionary of keyword arguments for the output function.

    The convert function is a shortcut function for Bio.SearchIO.parse(...)
    and Bio.SearchIO.write(...). It has the same return type as the write
    function. Format-specific arguments may be passed to the convert function,
    but only as dictionaries.

    Here is an example of using convert to convert from a BLAST+ XML file into a
    tabular file with comments:

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
    HMMER table file (e.g. conditional e-values, envelope coordinates). However,
    these values are likely to be meaningless as they are not true
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


def _test():
    """Run the Bio.SearchIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = 'Tests'
    outfiles = ['results.tab', 'search.idx']

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        # check that we're not overwriting any file, as the doctest
        # writes several output files
        for ofile in outfiles:
            assert not os.path.exists(ofile), ofile
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod(optionflags=doctest.ELLIPSIS)
        # delete example from SearchIO.write
        for ofile in outfiles:
            if os.path.exists(ofile):
                os.remove(ofile)
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
