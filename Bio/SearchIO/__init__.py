# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This is a work in progress, some useful links:
# http://bit.ly/searchio-terms (Spreadsheet of proposed obj/attr name scheme)

"""Biopython representation of sequence homology search program outputs.


Input
=====
The main interface for this module is the Bio.SearchIO.parse(...) function.
It takes two arguments: 1) a file handle or a filename of the input file
(the search output) and 2) a string of one of the supported formats.

    >>> # TODO: example


Input - Search Result with a Single Query
=========================================

    >>> # TODO: example


Input - Search Result with Multiple Queries
===========================================

    >>> # TODO: example


Format Conversion
=================

    >>> # TODO: example


Supported Homology Search Programs
==================================
Below is a list of supported search program output formats. Note that it is
not always enough to specify the program name alone; its format must also be
specified as well (not all output formats from a given program is supported).

 - blast-tab    - BLAST+ tabular output. Variants with or without header
                  comments are both supported, but for the variant without
                  header comments, only the default column ordering is
                  supported.
 - blast-xml    - BLAST+ XML output. 
 - blast-text   - BLAST+ plain text output.
 - blat-psl     - The default output of BLAT (PSL format). Variants with or
                  without header are both supported.
 - fasta        - Bill Pearson's FASTA -m 10 output.
 - hmmer-text   - HMMER regular text output format. Supported HMMER
                  subprograms are hmmerscan and hmmersearch.

"""

# For using with statement in Python 2.5 or Jython
from __future__ import with_statement

__docformat__ = 'epytext en'


from Bio.File import as_handle


# dictionary of supported formats for parse() and read()
_ITERATOR_MAP = {
        'blast-tab': ('BlastIO', 'blast_tabular_iterator'),
        'blast-text': ('BlastIO', 'blast_text_iterator'),
        'blast-xml': ('BlastIO', 'blast_xml_iterator'),
        'blat-psl': ('BlatIO', 'blast_psl_iterator'),
        'fasta': ('FastaIO', 'fasta_m10_iterator'),
        'hmmer-text': ('HmmerIO', 'hmmer_text_iterator'),
}

# dictionary of supported formats for index()
_INDEXER_MAP = {
        'blast-tab': ('BlastIO', 'BlastTabularIndexer'),
        'blast-text': ('BlastIO', 'BlastTextIndexer'),
        'blast-xml': ('BlastIO', 'BlastXmlIndexer'),
        'blat-psl': ('BlatIO', 'BlatPslIndexer'),
        'fasta': ('FastaIO', 'FastaM10Indexer'),
        'hmmer-text': ('HmmerIO', 'HmmerTextIndexer'),
}

# dictionary of supported conversions for convert()
_CONVERSION_MAP = {
        ('blast-xml', 'blast-tab'): ('_convert', '_blastxml_to_blasttab'),
        ('blast-xml', 'blast-text'): ('_convert', '_blastxml_to_blasttext'),
        # ...
}


def _get_handler(format, mapping):
    """Returns the object to handle the given format according to the mapping.

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
            raise ValueError("Format string '%s' should be lower case" % \
                    format)
        else:
            raise ValueError("Unknown format '%s'. Supported formats are "
                    "'%s'" % (format, "', '".join(mapping.keys())))

    mod_name, obj_name = obj_info
    mod = __import__('Bio.SearchIO.%s' % mod_name, fromlist=[1])

    return getattr(mod, obj_name)


def parse(handle, format=None):
    """Turns a search output file into an iterator returning Result objects.

    handle -- Handle to the file, or the filename as a string.
    format -- Lower case string denoting one of the supported formats.

    """
    # check if handle type is correct
    if not isinstance(handle, (basestring, file)):
        raise TypeError("Handle must either be a handle to a file or its "
                "name as string")

    # get the iterator object and do error checking
    iterator = _get_handler(format, _ITERATOR_MAP)

    # and start iterating
    with as_handle(handle) as source_file:
        generator = iterator(source_file)

        for result in generator:
            yield result


def read(handle, format):
    """Turns a search output file into a single Result.

    handle -- Handle to the file, or the filename as a string.
    format -- Lower case string denoting one of the supported formats.

    """
    generator = parse(handle, format)

    try:
        return generator.next()
    except StopIteration:
        raise ValueError("No results found in handle")


def to_dict(queries, key_function=lambda rec: rec.id):
    """Turns a Result iterator or list into a dictionary.

    queries -- Iterator returning Result objects or a list containing
               Result objects.
    key_function -- Optional callback function which when given a
                    Result should return a unique key for the dictionary.

    e.g. key_function = lambda rec : rec.id
    or,  key_function = lambda rec : rec.id.split('|')[0]

    If key_function is ommitted then query.id is used, on the assumption
    that the objects returned are Queries with a unique id. If duplicate
    keys are present, an error will be raised.

    """
    results = {}
    for query in queries:
        key = key_function(query)
        if key in results:
            raise ValueError("Duplicate key '%s'" % key)
        results[key] = query
    return results


def index(handle, format, key_function=lambda rec: rec.id):
    """Indexes a search output file and returns a dictionary-like object.

    - handle    - Handle to the file, or the filename as a string.
    - format    - Lower case string denoting one of the supported formats.
    - key_function - Optional callback function which when given a
                     Result should return a unique key for the dictionary.

    """
    # check if handle type is correct
    if not isinstance(handle, basestring):
        raise TypeError("Handle must be a string of filename")

    # get the indexer object and do error checking
    indexer = _get_handler(format, _INDEXER_MAP)

    from Bio.SearchIO._index import IndexedSearch
    return IndexedSearch(handle, format, indexer, key_function)


def index_db(index_filename, filenames=None, format=None, \
        key_function=lambda rec: rec.id):
    """Indexes several search output files into an SQLite database.

    - index_filename - The SQLite filename.
    - filenames - List of strings specifying file(s) to be indexed, or when
                  indexing a single file this can be given as a string.
                  (optional if reloading an existing index, but must match)
    - format    - Lower case string denoting one of the supported formats.
                  (optional if reloading an existing index, but must match)
    - key_function - Optional callback function which when given a
                     Result identifier string should return a unique
                     key for the dictionary.

    """
    # check index_filename and filenames
    if not isinstance(index_filename, basestring):
        raise TypeError("Need a string as index filename")
    if not isinstance(filenames, (basestring, list, tuple)):
        raise TypeError("Handle must either be a string as filename or a "
                "list of filenames")

    # cast filenames to list if it's a string
    if isinstance(filenames, basestring):
        filenames = [filenames]

    # get the indexer object and do error checking
    indexer = _get_handler(format, _INDEXER_MAP)

    from Bio.SearchIO._index import DbIndexedSearch
    return DbIndexedSearch(index_filename, filenames, format, indexer, \
            key_function)


def convert(in_file, in_format, out_file, out_format):
    """Convert between two search output formats, return number of records.

     - in_file      - Handle to the input file, or the filename as string.
     - in_format    - Lower case string denoting the format of the input file.
     - out_file     - Handle to the output file, or the filename as string.
     - out_format   - Lower case string denoting the format of the output file.

    Conversion is currently limited to operations between search output files
    from one homology search program. An error will be raised if the attempted
    conversion is not supported yet.

    Note that some conversion are lossy, so it can only go one way. For
    example conversion from blast-xml to blast-tab is possible as blast-xml
    has richer information than blast-tab. However the reverse is not possible,
    as blast-tab does not always contain the same information content as
    blast-xml (e.g. the HSP alignment is not always present in blast-tab). 

    """
    convert_pair = (in_format, out_format)

    try:
        converter = _CONVERSION_MAP[convert_pair]
    except KeyError:
        raise ValueError("Conversion from '%s' to '%s' is not supported" % \
                (in_format, out_format))

    converter(in_format, out_format)


def _test():
    """Run the Bio.SearchIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = 'Tests'

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod()
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
