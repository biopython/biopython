# Copyright 2014-2016 by Marco Galardini.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
r"""phenotype data input/output.

Input
=====
The main function is Bio.phenotype.parse(...) which takes an input file,
and format string.  This returns an iterator giving PlateRecord objects:

    >>> from Bio import phenotype
    >>> for record in phenotype.parse("phenotype/Plates.csv", "pm-csv"):
    ...     print("%s %i" % (record.id, len(record)))
    ...
    PM01 96
    PM09 96

Note that the parse() function will invoke the relevant parser for the
format with its default settings.  You may want more control, in which case
you need to create a format specific sequence iterator directly.

Input - Single Records
======================
If you expect your file to contain one-and-only-one record, then we provide
the following 'helper' function which will return a single PlateRecord, or
raise an exception if there are no records or more than one record:

    >>> from Bio import phenotype
    >>> record = phenotype.read("phenotype/Plate.json", "pm-json")
    >>> print("%s %i" % (record.id, len(record)))
    PM01 96

This style is useful when you expect a single record only (and would
consider multiple records an error).  For example, when dealing with PM
JSON files saved by the opm library.

However, if you just want the first record from a file containing multiple
record, use the next() function on the iterator:

    >>> from Bio import phenotype
    >>> record = next(phenotype.parse("phenotype/Plates.csv", "pm-csv"))
    >>> print("%s %i" % (record.id, len(record)))
    PM01 96

The above code will work as long as the file contains at least one record.
Note that if there is more than one record, the remaining records will be
silently ignored.

Output
======
Use the function Bio.phenotype.write(...), which takes a complete set of
PlateRecord objects (either as a list, or an iterator), an output file handle
(or in recent versions of Biopython an output filename as a string) and of
course the file format::

        from Bio import phenotype
        records = ...
        phenotype.write(records, "example.json", "pm-json")

Or, using a handle::

        from Bio import phenotype
        records = ...
        with open("example.json", "w") as handle:
           phenotype.write(records, handle, "pm-json")

You are expected to call this function once (with all your records) and if
using a handle, make sure you close it to flush the data to the hard disk.


File Formats
============
When specifying the file format, use lowercase strings.

 - pm-json - Phenotype Microarray plates in JSON format.
 - pm-csv  - Phenotype Microarray plates in CSV format, which is the
             machine vendor format

Note that while Bio.phenotype can read the above file formats, it can only
write in JSON format.
"""

try:
    # Both phen_micro.py and pm_fitting require NumPy, so require NumPy here
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install NumPy if you want to use Bio.phenotype. "
        "See http://www.numpy.org/"
    ) from None

from Bio.File import as_handle
from Bio.phenotype import phen_micro


# Convention for format names is "mainname-format" in lower case.

_FormatToIterator = {
    "pm-csv": phen_micro.CsvIterator,
    "pm-json": phen_micro.JsonIterator,
}

_FormatToWriter = {"pm-json": phen_micro.JsonWriter}


def write(plates, handle, format):
    """Write complete set of PlateRecords to a file.

     - plates    - A list (or iterator) of PlateRecord objects.
     - handle    - File handle object to write to, or filename as string
                   (note older versions of Biopython only took a handle).
     - format    - lower case string describing the file format to write.

    You should close the handle after calling this function.

    Returns the number of records written (as an integer).
    """
    # Try and give helpful error messages:
    if not isinstance(format, str):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError(f"Format string '{format}' should be lower case")

    if isinstance(plates, phen_micro.PlateRecord):
        plates = [plates]

    with as_handle(handle, "w") as fp:
        # Map the file format to a writer class
        if format in _FormatToWriter:
            writer_class = _FormatToWriter[format]
            count = writer_class(plates).write(fp)
        else:
            raise ValueError(f"Unknown format '{format}'")

        if not isinstance(count, int):
            raise TypeError(
                "Internal error - the underlying %s "
                "writer should have returned the record count, not %r" % (format, count)
            )

    return count


def parse(handle, format):
    """Turn a phenotype file into an iterator returning PlateRecords.

     - handle   - handle to the file, or the filename as a string
                  (note older versions of Biopython only took a handle).
     - format   - lower case string describing the file format.

    Typical usage, opening a file to read in, and looping over the record(s):

    >>> from Bio import phenotype
    >>> filename = "phenotype/Plates.csv"
    >>> for record in phenotype.parse(filename, "pm-csv"):
    ...    print("ID %s" % record.id)
    ...    print("Number of wells %i" % len(record))
    ...
    ID PM01
    Number of wells 96
    ID PM09
    Number of wells 96

    Use the Bio.phenotype.read(...) function when you expect a single record
    only.
    """
    # Try and give helpful error messages:
    if not isinstance(format, str):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError(f"Format string '{format}' should be lower case")
    if format not in _FormatToIterator:
        raise ValueError(f"Unknown format '{format}'")

    with as_handle(handle) as fp:
        yield from _FormatToIterator[format](fp)


def read(handle, format):
    """Turn a phenotype file into a single PlateRecord.

     - handle   - handle to the file, or the filename as a string
                  (note older versions of Biopython only took a handle).
     - format   - string describing the file format.

    This function is for use parsing phenotype files containing
    exactly one record.  For example, reading a PM JSON file:

    >>> from Bio import phenotype
    >>> record = phenotype.read("phenotype/Plate.json", "pm-json")
    >>> print("ID %s" % record.id)
    ID PM01
    >>> print("Number of wells %i" % len(record))
    Number of wells 96

    If the handle contains no records, or more than one record,
    an exception is raised.  For example::

        from Bio import phenotype
        record = phenotype.read("plates.csv", "pm-csv")
        Traceback (most recent call last):
        ...
        ValueError: More than one record found in handle

    If however you want the first record from a file containing
    multiple records this function would raise an exception (as
    shown in the example above).  Instead use:

    >>> from Bio import phenotype
    >>> record = next(phenotype.parse("phenotype/Plates.csv", "pm-csv"))
    >>> print("First record's ID %s" % record.id)
    First record's ID PM01

    Use the Bio.phenotype.parse(handle, format) function if you want
    to read multiple records from the handle.
    """
    iterator = parse(handle, format)
    try:
        first = next(iterator)
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No records found in handle")
    try:
        second = next(iterator)
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one record found in handle")
    return first


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
