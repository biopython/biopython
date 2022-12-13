# Copyright 2016 by Stephen Marshall.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parser for the cellosaurus.txt file from ExPASy.

See https://web.expasy.org/cellosaurus/

Tested with the release of Version 18 (July 2016).

Functions:
 - read       Reads a file containing one cell line entry
 - parse      Reads a file containing multiple cell line entries

Classes:
 - Record     Holds cell line data.

Examples
--------
This example downloads the Cellosaurus database and parses it. Note that
urlopen returns a stream of bytes, while the parser expects a stream of plain
string, so we use TextIOWrapper to convert bytes to string using the UTF-8
encoding. This is not needed if you download the cellosaurus.txt file in
advance and open it (see the comment below).

    >>> from urllib.request import urlopen
    >>> from io import TextIOWrapper
    >>> from Bio.ExPASy import cellosaurus
    >>> url = "ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt"
    >>> bytestream = urlopen(url)
    >>> textstream = TextIOWrapper(bytestream, "UTF-8")
    >>> # alternatively, use
    >>> # textstream = open("cellosaurus.txt")
    >>> # if you downloaded the cellosaurus.txt file in advance.
    >>> records = cellosaurus.parse(textstream)
    >>> for record in records:
    ...     if 'Homo sapiens' in record['OX'][0]:
    ...         print(record['ID'])  # doctest:+ELLIPSIS
    ...
    #15310-LN
    #W7079
    (L)PC6
    0.5alpha
    ...

"""


def parse(handle):
    """Parse cell line records.

    This function is for parsing cell line files containing multiple
    records.

    Arguments:
     - handle   - handle to the file.

    """
    while True:
        record = __read(handle)
        if not record:
            break
        yield record


def read(handle):
    """Read one cell line record.

    This function is for parsing cell line files containing
    exactly one record.

    Arguments:
     - handle   - handle to the file.

    """
    record = __read(handle)
    # We should have reached the end of the record by now
    remainder = handle.read()
    if remainder:
        raise ValueError("More than one cell line record found")
    return record


class Record(dict):
    """Holds information from an ExPASy Cellosaurus record as a Python dictionary.

    Each record contains the following keys:

     ---------  ---------------------------     ----------------------
     Line code  Content                         Occurrence in an entry
     ---------  ---------------------------     ----------------------
     ID         Identifier (cell line name)     Once; starts an entry
     AC         Accession (CVCL_xxxx)           Once
     AS         Secondary accession number(s)   Optional; once
     SY         Synonyms                        Optional; once
     DR         Cross-references                Optional; once or more
     RX         References identifiers          Optional: once or more
     WW         Web pages                       Optional; once or more
     CC         Comments                        Optional; once or more
     ST         STR profile data                Optional; once or more
     DI         Diseases                        Optional; once or more
     OX         Species of origin               Once or more
     HI         Hierarchy                       Optional; once or more
     OI         Originate from same individual  Optional; once or more
     SX         Sex (gender) of cell            Optional; once
     CA         Category                        Once
     //         Terminator                      Once; ends an entry

    """

    def __init__(self):
        """Initialize the class."""
        dict.__init__(self)
        self["ID"] = ""
        self["AC"] = ""
        self["AS"] = ""
        self["SY"] = ""
        self["DR"] = []
        self["RX"] = []
        self["WW"] = []
        self["CC"] = []
        self["ST"] = []
        self["DI"] = []
        self["OX"] = []
        self["HI"] = []
        self["OI"] = []
        self["SX"] = ""
        self["CA"] = ""

    def __repr__(self):
        """Return the canonical string representation of the Record object."""
        if self["ID"]:
            if self["AC"]:
                return f"{self.__class__.__name__} ({self['ID']}, {self['AC']})"
            else:
                return f"{self.__class__.__name__} ({self['ID']})"
        else:
            return f"{self.__class__.__name__} ( )"

    def __str__(self):
        """Return a readable string representation of the Record object."""
        output = "ID: " + self["ID"]
        output += " AC: " + self["AC"]
        output += " AS: " + self["AS"]
        output += " SY: " + self["SY"]
        output += " DR: " + repr(self["DR"])
        output += " RX: " + repr(self["RX"])
        output += " WW: " + repr(self["WW"])
        output += " CC: " + repr(self["CC"])
        output += " ST: " + repr(self["ST"])
        output += " DI: " + repr(self["DI"])
        output += " OX: " + repr(self["OX"])
        output += " HI: " + repr(self["HI"])
        output += " OI: " + repr(self["OI"])
        output += " SX: " + self["SX"]
        output += " CA: " + self["CA"]
        return output


# Everything below is private


def __read(handle):
    record = None

    for line in handle:

        key, value = line[:2], line[5:].rstrip()
        if key == "ID":
            record = Record()
            record["ID"] = value
        elif key in ["AC", "AS", "SY", "SX", "CA"]:
            record[key] += value
        elif key in [
            "AC",
            "AS",
            "SY",
            "RX",
            "WW",
            "CC",
            "ST",
            "DI",
            "OX",
            "HI",
            "OI",
            "SX",
            "CA",
        ]:
            record[key].append(value)
        elif key == "DR":
            k, v = value.split(";")
            record["DR"].append((k.strip(), v.strip()))
        elif key == "//":
            if record:
                return record
            else:
                continue
    if record:
        raise ValueError("Unexpected end of stream")


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
