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
You need to download the Cellosaurus database for this examples to
run, e.g. from ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt

    >> from Bio.ExPASy import cellosaurus
    >> with open('cellosaurus.txt') as handle:
    ...    records = cellosaurus.parse(handle)
    ...    for record in records:
    ...        if 'Homo sapiens' in record['OX'][0]:
    ...            print(record['ID'])
    ...
    #15310-LN
    #W7079
    (L)PC6
    00136
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
        if self["ID"]:
            if self["AC"]:
                return "%s (%s, %s)" % (self.__class__.__name__, self["ID"], self["AC"])
            else:
                return "%s (%s)" % (self.__class__.__name__, self["ID"])
        else:
            return "%s ( )" % (self.__class__.__name__)

    def __str__(self):
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
