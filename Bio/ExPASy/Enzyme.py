# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# Copyright 2009 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with the enzyme.dat file from
Enzyme.
http://www.expasy.ch/enzyme/

Tested with the release of 03-Mar-2009.

Functions:
read       Reads a file containing one ENZYME entry
parse      Reads a file containing multiple ENZYME entries

Classes:
Record     Holds ENZYME data.

"""

def parse(handle):
    """Parse ENZYME records.

    This function is for parsing ENZYME files containing multiple
    records.

    handle   - handle to the file."""

    while True:
        record = __read(handle)
        if not record:
            break
        yield record

def read(handle):
    """Read one ENZYME record.

    This function is for parsing ENZYME files containing
    exactly one record.

    handle   - handle to the file."""

    record = __read(handle)
    # We should have reached the end of the record by now
    remainder = handle.read()
    if remainder:
        raise ValueError("More than one ENZYME record found")
    return record


class Record(dict):
    """\
Holds information from an ExPASy ENZYME record as a Python dictionary.

Each record contains the following keys:
    ID: EC number
    DE: Recommended name
    AN: Alternative names (if any)
    CA: Catalytic activity
    CF: Cofactors (if any)
    PR: Pointers to the Prosite documentation entrie(s) that
        correspond to the enzyme (if any)
    DR: Pointers to the Swiss-Prot protein sequence entrie(s)
        that correspond to the enzyme (if any)
    CC: Comments
"""

    def __init__(self):
        dict.__init__(self)
        self["ID"] = ''
        self["DE"] = ''
        self["AN"] = []
        self["CA"] = ''
        self["CF"] = ''
        self["CC"] = []   # one comment per line
        self["PR"] = []
        self["DR"] = []
    
    def __repr__(self):
        if self["ID"]:
            if self["DE"]:
                return "%s (%s, %s)" % (self.__class__.__name__, 
                                        self["ID"], self["DE"])
            else:
                return "%s (%s)" % (self.__class__.__name__, 
                                       self["ID"])
        else:
            return "%s ( )" % (self.__class__.__name__)
            
    def __str__(self):
        output = "ID: " + self["ID"]
        output += " DE: " + self["DE"]
        output += " AN: " + repr(self["AN"])
        output += " CA: '" + self["CA"] + "'"
        output += " CF: " + self["CF"]
        output += " CC: " + repr(self["CC"])
        output += " PR: " + repr(self["PR"])
        output += " DR: %d Records" % len(self["DR"])
        return output

# Everything below is private

def __read(handle):
    record = None
    for line in handle:
        key, value = line[:2], line[5:].rstrip()
        if key=="ID":
            record = Record()
            record["ID"] = value
        elif key=="DE":
            record["DE"]+=value
        elif key=="AN":
            if record["AN"] and not record["AN"][-1].endswith("."):
                record["AN"][-1] += " " + value
            else:
                record["AN"].append(value)
        elif key=="CA":
            record["CA"] += value
        elif key=="DR":
            pair_data = value.rstrip(";").split(';')
            for pair in pair_data:
                t1, t2 = pair.split(',')
                row = [t1.strip(), t2.strip()]
                record["DR"].append(row)
        elif key=="CF":
            if record["CF"]:
                record["CF"] += " " + value
            else:
                record["CF"] = value
        elif key=="PR":
            assert value.startswith("PROSITE; ")
            value = value[9:].rstrip(";")
            record["PR"].append(value)
        elif key=='CC':
            if value.startswith("-!- "):
                record["CC"].append(value[4:])
            elif value.startswith("    ") and record["CC"]:
                record["CC"][-1] += value[3:]
            # copyright notice is silently skipped
        elif key=="//":
            if record:
                return record
            else: # This was the copyright notice
                continue
    if record:
        raise ValueError("Unexpected end of stream")
