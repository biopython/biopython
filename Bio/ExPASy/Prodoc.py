# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to work with the prosite.doc file from Prosite.

See https://www.expasy.org/prosite/

Tested with:
 - Release 15.0, July 1998
 - Release 16.0, July 1999
 - Release 20.22, 13 November 2007
 - Release 20.43, 10 February 2009

Functions:
 - read               Read a Prodoc file containing exactly one Prodoc entry.
 - parse              Iterates over entries in a Prodoc file.

Classes:
 - Record             Holds Prodoc data.
 - Reference          Holds data from a Prodoc reference.

"""


def read(handle):
    """Read in a record from a file with exactly one Prodoc record."""
    record = __read(handle)
    # We should have reached the end of the record by now
    line = handle.readline()
    if line:
        raise ValueError("More than one Prodoc record found")
    return record


def parse(handle):
    """Iterate over the records in a Prodoc file."""
    while True:
        record = __read(handle)
        if not record:
            return
        yield record


class Record:
    """Holds information from a Prodoc record.

    Attributes:
     - accession      Accession number of the record.
     - prosite_refs   List of tuples (prosite accession, prosite name).
     - text           Free format text.
     - references     List of reference objects.

    """

    def __init__(self):
        """Initialize the class."""
        self.accession = ""
        self.prosite_refs = []
        self.text = ""
        self.references = []


class Reference:
    """Holds information from a Prodoc citation.

    Attributes:
     - number     Number of the reference. (string)
     - authors    Names of the authors.
     - citation   Describes the citation.

    """

    def __init__(self):
        """Initialize the class."""
        self.number = ""
        self.authors = ""
        self.citation = ""


# Below are private functions


def __read_prosite_reference_line(record, line):
    line = line.rstrip()
    if line[-1] != "}":
        raise ValueError("I don't understand the Prosite reference on line\n%s" % line)
    acc, name = line[1:-1].split("; ")
    record.prosite_refs.append((acc, name))


def __read_text_line(record, line):
    record.text += line
    return True


def __read_reference_start(record, line):
    # Read the references
    reference = Reference()
    reference.number = line[1:3].strip()
    if line[1] == "E":
        # If it's an electronic reference, then the URL is on the
        # line, instead of the author.
        reference.citation = line[4:].strip()
    else:
        reference.authors = line[4:].strip()
    record.references.append(reference)


def __read_reference_line(record, line):
    if not line.strip():
        return False
    reference = record.references[-1]
    if line.startswith("     "):
        if reference.authors[-1] == ",":
            reference.authors += line[4:].rstrip()
        else:
            reference.citation += line[5:]
        return True
    raise Exception("I don't understand the reference line\n%s" % line)


def __read_copyright_line(record, line):
    # Skip the copyright statement
    if line.startswith("+----"):
        return False
    return True


def __read(handle):
    # Skip blank lines between records
    for line in handle:
        line = line.rstrip()
        if line and not line.startswith("//"):
            break
    else:
        return None
    record = Record()
    # Read the accession number
    if not line.startswith("{PDOC"):
        raise ValueError("Line does not start with '{PDOC':\n%s" % line)
    if line[-1] != "}":
        raise ValueError("I don't understand accession line\n%s" % line)
    record.accession = line[1:-1]
    # Read the Prosite references
    for line in handle:
        if line.startswith("{PS"):
            __read_prosite_reference_line(record, line)
        else:
            break
    else:
        raise ValueError("Unexpected end of stream.")
    # Read the actual text
    if not line.startswith("{BEGIN"):
        raise ValueError("Line does not start with '{BEGIN':\n%s" % line)
    read_line = __read_text_line
    for line in handle:
        if line.startswith("{END}"):
            # Clean up the record and return
            for reference in record.references:
                reference.citation = reference.citation.rstrip()
                reference.authors = reference.authors.rstrip()
            return record
        elif line[0] == "[" and line[3] == "]" and line[4] == " ":
            __read_reference_start(record, line)
            read_line = __read_reference_line
        elif line.startswith("+----"):
            read_line = __read_copyright_line
        elif read_line:
            if not read_line(record, line):
                read_line = None
    raise ValueError("Unexpected end of stream.")
