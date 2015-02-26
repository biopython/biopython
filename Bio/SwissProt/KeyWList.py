# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to parse the keywlist.txt file from SwissProt/UniProt

See:
http://www.expasy.ch/sprot/sprot-top.html
ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/docs/keywlist.txt

Classes:

    - Record            Stores the information about one keyword or one category
      in the keywlist.txt file.

Functions:

    - parse             Parses the keywlist.txt file and returns an iterator to
      the records it contains.
"""

from __future__ import print_function

__docformat__ = "restructuredtext en"


class Record(dict):
    """
    This record stores the information of one keyword or category in the
    keywlist.txt as a Python dictionary. The keys in this dictionary are
    the line codes that can appear in the keywlist.txt file::

        ---------  ---------------------------     ----------------------
        Line code  Content                         Occurrence in an entry
        ---------  ---------------------------     ----------------------
        ID         Identifier (keyword)            Once; starts a keyword entry
        IC         Identifier (category)           Once; starts a category entry
        AC         Accession (KW-xxxx)             Once
        DE         Definition                      Once or more
        SY         Synonyms                        Optional; once or more
        GO         Gene ontology (GO) mapping      Optional; once or more
        HI         Hierarchy                       Optional; once or more
        WW         Relevant WWW site               Optional; once or more
        CA         Category                        Once per keyword entry; absent
                                                   in category entries
    """
    def __init__(self):
        dict.__init__(self)
        for keyword in ("DE", "SY", "GO", "HI", "WW"):
            self[keyword] = []


def parse(handle):
    record = Record()
    # First, skip the header - look for start of a record
    for line in handle:
        if line.startswith("ID   "):
            # Looks like there was no header
            record["ID"] = line[5:].strip()
            break
        if line.startswith("IC   "):
            # Looks like there was no header
            record["IC"] = line[5:].strip()
            break
    # Now parse the records
    for line in handle:
        if line.startswith("-------------------------------------"):
            # We have reached the footer
            break
        key = line[:2]
        if key == "//":
            record["DE"] = " ".join(record["DE"])
            record["SY"] = " ".join(record["SY"])
            yield record
            record = Record()
        elif line[2:5] == "   ":
            value = line[5:].strip()
            if key in ("ID", "IC", "AC", "CA"):
                record[key] = value
            elif key in ("DE", "SY", "GO", "HI", "WW"):
                record[key].append(value)
            else:
                print("Ignoring: %s" % line.strip())
    # Read the footer and throw it away
    for line in handle:
        pass
