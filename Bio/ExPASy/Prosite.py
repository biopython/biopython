# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# Revisions Copyright 2007 by Peter Cock.  All rights reserved.
# Revisions Copyright 2009 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for the prosite dat file from Prosite at ExPASy.

See https://www.expasy.org/prosite/

Tested with:
 - Release 20.43, 10-Feb-2009
 - Release 2017_03 of 15-Mar-2017.

Functions:
 - read                  Reads a Prosite file containing one Prosite record
 - parse                 Iterates over records in a Prosite file.

Classes:
 - Record                Holds Prosite data.

"""


def parse(handle):
    """Parse Prosite records.

    This function is for parsing Prosite files containing multiple
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
    """Read one Prosite record.

    This function is for parsing Prosite files containing
    exactly one record.

    Arguments:
     - handle   - handle to the file.

    """
    record = __read(handle)
    # We should have reached the end of the record by now
    remainder = handle.read()
    if remainder:
        raise ValueError("More than one Prosite record found")
    return record


class Record:
    """Holds information from a Prosite record.

    Main attributes:
     - name           ID of the record.  e.g. ADH_ZINC
     - type           Type of entry.  e.g. PATTERN, MATRIX, or RULE
     - accession      e.g. PS00387
     - created        Date the entry was created.  (MMM-YYYY for releases
       before January 2017, DD-MMM-YYYY since January 2017)
     - data_update    Date the 'primary' data was last updated.
     - info_update    Date data other than 'primary' data was last updated.
     - pdoc           ID of the PROSITE DOCumentation.
     - description    Free-format description.
     - pattern        The PROSITE pattern.  See docs.
     - matrix         List of strings that describes a matrix entry.
     - rules          List of rule definitions (from RU lines).  (strings)
     - prorules       List of prorules (from PR lines). (strings)

    NUMERICAL RESULTS:
     - nr_sp_release  SwissProt release.
     - nr_sp_seqs     Number of seqs in that release of Swiss-Prot. (int)
     - nr_total       Number of hits in Swiss-Prot.  tuple of (hits, seqs)
     - nr_positive    True positives.  tuple of (hits, seqs)
     - nr_unknown     Could be positives.  tuple of (hits, seqs)
     - nr_false_pos   False positives.  tuple of (hits, seqs)
     - nr_false_neg   False negatives.  (int)
     - nr_partial     False negatives, because they are fragments. (int)

    COMMENTS:
     - cc_taxo_range  Taxonomic range.  See docs for format
     - cc_max_repeat  Maximum number of repetitions in a protein
     - cc_site        Interesting site.  list of tuples (pattern pos, desc.)
     - cc_skip_flag   Can this entry be ignored?
     - cc_matrix_type
     - cc_scaling_db
     - cc_author
     - cc_ft_key
     - cc_ft_desc
     - cc_version     version number (introduced in release 19.0)

    The following are all lists if tuples (swiss-prot accession, swiss-prot name).

    DATA BANK REFERENCES:
     - dr_positive
     - dr_false_neg
     - dr_false_pos
     - dr_potential   Potential hits, but fingerprint region not yet available.
     - dr_unknown     Could possibly belong
     - pdb_structs    List of PDB entries.

    """

    def __init__(self):
        """Initialize the class."""
        self.name = ""
        self.type = ""
        self.accession = ""
        self.created = ""
        self.data_update = ""
        self.info_update = ""
        self.pdoc = ""

        self.description = ""
        self.pattern = ""
        self.matrix = []
        self.rules = []
        self.prorules = []
        self.postprocessing = []

        self.nr_sp_release = ""
        self.nr_sp_seqs = ""
        self.nr_total = (None, None)
        self.nr_positive = (None, None)
        self.nr_unknown = (None, None)
        self.nr_false_pos = (None, None)
        self.nr_false_neg = None
        self.nr_partial = None

        self.cc_taxo_range = ""
        self.cc_max_repeat = ""
        self.cc_site = []
        self.cc_skip_flag = ""

        self.dr_positive = []
        self.dr_false_neg = []
        self.dr_false_pos = []
        self.dr_potential = []
        self.dr_unknown = []

        self.pdb_structs = []


# Everything below are private functions


def __read(handle):
    import re

    record = None
    for line in handle:
        keyword, value = line[:2], line[5:].rstrip()
        if keyword == "ID":
            record = Record()
            cols = value.split("; ")
            if len(cols) != 2:
                raise ValueError(f"I don't understand identification line\n{line}")
            record.name = cols[0]
            record.type = cols[1].rstrip(".")  # don't want '.'
        elif keyword == "AC":
            record.accession = value.rstrip(";")
        elif keyword == "DT":
            # e.g. from January 2017,
            # DT   01-APR-1990 CREATED; 01-APR-1990 DATA UPDATE; 01-APR-1990 INFO UPDATE.
            # Older files had brackets round the date descriptions and used MMM-YYYY
            dates = value.rstrip(".").split("; ")
            if dates[0].endswith((" (CREATED)", " CREATED")):
                # Remove last word
                record.created = dates[0].rsplit(" ", 1)[0]
            else:
                raise ValueError(f"I don't understand date line\n{line}")
            if dates[1].endswith((" (DATA UPDATE)", " DATA UPDATE")):
                # Remove last two words
                record.data_update = dates[1].rsplit(" ", 2)[0]
            else:
                raise ValueError(f"I don't understand date line\n{line}")
            if dates[2].endswith((" (INFO UPDATE)", " INFO UPDATE")):
                # Remove last two words
                record.info_update = dates[2].rsplit(" ", 2)[0]
            else:
                raise ValueError(f"I don't understand date line\n{line}")
        elif keyword == "DE":
            record.description = value
        elif keyword == "PA":
            record.pattern += value
        elif keyword == "MA":
            record.matrix.append(value)
        elif keyword == "PP":
            record.postprocessing.extend(value.split(";"))
        elif keyword == "RU":
            record.rules.append(value)
        elif keyword == "NR":
            cols = value.split(";")
            for col in cols:
                if not col:
                    continue
                qual, data = (word.lstrip() for word in col.split("="))
                if qual == "/RELEASE":
                    release, seqs = data.split(",")
                    record.nr_sp_release = release
                    record.nr_sp_seqs = int(seqs)
                elif qual == "/FALSE_NEG":
                    record.nr_false_neg = int(data)
                elif qual == "/PARTIAL":
                    record.nr_partial = int(data)
                elif qual in ["/TOTAL", "/POSITIVE", "/UNKNOWN", "/FALSE_POS"]:
                    m = re.match(r"(\d+)\((\d+)\)", data)
                    if not m:
                        raise Exception(f"Broken data {data} in comment line\n{line!r}")
                    hits = tuple(map(int, m.groups()))
                    if qual == "/TOTAL":
                        record.nr_total = hits
                    elif qual == "/POSITIVE":
                        record.nr_positive = hits
                    elif qual == "/UNKNOWN":
                        record.nr_unknown = hits
                    elif qual == "/FALSE_POS":
                        record.nr_false_pos = hits
                else:
                    raise ValueError(f"Unknown qual {qual} in comment line\n{line!r}")
        elif keyword == "CC":
            # Expect CC lines like this:
            # CC   /TAXO-RANGE=??EPV; /MAX-REPEAT=2;
            # Can (normally) split on ";" and then on "="
            cols = value.split(";")
            for col in cols:
                if not col or col[:17] == "Automatic scaling":
                    # DNAJ_2 in Release 15 has a non-standard comment line:
                    # CC   Automatic scaling using reversed database
                    # Throw it away.  (Should I keep it?)
                    continue
                if col.count("=") == 0:
                    # Missing qualifier!  Can we recover gracefully?
                    # For example, from Bug 2403, in PS50293 have:
                    # CC /AUTHOR=K_Hofmann; N_Hulo
                    continue
                qual, data = (word.lstrip() for word in col.split("="))
                if qual == "/TAXO-RANGE":
                    record.cc_taxo_range = data
                elif qual == "/MAX-REPEAT":
                    record.cc_max_repeat = data
                elif qual == "/SITE":
                    pos, desc = data.split(",")
                    record.cc_site.append((int(pos), desc))
                elif qual == "/SKIP-FLAG":
                    record.cc_skip_flag = data
                elif qual == "/MATRIX_TYPE":
                    record.cc_matrix_type = data
                elif qual == "/SCALING_DB":
                    record.cc_scaling_db = data
                elif qual == "/AUTHOR":
                    record.cc_author = data
                elif qual == "/FT_KEY":
                    record.cc_ft_key = data
                elif qual == "/FT_DESC":
                    record.cc_ft_desc = data
                elif qual == "/VERSION":
                    record.cc_version = data
                else:
                    raise ValueError(f"Unknown qual {qual} in comment line\n{line!r}")
        elif keyword == "DR":
            refs = value.split(";")
            for ref in refs:
                if not ref:
                    continue
                acc, name, type = (word.strip() for word in ref.split(","))
                if type == "T":
                    record.dr_positive.append((acc, name))
                elif type == "F":
                    record.dr_false_pos.append((acc, name))
                elif type == "N":
                    record.dr_false_neg.append((acc, name))
                elif type == "P":
                    record.dr_potential.append((acc, name))
                elif type == "?":
                    record.dr_unknown.append((acc, name))
                else:
                    raise ValueError(f"I don't understand type flag {type}")
        elif keyword == "3D":
            cols = value.split()
            for id in cols:
                record.pdb_structs.append(id.rstrip(";"))
        elif keyword == "PR":
            rules = value.split(";")
            record.prorules.extend(rules)
        elif keyword == "DO":
            record.pdoc = value.rstrip(";")
        elif keyword == "//":
            if not record:
                # Then this was the copyright statement
                continue
            break
        else:
            raise ValueError(f"Unknown keyword {keyword} found")
    else:
        return
    if not record:
        raise ValueError("Unexpected end of stream.")
    return record
