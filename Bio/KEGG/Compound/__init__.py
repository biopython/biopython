# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# Copyright 2007 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Code to work with the KEGG Ligand/Compound database.

Functions:
 - parse - Returns an iterator giving Record objects.

Classes:
 - Record - A representation of a KEGG Ligand/Compound.
"""


from Bio.KEGG import _default_wrap, _struct_wrap, _wrap_kegg, _write_kegg


# Set up line wrapping rules (see Bio.KEGG._wrap_kegg)
name_wrap = [0, "", (" ", "$", 1, 1), ("-", "$", 1, 1)]
id_wrap = _default_wrap
struct_wrap = _struct_wrap


class Record:
    """Holds info from a KEGG Ligand/Compound record.

    Attributes:
     - entry       The entry identifier.
     - name        A list of the compound names.
     - formula     The chemical formula for the compound
     - mass        The molecular weight for the compound
     - pathway     A list of 3-tuples: ('PATH', pathway id, pathway)
     - enzyme      A list of the EC numbers.
     - structures  A list of 2-tuples: (database, list of struct ids)
     - dblinks     A list of 2-tuples: (database, list of link ids)

    """

    def __init__(self):
        """Initialize as new record."""
        self.entry = ""
        self.name = []
        self.formula = ""
        self.mass = ""
        self.pathway = []
        self.enzyme = []
        self.structures = []
        self.dblinks = []

    def __str__(self):
        """Return a string representation of this Record."""
        return (
            self._entry()
            + self._name()
            + self._formula()
            + self._mass()
            + self._pathway()
            + self._enzyme()
            + self._structures()
            + self._dblinks()
            + "///"
        )

    def _entry(self):
        return _write_kegg("ENTRY", [self.entry])

    def _name(self):
        return _write_kegg(
            "NAME", [_wrap_kegg(line, wrap_rule=name_wrap) for line in self.name]
        )

    def _formula(self):
        return _write_kegg("FORMULA", [self.formula])

    def _mass(self):
        return _write_kegg("MASS", [self.mass])

    def _pathway(self):
        s = []
        for entry in self.pathway:
            s.append(entry[0] + "  " + entry[1])
        return _write_kegg(
            "PATHWAY", [_wrap_kegg(line, wrap_rule=id_wrap(16)) for line in s]
        )

    def _enzyme(self):
        return _write_kegg(
            "ENZYME", [_wrap_kegg(line, wrap_rule=name_wrap) for line in self.enzyme]
        )

    def _structures(self):
        s = []
        for entry in self.structures:
            s.append(entry[0] + ": " + "  ".join(entry[1]) + "  ")
        return _write_kegg(
            "STRUCTURES", [_wrap_kegg(line, wrap_rule=struct_wrap(5)) for line in s]
        )

    def _dblinks(self):
        s = []
        for entry in self.dblinks:
            s.append(entry[0] + ": " + " ".join(entry[1]))
        return _write_kegg(
            "DBLINKS", [_wrap_kegg(line, wrap_rule=id_wrap(9)) for line in s]
        )


def parse(handle):
    """Parse a KEGG Ligan/Compound file, returning Record objects.

    This is an iterator function, typically used in a for loop.  For
    example, using one of the example KEGG files in the Biopython
    test suite,

    >>> with open("KEGG/compound.sample") as handle:
    ...     for record in parse(handle):
    ...         print("%s %s" % (record.entry, record.name[0]))
    ...
    C00023 Iron
    C00017 Protein
    C00099 beta-Alanine
    C00294 Inosine
    C00298 Trypsin
    C00348 all-trans-Undecaprenyl phosphate
    C00349 2-Methyl-3-oxopropanoate
    C01386 NH2Mec

    """
    record = Record()
    for line in handle:
        if line[:3] == "///":
            yield record
            record = Record()
            continue
        if line[:12] != "            ":
            keyword = line[:12]
        data = line[12:].strip()
        if keyword == "ENTRY       ":
            words = data.split()
            record.entry = words[0]
        elif keyword == "NAME        ":
            data = data.strip(";")
            record.name.append(data)
        elif keyword == "ENZYME      ":
            while data:
                column = data[:16]
                data = data[16:]
                enzyme = column.strip()
                record.enzyme.append(enzyme)
        elif keyword == "PATHWAY     ":
            map, name = data.split("  ")
            pathway = ("PATH", map, name)
            record.pathway.append(pathway)
        elif keyword == "FORMULA     ":
            record.formula = data
        elif keyword in ("MASS        ", "EXACT_MASS  "):
            record.mass = data
        elif keyword == "DBLINKS     ":
            if ":" in data:
                key, values = data.split(":")
                values = values.split()
                row = (key, values)
                record.dblinks.append(row)
            else:
                row = record.dblinks[-1]
                key, values = row
                values.extend(data.split())
                row = key, values
                record.dblinks[-1] = row


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
