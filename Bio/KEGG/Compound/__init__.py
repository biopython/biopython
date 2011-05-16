# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# Copyright 2007 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with the KEGG Ligand/Compound database.

Functions:
parse - Returns an iterator giving Record objects.

Classes:
Record - A representation of a KEGG Ligand/Compound.
"""

# other Biopython stuff
from Bio.KEGG import _write_kegg
from Bio.KEGG import _wrap_kegg


# Set up line wrapping rules (see Bio.KEGG._wrap_kegg)
name_wrap = [0, "",
             (" ","$",1,1),
             ("-","$",1,1)]
id_wrap = lambda indent : [indent, "",
                           (" ","",1,0)]
struct_wrap = lambda indent : [indent, "",
                               ("  ","",1,1)]

class Record(object):
    """Holds info from a KEGG Ligand/Compound record.

    Members:
    entry       The entry identifier.
    name        A list of the compund names.
    formula     The chemical formula for the compound 
    mass        The molecular weight for the compound
    pathway     A list of 3-tuples: (database, id, pathway)
    enzyme      A list of 2-tuples: (enzyme id, role)
    structures  A list of 2-tuples: (database, list of struct ids)
    dblinks     A list of 2-tuples: (database, list of link ids)

    """
    def __init__(self):
        """__init___(self)

        Create a new Record.
        """
        self.entry      = ""
        self.name       = []
        self.formula    = ""
        self.mass       = ""
        self.pathway    = []
        self.enzyme     = []
        self.structures = []
        self.dblinks    = []
    def __str__(self):
        """__str__(self)

        Returns a string representation of this Record.
        """
        return self._entry() + \
               self._name()  + \
               self._formula() + \
               self._mass() + \
               self._pathway() + \
               self._enzyme() + \
               self._structures() + \
               self._dblinks() + \
               "///"
    def _entry(self):
        return _write_kegg("ENTRY",
                           [self.entry])
    def _name(self):
        return _write_kegg("NAME",
                           [_wrap_kegg(l, wrap_rule = name_wrap) \
                            for l in self.name])
    def _formula(self):
        return _write_kegg("FORMULA",
                           [self.formula])

    def _mass(self):
        return _write_kegg("MASS",
                           [self.mass])
    
    def _pathway(self):
        s = []
        for entry in self.pathway:
            s.append(entry[0] + ": " + entry[1] + "  " + entry[2])
        return _write_kegg("PATHWAY",
                           [_wrap_kegg(l, wrap_rule = id_wrap(16)) \
                            for l in s])
    def _enzyme(self):
        s = ""
        for entry in self.enzyme:
            if entry[1]:
                t = entry[0] + " (" + entry[1] + ")"
            else:
                t = entry[0]
            s = s + t.ljust(16)
        return _write_kegg("ENZYME",
                            [_wrap_kegg(s, wrap_rule = id_wrap(0))])
    def _structures(self):
        s = []
        for entry in self.structures:
            s.append(entry[0] + ": " + "  ".join(entry[1]) + "  ")
        return _write_kegg("STRUCTURES",
                           [_wrap_kegg(l, wrap_rule = struct_wrap(5)) \
                            for l in s])
    def _dblinks(self):
        s = []
        for entry in self.dblinks:
            s.append(entry[0] + ": " + " ".join(entry[1]))
        return _write_kegg("DBLINKS",
                           [_wrap_kegg(l, wrap_rule = id_wrap(9)) \
                            for l in s])


def parse(handle):
    """Parse a KEGG Ligan/Compound file, returning Record objects.

    This is an iterator function, typically used in a for loop.  For
    example, using one of the example KEGG files in the Biopython
    test suite,

    >>> handle = open("KEGG/compound.sample")
    >>> for record in parse(handle):
    ...     print record.entry, record.name[0]
    ...
    C00023 Iron
    C00017 Protein
    C00099 beta-Alanine
    C00294 Inosine
    C00298 Trypsin
    C00348 Undecaprenyl phosphate
    C00349 2-Methyl-3-oxopropanoate
    C01386 NH2Mec
    """
    record = Record()
    for line in handle:
        if line[:3]=="///":
            yield record
            record = Record()
            continue
        if line[:12]!="            ":
            keyword = line[:12]
        data = line[12:].strip()
        if keyword=="ENTRY       ":
            words = data.split()
            record.entry = words[0]
        elif keyword=="NAME        ":
            data = data.strip(";")
            record.name.append(data)
        elif keyword=="ENZYME      ":
            while data:
                column = data[:16]
                data = data[16:]
                if '(' in column:
                    entry = column.split()
                    enzyme = (entry[0], entry[1][1:-1])
                else:
                    enzyme = (column.strip(), "")
                record.enzyme.append(enzyme)
        elif keyword=="PATHWAY     ":
            if data[:5]=='PATH:':
                path, map, name = data.split(None,2)
                pathway = (path[:-1], map, name)
                record.pathway.append(pathway)
            else:
                pathway = record.pathway[-1]
                path, map, name = pathway
                name = name + " " + data
                pathway = path, map, name
                record.pathway[-1] = pathway
        elif keyword=="FORMULA     ":
            record.formula = data
        elif keyword=="MASS        ":
            record.mass = data
        elif keyword=="DBLINKS     ":
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

def _test():
    """Run the Bio.KEGG.Compound module's doctests.
    
    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","..","Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","..","Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    _test()
