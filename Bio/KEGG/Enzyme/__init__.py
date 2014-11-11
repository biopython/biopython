# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# Copyright 2007 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to work with the KEGG Enzyme database.

Functions:
parse - Returns an iterator giving Record objects.

Classes:
Record               -- Holds the information from a KEGG Enzyme record.
"""

from __future__ import print_function

from Bio.KEGG import _write_kegg
from Bio.KEGG import _wrap_kegg


# Set up line wrapping rules (see Bio.KEGG._wrap_kegg)
rxn_wrap = [0, "",
            (" + ", "", 1, 1),
            (" = ", "", 1, 1),
            (" ", "$", 1, 1),
            ("-", "$", 1, 1)]
name_wrap = [0, "",
             (" ", "$", 1, 1),
             ("-", "$", 1, 1)]
id_wrap = lambda indent: [indent, "", (" ", "", 1, 0)]
struct_wrap = lambda indent: [indent, "", ("  ", "", 1, 1)]


class Record(object):
    """Holds info from a KEGG Enzyme record.

    Members:
    entry       The EC number (withou the 'EC ').
    name        A list of the enzyme names.
    classname   A list of the classification terms.
    sysname     The systematic name of the enzyme.
    reaction    A list of the reaction description strings.
    substrate   A list of the substrates.
    product     A list of the products.
    inhibitor   A list of the inhibitors.
    cofactor    A list of the cofactors.
    effector    A list of the effectors.
    comment     A list of the comment strings.
    pathway     A list of 3-tuples: (database, id, pathway)
    genes       A list of 2-tuples: (organism, list of gene ids)
    disease     A list of 3-tuples: (database, id, disease)
    structures  A list of 2-tuples: (database, list of struct ids)
    dblinks     A list of 2-tuples: (database, list of db ids)
    """
    def __init__(self):
        """__init___(self)

        Create a new Record.
        """
        self.entry = ""
        self.name = []
        self.classname = []
        self.sysname = []
        self.reaction = []
        self.substrate = []
        self.product = []
        self.inhibitor = []
        self.cofactor = []
        self.effector = []
        self.comment = []
        self.pathway = []
        self.genes = []
        self.disease = []
        self.structures = []
        self.dblinks = []

    def __str__(self):
        """__str__(self)

        Returns a string representation of this Record.
        """
        return self._entry() + \
               self._name() + \
               self._classname() + \
               self._sysname() + \
               self._reaction() + \
               self._substrate() + \
               self._product() + \
               self._inhibitor() + \
               self._cofactor() + \
               self._effector() + \
               self._comment() + \
               self._pathway() + \
               self._genes() + \
               self._disease() + \
               self._structures() + \
               self._dblinks() + \
               "///"

    def _entry(self):
        return _write_kegg("ENTRY",
                           ["EC " + self.entry])

    def _name(self):
        return _write_kegg("NAME",
                           [_wrap_kegg(l, wrap_rule=name_wrap)
                            for l in self.name])

    def _classname(self):
        return _write_kegg("CLASS",
                           self.classname)

    def _sysname(self):
        return _write_kegg("SYSNAME",
                           [_wrap_kegg(l, wrap_rule=name_wrap)
                            for l in self.sysname])

    def _reaction(self):
        return _write_kegg("REACTION",
                           [_wrap_kegg(l, wrap_rule=rxn_wrap)
                            for l in self.reaction])

    def _substrate(self):
        return _write_kegg("SUBSTRATE",
                           [_wrap_kegg(l, wrap_rule=name_wrap)
                            for l in self.substrate])

    def _product(self):
        return _write_kegg("PRODUCT",
                           [_wrap_kegg(l, wrap_rule=name_wrap)
                            for l in self.product])

    def _inhibitor(self):
        return _write_kegg("INHIBITOR",
                           [_wrap_kegg(l, wrap_rule=name_wrap)
                            for l in self.inhibitor])

    def _cofactor(self):
        return _write_kegg("COFACTOR",
                           [_wrap_kegg(l, wrap_rule=name_wrap)
                            for l in self.cofactor])

    def _effector(self):
        return _write_kegg("EFFECTOR",
                           [_wrap_kegg(l, wrap_rule=name_wrap)
                            for l in self.effector])

    def _comment(self):
        return _write_kegg("COMMENT",
                           [_wrap_kegg(l, wrap_rule=id_wrap(0))
                            for l in self.comment])

    def _pathway(self):
        s = []
        for entry in self.pathway:
            s.append(entry[0] + ": " + entry[1] + "  " + entry[2])
        return _write_kegg("PATHWAY",
                           [_wrap_kegg(l, wrap_rule=id_wrap(16))
                            for l in s])

    def _genes(self):
        s = []
        for entry in self.genes:
            s.append(entry[0] + ": " + " ".join(entry[1]))
        return _write_kegg("GENES",
                           [_wrap_kegg(l, wrap_rule=id_wrap(5))
                            for l in s])

    def _disease(self):
        s = []
        for entry in self.disease:
            s.append(entry[0] + ": " + entry[1] + "  " + entry[2])
        return _write_kegg("DISEASE",
                           [_wrap_kegg(l, wrap_rule=id_wrap(13))
                            for l in s])

    def _structures(self):
        s = []
        for entry in self.structures:
            s.append(entry[0] + ": " + "  ".join(entry[1]) + "  ")
        return _write_kegg("STRUCTURES",
                           [_wrap_kegg(l, wrap_rule=struct_wrap(5))
                            for l in s])

    def _dblinks(self):
        # This is a bit of a cheat that won't work if enzyme entries
        # have more than one link id per db id. For now, that's not
        # the case - storing links ids in a list is only to make
        # this class similar to the Compound.Record class.
        s = []
        for entry in self.dblinks:
            s.append(entry[0] + ": " + "  ".join(entry[1]))
        return _write_kegg("DBLINKS", s)


def parse(handle):
    """Parse a KEGG Enzyme file, returning Record objects.

    This is an iterator function, typically used in a for loop.  For
    example, using one of the example KEGG files in the Biopython
    test suite,

    >>> with open("KEGG/enzyme.sample") as handle:
    ...     for record in parse(handle):
    ...         print("%s %s" % (record.entry, record.name[0]))
    ...
    1.1.1.1 Alcohol dehydrogenase
    1.1.1.62 Estradiol 17beta-dehydrogenase
    1.1.1.68 Transferred to EC 1.7.99.5
    1.6.5.3 NADH dehydrogenase (ubiquinone)
    1.14.13.28 3,9-Dihydroxypterocarpan 6a-monooxygenase
    2.4.1.68 Glycoprotein 6-alpha-L-fucosyltransferase
    3.1.1.6 Acetylesterase
    2.7.2.1 Acetate kinase

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
            record.entry = words[1]
        elif keyword=="CLASS       ":
            record.classname.append(data)
        elif keyword=="COFACTOR    ":
            record.cofactor.append(data)
        elif keyword=="COMMENT     ":
            record.comment.append(data)
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
        elif keyword=="DISEASE     ":
            if ":" in data:
                database, data = data.split(":")
                number, name = data.split(None, 1)
                row = (database, number, name)
                record.disease.append(row)
            else:
                row = record.disease[-1]
                database, number, name = row
                name = name + " " + data
                row = database, number, name
                record.disease[-1] = row
        elif keyword=="EFFECTOR    ":
            record.effector.append(data.strip(";"))
        elif keyword=="GENES       ":
            if data[3:5]==': ':
                key, values = data.split(":", 1)
                values = [value.split("(")[0] for value in values.split()]
                row = (key, values)
                record.genes.append(row)
            else:
                row = record.genes[-1]
                key, values = row
                for value in data.split():
                    value = value.split("(")[0]
                    values.append(value)
                row = key, values
                record.genes[-1] = row
        elif keyword=="INHIBITOR   ":
            record.inhibitor.append(data.strip(";"))
        elif keyword=="NAME        ":
            record.name.append(data.strip(";"))
        elif keyword=="PATHWAY     ":
            if data[:5]=='PATH:':
                _, map_num, name = data.split(None, 2)
                pathway = ('PATH', map_num, name)
                record.pathway.append(pathway)
            else:
                ec_num, name = data.split(None, 1)
                pathway = 'PATH', ec_num, name
                record.pathway.append(pathway)
        elif keyword=="PRODUCT     ":
            record.product.append(data.strip(";"))
        elif keyword=="REACTION    ":
            record.reaction.append(data.strip(";"))
        elif keyword=="STRUCTURES  ":
            if data[:4]=='PDB:':
                database = data[:3]
                accessions = data[4:].split()
                row = (database, accessions)
                record.structures.append(row)
            else:
                row = record.structures[-1]
                database, accessions = row
                accessions.extend(data.split())
                row = (database, accessions)
                record.structures[-1] = row
        elif keyword=="SUBSTRATE   ":
            record.substrate.append(data.strip(";"))
        elif keyword=="SYSNAME     ":
            record.sysname.append(data.strip(";"))


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
