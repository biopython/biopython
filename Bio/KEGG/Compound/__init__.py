# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""KEGG Compund package

This module provides code to work with the KEGG Ligand/Compund database.


Classes:
Record
Iterator
Parser

_Scanner
_Consumer

"""
__all__ = ['compound_format']

# XML from python
from xml.sax import handler

# Martel
import Martel
from Martel import RecordReader

# other Biopython stuff
from Bio import File
from Bio.KEGG import _write_kegg
from Bio.KEGG import _wrap_kegg
from Bio.ParserSupport import AbstractConsumer
from Bio.ParserSupport import EventGenerator

import compound_format

# Set up line wrapping rules (see Bio.KEGG._wrap_kegg)
name_wrap = [0, "",
             (" ","$",1,1),
             ("-","$",1,1)]
id_wrap = lambda indent : [indent, "",
                           (" ","",1,0)]
struct_wrap = lambda indent : [indent, "",
                               ("  ","",1,1)]

class Record:
    """Holds info from a KEGG Ligand/Compound record.

    Members:
    entry       The entry identifier.
    name        A list of the compund names.
    formula     The chemical formula for the compound 
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
                           map(lambda l:
                               _wrap_kegg(l, wrap_rule = name_wrap),
                               self.name))
    def _formula(self):
        return _write_kegg("FORMULA",
                           [self.formula])
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
            t = entry[0] + " (" + entry[1] + ")"
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


class Iterator:
    """Iterator to read a file of KEGG Ligand/Compund entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with Compound entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self._reader = RecordReader.EndsWith(handle, "///")          
        self._parser = parser

    def next(self):
        """Return the next KEGG Ligand/Compound record from the handle.

        Will return None if we ran out of records.
        """
        data = self._reader.next()
        if self._parser is not None:
            if data:
                return self._parser.parse(File.StringHandle(data))
        return data


class Parser:
    """Parse KEGG Ligand/Compound files into Record objects
    """
    def __init__(self, debug_level = 0):
        """Initialize the parser.

        Arguments:
        o debug_level - An optional argument that species the amount of
        debugging information Martel should spit out. By default we have
        no debugging info (the fastest way to do things), but if you want
        you can set this as high as two and see exactly where a parse fails.
        """
        self._scanner = _Scanner(debug_level)

    def parse(self, handle):
        """Parse the specified handle into a KEGG Ligand/Coumpund record.
        """
        self._consumer = _Consumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data


class _Consumer(AbstractConsumer):
    """Create a KEGG Ligand/Coumpund Record from scanner events.

    """
    def __init__(self):
        self.data = Record()
        self._current_enzyme_id    = ""
        self._current_path         = []
        self._current_dblinks      = []
        self._current_structure_db = ""
    def _unwrap(self, data, add_space = 0):
        lines = data.split("\n")
        if len(lines) == 1:
            return data
        else:
            s = ""
            for l in lines:
                l = l.lstrip()
                if add_space and l[0] != "$" and s[-1] != " ":
                    l = " " + l
                elif l[0] == "$":
                    l = l[1:]
                s = s + l
        return s
    def entry(self, entry):
        self.data.entry = entry[0]
    def name(self, name):
        self.data.name = map(self._unwrap, name)
    def formula(self, formula):
        self.data.formula = formula[0]
    def enzyme_id(self, enzyme_id):
        self._current_enzyme_id = enzyme_id[0]
    def enzyme_role(self, enzyme_role):
        self.data.enzyme.append((self._current_enzyme_id, enzyme_role[0]))
        self._current_enzyme_id = ""
    def pathway_db(self, pathway_db):
        self._current_path.append(pathway_db[0])
    def pathway_id(self, pathway_id):
        self._current_path.append(pathway_id[0])
    def pathway_desc(self, pathway_desc):
        self._current_path.append(" ".join(pathway_desc))
        self.data.pathway.append(tuple(self._current_path))
        self._current_path = []
    def structure_db(self, structure_db):
        self._current_structure_db = structure_db[0]
    def structure_id(self, structure_id):
        self.data.structures.append((self._current_structure_db, structure_id))
        self._current_structure_db = ""
    def dblinks_db(self, dblinks_db):
        self._current_dblinks.append(dblinks_db[0])
    def dblinks_id(self, dblinks_id):
        if not self._current_dblinks == []:
            self._current_dblinks.append(dblinks_id)
            self.data.dblinks.append(tuple(self._current_dblinks))
            self._current_dblinks = []
    def record_end(self, end):
        pass


def _strip(line_list):
    """Combine multiple lines of content separated by spaces.

    This function is used by the EventGenerator callback function to
    combine multiple lines of information. The lines are stripped to
    remove whitespace.
    """
    # first strip out extra whitespace
    return [x.strip() for x in line_list]
      

class _Scanner:
    """Start up Martel to do the scanning of the file.

    This initialzes the Martel based parser and connects it to a handler
    that will generate events for a Consumer.
    """
    def __init__(self, debug = 0):
        """Initialize the scanner by setting up our caches.

        Arguments:
        o debug - The level of debugging that the parser should display.
                  Level 0 is no debugging, Level 2 displays the most
                  debugging info (but is much slower).
                  See Martel documentation for more info on this.
        """
        # a listing of all tags we are interested in scanning for
        # in the Martel parser
        self.interest_tags = ["entry", "name", "formula",
                              "pathway_db", "pathway_id", "pathway_desc",
                              "enzyme_id", "enzyme_role",
                              "structure_db", "structure_id",
                              "dblinks_db", "dblinks_id",
                              "record_end"]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names(compound_format.record,
                                         self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        """Feeed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        self._parser.setContentHandler(EventGenerator(consumer,
                                                      self.interest_tags,
                                                      _strip))
        self._parser.setErrorHandler(handler.ErrorHandler())
        self._parser.parseFile(handle)


