# Copyright 2001 by Tarjei Mikkelsen. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to import KEGG Pathway maps for use with
the Biopython Pathway module.

The pathway maps are in the format:

RXXXXX:[X.X.X.X:] A + 2 B <=> C
RXXXXX:[X.X.X.X:] 3C <=> 2 D + E
...

where RXXXXX is a five-digit reaction id, and X.X.X.X is the optional
EC number of the enzyme that catalyze the reaction.


Classes:
Iterator             -- Iterates through a file of map file.
Parser               -- Parses KEGG map lines into Pathway.Reaction objects.

_Scanner
_Consumer

"""

# XML from python
from xml.sax import handler

# Martel
import Martel
from Martel import RecordReader

# other Biopython stuff
from Bio import File
from Bio.Pathway import Reaction
from Bio.ParserSupport import AbstractConsumer
from Bio.ParserSupport import EventGenerator

import map_format

class Iterator:
    """Iterator to read a file of KEGG reactions one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle to a map file to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self._reader = RecordReader.CountLines(handle, 1)          
        self._parser = parser

    def next(self):
        """Return the next Pathway.Reaction object from the handle.

        Will return None if we ran out of records.
        """
        data = self._reader.next()
        if self._parser is not None:
            if data:
                return self._parser.parse(File.StringHandle(data))
        return data
    
    def __iter__(self):
        return iter(self.next, None)


class Parser:
    """Parse KEGG Map files into Pathway.Reaction objects
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
        """Parse the specified handle into a Pathway.Reaction object.
        """
        self._consumer = _Consumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data


class _Consumer(AbstractConsumer):
    """Create a Pathway.Reaction object from scanner events.

    """
    def __init__(self):
        self.data = None
        self._reaction_id  = ""
        self._ec_id        = ""
        self._stoch_in     = []
        self._substrates   = []
        self._stoch_out    = []
        self._products     = []
        self._arrow        = 0
        self._stoch        = 0
    def stoch(self, stoch):
        if not self._arrow:
            self._stoch_in.append(int(stoch[0]))
        else:
            self._stoch_out.append(int(stoch[0]))
        self._stoch = 1
    def reactant(self, reactants):
        if not self._arrow:
            self._substrates.extend(reactants)
            self._stoch_in.extend([1] * (len(reactants) - self._stoch))
        else:
            self._products.extend(reactants)
            self._stoch_out.extend([1] * (len(reactants) - self._stoch))
        self._stoch = 0
    def arrows(self, arrows):
        self._arrow = 1
    def reaction_id(self, reaction_id):
        self._reaction_id = reaction_id[0]
    def ec_id(self, ec_id):
        self._ec_id = ec_id[0]
    def end(self, end):
        reactants = {}
        for i in range(len(self._products)):
            reactants[self._products[i]] = self._stoch_out[i]
        for i in range(len(self._substrates)):
            reactants[self._substrates[i]] = - self._stoch_in[i]
        self.data  = Reaction(reactants, [(self._ec_id,)], 1, self._reaction_id)


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
        self.interest_tags = ["stoch", "reactant", "arrows", "reaction_id",
                              "ec_id", "end"]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names(map_format.reaction,
                                         self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        """Feed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        self._parser.setContentHandler(EventGenerator(consumer,
                                                      self.interest_tags))
        self._parser.setErrorHandler(handler.ErrorHandler())
        self._parser.parseFile(handle)





