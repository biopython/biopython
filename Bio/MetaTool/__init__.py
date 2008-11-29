# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for output from MetaTool 3.5 (DEPRECATED).

MetaTool is a program which defines metabolic routes within networks.
This parser does not support the current version, MetaTool 5.0.

http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&
list_uids=10222413&dopt=Abstract
"""

import warnings
warnings.warn("Bio.MetalTool was deprecated, as it only supported the obsolete"
              +"MetaTool 3.5 output.")

# standard library
import string

import numpy.oldnumeric.matrix as Matrix

# XML from python 2.0
from xml.sax import handler

# Martel
import Martel
from Martel import RecordReader

from Bio.ParserSupport import EventGenerator
from Bio import File
import metatool_format
import Record

class Iterator:
    """Iterator interface to move over a file of MetaTool entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with Kabat entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self._reader = RecordReader.StartsWith(handle, "METATOOL")
#        self._reader = RecordReader.EndsWith(handle, "RECEND|\n")
        self._parser = parser

    def next(self):
        """Return the next MetaTool record from the handle.

        Will return None if we ran out of records.
        """
        data = self._reader.next()

        if self._parser is not None:
            if data:
                return self._parser.parse(File.StringHandle(data))

        return data
    
    def __iter__(self):
        return iter(self.next, None)

class _RecordConsumer:
    """Create a MetaTool Record object from scanner generated information.
    """
    def __init__(self):
        self.data = Record.Record()
        self.data.internal_metabolites = []
        self.data.external_metabolites = []

    def input_file_name( self, content ):
        self.data.input_file_name = content[ 0 ]

    def input_file_tag( self, content ):
        self.state = "input_file_state"

    def metabolite_count_tag( self, content ):
        self.state = "metabolite_count_state"

    def reaction_count_tag( self, content ):
        self.state = "reaction_count_state"

    def matrix_row( self, matrix_rows ):
        for matrix_row in matrix_rows:
            elements = matrix_row.split()
            vector = []
            for element in elements:
                vector.append( int( element ) )
            self._vectors.append( vector )

    def unbalanced_metabolite( self, content ):
        for metabolite in content:
            self.data.unbalanced_metabolites.append( metabolite )

    def branch_metabolite( self, content ):
        for metabolite in content:
            items = metabolite.split()
            name = items[ 0 ]
            consumed = int( items[ 1 ] )
            built = int( items[ 2 ] )
            vector = items[ 4 ].replace( 'r', '0' )

            vector = vector.replace( 'i', '1' )
            vector = list( vector )
            map( int, vector )
            entry = Record.MetaboliteRole( name, consumed, built, vector )
            self.data.branch_metabolites.append( entry )

    def non_branch_metabolite( self, content ):
        for metabolite in content:
            items = metabolite.split()
            name = items[ 0 ]
            consumed = int( items[ 1 ] )
            built = int( items[ 2 ] )
            vector = items[ 4 ].replace( 'r', '0' )
            vector = vector.replace( 'i', '1' )
            vector = list( vector )
            entry = Record.MetaboliteRole( name, consumed, built, vector )
            self.data.non_branch_metabolites.append( entry )

    def stoichiometric_tag( self, content ):
        self.state = "stoichiometry_state"
        self._vectors = []
        self._enzymes = []
        self._reactions = []

    def kernel_tag( self, kernel_tag ):
        self.state = "kernel_state"
        self._vectors = []
        self._enzymes = []
        self._reactions = []

    def subsets_tag( self, content ):
        self.state = "subsets_state"
        self._vectors = []
        self._enzymes = []
        self._reactions = []

    def reduced_system_tag( self, content ):
        self.state = "reduced_system_state"
        self._vectors = []
        self._enzymes = []
        self._reactions = []

    def convex_basis_tag( self, content ):
        self.state = "convex_basis_state"
        self._vectors = []
        self._enzymes = []
        self._reactions = []

    def conservation_relations_tag( self, content ):
        self.state = "conservation_relations_state"
        self._vectors = []
        self._enzymes = []
        self._reactions = []

    def elementary_modes_tag( self, content ):
        self.state = "elementary_modes_state"
        self._vectors = []
        self._enzymes = []
        self._reactions = []

    def metabolite_line( self, content ):
        self.data.external_metabolites = []
        self.data.internal_metabolites = []
        for metabolite in content:
            items = metabolite.split()
            entry = Record.Metabolite( int( items[ 0 ] ), items[ 2 ] )

            if( items[ 1 ] == "external" ):
                self.data.external_metabolites.append( entry )
            else:
                self.data.internal_metabolites.append( entry )


    def num_int_metabolites( self, content ):
        num_int_metabolites = content[ 0 ]
        self.data.num_int_metabolites = int( num_int_metabolites )

    def num_reactions( self, content ):
        num_reactions = content[ 0 ]
        self.data.num_reactions = int( num_reactions )

    def irreversible_vector( self, content ):
        self._irreversible_vector = content[ 0 ].split()

    def reaction( self, reactions ):
        for reaction in reactions:
            items = reaction.split()
            item = ' '.join( items[ 1: ] )
            self._reactions.append(item)

    def enzyme( self, enzymes ):
        for enzyme in enzymes:
            items = enzyme.split()
            item = ' '.join( items[ 1: ] )
            self._enzymes.append(item)

    def sum_is_constant_line( self, lines ):
        for line in lines:
            items = line.split( ':')
            items = items[ 1 ].split( '=' )
            self.data.sum_is_constant_lines.append( items[ 0 ] )

    def num_rows( self, num_rows ):
        pass

    def num_cols( self, num_cols ):
        pass

    def metabolite_roles( self, content ):
        for metabolite_role in content:
            cols = metabolite_role.split()

    def end_stochiometric( self, content ):
        if( self._vectors != [] ):
            self.data.stochiometric.matrix = Matrix.Matrix( self._vectors  )
        self.data.stochiometric.enzymes = []
        for enzyme in self._enzymes:
            self.data.stochiometric.enzymes.append( enzyme )
        self.data.stochiometric.enzymes = []
        for reaction in self._reactions:
            self.data.stochiometric.reactions.append( reaction )
        for col in self._irreversible_vector:
            self.data.stochiometric.irreversible_vector.append( col )

    def end_kernel( self, content ):
        if( self._vectors != [] ):
            self.data.kernel.matrix = Matrix.Matrix( self._vectors )
        self.data.kernel.enzymes = []
        for enzyme in self._enzymes:
            self.data.kernel.enzymes.append( enzyme )
        for reaction in self._reactions:
            self.data.kernel.reactions.append( reaction )

    def end_subsets( self, content ):
        if( self._vectors != [] ):
            self.data.subsets.matrix = Matrix.Matrix( self._vectors )
        self.data.subsets.enzymes = []
        for enzyme in self._enzymes:
            self.data.subsets.enzymes.append( enzyme )
        for reaction in self._reactions:
            self.data.subsets.reactions.append( reaction )


    def end_reduced_system( self, content ):
        if( self._vectors != [] ):
            self.data.reduced_system.matrix = Matrix.Matrix( self._vectors[:14] )
        self.data.reduced_system.enzymes = []
        for enzyme in self._enzymes:
            self.data.reduced_system.enzymes.append( enzyme )
        for reaction in self._reactions:
            self.data.reduced_system.reactions.append( reaction )
        for col in self._irreversible_vector:
            self.data.reduced_system.irreversible_vector.append( col )


    def end_convex_basis( self, content ):
        if( self._vectors != [] ):
            self.data.convex_basis.matrix = Matrix.Matrix( self._vectors )
        self.data.convex_basis.enzymes = []
        for enzyme in self._enzymes:
            self.data.convex_basis.enzymes.append( enzyme )
        for reaction in self._reactions:
            self.data.convex_basis.reactions.append( reaction )

    def end_conservation_relations( self, content ):
        if( self._vectors != [] ):
            self.data.conservation_relations.matrix = Matrix.Matrix( self._vectors )
        self.data.conservation_relations.enzymes = []
        for enzyme in self._enzymes:
            self.data.conservation_relations.enzymes.append( enzyme )
        for reaction in self._reactions:
            self.data.conservation_relations.reactions.append( reaction )


    def end_elementary_modes( self, content ):
        if( self._vectors != [] ):
            self.data.elementary_modes.matrix = Matrix.Matrix( self._vectors )
        self.data.elementary_modes.enzymes = []
        for enzyme in self._enzymes:
            self.data.elementary_modes.enzymes.append( enzyme )
        for reaction in self._reactions:
            self.data.elementary_modes.reactions.append( reaction )

class _Scanner:
    """Start up Martel to do the scanning of the file.

    This initialzes the Martel based parser and connects it to a handler
    that will generate events for a Feature Consumer.
    """
    def __init__(self, debug = 0):
        """Initialize the scanner by setting up our caches.

        Creating the parser takes a long time, so we want to cache it
        to reduce parsing time.

        Arguments:
        o debug - The level of debugging that the parser should
        display. Level 0 is no debugging, Level 2 displays the most
        debugging info (but is much slower). See Martel documentation
        for more info on this.
        """
        # a listing of all tags we are interested in scanning for
        # in the MartelParser
        self.interest_tags = [ "input_file_name", "num_int_metabolites", \
            "num_reactions", "metabolite_line",  "unbalanced_metabolite", \
            "num_rows", "num_cols", "irreversible_vector", \
            "branch_metabolite", "non_branch_metabolite", \
            "stoichiometric_tag", "kernel_tag", "subsets_tag", \
            "reduced_system_tag", "convex_basis_tag", \
            "conservation_relations_tag", "elementary_modes_tag", \
            "reaction", "enzyme", "matrix_row", "sum_is_constant_line", \
            "end_stochiometric", "end_kernel", "end_subsets", \
            "end_reduced_system", "end_convex_basis", \
            "end_conservation_relations", "end_elementary_modes" ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names( metatool_format.metatool_record,
                                         self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        """Feeed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        self._parser.setContentHandler(EventGenerator(consumer,
                                                      self.interest_tags ))
#                                                      _strip_and_combine ))
        self._parser.setErrorHandler(handler.ErrorHandler())

        self._parser.parseFile(handle)

class RecordParser:
    """Parse MetaTool files into Record objects
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
        """Parse the specified handle into a MetaTool record.
        """
        self._consumer = _RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

def _strip_and_combine(line_list):
    """Combine multiple lines of content separated by spaces.

    This function is used by the EventGenerator callback function to
    combine multiple lines of information. The lines are first
    stripped to remove whitepsace, and then combined so they are separated
    by a space. This is a simple minded way to combine lines, but should
    work for most cases.
    """
    # first strip out extra whitespace
    stripped_line_list = map(string.strip, line_list)

    # now combine everything with spaces
    return ' '.join(stripped_line_list)





