# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__all__ = [
    'Record',
    'kabat_format'
    ]

# standard library
import string
import array
import os
import re
import sgmllib
import urlparse

# XML from python 2.0
from xml.sax import handler

# Martel
import Martel
from Martel import RecordReader

from Bio.ParserSupport import EventGenerator
from Bio.ParserSupport import AbstractConsumer
from Bio.FilteredReader import FilteredReader
from Bio.FilteredReader import remove_empty_line
from Bio.SeqFeature import Reference
from Bio import File
from Martel.Dispatch import Dispatcher
import kabat_format
import Record


def remove_kabat_header_line( line ):
    stripped_line = line.strip()
    if( stripped_line.startswith( '<<' ) ):
        if( stripped_line.endswith( '>>' ) ):
            return ''
        else:
            return line[ : ]
    else:
        return line[ : ]


class Iterator:
    """Iterator interface to move over a file of Kabat entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with Kabat entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.

        Instructions:
        Browse http://immuno.bme.nwu.edu/seqhunt.html
        Choose a search criterion
        Enter a search pattern
        Select ASCII
        Submit search
        Cut and paste search results into a text file
        Use code fragment as a template

        Code Fragment:
        src_handle = open( datafile )
        iterator = Kabat.Iterator(src_handle, record_parser)
        data = iterator.next()
        data.print_kabat()
        print '\n'

        """
        filtered_handle = FilteredReader( handle )
        filtered_handle.filter_chain = [ remove_empty_line, remove_kabat_header_line ]
        self._reader = RecordReader.StartsWith( filtered_handle, "KADBID")
        self._parser = parser

    def next(self):
        """Return the next Kabat record from the handle.

        Will return None if we ran out of records.
        """
        data = self._reader.next()

        if self._parser is not None:
            if data:
                dumpfile = open( 'dump', 'w' )
                dumpfile.write( data )
                dumpfile.close()
                return self._parser.parse(File.StringHandle(data))

        return data

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
        self.interest_tags = ["kabatid", "creation_date", "last_mod_date",
                              "definition", "species", "nucleotide_sequence_name",
                              "amino_acid_sequence_name", "nucleotide_ref_author",
                              "nucleotide_ref_pubmed", "nucleotide_ref_journal",
                              "amino_acid_ref_author", "amino_acid_ref_pubmed",
                              "amino_acid_ref_journal", "annotation_key",
                              "annotation_val", "codon",
                               "amino_1_letter_code" ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names(kabat_format.kabat_record, self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        """Feeed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        consumer.set_interest_tags( self.interest_tags )
        self._parser.setContentHandler( consumer )
#        self._parser.setErrorHandler(handle.ErrorHandler())

        self._parser.parseFile(handle)

class _RecordConsumer( Dispatcher ):
    """Create a Kabat Record object from scanner generated information.
    """
    def __init__(self):
        Dispatcher.__init__( self )
        self.data = Record.Record()

        self._cur_nucleotide_ref = None
        self._cur_amino_acid_ref = None

    def set_interest_tags( self, interest_tags ):
        self.interest_tags = interest_tags

    def start_kabatid(self, content, attrs ):
        self.save_characters()

    def end_kabatid(self, content ):
        next_line = self.get_characters()
        self.data.kabatid = next_line

    def start_creation_date(self, content, attrs ):
        self.save_characters()

    def end_creation_date(self, content):
        next_line = self.get_characters()
        self.data.creation_date = next_line

    def start_last_mod_date(self, content, attrs ):
        self.save_characters()

    def end_last_mod_date(self, content):
        next_line = self.get_characters()
        self.data.date_last_mod = next_line

    def start_definition(self, content, attrs ):
        self.save_characters()

    def end_definition(self, content):
        next_line = self.get_characters()
        self.data.definition = next_line

    def start_species(self, content, attrs ):
        self.save_characters()

    def end_species(self, content):
        next_line = self.get_characters()
        self.data.species = next_line

    def start_nucleotide_sequence_name(self, content, attrs ):
        self.save_characters()

    def end_nucleotide_sequence_name(self, content):
        self.data.nucleotide_sequence_name = self.get_characters()

    def start_amino_acid_sequence_name(self, content, attrs ):
        self.save_characters()

    def end_amino_acid_sequence_name(self, content):
        self.data.nucleotide_sequence_name = self.get_characters()

    def start_nucleotide_ref_author(self, content, attrs ):
        self.save_characters()

    def end_nucleotide_ref_author(self, content):
        next_line = self.get_characters()
        author_info = next_line[ 6: ]
        self._cur_nucleotide_ref = Record.KabatReference()
        self._cur_nucleotide_ref.authors = author_info
        continuation = next_line[ :6 ]
        continuation = continuation.strip()
        if( continuation != '1' ):
            self.data.nucleotide_refs.append( self._cur_nucleotide_ref )
            self._cur_nucleotide_ref = None

    def start_nucleotide_ref_journal(self, content, attrs ):
        self.save_characters()

    def end_nucleotide_ref_journal(self, content):
        next_line = self.get_characters()
        journal_info = next_line[ 6: ]
        self._cur_nucleotide_ref.journal = journal_info
        continuation = next_line[ :6 ]
        continuation = continuation.strip()
        if( continuation != '1' ):
            self.data.nucleotide_refs.append( self._cur_nucleotide_ref )
            self._cur_nucleotide_ref = None

    def start_nucleotide_ref_pubmed(self, content, attrs ):
        self.save_characters()

    def end_nucleotide_ref_pubmed(self, content):
        next_line = self.get_characters()
        pubmed_info = next_line[ 6: ]
        self._cur_nucleotide_ref.pubmed_id = pubmed_info.strip()
        continuation = next_line[ :6 ]
        continuation = continuation.strip()
        if( continuation != '1' ):
            self.data.nucleotide_refs.append( self._cur_nucleotide_ref )
            self._cur_nucleotide_ref = None

    def start_amino_acid_ref_author(self, content, attrs ):
        self.save_characters()

    def end_amino_acid_ref_author(self, content):
        next_line = self.get_characters()
        author_info = next_line[ 6: ]
        self._cur_amino_acid_ref = Record.KabatReference()
        self._cur_amino_acid_ref.authors = author_info
        continuation = next_line[ :6 ]
        try:
            continuation = continuation.strip()
        except:
            print 'continuation is %s' % continuation
        if( continuation != '1' ):
            self.data.amino_acid_refs.append( self._cur_amino_acid_ref )
            self._cur_amino_acid_ref = None

    def start_amino_acid_ref_journal(self, content, attrs ):
        self.save_characters()

    def end_amino_acid_ref_journal(self, content):
        next_line = self.get_characters()
        journal_info = next_line[ 6: ]
        self._cur_amino_acid_ref.journal = journal_info
        continuation = next_line[ :6 ]
        continuation = continuation.strip()
        if( continuation != '1' ):
            self.data.amino_acid_refs.append( self._cur_amino_acid_ref )
            self._cur_amino_acid_ref = None

    def start_amino_acid_ref_pubmed(self, content, attrs ):
        self.save_characters()

    def end_amino_acid_ref_pubmed(self, content):
        next_line = self.get_characters()
        pubmed_info = next_line[ 6: ]
        self._cur_amino_acid_ref.pubmed_id = pubmed_info.strip()
        continuation = next_line[ :6 ]
        continuation = continuation.strip()
        if( continuation != '1' ):
            self.data.amino_acid_refs.append( self._cur_amino_acid_ref )
            self._cur_amino_acid_ref = None

    def start_codon(self, content, attrs ):
        self.save_characters()

    def end_codon(self, content):
        nucleotides = self.get_characters()
        nucleotides = nucleotides.strip()
        for item in array.array( 'c', nucleotides ):
            self.data.nucleotide_sequence.append( item )
        num_dashes = 3 - len( nucleotides )
        for i in range( 0, num_dashes ):
            self.data.nucleotide_sequence.append( '-' )
        self.data.amino_acid_sequence.append( '-' )

    def start_amino_1_letter_code(self, content, attrs ):
        self.save_characters()

    def end_amino_1_letter_code(self, content):
        next_line = self.get_characters()
        self.data.amino_acid_sequence[ -1 ] = next_line[ 0 ]

    def start_annotation_key(self, content, attrs ):
        self.save_characters()

    def end_annotation_key( self, content ):
        text = self.get_characters()
        self.pending_key = text

    def start_annotation_val(self, content, attrs ):
        self.save_characters()

    def end_annotation_val( self, content ):
        val = self.get_characters()
        val = val.strip()

        if( self.pending_key != None ):
            key = self.pending_key
            if( self.data.annotation.has_key( key ) ):
                val = self.data.annotation[ key ] + val
            self.data.annotation[ key ] = val
            self.pending_key = None


class RecordParser:
    """Parse Kabat files into Record objects
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
        """Parse the specified handle into a GenBank record.
        """
        self._consumer = _RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

