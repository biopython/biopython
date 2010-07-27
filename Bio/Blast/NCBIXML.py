# Copyright 2000 by Bertrand Frottier .  All rights reserved.
# Revisions 2005-2006 copyright Michiel de Hoon
# Revisions 2006-2009 copyright Peter Cock
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""This module provides code to work with the BLAST XML output
following the DTD available on the NCBI FTP
ftp://ftp.ncbi.nlm.nih.gov/blast/documents/xml/NCBI_BlastOutput.dtd

Classes:
BlastParser         Parses XML output from BLAST (direct use discouraged).
                    This (now) returns a list of Blast records.
                    Historically it returned a single Blast record.
                    You are expected to use this via the parse or read functions.

_XMLParser          Generic SAX parser (private).

Functions:
parse               Incremental parser, this is an iterator that returns
                    Blast records.  It uses the BlastParser internally.
read                Returns a single Blast record. Uses the BlastParser internally.
"""
from Bio.Blast import Record
import xml.sax
from xml.sax.handler import ContentHandler

class _XMLparser(ContentHandler):
    """Generic SAX Parser

    Just a very basic SAX parser.

    Redefine the methods startElement, characters and endElement.
    """
    def __init__(self, debug=0):
        """Constructor

        debug - integer, amount of debug information to print
        """
        self._tag = []
        self._value = ''
        self._debug = debug
        self._debug_ignore_list = []

    def _secure_name(self, name):
        """Removes 'dangerous' from tag names

        name -- name to be 'secured'
        """
        # Replace '-' with '_' in XML tag names
        return name.replace('-', '_')
    
    def startElement(self, name, attr):
        """Found XML start tag

        No real need of attr, BLAST DTD doesn't use them

        name -- name of the tag

        attr -- tag attributes
        """
        self._tag.append(name)
        
        # Try to call a method (defined in subclasses)
        method = self._secure_name('_start_' + name)

        #Note could use try / except AttributeError
        #BUT I found often triggered by nested errors...
        if hasattr(self, method):
            eval("self.%s()" % method)
            if self._debug > 4:
                print "NCBIXML: Parsed:  " + method
        else:
            # Doesn't exist (yet)
            if method not in self._debug_ignore_list:
                if self._debug > 3:
                    print "NCBIXML: Ignored: " + method
                self._debug_ignore_list.append(method)

        #We don't care about white space in parent tags like Hsp,
        #but that white space doesn't belong to child tags like Hsp_midline
        if self._value.strip():
            raise ValueError("What should we do with %s before the %s tag?" \
                             % (repr(self._value), name))
        self._value = ""

    def characters(self, ch):
        """Found some text

        ch -- characters read
        """
        self._value += ch # You don't ever get the whole string

    def endElement(self, name):
        """Found XML end tag

        name -- tag name
        """
        # DON'T strip any white space, we may need it e.g. the hsp-midline
        
        # Try to call a method (defined in subclasses)
        method = self._secure_name('_end_' + name)
        #Note could use try / except AttributeError
        #BUT I found often triggered by nested errors...
        if hasattr(self, method):
            eval("self.%s()" % method)
            if self._debug > 2:
                print "NCBIXML: Parsed:  " + method, self._value
        else:
            # Doesn't exist (yet)
            if method not in self._debug_ignore_list:
                if self._debug > 1:
                    print "NCBIXML: Ignored: " + method, self._value
                self._debug_ignore_list.append(method)
        
        # Reset character buffer
        self._value = ''
        
class BlastParser(_XMLparser):
    """Parse XML BLAST data into a Record.Blast object

    All XML 'action' methods are private methods and may be:
    _start_TAG      called when the start tag is found
    _end_TAG        called when the end tag is found
    """

    def __init__(self, debug=0):
        """Constructor

        debug - integer, amount of debug information to print
        """
        # Calling superclass method
        _XMLparser.__init__(self, debug)
        
        self._parser = xml.sax.make_parser()
        self._parser.setContentHandler(self)
        
        # To avoid ValueError: unknown url type: NCBI_BlastOutput.dtd
        self._parser.setFeature(xml.sax.handler.feature_validation, 0)
        self._parser.setFeature(xml.sax.handler.feature_namespaces, 0)
        self._parser.setFeature(xml.sax.handler.feature_external_pes, 0)
        self._parser.setFeature(xml.sax.handler.feature_external_ges, 0)

        self.reset()

    def reset(self):
        """Reset all the data allowing reuse of the BlastParser() object"""
        self._records = []
        self._header = Record.Header()
        self._parameters = Record.Parameters()
        self._parameters.filter = None #Maybe I should update the class?

    def _start_Iteration(self):
        self._blast = Record.Blast()
        pass

    def _end_Iteration(self):
        # We stored a lot of generic "top level" information
        # in self._header (an object of type Record.Header)
        self._blast.reference = self._header.reference
        self._blast.date = self._header.date
        self._blast.version = self._header.version
        self._blast.database = self._header.database
        self._blast.application = self._header.application

        # These are required for "old" pre 2.2.14 files
        # where only <BlastOutput_query-ID>, <BlastOutput_query-def>
        # and <BlastOutput_query-len> were used.  Now they
        # are suplemented/replaced by <Iteration_query-ID>,
        # <Iteration_query-def> and <Iteration_query-len>
        if not hasattr(self._blast, "query") \
        or not self._blast.query:
            self._blast.query = self._header.query
        if not hasattr(self._blast, "query_id") \
        or not self._blast.query_id:
            self._blast.query_id = self._header.query_id
        if not hasattr(self._blast, "query_letters") \
        or not self._blast.query_letters:
            self._blast.query_letters = self._header.query_letters

        # Hack to record the query length as both the query_letters and
        # query_length properties (as in the plain text parser, see
        # Bug 2176 comment 12):
        self._blast.query_length = self._blast.query_letters
        # Perhaps in the long term we should deprecate one, but I would
        # prefer to drop query_letters - so we need a transition period
        # with both.

        # Hack to record the claimed database size as database_length
        # (as well as in num_letters_in_database, see Bug 2176 comment 13):
        self._blast.database_length = self._blast.num_letters_in_database
        # TODO? Deprecate database_letters next?

        # Hack to record the claimed database sequence count as database_sequences
        self._blast.database_sequences = self._blast.num_sequences_in_database

        # Apply the "top level" parameter information
        self._blast.matrix = self._parameters.matrix
        self._blast.num_seqs_better_e = self._parameters.num_seqs_better_e
        self._blast.gap_penalties = self._parameters.gap_penalties
        self._blast.filter = self._parameters.filter
        self._blast.expect = self._parameters.expect
        self._blast.sc_match = self._parameters.sc_match
        self._blast.sc_mismatch = self._parameters.sc_mismatch

        #Add to the list
        self._records.append(self._blast)
        #Clear the object (a new empty one is create in _start_Iteration)
        self._blast = None

        if self._debug : "NCBIXML: Added Blast record to results"

    # Header
    def _end_BlastOutput_program(self):
        """BLAST program, e.g., blastp, blastn, etc.

        Save this to put on each blast record object
        """
        self._header.application = self._value.upper()

    def _end_BlastOutput_version(self):
        """version number and date of the BLAST engine.

        e.g. "BLASTX 2.2.12 [Aug-07-2005]" but there can also be
        variants like "BLASTP 2.2.18+" without the date.

        Save this to put on each blast record object
        """
        parts = self._value.split()
        #TODO - Check the first word starts with BLAST?

        #The version is the second word (field one)
        self._header.version = parts[1]
        
        #Check there is a third word (the date)
        if len(parts) >= 3:
            if parts[2][0] == "[" and parts[2][-1] == "]":
                self._header.date = parts[2][1:-1]
            else:
                #Assume this is still a date, but without the
                #square brackets
                self._header.date = parts[2]

    def _end_BlastOutput_reference(self):
        """a reference to the article describing the algorithm

        Save this to put on each blast record object
        """
        self._header.reference = self._value

    def _end_BlastOutput_db(self):
        """the database(s) searched

        Save this to put on each blast record object
        """
        self._header.database = self._value

    def _end_BlastOutput_query_ID(self):
        """the identifier of the query

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-ID> is enough
        """
        self._header.query_id = self._value

    def _end_BlastOutput_query_def(self):
        """the definition line of the query

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-def> is enough
        """
        self._header.query = self._value

    def _end_BlastOutput_query_len(self):
        """the length of the query

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-len> is enough
        """
        self._header.query_letters = int(self._value)

    def _end_Iteration_query_ID(self):
        """the identifier of the query
        """
        self._blast.query_id = self._value

    def _end_Iteration_query_def(self):
        """the definition line of the query
        """
        self._blast.query = self._value

    def _end_Iteration_query_len(self):
        """the length of the query
        """
        self._blast.query_letters = int(self._value)

##     def _end_BlastOutput_query_seq(self):
##         """the query sequence
##         """
##         pass # XXX Missing in Record.Blast ?

##     def _end_BlastOutput_iter_num(self):
##         """the psi-blast iteration number
##         """
##         pass # XXX TODO PSI

    def _end_BlastOutput_hits(self):
        """hits to the database sequences, one for every sequence
        """
        self._blast.num_hits = int(self._value)

##     def _end_BlastOutput_message(self):
##         """error messages
##         """
##         pass # XXX What to do ?

    # Parameters
    def _end_Parameters_matrix(self):
        """matrix used (-M)
        """
        self._parameters.matrix = self._value
        
    def _end_Parameters_expect(self):
        """expect values cutoff (-e)
        """
        # NOTE: In old text output there was a line:
        # Number of sequences better than 1.0e-004: 1
        # As far as I can see, parameters.num_seqs_better_e
        # would take the value of 1, and the expectation
        # value was not recorded.
        #
        # Anyway we should NOT record this against num_seqs_better_e
        self._parameters.expect = self._value

##     def _end_Parameters_include(self):
##         """inclusion threshold for a psi-blast iteration (-h)
##         """
##         pass # XXX TODO PSI
    
    def _end_Parameters_sc_match(self):
        """match score for nucleotide-nucleotide comparaison (-r)
        """
        self._parameters.sc_match = int(self._value)

    def _end_Parameters_sc_mismatch(self):
        """mismatch penalty for nucleotide-nucleotide comparaison (-r)
        """
        self._parameters.sc_mismatch = int(self._value)

    def _end_Parameters_gap_open(self):
        """gap existence cost (-G)
        """
        self._parameters.gap_penalties = int(self._value)

    def _end_Parameters_gap_extend(self):
        """gap extension cose (-E)
        """
        self._parameters.gap_penalties = (self._parameters.gap_penalties,
                                         int(self._value))

    def _end_Parameters_filter(self):
        """filtering options (-F)
        """
        self._parameters.filter = self._value

##     def _end_Parameters_pattern(self):
##         """pattern used for phi-blast search
##         """
##         pass # XXX TODO PSI

##     def _end_Parameters_entrez_query(self):
##         """entrez query used to limit search
##         """
##         pass # XXX TODO PSI

    # Hits
    def _start_Hit(self):
        self._blast.alignments.append(Record.Alignment())
        self._blast.descriptions.append(Record.Description())
        self._blast.multiple_alignment = []
        self._hit = self._blast.alignments[-1]
        self._descr = self._blast.descriptions[-1]
        self._descr.num_alignments = 0

    def _end_Hit(self):
        #Cleanup
        self._blast.multiple_alignment = None
        self._hit = None
        self._descr = None

    def _end_Hit_id(self):
        """identifier of the database sequence
        """
        self._hit.hit_id = self._value
        self._hit.title = self._value + ' '

    def _end_Hit_def(self):
        """definition line of the database sequence
        """
        self._hit.hit_def = self._value
        self._hit.title += self._value
        self._descr.title = self._hit.title

    def _end_Hit_accession(self):
        """accession of the database sequence
        """
        self._hit.accession = self._value
        self._descr.accession = self._value

    def _end_Hit_len(self):
        self._hit.length = int(self._value)

    # HSPs
    def _start_Hsp(self):
        #Note that self._start_Hit() should have been called
        #to setup things like self._blast.multiple_alignment
        self._hit.hsps.append(Record.HSP())
        self._hsp = self._hit.hsps[-1]
        self._descr.num_alignments += 1
        self._blast.multiple_alignment.append(Record.MultipleAlignment())
        self._mult_al = self._blast.multiple_alignment[-1]

    # Hsp_num is useless
    def _end_Hsp_score(self):
        """raw score of HSP
        """
        self._hsp.score = float(self._value)
        if self._descr.score == None:
            self._descr.score = float(self._value)

    def _end_Hsp_bit_score(self):
        """bit score of HSP
        """
        self._hsp.bits = float(self._value)
        if self._descr.bits == None:
            self._descr.bits = float(self._value)

    def _end_Hsp_evalue(self):
        """expect value value of the HSP
        """
        self._hsp.expect = float(self._value)
        if self._descr.e == None:
            self._descr.e = float(self._value)

    def _end_Hsp_query_from(self):
        """offset of query at the start of the alignment (one-offset)
        """
        self._hsp.query_start = int(self._value)

    def _end_Hsp_query_to(self):
        """offset of query at the end of the alignment (one-offset)
        """
        self._hsp.query_end = int(self._value)

    def _end_Hsp_hit_from(self):
        """offset of the database at the start of the alignment (one-offset)
        """
        self._hsp.sbjct_start = int(self._value)

    def _end_Hsp_hit_to(self):
        """offset of the database at the end of the alignment (one-offset)
        """
        self._hsp.sbjct_end = int(self._value)

##     def _end_Hsp_pattern_from(self):
##         """start of phi-blast pattern on the query (one-offset)
##         """
##         pass # XXX TODO PSI

##     def _end_Hsp_pattern_to(self):
##         """end of phi-blast pattern on the query (one-offset)
##         """
##         pass # XXX TODO PSI

    def _end_Hsp_query_frame(self):
        """frame of the query if applicable
        """
        self._hsp.frame = (int(self._value),)

    def _end_Hsp_hit_frame(self):
        """frame of the database sequence if applicable
        """
        self._hsp.frame += (int(self._value),)

    def _end_Hsp_identity(self):
        """number of identities in the alignment
        """
        self._hsp.identities = int(self._value)

    def _end_Hsp_positive(self):
        """number of positive (conservative) substitutions in the alignment
        """
        self._hsp.positives = int(self._value)

    def _end_Hsp_gaps(self):
        """number of gaps in the alignment
        """
        self._hsp.gaps = int(self._value)

    def _end_Hsp_align_len(self):
        """length of the alignment
        """
        self._hsp.align_length = int(self._value)

##     def _en_Hsp_density(self):
##         """score density
##         """
##         pass # XXX ???

    def _end_Hsp_qseq(self):
        """alignment string for the query
        """
        self._hsp.query = self._value

    def _end_Hsp_hseq(self):
        """alignment string for the database
        """
        self._hsp.sbjct = self._value

    def _end_Hsp_midline(self):
        """Formatting middle line as normally seen in BLAST report
        """
        self._hsp.match = self._value # do NOT strip spaces!
        assert len(self._hsp.match)==len(self._hsp.query)
        assert len(self._hsp.match)==len(self._hsp.sbjct)

    # Statistics
    def _end_Statistics_db_num(self):
        """number of sequences in the database
        """
        self._blast.num_sequences_in_database = int(self._value)

    def _end_Statistics_db_len(self):
        """number of letters in the database
        """
        self._blast.num_letters_in_database = int(self._value)

    def _end_Statistics_hsp_len(self):
        """the effective HSP length
        """
        self._blast.effective_hsp_length = int(self._value)

    def _end_Statistics_eff_space(self):
        """the effective search space
        """
        self._blast.effective_search_space = float(self._value)

    def _end_Statistics_kappa(self):
        """Karlin-Altschul parameter K
        """
        self._blast.ka_params = float(self._value)

    def _end_Statistics_lambda(self):
        """Karlin-Altschul parameter Lambda
        """
        self._blast.ka_params = (float(self._value),
                                 self._blast.ka_params)

    def _end_Statistics_entropy(self):
        """Karlin-Altschul parameter H
        """
        self._blast.ka_params = self._blast.ka_params + (float(self._value),)
    
def read(handle, debug=0):
   """Returns a single Blast record (assumes just one query).

   This function is for use when there is one and only one BLAST
   result in your XML file.

   Use the Bio.Blast.NCBIXML.parse() function if you expect more than
   one BLAST record (i.e. if you have more than one query sequence).

   """
   iterator = parse(handle, debug)
   try:
       first = iterator.next()
   except StopIteration:
       first = None
   if first is None:
       raise ValueError("No records found in handle")
   try:
       second = iterator.next()
   except StopIteration:
       second = None
   if second is not None:
       raise ValueError("More than one record found in handle")
   return first


def parse(handle, debug=0):
    """Returns an iterator a Blast record for each query.

    handle - file handle to and XML file to parse
    debug - integer, amount of debug information to print

    This is a generator function that returns multiple Blast records
    objects - one for each query sequence given to blast.  The file
    is read incrementally, returning complete records as they are read
    in.

    Should cope with new BLAST 2.2.14+ which gives a single XML file
    for mutliple query records.

    Should also cope with XML output from older versions BLAST which
    gave multiple XML files concatenated together (giving a single file
    which strictly speaking wasn't valid XML)."""
    from xml.parsers import expat
    BLOCK = 1024
    MARGIN = 10 # must be at least length of newline + XML start
    XML_START = "<?xml"

    text = handle.read(BLOCK)
    pending = ""

    if not text:
        #NO DATA FOUND!
        raise ValueError("Your XML file was empty")
    
    while text:
        #We are now starting a new XML file
        if not text.startswith(XML_START):
            raise ValueError("Your XML file did not start with %s... "
                             "but instead %s" \
                             % (XML_START, repr(text[:20])))

        expat_parser = expat.ParserCreate()
        blast_parser = BlastParser(debug)
        expat_parser.StartElementHandler = blast_parser.startElement
        expat_parser.EndElementHandler = blast_parser.endElement
        expat_parser.CharacterDataHandler = blast_parser.characters

        expat_parser.Parse(text, False)
        while blast_parser._records:
            record = blast_parser._records[0]
            blast_parser._records = blast_parser._records[1:]
            yield record

        while True:
            #Read in another block of the file...
            text, pending = pending + handle.read(BLOCK), ""
            if not text:
                #End of the file!
                expat_parser.Parse("", True) # End of XML record
                break

            #Now read a little bit more so we can check for the
            #start of another XML file...
            pending = handle.read(MARGIN)

            if (text+pending).find("\n" + XML_START) == -1:
                # Good - still dealing with the same XML file
                expat_parser.Parse(text, False)        
                while blast_parser._records:
                    yield blast_parser._records.pop(0)
            else:
                # This is output from pre 2.2.14 BLAST,
                # one XML file for each query!
                
                # Finish the old file:
                text, pending = (text+pending).split("\n" + XML_START,1)
                pending = XML_START + pending

                expat_parser.Parse(text, True) # End of XML record
                while blast_parser._records:
                    yield blast_parser._records.pop(0)
               
                #Now we are going to re-loop, reset the
                #parsers and start reading the next XML file
                text, pending = pending, ""
                break

        #this was added because it seems that the Jython expat parser
        #was adding records later then the Python one
        while blast_parser._records:
            yield blast_parser._records.pop(0)
            
        #At this point we have finished the first XML record.
        #If the file is from an old version of blast, it may
        #contain more XML records (check if text=="").
        assert pending==""
        assert len(blast_parser._records) == 0
        
    #We should have finished the file!
    assert text==""
    assert pending==""
    assert len(blast_parser._records) == 0

if __name__ == '__main__':
    import sys
    import os
    handle = open(sys.argv[1])
    r_list = parse(handle)

    for r in r_list:
        # Small test
        print 'Blast of', r.query
        print 'Found %s alignments with a total of %s HSPs' % (len(r.alignments),
                  reduce(lambda a,b: a+b,
                         [len(a.hsps) for a in r.alignments]))

        for al in r.alignments:
            print al.title[:50], al.length, 'bp', len(al.hsps), 'HSPs'

        # Cookbook example
        E_VALUE_THRESH = 0.04
        for alignment in r.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    print '*****'
                    print 'sequence', alignment.title
                    print 'length', alignment.length
                    print 'e value', hsp.expect
                    print hsp.query[:75] + '...'
                    print hsp.match[:75] + '...'
                    print hsp.sbjct[:75] + '...'
