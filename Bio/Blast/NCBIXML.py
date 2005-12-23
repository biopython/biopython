# BLAST XML parsing
"""This module provides code to work with the BLAST XML output
following the DTD available on the NCBI FTP
ftp://ftp.ncbi.nlm.nih.gov/blast/documents/xml/NCBI_BlastOutput.dtd

Classes:
BlastParser         Parses XML output from BLAST.

_XMLParser          Generic SAX parser.
"""
from Bio.Blast import Record
import xml.sax
from xml.sax.handler import ContentHandler

class _XMLparser(ContentHandler):
    """Generic SAX Parser

    Just a very basic SAX parser.

    Redefine the methods startElement, characters and endElement.
    """
    def __init__(self):
        """Constructor
        """
        self._tag = []
        self._value = ''

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
        
        # Try to call a method
        try:
            eval(self._secure_name('self._start_' + name))()
        except AttributeError:
            # Doesn't exist (yet)
            pass

    def characters(self, ch):
        """Found some text

        ch -- characters read
        """
        self._value += ch # You don't ever get the whole string

    def endElement(self, name):
        """Found XML end tag

        name -- tag name
        """
        # Strip character buffer
        self._value = self._value.strip()
        
        # Try to call a method (defined in subclasses)
        try:
            eval(self._secure_name('self._end_' + name))()
        except AttributeError: # Method doesn't exist (yet ?)
            pass    

        # Reset character buffer
        self._value = ''
        
class BlastParser(_XMLparser):
    """Parse XML BLAST data into a Record.Blast object

    Methods:
    parse           Parses BLAST XML data.

    All XML 'action' methods are private methods and may be:
    _start_TAG      called when the start tag is found
    _end_TAG        called when the end tag is found
    """

    def __init__(self):
        """Constructor
        """
        # Calling superclass method
        _XMLparser.__init__(self)
        
        self._parser = xml.sax.make_parser()
        self._parser.setContentHandler(self)
        
        # To avoid ValueError: unknown url type: NCBI_BlastOutput.dtd
        self._parser.setFeature(xml.sax.handler.feature_validation, 0)
        self._parser.setFeature(xml.sax.handler.feature_namespaces, 0)
        self._parser.setFeature(xml.sax.handler.feature_external_pes, 0)
        self._parser.setFeature(xml.sax.handler.feature_external_ges, 0)

        self._blast = Record.Blast()

    def parse(self, handler):
        """Parses the XML data

        handler -- file handler or StringIO
        """
        # bugfix: changed `filename` to `handler`. Iddo 12/20/2004
        self._parser.parse(handler)
        return self._blast

    # Header
    def _end_BlastOutput_program(self):
        """BLAST program, e.g., blastp, blastn, etc.
        """
        self._blast.application = self._value.uuper()

    def _end_BlastOutput_version(self):
        """version number of the BLAST engine (e.g., 2.1.2)
        """
        self._blast.version = self._value.split()[1]
        self._blast.date = self._value.split()[2][1:-1]

    def _end_BlastOutput_reference(self):
        """a reference to the article describing the algorithm
        """
        self._blast.reference = self._value

    def _end_BlastOutput_db(self):
        """the database(s) searched
        """
        self._blast.database = self._value

##     def _end_BlastOutput_query_ID(self):
##         """the identifier of the query
##         """
##         pass # XXX MISSING in Record.Blast ? Useless ?

    def _end_BlastOutput_query_def(self):
        """the definition line of the query
        """
        self._blast.query = self._value

    def _end_BlastOutput_query_len(self):
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
        self._blast.matrix = self._value
        
    def _end_Parameters_expect(self):
        """expect values cutoff (-e)
        """
        self._blast.num_seqs_better_e = self._value

##     def _end_Parameters_include(self):
##         """inclusion threshold for a psi-blast iteration (-h)
##         """
##         pass # XXX TODO PSI
    
    def _end_Parameters_sc_match(self):
        """match score for nucleotide-nucleotide comparaison (-r)
        """
        self._blast.sc_match = int(self._value)

    def _end_Parameters_sc_mismatch(self):
        """mismatch penalty for nucleotide-nucleotide comparaison (-r)
        """
        self._blast.sc_mismatch = int(self._value)

    def _end_Parameters_gap_open(self):
        """gap existence cost (-G)
        """
        self._blast.gap_penalties = int(self._value)

    def _end_Parameters_gap_extend(self):
        """gap extension cose (-E)
        """
        self._blast.gap_penalties = (self._blast.gap_penalties,
                                     int(self._value))

    def _end_Parameters_filter(self):
        """filtering options (-F)
        """
        self._blast.filter = self._value

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
        
    # Hit_num is useless

    def _end_Hit_id(self):
        """identifier of the database sequence
        """
        self._hit.title = self._value + ' '

    def _end_Hit_def(self):
        """definition line of the database sequence
        """
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
        self._hit.hsps.append(Record.HSP())
        self._hsp = self._hit.hsps[-1]
        self._descr.num_alignments += 1
        self._blast.multiple_alignment.append(Record.MultipleAligment())
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
    # No need for Hsp_query_to

    def _end_Hsp_hit_from(self):
        """offset of the database at the start of the alignment (one-offset)
        """
        self._hsp.sbjct_start = int(self._value)
    # No need for Hsp_hit_to

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
        self._hsp.match = self._value

    # Statistics
    def _end_Statistics_db_num(self):
        """number of sequences in the database
        """
        self._blast.num_sequences_in_database = int(self._value)

    def _end_Statistics_db_len(self):
        """number of letters in the database
        """
        self._blast._num_letters_in_database = int(self._value)

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
    

    
if __name__ == '__main__':
    import sys
    p = BlastParser()
    r = p.parse(sys.argv[1])

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
