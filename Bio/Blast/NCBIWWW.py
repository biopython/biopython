# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""NCBIWWW.py

This module provides code to work with the WWW version of BLAST
provided by the NCBI.
http://www.ncbi.nlm.nih.gov/BLAST/

Classes:
BlastParser   Parses output from WWW blast.
_Scanner      Scans output from NCBI's BLAST WWW server.

Functions:
blast        Do a BLAST search.

"""
import string

from Bio import File
from Bio.WWW import NCBI
from Bio.ParserSupport import *
import NCBIStandalone

import time
import re
import sgmllib
import urlparse

class BlastParser:
    """Parses WWW BLAST data into a Record.Blast object.

    """
    def __init__(self):
        """__init__(self)"""
        self._scanner = _Scanner()
        self._consumer = SGMLStrippingConsumer(NCBIStandalone._BlastConsumer())

    def parse(self, handle):
        """parse(self, handle)"""
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data
    
class _Scanner:
    """Scan BLAST output from NCBI's web server at:
    http://www.ncbi.nlm.nih.gov/BLAST/
    
    Tested with BLAST v2.0.10

    Methods:
    feed     Feed data into the scanner.
    
    """
    def feed(self, handle, consumer):
        """S.feed(handle, consumer)

        Feed in a BLAST report for scanning.  handle is a file-like
        object that contains the BLAST report.  consumer is a Consumer
        object that will receive events as the report is scanned.

        """
        # This stuff appears in 2.0.12.
        # <p><!--
        # QBlastInfoBegin
        #         Status=READY
        # QBlastInfoEnd
        # --><p>

        # <HTML>
        # <HEAD>
        # <TITLE>BLAST Search Results </TITLE>
        # </HEAD>
        # <BODY BGCOLOR="#FFFFFF" LINK="#0000FF" VLINK="#660099" ALINK="#660099
        # <A HREF="http://www.ncbi.nlm.nih.gov/BLAST/blast_form.map"> <IMG SRC=
        # <BR><BR><PRE>

        # BLAST Formatted information
        
        # 
        # </BODY>
        # </HTML>
        # </BODY>
        # </HTML>
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        # Read HTML formatting up to the "BLAST" version line.
        read_and_call_until(uhandle, consumer.noevent, start='<b>',
                            contains='BLAST')

        self._scan_header(uhandle, consumer)
	self._scan_rounds(uhandle, consumer)
        self._scan_database_report(uhandle, consumer)
        self._scan_parameters(uhandle, consumer)

        # Read HTML footer information.
        while uhandle.peekline():
            read_and_call(uhandle, consumer.noevent)

    def _scan_header(self, uhandle, consumer):
        # <b>BLASTP 2.0.10 [Aug-26-1999]</b>
        # 
        # 
        # <b><a href="http://www.ncbi.nlm.nih.gov/htbin-
        # post/Entrez/query?uid=9254694&form=6&db=m&Dopt=r">Reference</a>:</b>
        # Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch&auml;ffer, 
        # Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), 
        # "Gapped BLAST and PSI-BLAST: a new generation of protein database sea
        # programs",  Nucleic Acids Res. 25:3389-3402.
        # <p>
        # <b>Query=</b> gi|120291|sp|P21297|FLBT_CAUCR FLBT PROTEIN.
        #          (141 letters)
        # 
        # <b>Database:</b> Non-redundant SwissProt sequences
        #            82,258 sequences; 29,652,561 total letters
        # 
        # <p> <p>If you have any problems or questions with the results of this

        # If there are hits, and Graphical Overview was selected:
        # <FORM NAME="BLASTFORM">
        # </PRE>
        # <CENTER>
        # <H3><a href="/BLAST/newoptions.html#graphical-overview"> Distribution
        # <input name=defline size=80 value="Mouse-over to show defline and sco
        # </CENTER>
        # <map name=img_map>
        # <area shape=rect coords=69,101,476,106 href="#120291" ONMOUSEOVER='do
        # <area shape=rect coords=156,108,305,113 href="#3024946" ONMOUSEOVER='
        # </map>
        # <CENTER>
        # <IMG WIDTH=529 HEIGHT=115 USEMAP=#img_map BORDER=1 SRC="nph-getgif.cg
        # <HR>
        # <PRE>

        consumer.start_header()

        # Read the "BLAST" version line and the following blanks.
        read_and_call(uhandle, consumer.version, contains='BLAST')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        # Read the reference lines and the '<p>' line.
        read_and_call_until(uhandle, consumer.reference, start='<p>')
        read_and_call(uhandle, consumer.noevent)

        # Read the RID line, for version 2.0.12 (2.0.11?) and above.
        attempt_read_and_call(uhandle, consumer.noevent, start='RID')

        # Read the Query lines and the following blank line.
        read_and_call(uhandle, consumer.query_info, contains='Query=')
        read_and_call_until(uhandle, consumer.query_info, blank=1)
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        # Read the database lines and the following blank line.
        attempt_read_and_call(uhandle, consumer.noevent, start='<p>')
        read_and_call(uhandle, consumer.database_info, contains='Database')
        read_and_call(uhandle, consumer.database_info, contains='sequences')
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent,
                      contains='problems or questions')

        # Read the blast form, if it exists. 
        if attempt_read_and_call(uhandle, consumer.noevent,
                                 contains='BLASTFORM'):
            read_and_call_until(uhandle, consumer.noevent, blank=1)
        elif attempt_read_and_call(uhandle, consumer.noevent, start='<PRE>'):
            read_and_call_until(uhandle, consumer.noevent, blank=1)

        # Read the blank lines until the next section.
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        consumer.end_header()

    def _scan_rounds(self, uhandle, consumer):
        self._scan_descriptions(uhandle, consumer)
        self._scan_alignments(uhandle, consumer)

    def _scan_descriptions(self, uhandle, consumer):
        consumer.start_descriptions()

        # Three things can happen here:
        # 1.  line contains 'Score     E'
        # 2.  line contains "No significant similarity"
        # 3.  no descriptions
        if not attempt_read_and_call(
            uhandle, consumer.noevent, contains='Score     E'):
            # Either case 2 or 3.  Look for "No hits found".
            attempt_read_and_call(uhandle, consumer.no_hits,
                                  contains='No significant similarity')
            read_and_call_while(uhandle, consumer.noevent, blank=1)
            consumer.end_descriptions()
            # Stop processing.
            return
        # Sequences producing significant alignments:                      
        # 
        # <a href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Ret
        # <a href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Ret
        # 
        # Read the score header lines and a blank line.
        read_and_call(uhandle, consumer.noevent, start='Sequences producing')
        read_and_call(uhandle, consumer.noevent, blank=1)

        # Read the descriptions and the following blank line.
        read_and_call_until(uhandle, consumer.description, blank=1)
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        consumer.end_descriptions()

    def _scan_alignments(self, uhandle, consumer):
        # Check to see whether I'm at an alignment or database report.
        # Possibilities:
        # 1) BLASTP 2.0.14, pairwise alignment
        #   <CENTER><b><FONT color="green">Alignments</FONT></b></CENTER>
        #   <PRE>
        #   ><a name = 121837></a><a href="http://www.ncbi.nlm.nih.gov:80/entre
        # 2) BLASTP 2.0.10, pairwise alignment
        #   <PRE>
        #   <a name = 120291> </a><a href="http://www.ncbi.nlm.nih.gov:80/entre
        # 3) BLASTP 2.0.10, master-slave
        #   <PRE>
        #   blast_tmp 1    MFQQIGAVQAKSGTDEPAHPCEKFPPERKCEAVFWKPLPRHEAREILLAARK
        # 4) BLASTP 2.0.10, 2.0.14, database
        #   <PRE>
        #     Database: Non-redundant SwissProt sequences

        # Get the first two lines and examine them.
        line1 = safe_readline(uhandle)
        line2 = safe_readline(uhandle)
        uhandle.saveline(line2)
        uhandle.saveline(line1)

        is_pairwise = is_masterslave = 0
        if string.find(line1, 'Alignments') >= 0:
            is_pairwise = 1
        elif line2[:10] == '  Database':
            pass
        elif line2[:9] == 'blast_tmp':
            is_masterslave = 1
        elif line1[:5] == '<PRE>':
            is_pairwise = 1
        else:
            raise SyntaxError, "Cannot resolve location at line:\n%s" % line1

        if is_pairwise:
            self._scan_pairwise_alignments(uhandle, consumer)
        elif is_masterslave:
            self._scan_masterslave_alignment(uhandle, consumer)

    def _scan_pairwise_alignments(self, uhandle, consumer):
        while 1:
            # The first line is <PRE>.  Check the second line to see if
            # I'm still at an alignment.
            line1 = safe_readline(uhandle)
            line2 = safe_readline(uhandle)
            uhandle.saveline(line2)
            uhandle.saveline(line1)
            if line2[:10] == '  Database':
                break

            # Occasionally, there's a bug where the alignment_header and
            # hsp_header are skipped, leaving only the hsp_alignment.
            # Detect this and handle it accordingly.
            if line2[:6] == 'Query:':
                self._scan_abbreviated_pairwise_alignment(uhandle, consumer)
            else:
                self._scan_one_pairwise_alignment(uhandle, consumer)

    def _scan_abbreviated_pairwise_alignment(self, uhandle, consumer):
        # Sometimes all header information is skipped, leaving
        # only the raw alignments.  I believe this is a bug because
        # without the header information, you lose vital information such
        # as score, target sequence id, etc.
        # Format:
        # <PRE>
        # hsp_alignment

        consumer.start_alignment()
        consumer.start_hsp()
        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        self._scan_hsp_alignment(uhandle, consumer)
        consumer.end_hsp()
        consumer.end_alignment()
        
    def _scan_one_pairwise_alignment(self, uhandle, consumer):
        # Alignment format:
        # <CENTER><b><FONT color="green">Alignments</FONT></b></CENTER>
        #       (BLAST 2.0.14)
        # <PRE>
        # alignment_header
        #   hsp_header
        #   hsp_alignment
        #   [...]
        # The hsp_header and hsp_alignment blocks can be repeated.

        consumer.start_alignment()
        attempt_read_and_call(uhandle, consumer.noevent, contains='Alignments')
        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        self._scan_alignment_header(uhandle, consumer)

        # Scan a bunch of score/alignment's.
        while 1:
            # An HSP header starts with ' Score'.
            # However, if the HSP header is not the first one in the
            # alignment, there will be a '<PRE>' line first.  Therefore,
            # I will need to check either of the first two lines to
            # see if I'm at an HSP header.
            line1 = safe_readline(uhandle)
            line2 = safe_readline(uhandle)
            uhandle.saveline(line2)
            uhandle.saveline(line1)
            if line1[:6] != ' Score' and line2[:6] != ' Score':
                break
            self._scan_hsp(uhandle, consumer)
                
        consumer.end_alignment()

    def _scan_alignment_header(self, uhandle, consumer):
        # <a name = 120291> </a><a href="http://www.ncbi.nlm.nih.gov:80/entrez/
        #            Length = 141
        #            
        while 1:
            line = safe_readline(uhandle)
            if string.lstrip(line)[:8] == 'Length =':
                consumer.length(line)
                break
            elif is_blank_line(line):
                # Check to make sure I haven't missed the Length line
                raise SyntaxError, "I missed the Length in an alignment header"
            consumer.title(line)

        if not attempt_read_and_call(uhandle, consumer.noevent,
                                     start='          '):
            read_and_call(uhandle, consumer.noevent, blank=1)

    def _scan_hsp(self, uhandle, consumer):
        consumer.start_hsp()
        self._scan_hsp_header(uhandle, consumer)
        self._scan_hsp_alignment(uhandle, consumer)
        consumer.end_hsp()

    def _scan_hsp_header(self, uhandle, consumer):
        # If the HSP is not the first one within an alignment, includes:
        # <PRE>
        
        #  Score = 22.7 bits (47), Expect = 2.5
        #  Identities = 10/36 (27%), Positives = 18/36 (49%)
        #  Strand = Plus / Plus
        #  Frame = +3
        #

        attempt_read_and_call(uhandle, consumer.noevent, start='<PRE>')
        read_and_call(uhandle, consumer.score, start=' Score')
        read_and_call(uhandle, consumer.identities, start=' Identities')
        # BLASTN
        attempt_read_and_call(uhandle, consumer.strand, start = ' Strand')
        # BLASTX, TBLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.frame, start = ' Frame')
        read_and_call(uhandle, consumer.noevent, blank=1)

    def _scan_hsp_alignment(self, uhandle, consumer):
        # Query: 11 GRGVSACA-------TCDGFFYRNQKVAVIGGGNTAVEEALYLSNIASEVHLIHRRDGF
        #           GRGVS+         TC    Y  + + V GGG+ + EE   L     +   I R+
        # Sbjct: 12 GRGVSSVVRRCIHKPTCKE--YAVKIIDVTGGGSFSAEEVQELREATLKEVDILRKVSG
        # 
        # Query: 64 AEKILIKR 71
        #              I +K 
        # Sbjct: 70 PNIIQLKD 77
        # </PRE>
        #
        # 

        while 1:
            # Blastn adds an extra line filled with spaces before Query
            attempt_read_and_call(uhandle, consumer.noevent, start='     ')
            read_and_call(uhandle, consumer.query, start='Query')
            read_and_call(uhandle, consumer.align, start='     ')
            read_and_call(uhandle, consumer.sbjct, start='Sbjct')
            if not attempt_read_and_call(uhandle, consumer.noevent, blank=1):
                break
        read_and_call(uhandle, consumer.noevent, start='</PRE>')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

    def _scan_masterslave_alignment(self, uhandle, consumer):
        consumer.start_alignment()
        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        while 1:
            line = safe_readline(uhandle)
            if is_blank_line(line):
                consumer.noevent(line)
            elif line[:6] == '</PRE>':
                consumer.noevent(line)
                break
            else:
                consumer.multalign(line)
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        consumer.end_alignment()

    def _scan_database_report(self, uhandle, consumer):
        # <PRE>
        #   Database: Non-redundant SwissProt sequences
        #     Posted date:  Dec 18, 1999  8:26 PM
        #   Number of letters in database: 29,652,561
        #   Number of sequences in database:  82,258
        #   
        # Lambda     K      H
        #    0.317    0.133    0.395 
        # 
        # Gapped
        # Lambda     K      H
        #    0.270   0.0470    0.230 
        # 

        consumer.start_database_report()

        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        read_and_call(uhandle, consumer.database, start='  Database')
        read_and_call(uhandle, consumer.posted_date, start='    Posted')
        read_and_call(uhandle, consumer.num_letters_in_database,
                      start='  Number of letters')
        read_and_call(uhandle, consumer.num_sequences_in_database,
                      start='  Number of sequences')
        read_and_call(uhandle, consumer.noevent, start='  ')

        read_and_call(uhandle, consumer.noevent, start='Lambda')
        read_and_call(uhandle, consumer.ka_params)
        read_and_call(uhandle, consumer.noevent, blank=1)

        # not BLASTP
        attempt_read_and_call(uhandle, consumer.gapped, start='Gapped')
        # not TBLASTX
        if attempt_read_and_call(uhandle, consumer.noevent, start='Lambda'):
            read_and_call(uhandle, consumer.ka_params_gap)
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        consumer.end_database_report()

    def _scan_parameters(self, uhandle, consumer):
        # Matrix: BLOSUM62
        # Number of Hits to DB: 1st pass: 41542626, 2nd pass: 9765
        # Number of Sequences: 1st pass: 89405, 2nd pass: 84
        # Number of extensions: 1st pass: 500847, 2nd pass: 6747
        # Number of successful extensions: 1st pass: 14, 2nd pass: 49
        # Number of sequences better than 10.0: 20
        # length of query: 205
        # length of database: 10,955,950
        # effective HSP length: 46
        # effective length of query: 158
        # effective length of database: 6,843,320
        # effective search space: 1081244560
        # effective search space used: 1081244560
        # frameshift window, decay const: 50,  0.5
        # T: 13
        # A: 40
        # X1: 16 ( 7.3 bits)
        # X2: 0 ( 0.0 bits)
        # S1: 41 (21.7 bits)
        # S2: 52 (26.7 bits)
        # 
        # </PRE>

        consumer.start_parameters()

        read_and_call(uhandle, consumer.matrix, start='Matrix')
        # not TBLASTX
        attempt_read_and_call(uhandle, consumer.gap_penalties, start='Gap')
        read_and_call(uhandle, consumer.num_hits,
                      start='Number of Hits')
        read_and_call(uhandle, consumer.num_sequences,
                      start='Number of Sequences')
        read_and_call(uhandle, consumer.num_extends,
                      start='Number of extensions')
        read_and_call(uhandle, consumer.num_good_extends,
                      start='Number of successful')

        read_and_call(uhandle, consumer.num_seqs_better_e,
                      start='Number of sequences')

        # not BLASTN, TBLASTX
        if attempt_read_and_call(uhandle, consumer.hsps_no_gap,
                                 start="Number of HSP's better"):
            read_and_call(uhandle, consumer.hsps_prelim_gapped,
                          start="Number of HSP's successfully")
            read_and_call(uhandle, consumer.hsps_prelim_gap_attempted,
                          start="Number of HSP's that")
            read_and_call(uhandle, consumer.hsps_gapped,
                          start="Number of HSP's gapped")

        read_and_call(uhandle, consumer.query_length,
                      start='length of query')
        read_and_call(uhandle, consumer.database_length,
                      start='length of database')

        read_and_call(uhandle, consumer.effective_hsp_length,
                      start='effective HSP')
        read_and_call(uhandle, consumer.effective_query_length,
                      start='effective length of query')
        read_and_call(uhandle, consumer.effective_database_length,
                      start='effective length of database')
        read_and_call(uhandle, consumer.effective_search_space,
                      start='effective search space')
        read_and_call(uhandle, consumer.effective_search_space_used,
                      start='effective search space used')

        # BLASTX, TBLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.frameshift, start='frameshift')
        read_and_call(uhandle, consumer.threshold, start='T')
        read_and_call(uhandle, consumer.window_size, start='A')
        read_and_call(uhandle, consumer.dropoff_1st_pass, start='X1')
        read_and_call(uhandle, consumer.gap_x_dropoff, start='X2')
        # not BLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.gap_x_dropoff_final,
                              start='X3')
        read_and_call(uhandle, consumer.gap_trigger, start='S1')
        read_and_call(uhandle, consumer.blast_cutoff, start='S2')

        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start='</PRE>')

        consumer.end_parameters()


def blast(program, datalib, sequence,
          input_type='Sequence in FASTA format',
          double_window=None, gi_list='(None)', expect='10',
          filter='L', genetic_code='Standard (1)',
          mat_param='PAM30     9       1',
          other_advanced=None, ncbi_gi=None, overview=None,
          alignment_view='0', descriptions=None, alignments=None,
          email=None, path=None, html=None, 
          cgi='http://www.ncbi.nlm.nih.gov/blast/blast.cgi',
          timeout=20
          ):
    """blast(program, datalib, sequence,
    input_type='Sequence in FASTA format',
    double_window=None, gi_list='(None)', expect='10',
    filter='L', genetic_code='Standard (1)',
    mat_param='PAM30     9       1',
    other_advanced=None, ncbi_gi=None, overview=None,
    alignment_view='0', descriptions=None, alignments=None,
    email=None, path=None, html=None, 
    cgi='http://www.ncbi.nlm.nih.gov/blast/blast.cgi',
    timeout=20) -> handle

    Do a BLAST search against NCBI.  Returns a handle to the results.
    timeout is the number of seconds to wait for the results before timing
    out.  The other parameters are provided to BLAST.  A description
    can be found online at:
    http://www.ncbi.nlm.nih.gov/BLAST/newoptions.html

    """
    # NCBI Blast is hard to work with.  The user enters a query, and then
    # it returns a "reference" page which contains a button that the user
    # clicks to retrieve the results.  This will retrieve the "results"
    # page.  However, this page may not contain BLAST results if the
    # search isn't done.
    # This function will send off the query and parse the reference
    # page to figure out how to retrieve the results.  Then, it needs to
    # check the results to see if the search has been finished.
    params = {'PROGRAM' : program,
              'DATALIB' : datalib,
              'SEQUENCE' : sequence,
              'DOUBLE_WINDOW' : double_window,
              'INPUT_TYPE' : input_type,
              'EXPECT' : expect,
              'FILTER' : filter,
              'GENETIC_CODE' : genetic_code,
              'MAT_PARAM' : mat_param,
              'OTHER_ADVANCED' : other_advanced,
              'NCBI_GI' : ncbi_gi,
              'OVERVIEW' : overview,
              'ALIGNMENT_VIEW' : alignment_view,
              'DESCRIPTIONS' : descriptions,
              'ALIGNMENTS' : alignments,
              'EMAIL' : email,
              'PATH' : path,
              'HTML' : html
              }
    variables = {}
    for k in params.keys():
        if params[k] is not None:
            variables[k] = str(params[k])
    # This returns a handle to the HTML file that points to the results.
    handle = NCBI._open(cgi, variables, get=0)
    # Now parse the HTML from the handle and figure out how to retrieve
    # the results.
    refcgi, params = _parse_blast_ref_page(handle, cgi)

    start = time.time()
    while 1:
        # Sometimes the BLAST results aren't done yet.  Look at the page
        # to see if the results are there.  If not, then try again later.
        handle = NCBI._open(cgi, params, get=0)
        ready, results, refresh_delay = _parse_blast_results_page(handle)
        if ready:
            break
        # Time out if it's not done after timeout minutes.
        if time.time() - start > timeout*60:
            raise IOError, "timed out after %d minutes" % timeout
        # pause and try again.
        time.sleep(refresh_delay)
    return File.UndoHandle(File.StringHandle(results))

def _parse_blast_ref_page(handle, base_cgi):
    """_parse_blast_ref_page(handle, base_cgi) -> cgi, parameters"""
    # I can speed things up by putting the class declarations into the
    # module scope, instead of recreating them in every function call.
    # However, since the running time for the blast call will be dominated
    # by NCBI's BLAST, it probably won't make much of a difference.
    # This way, the implementation details are hidden in the function.
    class RefPageParser(sgmllib.SGMLParser):
        def __init__(self, cgi):
            sgmllib.SGMLParser.__init__(self)
            self.cgi = cgi
            self.params = {}
        def do_form(self, attributes):
            # parse the "FORM" tag to see where the CGI script should be.
            for attr, value in attributes:
                attr = string.upper(attr)
                if attr == 'ACTION':
                    self.cgi = urlparse.urljoin(self.cgi, value)
        def do_input(self, attributes):
            # parse the "INPUT" tags to try and find the reference ID (RID)
            is_rid = 0
            rid = None
            for attr, value in attributes:
                attr, value = string.upper(attr), string.upper(value)
                if attr == 'NAME' and value == 'RID':
                    is_rid = 1
                elif attr == 'VALUE':
                    rid = value
            if is_rid and rid:
                self.params['RID'] = rid
                
    parser = RefPageParser(base_cgi)
    parser.feed(handle.read())
    if not parser.params.has_key('RID'):
        raise SyntaxError, "Error getting BLAST results: RID not found"
    return parser.cgi, parser.params
    
def _parse_blast_results_page(handle):
    """_parse_blast_results_page(handle) -> ready, results, refresh_delay"""
    class ResultsParser(sgmllib.SGMLParser):
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self.ready = 0
            self.refresh = 5
        def handle_comment(self, comment):
            comment = string.lower(comment)
            if string.find(comment, 'status=ready') >= 0:
                self.ready = 1
        _refresh_re = re.compile('REFRESH_DELAY=(\d+)', re.IGNORECASE)
        def do_meta(self, attributes):
            for attr, value in attributes:
                m = self._refresh_re.search(value)
                if m:
                    self.refresh = int(m.group(1))
    results = handle.read()
    parser = ResultsParser()
    parser.feed(results)
    return parser.ready, results, parser.refresh

