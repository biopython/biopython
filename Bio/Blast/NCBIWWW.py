# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Patched by Brad Chapman.

"""NCBIWWW.py

This module provides code to work with the WWW version of BLAST
provided by the NCBI.
http://www.ncbi.nlm.nih.gov/BLAST/

Classes:
BlastParser   Parses output from WWW blast.
_Scanner      Scans output from NCBI's BLAST WWW server.

Functions:
blast         Do a BLAST search against the WWW page.
blasturl      Do a BLAST search against the stable blasturl.

"""
import string
import time
import re
import sgmllib
import urlparse
import socket
import cStringIO

from Bio import File
from Bio.WWW import NCBI
from Bio.ParserSupport import *
import NCBIStandalone

class BlastParser(AbstractParser):
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
        read_and_call_until(uhandle, consumer.noevent,
                            has_re=re.compile(r'<b>.?BLAST'))

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
        # <PRE>  XXX
        consumer.start_header()

        # Read the "BLAST" version line and the following blanks.
        read_and_call(uhandle, consumer.version, contains='BLAST')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        # Read the reference lines and the '<p>' line.
        read_and_call_until(uhandle, consumer.reference, start='<p>')
        read_and_call(uhandle, consumer.noevent)

        # Read the RID line, for version 2.0.12 (2.0.11?) and above.
        attempt_read_and_call(uhandle, consumer.noevent, start='RID')
        # Brad Chapman noticed a '<p>' line in BLASTN 2.1.1
        attempt_read_and_call(uhandle, consumer.noevent, start='<p>')

        # Apparently, there's some discrepancy between whether the
        # Query or database comes first.  Usually the Query does, but
        # Brad noticed a case where the database came first.
        if uhandle.peekline().find("Query=") >= 0:
            self._scan_query_info(uhandle, consumer)
            self._scan_database_info(uhandle, consumer)
        else:
            self._scan_database_info(uhandle, consumer)
            self._scan_query_info(uhandle, consumer)
        consumer.end_header()

    def _scan_database_info(self, uhandle, consumer):
        attempt_read_and_call(uhandle, consumer.noevent, start='<p>')
        read_and_call(uhandle, consumer.database_info, contains='Database')
        read_and_call(uhandle, consumer.database_info, contains='sequences')
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent,
                      contains='problems or questions')
        if attempt_read_and_call(uhandle, consumer.noevent,
                                 contains="BLASTFORM"):
            while 1:
                line = uhandle.peekline()
                if is_blank_line(line):
                    break
                elif string.find(line, "Query=") >= 0:
                    break
                consumer.noevent(uhandle.readline())
        if attempt_read_and_call(uhandle, consumer.noevent,
                                 contains="Taxonomy reports"):
            read_and_call(uhandle, consumer.noevent, start="<BR>")
        attempt_read_and_call(uhandle, consumer.noevent, start="<PRE>")
        read_and_call_while(uhandle, consumer.noevent, blank=1)

    def _scan_query_info(self, uhandle, consumer):
        # Read the Query lines and the following blank line.
        read_and_call(uhandle, consumer.query_info, contains='Query=')
        read_and_call_until(uhandle, consumer.query_info, blank=1)
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        if attempt_read_and_call(uhandle, consumer.noevent, start="<PRE>"):
            read_and_call_while(uhandle, consumer.noevent, blank=1)
            
        
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
            uhandle, consumer.description_header, contains='Score     E'):
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
        read_and_call(uhandle, consumer.description_header,
                      start='Sequences producing')
        read_and_call(uhandle, consumer.noevent, blank=1)

        # Read the descriptions
        read_and_call_while(uhandle, consumer.description, blank=0, start='<a')

        # two choices here, either blank lines or a </PRE>
        if not attempt_read_and_call(uhandle, consumer.noevent,
                                     contains='</PRE>'):
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
        
        # 6/3/2001, </PRE> is gone, replaced by </form>
        

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

        attempt_read_and_call(uhandle, consumer.query_length,
                              start='length of query')
        read_and_call(uhandle, consumer.database_length,
                      start='length of database')

        read_and_call(uhandle, consumer.effective_hsp_length,
                      start='effective HSP')
        attempt_read_and_call(uhandle, consumer.effective_query_length,
                              start='effective length of query')
        read_and_call(uhandle, consumer.effective_database_length,
                      start='effective length of database')
        attempt_read_and_call(uhandle, consumer.effective_search_space,
                              start='effective search space:')
        attempt_read_and_call(uhandle, consumer.effective_search_space_used,
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
        attempt_read_and_call(uhandle, consumer.blast_cutoff, start='S2')

        read_and_call(uhandle, consumer.noevent, blank=1)
        attempt_read_and_call(uhandle, consumer.noevent, start="</PRE>")
        attempt_read_and_call(uhandle, consumer.noevent, start="</form>")

        consumer.end_parameters()

def blast(program, database, query,
          query_from='', query_to='',
          entrez_query='(none)',
          filter='L',
          expect='10',
          word_size=None,
          ungapped_alignment='no',    # deprecated, not on webpage anymore
          other_advanced=None,
          cdd_search='on',
          composition_based_statistics=None,
          matrix_name=None,
          run_psiblast=None,
          i_thresh='0.001',
          genetic_code='1',
          show_overview='on',
          ncbi_gi='on',
          format_object='alignment',
          format_type='html',
          descriptions='100',
          alignments='50',
          alignment_view='Pairwise',
          auto_format='on',
          cgi='http://www.ncbi.nlm.nih.gov/blast/Blast.cgi',
          timeout=20, output_fn=None):
    
    """blast(program, database, query[, query_from][, query_to]
    [, entrez_query][, filter][, expect]
    [, word_size][, other_advanced][, cdd_search]
    [, composition_based_statistics][, matrix_name][, run_psiblast]
    [, i_thresh][, genetic_code][, show_overview][, ncbi_gi]
    [, format_object][, format_type][, descriptions][, alignments]
    [, alignment_view][, auto_format][, cgi][, timeout]) -> handle

    Blast against the NCBI Blast web page.  This uses the NCBI web
    page cgi script to BLAST, and returns a handle to the
    results. See:
    
    http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.html
    
    for more descriptions about the options.

    Required Inputs:
    o program - The name of the blast program to run (ie. blastn, blastx...)
    o database - The database to search against (ie. nr, dbest...)
    o query - The input for the search, which NCBI tries to autodetermine
    the type of. Ideally, this would be a sequence in FASTA format.

    General Options:
    filter, expect, word_size, other_advanced

    Formatting Options:
    show_overview, ncbi_gi, format_object, format_type, descriptions,
    alignments, alignment_view, auto_format

    Protein specific options:
    cdd_search, composition_based_statistics, matrix_name, run_psiblast,
    i_thresh

    Translated specific options:
    genetic code
    
    """
    # NCBI Blast is hard to work with.  The user enters a query, and then
    # it returns a "reference" page which contains a button that the user
    # clicks to retrieve the results.  This will retrieve the "results"
    # page.  However, this page may not contain BLAST results if the
    # search isn't done.
    # This function will send off the query and parse the reference
    # page to figure out how to retrieve the results.  Then, it needs to
    # periodically query the results to see if the search has finished.
    # When it has, then it can retrieve the actual blast results.
    params = {'PROGRAM' : program,
              'QUERY_FROM' : query_from,
              'QUERY_TO' : query_to,
              'DATABASE' : database,
              'QUERY' : query,
              'ENTREZ_QUERY' : entrez_query,
              'FILTER' : filter,
              'EXPECT' : expect,
              'WORD_SIZE' : word_size,
              'OTHER_ADVANCED': other_advanced,
              'CDD_SEARCH' : cdd_search,
              'COMPOSITION_BASED_STATISTICS' : composition_based_statistics,
              'MATRIX_NAME' : matrix_name,
              'RUN_PSIBLAST' : run_psiblast,
              'I_THRESH' : i_thresh,
              'GENETIC_CODE' : genetic_code,
              'SHOW_OVERVIEW' : show_overview,
              'NCBI_GI' : ncbi_gi,
              'FORMAT_OBJECT' : format_object,
              'FORMAT_TYPE' : format_type,
              'DESCRIPTIONS' : descriptions,
              'ALIGNMENTS' : alignments,
              'ALIGNMENT_VIEW' : alignment_view,
              'AUTO_FORMAT' : auto_format}

    default_word_sizes = {
        'blastp' : 3,
        'blastn' : 11,
        'blastx' : 3,
        'tblastn' : 3,
        'tblastx' : 3
        }
    if not params['WORD_SIZE']:
        params['WORD_SIZE'] = default_word_sizes.get(params['PROGRAM'], 3)
    
    variables = {}
    for k in params.keys():
        if params[k] is not None:
            variables[k] = str(params[k])
            
    variables['CLIENT'] = 'web'
    variables['SERVICE'] = 'plain'
    variables['CMD'] = 'Put'
    variables['LAYOUT'] = 'OneWindow'

    if program.upper() == 'BLASTN':
        variables['PAGE'] = 'Nucleotides'
    elif program.upper() == 'BLASTP':
        variables['PAGE'] = 'Proteins'
    elif program.upper() in ['BLASTX', 'TBLASTN','TBLASTX']:
        variables['PAGE'] = 'Translations'
    else:
        raise ValueError("Unexpected program name %s" % program)

    # These parameters are not yet implemented.
    # LCASE_MASK=''
    # GAPCOSTS=''
    # PSSM=''
    # PHI_PATTERN=''
    # FORMAT_BLOCK_ON_RESPAGE='None'
    # EMAIL_ADDRESS=''

    # This returns a handle to the HTML file that points to the results.
    handle = NCBI._open(cgi, variables, get=0)
    # Now parse the HTML from the handle and figure out how to retrieve
    # the results.
    if output_fn is not None:
        results = handle.read()
        output_fn(results)
        handle = File.StringHandle(results)
    ref_cgi, ref_params = _parse_blast_ref_page(handle)
    ref_cgi = urlparse.urljoin(cgi, ref_cgi)  # convert to absolute URL

    # Start with the initial recommended delay.
    refresh_delay = int(ref_params.get("RTOE", 5))

    start = time.time()
    while 1:
        # pause before trying to get the results
        time.sleep(refresh_delay)
        
        # Sometimes the BLAST results aren't done yet.  Look at the page
        # to see if the results are there.  If not, then try again later.
        handle = NCBI._open(ref_cgi, ref_params, get=0)
        if output_fn is not None:
            results = handle.read()
            output_fn(results)
            handle = File.StringHandle(results)
        ready, results_cgi, results_params = _parse_blast_results_page(handle)
        results_cgi = urlparse.urljoin(cgi, results_cgi)    # to absolute URL
        if ready:
            break
        # Time out if it's not done after timeout minutes.
        if time.time() - start > timeout*60:
            raise IOError, "timed out after %d minutes" % timeout

    # Now query for the actual results.  To do this, the CGI script
    # needs CMD="Get", which should already be in results_params.
    # Also, for some reason, this fails if FORMAT_OBJECT is in
    # results_params, so we need to get rid of it.
    if results_params.has_key("FORMAT_OBJECT"):
        del results_params["FORMAT_OBJECT"]
    return NCBI._open(results_cgi, results_params, get=0)

class _FormParser(sgmllib.SGMLParser):
    """Parse a form in an HTML page.

    Members:
    forms   List of forms in the page.
            Each form is a tuple of (action, params) where action
            is a string to the CGI script and params is a dict of
            keys and values to pass to the script.

    """
    def __init__(self):
        sgmllib.SGMLParser.__init__(self)
        self.forms = []
        self._current_form = '', {}
    def start_form(self, attributes):
        # Parse the "FORM" tag to see where the CGI script is.
        attr_dict = self._attr2dict(attributes)
        self._current_form = (attr_dict.get('ACTION', self._current_form[0]),
                              self._current_form[1])
    def end_form(self):
        self.forms.append(self._current_form)
        self._current_form = '', {}
    def do_input(self, attributes):
        params = self._current_form[1]
        attr_dict = self._attr2dict(attributes)
        if attr_dict.has_key('NAME'):
            # Get the value, handling check boxes.
            value = attr_dict.get('VALUE', attr_dict.get('CHECKED', ''))
            params[attr_dict['NAME']] = value
    def _attr2dict(self, attributes):
        attr_dict = {}
        for name, value in attributes:
            attr_dict[string.upper(name)] = value
        return attr_dict
    # XXX should handle SELECT, not implemented yet        

def _parse_blast_ref_page(handle):
    """_parse_blast_ref_page(handle, base_cgi) -> cgi, parameters"""
    parser = _FormParser()
    parser.feed(handle.read())
    if len(parser.forms) != 1:
        raise SyntaxError, "Form broken in BLAST reference page"
    cgi, params = parser.forms[0]
    if not params.has_key('RID'):
        raise SyntaxError, "Error getting BLAST results: RID not found"
    return cgi, params
    
def _parse_blast_results_page(handle):
    """_parse_blast_results_page(handle) -> ready, cgi, params"""
    class _ResultsParser(_FormParser):
        def __init__(self):
            _FormParser.__init__(self)
            self.ready = 0
        #_refresh_re = re.compile(r",\s*(\d+)\s*\);")
        def handle_comment(self, comment):
            # There is lots of information in the comments of the results
            # page:
            # <!--
            # QBlastInfoBegin
            #         Status=WAITING
            # QBlastInfoEnd
            # -->
            # <SCRIPT LANGUAGE="JavaScript"><!--
            # setTimeout('document.forms[0].submit();',15000);
            # //--></SCRIPT>
            comment = string.lower(comment)
            if string.find(comment, 'status=ready') >= 0:
                self.ready = 1
            #elif string.find(comment, 'settimeout') >= 0:
            #    # parse the refresh delay out of the comment
            #    m = self._refresh_re.search(comment)
            #    assert m, "Failed to parse refresh time from %s" % comment
            #    self.refresh = int(m.group(1))/1000  # give in milliseconds
    parser = _ResultsParser()
    parser.feed(handle.read())
    
    # The results page has 2 forms.  The first one is used if the
    # results are ready.  Otherwise, return the second one.
    if len(parser.forms) != 2:
        raise SyntaxError, "I expected 2 forms in the results page."
    if parser.ready:
        cgi, params = parser.forms[0]
    else:
        cgi, params = parser.forms[1]
    return parser.ready, cgi, params

def blasturl(program, datalib, sequence,
             ncbi_gi=None, descriptions=None, alignments=None,
             expect=None, matrix=None,
             gap_existence=None, gap_extend=None, gapped=None,
             filter=None, html=None, gcode=None, path=None
             ):
    """blasturl(program, datalib, sequence[, ncbi_gi][, descriptions]
    [, alignments][, expect][, matrix][, gap_existence][, gap_extend]
    [, gapped][, filter][, html][, gcode]) -> handle

    Do a BLAST search using the stable URL provided by NCBI.
    program        BLASTP, BLASTN, BLASTX, TBLASTN, or TBLASTX.
    datalib        Which database to search against.
    sequence       The sequence to search.
    ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.  Def FALSE.
    descriptions   Number of descriptions to show.  Def 100.
    alignments     Number of alignments to show.  Def 50.
    expect         An expect value cutoff.
    matrix         Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
    gap_existence  Give a gap open penalty.
    gap_extend     Give a gap extension penalty.
    gapped         TRUE/FALSE for giving gapped alignments.  Def TRUE.
    filter         "none" turns off filtering.  Default uses 'seg' or 'dust'.
    html           TRUE/FALSE for html output.  Def FALSE.
    gcode          Specify an alternate genetic code for (T)BLASTX.

    This function does no checking of the validity of the parameters
    and passes the values to the server as is.  More help is available at:
    http://www.ncbi.nlm.nih.gov/BLAST/blast_overview.html
    
    """
    lines = []
    lines.append('PROGRAM %s' % program)
    lines.append('DATALIB %s' % datalib)

    parameters = [('NCBI_GI', ncbi_gi),
                  ('DESCRIPTIONS', descriptions),
                  ('ALIGNMENTS', alignments),
                  ('EXPECT', expect),
                  ('MATRIX', matrix),
                  ('GAP_EXISTENCE', gap_existence),
                  ('GAP_EXTEND', gap_extend),
                  ('GAPPED', gapped),
                  ('FILTER', filter),
                  ('HTML', html),
                  ('GCODE', gcode),
                  ('PATH', path)
                  ]
    for name, value in parameters:
        if value is not None:
            lines.append("%s %s" % (name, value))

    lines.append('')
    lines.append('BEGIN')
    while sequence:
        lines.append(sequence[:60])
        sequence = sequence[60:]

    message = string.join(lines, '\n')

    outhandle = cStringIO.StringIO()
    _send_to_blasturl(message, outhandle)
    outhandle.seek(0)   # Reset the handle to the beginning.
    return outhandle

def _send_to_blasturl(query, outhandle):
    """_send_to_blasturl(query, outhandle)

    Send a BLAST request to the stable blasturl server at the NCBI.
    ftp://ncbi.nlm.nih.gov/blast/blasturl/
    The results are written to outhandle.
    
    """
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect(('www.ncbi.nlm.nih.gov', 80))
    
    sock.send('POST /cgi-bin/BLAST/nph-blast_report HTTP/1.0\n')
    sock.send('User-Agent: BiopythonClient\n')
    sock.send('Connection: Keep-Alive\n')
    sock.send('Content-type: application/x-www-form-urlencoded\n')
    sock.send('Content-Length: %d\n' % len(query))
    sock.send('\n')
    sock.send(query)

    while 1:
        data = sock.recv(1024)
        if not data:
            break
        outhandle.write(data)
    sock.close()
