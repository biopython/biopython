# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""NCBIWWW.py

This module provides code to work with the WWW version of BLAST
provided by the NCBI.
http://www.ncbi.nlm.nih.gov/BLAST/

Classes:
Scanner      Scans output from NCBI's BLAST WWW server.
"""


import string

from Bio import File
from Bio.ParserSupport import *


class Scanner:
    """Scan BLAST output from NCBI's web server at:
    http://www.ncbi.nlm.nih.gov/BLAST/
    
    Tested with BLAST v2.0.10

    Methods:
    feed     Feed data into the scanner.
    """
    
    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a BLAST report for scanning.  handle is a file-like
        object that contains the BLAST report.  consumer is a Consumer
        object that will receive events as the report is scanned.

        """
        # <HTML>
        # <HEAD>
        # <TITLE>BLAST Search Results </TITLE>
        # </HEAD>
        # <BODY BGCOLOR="#FFFFFF" LINK="#0000FF" VLINK="#660099" ALINK="#660099
        # <A HREF="http://www.ncbi.nlm.nih.gov/BLAST/blast_form.map"> <IMG SRC=

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
        read_and_call(uhandle, consumer.noevent, start='<HTML>')
        read_and_call(uhandle, consumer.noevent, start='<HEAD>')
        read_and_call(uhandle, consumer.noevent, start='<TITLE>')
        read_and_call(uhandle, consumer.noevent, start='</HEAD>')
        read_and_call(uhandle, consumer.noevent, start='<BODY')
        read_and_call(uhandle, consumer.noevent, start='<A HREF')

        self._scan_header(uhandle, consumer)
	self._scan_rounds(uhandle, consumer)
        self._scan_database_report(uhandle, consumer)
        read_and_call(uhandle, consumer.noevent, blank=1)
        self._scan_parameters(uhandle, consumer)

        # Read HTML footer information.
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start='</BODY>')
        read_and_call(uhandle, consumer.noevent, start='</HTML>')
        read_and_call(uhandle, consumer.noevent, start='</BODY>')
        read_and_call(uhandle, consumer.noevent, start='</HTML>')

    def _scan_header(self, uhandle, consumer):
        # <BR><BR><PRE>
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

        # Read the "BLAST" version line and the two following blanks.
        read_and_call(uhandle, consumer.noevent, contains='PRE')
        read_and_call(uhandle, consumer.version, contains='BLAST')
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, blank=1)

        # Read the reference lines and the '<p>' line.
        read_and_call(uhandle, consumer.reference, start='<b><a href=')
        while 1:
            line = safe_readline(uhandle)
            if line[:3] == '<p>':
                consumer.noevent(line)
                break
            consumer.reference(line)

        # Read the Query lines and the following blank line.
        read_and_call(uhandle, consumer.query_info, contains='Query=')
        while 1:
            line = safe_readline(uhandle)
            if is_blank_line(line):
                consumer.noevent(line)
                break
            consumer.query_info(line)

        # Read the database lines and the following blank line.
        read_and_call(uhandle, consumer.database_info, contains='Database')
        read_and_call(uhandle, consumer.database_info, contains='sequences')
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent,
                      contains='problems or questions')

        # Read the blast form, if it exists. 
        if attempt_read_and_call(uhandle, consumer.noevent,
                                 contains='BLASTFORM'):
            while 1:
                line = safe_readline(uhandle)
                consumer.noevent(line)
                if line[:5] == '<PRE>':
                    break

        consumer.end_header()

    def _scan_rounds(self, uhandle, consumer):
        self._scan_descriptions(uhandle, consumer)
        self._scan_alignments(uhandle, consumer)

    def _scan_descriptions(self, uhandle, consumer):
        consumer.start_descriptions()

        line = safe_peekline(uhandle)
        if string.find(line, 'No significant similarity') >= 0:
            # no hits found:
            # <b>No significant similarity found.</b> For reasons why, <A HREF
            read_and_call(uhandle, consumer.no_hits)
        elif is_blank_line(line):
            # no descriptions:
            # 
            # 
            read_and_call(uhandle, consumer.noevent, blank=1)
            read_and_call(uhandle, consumer.noevent, blank=1)
        else:
            # Normal descriptions:
            # <PRE>
            # 
            #                                                                  
            # Sequences producing significant alignments:                      
            # 
            # <a href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Ret
            # <a href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Ret
            # 
            read_and_call(uhandle, consumer.noevent, start='<PRE>')
            read_and_call(uhandle, consumer.noevent, blank=1)

            # Read the score header lines and a blank line.
            read_and_call(uhandle, consumer.noevent,
                          contains='Score     E')
            read_and_call(uhandle, consumer.noevent,
                          start='Sequences producing')
            read_and_call(uhandle, consumer.noevent, blank=1)

            # Read the descriptions and the following blank line.
            while 1:
                if not attempt_read_and_call(uhandle, consumer.description,
                                             start='<a href'):
                    read_and_call(uhandle, consumer.noevent, blank=1)
                    break

        consumer.end_descriptions()

    def _scan_alignments(self, uhandle, consumer):
        # An alignment starts at a <PRE>.
        # If I'm at a blank line, then there's no alignment and I'm
        # at a database report.
        line1 = safe_readline(uhandle)
        line2 = safe_readline(uhandle)
        uhandle.saveline(line2)
        uhandle.saveline(line1)
        if is_blank_line(line1) or line2[:10] == '  Database':
            return

        # It appears that first sequence in a masterslave alignment
        # is generated by BLAST and contains no link to the descriptions.
        if line2[:9] == 'blast_tmp':
            self._scan_masterslave_alignment(uhandle, consumer)
        else:
            self._scan_pairwise_alignments(uhandle, consumer)

    def _scan_pairwise_alignments(self, uhandle, consumer):
        while 1:
            # If I'm at an alignment header, the first line should be
            # <PRE>.  If I'm at the database report, then the
            # first line will be blank.
            line1 = safe_readline(uhandle)
            line2 = safe_readline(uhandle)
            uhandle.saveline(line2)
            uhandle.saveline(line1)
            if line1[:5] != '<PRE>':
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
        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        self._scan_hsp_alignment(uhandle, consumer)
        consumer.end_alignment()
            
    def _scan_one_pairwise_alignment(self, uhandle, consumer):
        # Alignment format:
        # <PRE>
        # alignment_header
        #   hsp_header
        #   hsp_alignment
        #   [...]
        # The hsp_header and hsp_alignment blocks can be repeated.

        consumer.start_alignment()
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
            index = string.find(line, 'Length =')
            # if index == 10 or index == 11 or index == 12:
            if index >= 10:
                consumer.length(line)
                break
            elif is_blank_line(line):
                # Check to make sure I haven't missed the Length line
                raise SyntaxError, "I missed the Length in an alignment header"
            consumer.title(line)

        read_and_call(uhandle, consumer.noevent, start='          ')

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
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, blank=1)

    def _scan_masterslave_alignment(self, uhandle, consumer):
        consumer.start_alignment()
        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        while 1:
            line = safe_readline(uhandle)
            if is_blank_line(line):
                consumer.noevent(line)
                # # If the blank line is followed by '<PRE>', then
                # # the blank line belongs to the database report.
                # line2 = safe_readline(uhandle)
                # if line2[:5] == '<PRE>':
                #    uhandle.saveline(line2)
                #    uhandle.saveline(line)
                #    break
                #consumer.noevent(line)
                #uhandle.saveline(line2)
            elif line[:6] == '</PRE>':
                consumer.noevent(line)
                break
            else:
                consumer.multalign(line)
        consumer.end_alignment()

    def _scan_database_report(self, uhandle, consumer):
        # 
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

        # There's no blank line if there were no hits!
        attempt_read_and_call(uhandle, consumer.noevent, blank=1)
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
            read_and_call(uhandle, consumer.noevent, blank=1)

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

