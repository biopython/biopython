# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Patched by Brad Chapman.
# Chris Wroe added modifications for work in myGrid

"""
This module provides code to work with the WWW version of BLAST
provided by the NCBI.
http://blast.ncbi.nlm.nih.gov/

Functions:
qblast        Do a BLAST search using the QBLAST API.
"""

import sys
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

from Bio._py3k import _as_string, _as_bytes

def qblast(program, database, sequence,
           auto_format=None,composition_based_statistics=None,
           db_genetic_code=None,endpoints=None,entrez_query='(none)',
           expect=10.0,filter=None,gapcosts=None,genetic_code=None,
           hitlist_size=50,i_thresh=None,layout=None,lcase_mask=None,
           matrix_name=None,nucl_penalty=None,nucl_reward=None,
           other_advanced=None,perc_ident=None,phi_pattern=None,
           query_file=None,query_believe_defline=None,query_from=None,
           query_to=None,searchsp_eff=None,service=None,threshold=None,
           ungapped_alignment=None,word_size=None,
           alignments=500,alignment_view=None,descriptions=500,
           entrez_links_new_window=None,expect_low=None,expect_high=None,
           format_entrez_query=None,format_object=None,format_type='XML',
           ncbi_gi=None,results_file=None,show_overview=None, megablast=None,
           ):
    """Do a BLAST search using the QBLAST server at NCBI.

    Supports all parameters of the qblast API for Put and Get.
    Some useful parameters:
    program        blastn, blastp, blastx, tblastn, or tblastx (lower case)
    database       Which database to search against (e.g. "nr").
    sequence       The sequence to search.
    ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
    descriptions   Number of descriptions to show.  Def 500.
    alignments     Number of alignments to show.  Def 500.
    expect         An expect value cutoff.  Def 10.0.
    matrix_name    Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
    filter         "none" turns off filtering.  Default no filtering
    format_type    "HTML", "Text", "ASN.1", or "XML".  Def. "XML".
    entrez_query   Entrez query to limit Blast search
    hitlist_size   Number of hits to return. Default 50
    megablast      TRUE/FALSE whether to use MEga BLAST algorithm (blastn only)
    service        plain, psi, phi, rpsblast, megablast (lower case)

    This function does no checking of the validity of the parameters
    and passes the values to the server as is.  More help is available at:
    http://www.ncbi.nlm.nih.gov/BLAST/blast_overview.html

    """
    import urllib, urllib2
    import time

    assert program in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    # Format the "Put" command, which sends search requests to qblast.
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node5.html on 9 July 2007
    # Additional parameters are taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node9.html on 8 Oct 2010
    # To perform a PSI-BLAST or PHI-BLAST search the service ("Put" and "Get" commands) must be specified
    # (e.g. psi_blast = NCBIWWW.qblast("blastp", "refseq_protein", input_sequence, service="psi"))
    parameters = [
        ('AUTO_FORMAT',auto_format),
        ('COMPOSITION_BASED_STATISTICS',composition_based_statistics),
        ('DATABASE',database),
        ('DB_GENETIC_CODE',db_genetic_code),
        ('ENDPOINTS',endpoints),
        ('ENTREZ_QUERY',entrez_query),
        ('EXPECT',expect),
        ('FILTER',filter),
        ('GAPCOSTS',gapcosts),
        ('GENETIC_CODE',genetic_code),
        ('HITLIST_SIZE',hitlist_size),
        ('I_THRESH',i_thresh),
        ('LAYOUT',layout),
        ('LCASE_MASK',lcase_mask),
        ('MEGABLAST',megablast),
        ('MATRIX_NAME',matrix_name),
        ('NUCL_PENALTY',nucl_penalty),
        ('NUCL_REWARD',nucl_reward),
        ('OTHER_ADVANCED',other_advanced),
        ('PERC_IDENT',perc_ident),
        ('PHI_PATTERN',phi_pattern),
        ('PROGRAM',program),
        #('PSSM',pssm), - It is possible to use PSI-BLAST via this API?
        ('QUERY',sequence),
        ('QUERY_FILE',query_file),
        ('QUERY_BELIEVE_DEFLINE',query_believe_defline),
        ('QUERY_FROM',query_from),
        ('QUERY_TO',query_to),
        #('RESULTS_FILE',...), - Can we use this parameter?
        ('SEARCHSP_EFF',searchsp_eff),
        ('SERVICE',service),
        ('THRESHOLD',threshold),
        ('UNGAPPED_ALIGNMENT',ungapped_alignment),
        ('WORD_SIZE',word_size),
        ('CMD', 'Put'),
        ]
    query = [x for x in parameters if x[1] is not None]
    message = _as_bytes(urllib.urlencode(query))

    # Send off the initial query to qblast.
    # Note the NCBI do not currently impose a rate limit here, other
    # than the request not to make say 50 queries at once using multiple
    # threads.
    request = urllib2.Request("http://blast.ncbi.nlm.nih.gov/Blast.cgi",
                              message,
                              {"User-Agent":"BiopythonClient"})
    handle = urllib2.urlopen(request)

    # Format the "Get" command, which gets the formatted results from qblast
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007	
    rid, rtoe = _parse_qblast_ref_page(handle)
    parameters = [
        ('ALIGNMENTS',alignments),
        ('ALIGNMENT_VIEW',alignment_view),
        ('DESCRIPTIONS',descriptions),
        ('ENTREZ_LINKS_NEW_WINDOW',entrez_links_new_window),
        ('EXPECT_LOW',expect_low),
        ('EXPECT_HIGH',expect_high),
        ('FORMAT_ENTREZ_QUERY',format_entrez_query),
        ('FORMAT_OBJECT',format_object),
        ('FORMAT_TYPE',format_type),
        ('NCBI_GI',ncbi_gi),
        ('RID',rid),
        ('RESULTS_FILE',results_file),
        ('SERVICE',service),
        ('SHOW_OVERVIEW',show_overview),
        ('CMD', 'Get'),
        ]
    query = [x for x in parameters if x[1] is not None]
    message = _as_bytes(urllib.urlencode(query))

    # Poll NCBI until the results are ready.  Use a 3 second wait
    delay = 3.0
    previous = time.time()
    while True:
        current = time.time()
        wait = previous + delay - current
        if wait > 0:
            time.sleep(wait)
            previous = current + wait
        else:
            previous = current

        request = urllib2.Request("http://blast.ncbi.nlm.nih.gov/Blast.cgi",
                                  message,
                                  {"User-Agent":"BiopythonClient"})
        handle = urllib2.urlopen(request)
        results = _as_string(handle.read())

        # Can see an "\n\n" page while results are in progress,
        # if so just wait a bit longer...
        if results=="\n\n":
            continue
        # XML results don't have the Status tag when finished
        if results.find("Status=") < 0:
            break
        i = results.index("Status=")
        j = results.index("\n", i)
        status = results[i+len("Status="):j].strip()
        if status.upper() == "READY":
            break

    return StringIO(results)

def _parse_qblast_ref_page(handle):
    """Extract a tuple of RID, RTOE from the 'please wait' page (PRIVATE).

    The NCBI FAQ pages use TOE for 'Time of Execution', so RTOE is proably
    'Request Time of Execution' and RID would be 'Request Identifier'.
    """
    s = _as_string(handle.read())
    i = s.find("RID =")
    if i == -1:
        rid = None
    else:
        j = s.find("\n", i)
        rid = s[i+len("RID ="):j].strip()

    i = s.find("RTOE =")
    if i == -1:
        rtoe = None
    else:
        j = s.find("\n", i)
        rtoe = s[i+len("RTOE ="):j].strip()

    if not rid and not rtoe:
        #Can we reliably extract the error message from the HTML page?
        #e.g.  "Message ID#24 Error: Failed to read the Blast query:
        #       Nucleotide FASTA provided for protein sequence"
        #or    "Message ID#32 Error: Query contains no data: Query
        #       contains no sequence data"
        #
        #This used to occur inside a <div class="error msInf"> entry:
        i = s.find('<div class="error msInf">')
        if i != -1:
            msg = s[i+len('<div class="error msInf">'):].strip()
            msg = msg.split("</div>",1)[0].split("\n",1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        #In spring 2010 the markup was like this:
        i = s.find('<p class="error">')
        if i != -1:
            msg = s[i+len('<p class="error">'):].strip()
            msg = msg.split("</p>",1)[0].split("\n",1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        #Generic search based on the way the error messages start:
        i = s.find('Message ID#')
        if i != -1:
            #Break the message at the first HTML tag
            msg = s[i:].split("<",1)[0].split("\n",1)[0].strip()
            raise ValueError("Error message from NCBI: %s" % msg)
        #We didn't recognise the error layout :(
        #print s
        raise ValueError("No RID and no RTOE found in the 'please wait' page, "
                         "there was probably an error in your request but we "
                         "could not extract a helpful error message.")
    elif not rid:
        #Can this happen?
        raise ValueError("No RID found in the 'please wait' page."
                         " (although RTOE = %s)" % repr(rtoe))
    elif not rtoe:
        #Can this happen?
        raise ValueError("No RTOE found in the 'please wait' page."
                         " (although RID = %s)" % repr(rid))

    try:
        return rid, int(rtoe)
    except ValueError:
        raise ValueError("A non-integer RTOE found in " \
                         +"the 'please wait' page, %s" % repr(rtoe))

 
