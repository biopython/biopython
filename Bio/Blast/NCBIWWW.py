# Copyright 1999 by Jeffrey Chang.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Patched by Brad Chapman.
# Chris Wroe added modifications for work in myGrid

"""Code to invoke the NCBI BLAST server over the internet.

This module provides code to work with the WWW version of BLAST
provided by the NCBI. https://blast.ncbi.nlm.nih.gov/
"""


import warnings

from io import StringIO
import time

from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.request import Request

from Bio import BiopythonWarning


NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def qblast(
    program,
    database,
    sequence,
    url_base=NCBI_BLAST_URL,
    auto_format=None,
    composition_based_statistics=None,
    db_genetic_code=None,
    endpoints=None,
    entrez_query="(none)",
    expect=10.0,
    filter=None,
    gapcosts=None,
    genetic_code=None,
    hitlist_size=50,
    i_thresh=None,
    layout=None,
    lcase_mask=None,
    matrix_name=None,
    nucl_penalty=None,
    nucl_reward=None,
    other_advanced=None,
    perc_ident=None,
    phi_pattern=None,
    query_file=None,
    query_believe_defline=None,
    query_from=None,
    query_to=None,
    searchsp_eff=None,
    service=None,
    threshold=None,
    ungapped_alignment=None,
    word_size=None,
    short_query=None,
    alignments=500,
    alignment_view=None,
    descriptions=500,
    entrez_links_new_window=None,
    expect_low=None,
    expect_high=None,
    format_entrez_query=None,
    format_object=None,
    format_type="XML",
    ncbi_gi=None,
    results_file=None,
    show_overview=None,
    megablast=None,
    template_type=None,
    template_length=None,
):
    """BLAST search using NCBI's QBLAST server or a cloud service provider.

    Supports all parameters of the old qblast API for Put and Get.

    Please note that NCBI uses the new Common URL API for BLAST searches
    on the internet (http://ncbi.github.io/blast-cloud/dev/api.html). Thus,
    some of the parameters used by this function are not (or are no longer)
    officially supported by NCBI. Although they are still functioning, this
    may change in the future.

    The Common URL API (http://ncbi.github.io/blast-cloud/dev/api.html) allows
    doing BLAST searches on cloud servers. To use this feature, please set
    ``url_base='http://host.my.cloud.service.provider.com/cgi-bin/blast.cgi'``
    and ``format_object='Alignment'``. For more details, please see
    https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=CloudBlast

    Some useful parameters:

     - program        blastn, blastp, blastx, tblastn, or tblastx (lower case)
     - database       Which database to search against (e.g. "nr").
     - sequence       The sequence to search.
     - ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
     - descriptions   Number of descriptions to show.  Def 500.
     - alignments     Number of alignments to show.  Def 500.
     - expect         An expect value cutoff.  Def 10.0.
     - matrix_name    Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
     - filter         "none" turns off filtering.  Default no filtering
     - format_type    "HTML", "Text", "ASN.1", or "XML".  Def. "XML".
     - entrez_query   Entrez query to limit Blast search
     - hitlist_size   Number of hits to return. Default 50
     - megablast      TRUE/FALSE whether to use MEga BLAST algorithm (blastn only)
     - short_query    TRUE/FALSE whether to adjust the search parameters for a
                      short query sequence. Note that this will override
                      manually set parameters like word size and e value. Turns
                      off when sequence length is > 30 residues. Default: None.
     - service        plain, psi, phi, rpsblast, megablast (lower case)

    This function does no checking of the validity of the parameters
    and passes the values to the server as is.  More help is available at:
    https://ncbi.github.io/blast-cloud/dev/api.html

    """
    programs = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]
    if program not in programs:
        raise ValueError(
            "Program specified is %s. Expected one of %s"
            % (program, ", ".join(programs))
        )

    # SHORT_QUERY_ADJUST throws an error when using blastn (wrong parameter
    # assignment from NCBIs side).
    # Thus we set the (known) parameters directly:
    if short_query and program == "blastn":
        short_query = None
        # We only use the 'short-query' parameters for short sequences:
        if len(sequence) < 31:
            expect = 1000
            word_size = 7
            nucl_reward = 1
            filter = None
            lcase_mask = None
            warnings.warn(
                '"SHORT_QUERY_ADJUST" is incorrectly implemented (by NCBI) for blastn.'
                " We bypass the problem by manually adjusting the search parameters."
                " Thus, results may slightly differ from web page searches.",
                BiopythonWarning,
            )

    # Format the "Put" command, which sends search requests to qblast.
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node5.html on 9 July 2007
    # Additional parameters are taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node9.html on 8 Oct 2010
    # To perform a PSI-BLAST or PHI-BLAST search the service ("Put" and "Get" commands) must be specified
    # (e.g. psi_blast = NCBIWWW.qblast("blastp", "refseq_protein", input_sequence, service="psi"))
    parameters = [
        ("AUTO_FORMAT", auto_format),
        ("COMPOSITION_BASED_STATISTICS", composition_based_statistics),
        ("DATABASE", database),
        ("DB_GENETIC_CODE", db_genetic_code),
        ("ENDPOINTS", endpoints),
        ("ENTREZ_QUERY", entrez_query),
        ("EXPECT", expect),
        ("FILTER", filter),
        ("GAPCOSTS", gapcosts),
        ("GENETIC_CODE", genetic_code),
        ("HITLIST_SIZE", hitlist_size),
        ("I_THRESH", i_thresh),
        ("LAYOUT", layout),
        ("LCASE_MASK", lcase_mask),
        ("MEGABLAST", megablast),
        ("MATRIX_NAME", matrix_name),
        ("NUCL_PENALTY", nucl_penalty),
        ("NUCL_REWARD", nucl_reward),
        ("OTHER_ADVANCED", other_advanced),
        ("PERC_IDENT", perc_ident),
        ("PHI_PATTERN", phi_pattern),
        ("PROGRAM", program),
        # ('PSSM',pssm), - It is possible to use PSI-BLAST via this API?
        ("QUERY", sequence),
        ("QUERY_FILE", query_file),
        ("QUERY_BELIEVE_DEFLINE", query_believe_defline),
        ("QUERY_FROM", query_from),
        ("QUERY_TO", query_to),
        # ('RESULTS_FILE',...), - Can we use this parameter?
        ("SEARCHSP_EFF", searchsp_eff),
        ("SERVICE", service),
        ("SHORT_QUERY_ADJUST", short_query),
        ("TEMPLATE_TYPE", template_type),
        ("TEMPLATE_LENGTH", template_length),
        ("THRESHOLD", threshold),
        ("UNGAPPED_ALIGNMENT", ungapped_alignment),
        ("WORD_SIZE", word_size),
        ("CMD", "Put"),
    ]
    query = [x for x in parameters if x[1] is not None]
    message = urlencode(query).encode()

    # Send off the initial query to qblast.
    # Note the NCBI do not currently impose a rate limit here, other
    # than the request not to make say 50 queries at once using multiple
    # threads.
    request = Request(url_base, message, {"User-Agent": "BiopythonClient"})
    handle = urlopen(request)

    # Format the "Get" command, which gets the formatted results from qblast
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007
    rid, rtoe = _parse_qblast_ref_page(handle)
    parameters = [
        ("ALIGNMENTS", alignments),
        ("ALIGNMENT_VIEW", alignment_view),
        ("DESCRIPTIONS", descriptions),
        ("ENTREZ_LINKS_NEW_WINDOW", entrez_links_new_window),
        ("EXPECT_LOW", expect_low),
        ("EXPECT_HIGH", expect_high),
        ("FORMAT_ENTREZ_QUERY", format_entrez_query),
        ("FORMAT_OBJECT", format_object),
        ("FORMAT_TYPE", format_type),
        ("NCBI_GI", ncbi_gi),
        ("RID", rid),
        ("RESULTS_FILE", results_file),
        ("SERVICE", service),
        ("SHOW_OVERVIEW", show_overview),
        ("CMD", "Get"),
    ]
    query = [x for x in parameters if x[1] is not None]
    message = urlencode(query).encode()

    # Poll NCBI until the results are ready.
    # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
    # 1. Do not contact the server more often than once every 10 seconds.
    # 2. Do not poll for any single RID more often than once a minute.
    # 3. Use the URL parameter email and tool, so that the NCBI
    #    can contact you if there is a problem.
    # 4. Run scripts weekends or between 9 pm and 5 am Eastern time
    #    on weekdays if more than 50 searches will be submitted.
    # --
    # Could start with a 10s delay, but expect most short queries
    # will take longer thus at least 70s with delay. Therefore,
    # start with 20s delay, thereafter once a minute.
    delay = 20  # seconds
    while True:
        current = time.time()
        wait = qblast._previous + delay - current
        if wait > 0:
            time.sleep(wait)
            qblast._previous = current + wait
        else:
            qblast._previous = current
        # delay by at least 60 seconds only if running the request against the public NCBI API
        if delay < 60 and url_base == NCBI_BLAST_URL:
            # Wasn't a quick return, must wait at least a minute
            delay = 60

        request = Request(url_base, message, {"User-Agent": "BiopythonClient"})
        handle = urlopen(request)
        results = handle.read().decode()

        # Can see an "\n\n" page while results are in progress,
        # if so just wait a bit longer...
        if results == "\n\n":
            continue
        # XML results don't have the Status tag when finished
        if "Status=" not in results:
            break
        i = results.index("Status=")
        j = results.index("\n", i)
        status = results[i + len("Status=") : j].strip()
        if status.upper() == "READY":
            break
    return StringIO(results)


qblast._previous = 0


def _parse_qblast_ref_page(handle):
    """Extract a tuple of RID, RTOE from the 'please wait' page (PRIVATE).

    The NCBI FAQ pages use TOE for 'Time of Execution', so RTOE is probably
    'Request Time of Execution' and RID would be 'Request Identifier'.
    """
    s = handle.read().decode()
    i = s.find("RID =")
    if i == -1:
        rid = None
    else:
        j = s.find("\n", i)
        rid = s[i + len("RID =") : j].strip()

    i = s.find("RTOE =")
    if i == -1:
        rtoe = None
    else:
        j = s.find("\n", i)
        rtoe = s[i + len("RTOE =") : j].strip()

    if not rid and not rtoe:
        # Can we reliably extract the error message from the HTML page?
        # e.g.  "Message ID#24 Error: Failed to read the Blast query:
        #       Nucleotide FASTA provided for protein sequence"
        # or    "Message ID#32 Error: Query contains no data: Query
        #       contains no sequence data"
        #
        # This used to occur inside a <div class="error msInf"> entry:
        i = s.find('<div class="error msInf">')
        if i != -1:
            msg = s[i + len('<div class="error msInf">') :].strip()
            msg = msg.split("</div>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        # In spring 2010 the markup was like this:
        i = s.find('<p class="error">')
        if i != -1:
            msg = s[i + len('<p class="error">') :].strip()
            msg = msg.split("</p>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        # Generic search based on the way the error messages start:
        i = s.find("Message ID#")
        if i != -1:
            # Break the message at the first HTML tag
            msg = s[i:].split("<", 1)[0].split("\n", 1)[0].strip()
            raise ValueError("Error message from NCBI: %s" % msg)
        # We didn't recognise the error layout :(
        # print s
        raise ValueError(
            "No RID and no RTOE found in the 'please wait' page, "
            "there was probably an error in your request but we "
            "could not extract a helpful error message."
        )
    elif not rid:
        # Can this happen?
        raise ValueError(
            "No RID found in the 'please wait' page. (although RTOE = %r)" % rtoe
        )
    elif not rtoe:
        # Can this happen?
        raise ValueError(
            "No RTOE found in the 'please wait' page. (although RID = %r)" % rid
        )

    try:
        return rid, int(rtoe)
    except ValueError:
        raise ValueError(
            "A non-integer RTOE found in the 'please wait' page, %r" % rtoe
        ) from None
