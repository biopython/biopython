# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# Revisions 2023 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Code to parse and store BLAST XML output, and to invoke the NCBI BLAST web server.

This module provides code to parse and store BLAST XML output, following its
definition in the associated BLAST XML DTD file:
https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd

This module also provides code to invoke the BLAST web server provided by NCBI.
https://blast.ncbi.nlm.nih.gov/

Variables:

    - email        Set the Blast email parameter (default is None).
    - tool         Set the Blast tool parameter (default is ``biopython``).

"""


import warnings

import time

from urllib.parse import urlencode
from urllib.request import build_opener, install_opener
from urllib.request import urlopen
from urllib.request import HTTPPasswordMgrWithDefaultRealm, HTTPBasicAuthHandler
from urllib.request import Request

from Bio import BiopythonWarning
from Bio import StreamModeError
from Bio._utils import function_with_previous


email = None
tool = "biopython"


NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


class Record(list):
    """Stores the BLAST results for a single query."""

    def __init__(self):
        """Initialize the Record object."""
        self.query = None


class Records:
    """Stores the BLAST results of a single BLAST run."""

    def __init__(self, source):
        """Initialize the Records object."""
        from Bio.Blast._parser import XMLHandler

        self.source = source
        try:
            stream = open(source, "rb")
        except TypeError:  # not a path, assume we received a stream
            if source.read(0) != b"":
                raise StreamModeError(
                    "BLAST output files must be opened in binary mode."
                ) from None
            stream = source
        self._stream = stream
        handler = XMLHandler(stream)
        handler.read_header(self)
        self._handler = handler

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        try:
            stream = self._stream
        except AttributeError:
            return
        if stream is not self.source:
            stream.close()
        del self._stream

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._handler)


def parse(source):
    """Parse an XML file containing BLAST output and return a Bio.Blast.Records object.

    This returns an iterator object; iterating over it returns Bio.Blast.Record
    objects one by one.

    The source can be a file stream or the path to an XML file containing the
    BLAST output. If a file stream, source  must be in binary mode. This allows
    the parser to detect the encoding from the XML file,and to use it to convert
    any text in the XML to the correct Unicode string. The qblast function in
    Bio.Blast returns a file stream in binary mode. For files, please use mode
    "rb" when opening the file, as in

    >>> from Bio import Blast
    >>> stream = open("Blast/wnts.xml", "rb")  # opened in binary mode
    >>> records = Blast.parse(stream)
    >>> for record in records:
    ...     print(record.query.id, record.query.description)
    ...
    Query_1 gi|195230749:301-1383 Homo sapiens wingless-type MMTV integration site family member 2 (WNT2), transcript variant 1, mRNA
    Query_2 gi|325053704:108-1166 Homo sapiens wingless-type MMTV integration site family, member 3A (WNT3A), mRNA
    Query_3 gi|156630997:105-1160 Homo sapiens wingless-type MMTV integration site family, member 4 (WNT4), mRNA
    Query_4 gi|371502086:108-1205 Homo sapiens wingless-type MMTV integration site family, member 5A (WNT5A), transcript variant 2, mRNA
    Query_5 gi|53729353:216-1313 Homo sapiens wingless-type MMTV integration site family, member 6 (WNT6), mRNA
    >>> stream.close()

    """
    return Records(source)


def read(source):
    """Parse an XML file containing BLAST output for a single query and return it.

    Internally, this function uses Bio.Blast.parse to obtain an iterator over
    BLAST records.  The function then reads one record from the iterator,
    ensures that there are no further records, and returns the record it found
    as a Bio.Blast.Record object. An exception is raised if no records are
    found, or more than one record is found.

    The source can be a file stream or the path to an XML file containing the
    BLAST output. If a file stream, source  must be in binary mode. This allows
    the parser to detect the encoding from the XML file,and to use it to convert
    any text in the XML to the correct Unicode string. The qblast function in
    Bio.Blast returns a file stream in binary mode. For files, please use mode
    "rb" when opening the file, as in

    >>> from Bio import Blast
    >>> stream = open("Blast/xml_2900_blastn_001.xml", "rb")  # opened in binary mode
    >>> record = Blast.read(stream)
    >>> record.query.id
    'G26684.1'
    >>> record.query.description
    'human STS STS_D11570, sequence tagged site'
    >>> len(record)
    10
    >>> stream.close()

    Use the Bio.Blast.parse function if you want to read a file containing
    BLAST output for more than one query.
    """
    with parse(source) as records:
        try:
            record = next(records)
        except StopIteration:
            raise ValueError("No BLAST output found.") from None
        try:
            next(records)
            raise ValueError("BLAST output for more than one query found.")
        except StopIteration:
            pass
    return record


@function_with_previous
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
    username="blast",
    password=None,
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
     - format_type    "XML" (default), "HTML", "Text", "XML2", "JSON2",
                      or "Tabular".
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
            f"Program specified is {program}. Expected one of {', '.join(programs)}"
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
    parameters = {
        "AUTO_FORMAT": auto_format,
        "COMPOSITION_BASED_STATISTICS": composition_based_statistics,
        "DATABASE": database,
        "DB_GENETIC_CODE": db_genetic_code,
        "ENDPOINTS": endpoints,
        "ENTREZ_QUERY": entrez_query,
        "EXPECT": expect,
        "FILTER": filter,
        "GAPCOSTS": gapcosts,
        "GENETIC_CODE": genetic_code,
        "HITLIST_SIZE": hitlist_size,
        "I_THRESH": i_thresh,
        "LAYOUT": layout,
        "LCASE_MASK": lcase_mask,
        "MEGABLAST": megablast,
        "MATRIX_NAME": matrix_name,
        "NUCL_PENALTY": nucl_penalty,
        "NUCL_REWARD": nucl_reward,
        "OTHER_ADVANCED": other_advanced,
        "PERC_IDENT": perc_ident,
        "PHI_PATTERN": phi_pattern,
        "PROGRAM": program,
        # ('PSSM': pssm: - It is possible to use PSI-BLAST via this API?
        "QUERY": sequence,
        "QUERY_FILE": query_file,
        "QUERY_BELIEVE_DEFLINE": query_believe_defline,
        "QUERY_FROM": query_from,
        "QUERY_TO": query_to,
        # 'RESULTS_FILE': ...: - Can we use this parameter?
        "SEARCHSP_EFF": searchsp_eff,
        "SERVICE": service,
        "SHORT_QUERY_ADJUST": short_query,
        "TEMPLATE_TYPE": template_type,
        "TEMPLATE_LENGTH": template_length,
        "THRESHOLD": threshold,
        "UNGAPPED_ALIGNMENT": ungapped_alignment,
        "WORD_SIZE": word_size,
        "CMD": "Put",
    }

    if password is not None:
        # handle authentication for BLAST cloud
        password_mgr = HTTPPasswordMgrWithDefaultRealm()
        password_mgr.add_password(None, url_base, username, password)
        handler = HTTPBasicAuthHandler(password_mgr)
        opener = build_opener(handler)
        install_opener(opener)

    if url_base == NCBI_BLAST_URL:
        parameters.update({"email": email, "tool": tool})
    parameters = {key: value for key, value in parameters.items() if value is not None}
    message = urlencode(parameters).encode()
    request = Request(url_base, message, {"User-Agent": "BiopythonClient"})
    # Send off the initial query to qblast.
    # Note the NCBI do not currently impose a rate limit here, other
    # than the request not to make say 50 queries at once using multiple
    # threads.
    stream = urlopen(request)

    # Format the "Get" command, which gets the formatted results from qblast
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007
    rid, rtoe = _parse_qblast_ref_page(stream)
    parameters = {
        "ALIGNMENTS": alignments,
        "ALIGNMENT_VIEW": alignment_view,
        "DESCRIPTIONS": descriptions,
        "ENTREZ_LINKS_NEW_WINDOW": entrez_links_new_window,
        "EXPECT_LOW": expect_low,
        "EXPECT_HIGH": expect_high,
        "FORMAT_ENTREZ_QUERY": format_entrez_query,
        "FORMAT_OBJECT": format_object,
        "FORMAT_TYPE": format_type,
        "NCBI_GI": ncbi_gi,
        "RID": rid,
        "RESULTS_FILE": results_file,
        "SERVICE": service,
        "SHOW_OVERVIEW": show_overview,
        "CMD": "Get",
    }
    parameters = {key: value for key, value in parameters.items() if value is not None}
    message = urlencode(parameters).encode()

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
        wait = qblast.previous + delay - current
        if wait > 0:
            time.sleep(wait)
            qblast.previous = current + wait
        else:
            qblast.previous = current
        # delay by at least 60 seconds only if running the request against the public NCBI API
        if delay < 60 and url_base == NCBI_BLAST_URL:
            # Wasn't a quick return, must wait at least a minute
            delay = 60

        request = Request(url_base, message, {"User-Agent": "BiopythonClient"})
        stream = urlopen(request)
        data = stream.peek()
        if format_type == "HTML" and b"<title>NCBI Blast:</title>" in data:
            continue
        elif data.startswith(b"<!DOCTYPE html"):
            continue
        else:
            break
    if format_type == "XML":
        assert data.startswith(b"<?xml ")
    elif format_type == "HTML":
        assert data.startswith(b"<!DOCTYPE html ")
    elif format_type in ("Text", "Tabular"):
        assert data.startswith(b"<p><!--\nQBlastInfoBegin")
    elif format_type in ("XML2", "JSON2"):
        assert data.startswith(b"PK\x03\x04")  # zipped file
    return stream


qblast.previous = 0


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
                raise ValueError(f"Error message from NCBI: {msg}")
        # In spring 2010 the markup was like this:
        i = s.find('<p class="error">')
        if i != -1:
            msg = s[i + len('<p class="error">') :].strip()
            msg = msg.split("</p>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError(f"Error message from NCBI: {msg}")
        # Generic search based on the way the error messages start:
        i = s.find("Message ID#")
        if i != -1:
            # Break the message at the first HTML tag
            msg = s[i:].split("<", 1)[0].split("\n", 1)[0].strip()
            raise ValueError(f"Error message from NCBI: {msg}")
        # We didn't recognise the error layout :(
        # print(s)
        raise ValueError(
            "No RID and no RTOE found in the 'please wait' page, "
            "there was probably an error in your request but we "
            "could not extract a helpful error message."
        )
    elif not rid:
        # Can this happen?
        raise ValueError(
            f"No RID found in the 'please wait' page. (although RTOE = {rtoe!r})"
        )
    elif not rtoe:
        # Can this happen?
        raise ValueError(
            f"No RTOE found in the 'please wait' page. (although RID = {rid!r})"
        )

    try:
        return rid, int(rtoe)
    except ValueError:
        raise ValueError(
            f"A non-integer RTOE found in the 'please wait' page, {rtoe!r}"
        ) from None
