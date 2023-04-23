# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to access resources at ExPASy over the WWW.

See https://www.expasy.org/


Functions:
 - get_prodoc_entry  Interface to the get-prodoc-entry CGI script.
 - get_prosite_entry Interface to the get-prosite-entry CGI script.
 - get_prosite_raw   Interface to the get-prosite-raw CGI script.
 - get_sprot_raw     Interface to the get-sprot-raw CGI script.

"""

import io
from urllib.request import urlopen
from urllib.error import HTTPError


def get_prodoc_entry(
    id, cgi="https://prosite.expasy.org/cgi-bin/prosite/get-prodoc-entry"
):
    """Get a text handle to a PRODOC entry at ExPASy in HTML format.

    >>> from Bio import ExPASy
    >>> import os
    >>> with ExPASy.get_prodoc_entry('PDOC00001') as in_handle:
    ...     html = in_handle.read()
    ...
    >>> with open("myprodocrecord.html", "w") as out_handle:
    ...     length = out_handle.write(html)
    ...
    >>> os.remove("myprodocrecord.html")  # tidy up

    For a non-existing key XXX, ExPASy returns an HTML-formatted page
    containing this text: 'There is currently no PROSITE entry for'
    """
    return _open(f"{cgi}?{id}")


def get_prosite_entry(
    id, cgi="https://prosite.expasy.org/cgi-bin/prosite/get-prosite-entry"
):
    """Get a text handle to a PROSITE entry at ExPASy in HTML format.

    >>> from Bio import ExPASy
    >>> import os
    >>> with ExPASy.get_prosite_entry('PS00001') as in_handle:
    ...     html = in_handle.read()
    ...
    >>> with open("myprositerecord.html", "w") as out_handle:
    ...     length = out_handle.write(html)
    ...
    >>> os.remove("myprositerecord.html")  # tidy up

    For a non-existing key XXX, ExPASy returns an HTML-formatted page
    containing this text: 'There is currently no PROSITE entry for'
    """
    return _open(f"{cgi}?{id}")


def get_prosite_raw(id, cgi=None):
    """Get a text handle to a raw PROSITE or PRODOC record at ExPASy.

    The cgi argument is deprecated due to changes in the ExPASy
    website.

    >>> from Bio import ExPASy
    >>> from Bio.ExPASy import Prosite
    >>> with ExPASy.get_prosite_raw('PS00001') as handle:
    ...     record = Prosite.read(handle)
    ...
    >>> print(record.accession)
    PS00001

    This function raises a ValueError if the identifier does not exist:

    >>> handle = ExPASy.get_prosite_raw("DOES_NOT_EXIST")
    Traceback (most recent call last):
        ...
    ValueError: Failed to find entry 'DOES_NOT_EXIST' on ExPASy

    """
    handle = _open(f"https://prosite.expasy.org/{id}.txt")
    if handle.url == "https://www.expasy.org/":
        raise ValueError(f"Failed to find entry '{id}' on ExPASy") from None
    return handle


def get_sprot_raw(id):
    """Get a text handle to a raw SwissProt entry at ExPASy.

    For an ID of XXX, fetches http://www.uniprot.org/uniprot/XXX.txt
    (as per the https://www.expasy.org/expasy_urls.html documentation).

    >>> from Bio import ExPASy
    >>> from Bio import SwissProt
    >>> with ExPASy.get_sprot_raw("O23729") as handle:
    ...     record = SwissProt.read(handle)
    ...
    >>> print(record.entry_name)
    CHS3_BROFI

    This function raises a ValueError if the identifier does not exist:

    >>> ExPASy.get_sprot_raw("DOES_NOT_EXIST")
    Traceback (most recent call last):
        ...
    ValueError: Failed to find SwissProt entry 'DOES_NOT_EXIST'

    """
    try:
        handle = _open(f"http://www.uniprot.org/uniprot/{id}.txt")
    except HTTPError as exception:
        if exception.code in (400, 404):
            raise ValueError(f"Failed to find SwissProt entry '{id}'") from None
        else:
            raise
    return handle


def _open(url):
    """Open URL and convert to text assuming UTF-8 encoding (PRIVATE)."""
    handle = urlopen(url)
    text_handle = io.TextIOWrapper(handle, encoding="UTF-8")
    text_handle.url = handle.url
    return text_handle
