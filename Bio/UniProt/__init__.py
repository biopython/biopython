# Copyright 2013 by Iddo Friedberg idoerg@gmail.com
# Revision copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with assorted UniProt file formats.

This currently include parsers for the GAF, GPA and GPI formats
from UniProt-GOA as the module Bio.UniProt.GOA.

See also Bio.SwissProt and the "swiss" support in Bio.SeqIO for
the legacy plain text sequence format still used in UniProt.
"""

from Bio._py3k import urlopen as _urlopen
from Bio._py3k import urlencode as _urlencode
from Bio._py3k import _binary_to_string_handle

def search(text, format="tab", sort="score", oragnism="", columns=(),
           isoform=False, compress=False, offset=0, limit=0):
    """
    Performs a query over the UniProt API.
    More at: https://www.uniprot.org/help/api_queries

    :param text: Text to be queried.
    :param format: The format to be retrieved. Possible:
                    - tab
                    - html
                    - xls
                    - fasta
                    - gff
                    - txt
                    - xml
                    - rdf
                    - list
                    - rss
    :param sort: How the query results should be sorted. Default is by the score
    of the search.
    :param oragnism: Specifies the organism for which the query should be done.
    :param columns: Which columns should be retrieved in a tuple format.
    Possible:
        citation, clusters, comments, domains, domain, ec, id, entry name,
        existence, families, features, genes, go, go-id, interactor, keywords,
        last-modified, length, organism, organism-id, pathway, protein names,
        reviewed, 3d, version, virus hosts.
    :param isoform: Include isoform sequences (used together with fasta output
    format).
    :param compress: Returns results gzipped.
    :param offset: Offset of the first result taken.
    :param limit: The amount of results taken.
    :return: Result of the query as TextIOWrapper.
    """
    cgi = "https://www.uniprot.org/uniprot/?"
    variables = {"query": text,
                 "format": format,
                 "sort": sort,
                 "offset": str(offset)}
    if oragnism:
        variables["organism"] = oragnism
    if columns:
        variables["columns"] = ",".join(columns)
    if isoform:
        variables["isoform"] = "Yes"
    if compress:
        variables["compress"] = "Yes"
    if limit:
        variables["limit"] = str(limit)

    fullcgi = "".join((cgi, _urlencode(variables)))
    return _binary_to_string_handle(_urlopen(fullcgi))
