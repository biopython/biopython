# Copyright 2013 by Iddo Friedberg idoerg@gmail.com
# Revision copyright 2013 by Peter Cock.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for dealing with assorted UniProt file formats and interacting with the UniProt database.

This currently include parsers for the GAF, GPA and GPI formats
from UniProt-GOA as the module Bio.UniProt.GOA.

See also Bio.SwissProt and the "swiss" support in Bio.SeqIO for
the legacy plain text sequence format still used in UniProt.

See also Bio.SeqIO.SwissIO for the "uniprot-xml" support in
Bio.SeqIO.
"""

import json
import re
import urllib.parse
from http.client import HTTPResponse
from typing import Optional
from urllib.request import urlopen

_re_next_link = re.compile(r'<(.+)>; rel="next"')


def _get_next_link(response: HTTPResponse) -> Optional[str]:
    headers = response.headers

    if "Link" in headers and headers["Link"]:
        match = _re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

    return None


def _get_results(response: HTTPResponse) -> list[dict]:
    return json.loads(response.read().decode())["results"]


def _get_search_result_count(response: HTTPResponse) -> int:
    headers = response.headers

    if "x-total-results" in headers and headers["x-total-results"]:
        return int(headers["x-total-results"])
    else:
        return 0


class _UniProtSearchResults:
    """A sequence over the results of a UniProt search.

    Do not use this class directly. Instead, use the :meth:`UniProt.search` method.
    """

    def _fetch_next_batch(self) -> HTTPResponse:
        assert self.next_url is not None  # type: ignore
        with urlopen(self.next_url) as response:  # type: ignore
            self.results_cache += _get_results(response)
            self.next_url = _get_next_link(response)
            return response

    def __init__(self, first_url: str):
        self.next_url = first_url
        self.results_cache: list[dict] = []
        self.next_result_index = 0
        response = self._fetch_next_batch()
        self.search_result_count = _get_search_result_count(response)

    def __iter__(self):
        return self

    def __len__(self) -> int:
        """Return the total number of search results, regardless of the batch size."""
        return self.search_result_count

    def _fetch_for(self, index: int) -> None:
        """Fetch batches until the given index is in the cache."""
        assert index in range(len(self))
        while index >= len(self.results_cache):
            self._fetch_next_batch()

    def __next__(self) -> dict:
        if self.next_result_index < len(self):
            self._fetch_for(self.next_result_index)
            try:
                next_result = self.results_cache[self.next_result_index]
                self.next_result_index += 1
                return next_result
            except IndexError:
                raise StopIteration
        else:
            raise StopIteration

    def __getitem__(self, index):
        if isinstance(index, slice):
            start, stop, step = index.indices(len(self))
            # The assertions below should be guaranteed by the indices method
            assert 0 <= start < len(self) and 0 <= stop <= len(self)
            if step > 0:
                if start <= stop and stop > 0:
                    self._fetch_for(stop - 1)
            else:
                # If start is 0 and step is negative, the slice is empty
                if start >= stop and start > 0:
                    self._fetch_for(start)
        elif isinstance(index, int):
            if index not in range(-len(self), len(self)):
                raise IndexError("Index out of bounds.")
            self._fetch_for(index % len(self))
        return self.results_cache[index]


def search(
    query: str, fields: Optional[list[str]] = None, batch_size: int = 500
) -> _UniProtSearchResults:
    """Search the UniProt database.

    Consider using `query syntax <https://www.uniprot.org/help/text-search>`_ and
    `query fields <https://www.uniprot.org/help/query-fields>`_ to refine your search.

    See the API details `here <https://www.uniprot.org/help/api_queries>`_.

    >>> from Bio import UniProt
    >>> from itertools import islice
    >>> # Get the first 10 results
    >>> results = UniProt.search("(organism_id:2697049) AND (reviewed:true)")[:10]

    :param query: The query string to search UniProt with
    :type query: str
    :param fields: The columns to retrieve in the results, defaults to all fields
    :type fields: List[str], optional
    :param batch_size: The number of results to retrieve in each batch, defaults to 500
    :type batch_size: int
    :return: An iterator over the search results
    :rtype: _UniProtSearchResults
    """
    parameters = {
        "query": query,
        "size": batch_size,
        "format": "json",
    }
    if fields:
        parameters["fields"] = ",".join(fields)
    url = f"https://rest.uniprot.org/uniprotkb/search?{urllib.parse.urlencode(parameters)}"

    return _UniProtSearchResults(url)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
