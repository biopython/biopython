"""Code to access CATHDB API.

See http://www.cathdb.info


Functions:
    - by_funfhmmer : interface to search by sequence
    - retrieve_results : interface to retrieve results
    - check_progress : interface to check search task progress
"""

import json
# Importing these functions with leading underscore as not intended for reuse
from Bio._py3k import _as_bytes
from Bio._py3k import Request as _request
from Bio._py3k import urlopen as _urlopen
from Bio._py3k import urlencode as _urlencode


def _make_json_request(url, data=None):
    """Make an application/json accepting request"""
    if data:
        req = _request(url, _as_bytes(data), {'Accept': 'application/json'})
    else:
        req = _request(url, headers={'Accept': 'application/json'})
    return req


def search_by_sequence(fasta, url='http://www.cathdb.info/search/by_funfhmmer'):
    """Search CATH by sequence fasta"""
    variables = {'fasta': fasta}
    data = _urlencode(variables)
    req = _make_json_request(url, data)
    response = _urlopen(req)
    try:
        json_response = json.load(response)
    except AttributeError as err:
        raise ValueError('Error receiving json response from {}. Received response: {}'.format(url, response.read()))
    return json_response['task_id']


def check_progress(task_id, url='http://www.cathdb.info/search/by_funfhmmer/check'):
    """Check the progress of search job."""
    url = '{}/{}'.format(url, task_id)
    response = _urlopen(_make_json_request(url))
    try:
        json_response = json.load(response)
    except AttributeError as err:
        raise ValueError('Error receiving json response from {}. Received response: {}'.format(url, response.read()))
    return json_response


def retrieve_results(task_id, url='http://www.cathdb.info/search/by_funfhmmer/results'):
    """Retrieve the results of search job."""
    url = '{}/{}'.format(url, task_id)
    response = _urlopen(_make_json_request(url))
    try:
        json_response = json.load(response)
    except AttributeError as err:
        raise ValueError('Error receiving json response from {}. Received response: {}'.format(url, response.read()))
    return json_response
