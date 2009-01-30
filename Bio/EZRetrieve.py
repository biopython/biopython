"""This module contains code to access EZRetrieve (OBSOLETE).

This is a very simple interface to the EZRetrieve website described in:

Zhang, H., Ramanathan, Y., Soteropoulos, P., Recce, M., and Tolias, P.P. (2002).
EZ-Retrieve: A web-server for batch retrieval of coordinate-specified human
DNA sequences and underscoring putative transcription factor-binding sites.
Nucl. Acids. Res. 2002 30: e121.
http://dx.doi.org/10.1093/nar/gnf120

Functions:
retrieve_single  Retrieve a single sequence from EZRetrieve.
parse_single     Parse the results from EZRetrieve into FASTA format.

This module is now considered to be obsolete, and is likely to be deprecated
in a future release of Biopython, and later removed.
"""

def retrieve_single(id, from_, to, retrieve_by=None, organism=None,
                    parse_results=1):
    import urllib
    
    CGI = "http://siriusb.umdnj.edu:18080/EZRetrieve/single_r_run.jsp"
    org2value = {"Hs" : "0", "Mm" : "1", "Rn" : 2}
    organism = organism or "Hs"
    assert organism in org2value

    acctype2value = {"genbank":0, "unigene":1, "locuslink":2, "image":3}
    retrieve_by = retrieve_by or "GenBank"
    retrieve_by = retrieve_by.lower()
    assert retrieve_by in acctype2value

    params = {
        "input" : str(id),
        "from" : str(from_),
        "to" : str(to),
        "org" : org2value[organism],
        "AccType" : acctype2value[retrieve_by],
        }
    options = urllib.urlencode(params)
    handle = urllib.urlopen(CGI, options)
    if parse_results:
        results = parse_single(handle)
    else:
        results = handle.read()
    return results

def parse_single(handle):
    """Return a FASTA-formatted string for the sequence.  May raise an
    AssertionError if there was a problem retrieving the sequence.

    """
    import re
    results = handle.read()
    lresults = results.lower()
    
    i = results.find("Error: ")
    if i >= 0:
        j = lresults.index("<br>", i)
        errmsg = results[i:j].strip()
        raise AssertionError(errmsg)

    i = lresults.find("<b>>")
    assert i >= 0, "Couldn't find sequence."
    j = lresults.find("<br><br>", i)
    seqdata = results[i:j]
    reobj = re.compile(r"<[^>]*>", re.IGNORECASE|re.DOTALL)
    seqdata = reobj.sub("", seqdata)
    seqdata = re.sub(r"\s+", r"\n", seqdata)
    seqdata = seqdata.strip() + "\n"
    return seqdata
