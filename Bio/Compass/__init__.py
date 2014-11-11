# Copyright 2004 by James Casbon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to deal with COMPASS output, a program for profile/profile comparison.

Compass is described in:

Sadreyev R, Grishin N. COMPASS: a tool for comparison of multiple protein
alignments with assessment of statistical significance. J Mol Biol. 2003 Feb
7;326(1):317-36.

Tested with COMPASS 1.24.
"""

import re

__docformat__ = "restructuredtext en"


def read(handle):
    """Reads a COMPASS file containing one COMPASS record."""
    record = None
    try:
        line = next(handle)
        record = Record()
        __read_names(record, line)
        line = next(handle)
        __read_threshold(record, line)
        line = next(handle)
        __read_lengths(record, line)
        line = next(handle)
        __read_profilewidth(record, line)
        line = next(handle)
        __read_scores(record, line)
    except StopIteration:
        if not record:
            raise ValueError("No record found in handle")
        else:
            raise ValueError("Unexpected end of stream.")
    for line in handle:
        if not line.strip():  # skip empty lines
            continue
        __read_query_alignment(record, line)
        try:
            line = next(handle)
            __read_positive_alignment(record, line)
            line = next(handle)
            __read_hit_alignment(record, line)
        except StopIteration:
            raise ValueError("Unexpected end of stream.")
    return record


def parse(handle):
    """Iterates over records in a COMPASS file."""
    record = None
    try:
        line = next(handle)
    except StopIteration:
        return
    while True:
        try:
            record = Record()
            __read_names(record, line)
            line = next(handle)
            __read_threshold(record, line)
            line = next(handle)
            __read_lengths(record, line)
            line = next(handle)
            __read_profilewidth(record, line)
            line = next(handle)
            __read_scores(record, line)
        except StopIteration:
            raise ValueError("Unexpected end of stream.")
        for line in handle:
            if not line.strip():
                continue
            if "Ali1:" in line:
                yield record
                break
            __read_query_alignment(record, line)
            try:
                line = next(handle)
                __read_positive_alignment(record, line)
                line = next(handle)
                __read_hit_alignment(record, line)
            except StopIteration:
                raise ValueError("Unexpected end of stream.")
        else:
            yield record
            break


class Record(object):
    """Hold information from one compass hit.

    Ali1 is the query, Ali2 the hit.
    """

    def __init__(self):
        self.query = ''
        self.hit = ''
        self.gap_threshold = 0
        self.query_length = 0
        self.query_filtered_length = 0
        self.query_nseqs = 0
        self.query_neffseqs = 0
        self.hit_length = 0
        self.hit_filtered_length = 0
        self.hit_nseqs = 0
        self.hit_neffseqs = 0
        self.sw_score = 0
        self.evalue = -1
        self.query_start = -1
        self.hit_start = -1
        self.query_aln = ''
        self.hit_aln = ''
        self.positives = ''

    def query_coverage(self):
        """Return the length of the query covered in the alignment."""
        s = self.query_aln.replace("=", "")
        return len(s)

    def hit_coverage(self):
        """Return the length of the hit covered in the alignment."""
        s = self.hit_aln.replace("=", "")
        return len(s)

# Everything below is private

__regex = {"names": re.compile("Ali1:\s+(\S+)\s+Ali2:\s+(\S+)\s+"),
           "threshold": re.compile("Threshold of effective gap content in columns: (\S+)"),
           "lengths": re.compile("length1=(\S+)\s+filtered_length1=(\S+)\s+length2=(\S+)\s+filtered_length2=(\S+)"),
           "profilewidth": re.compile("Nseqs1=(\S+)\s+Neff1=(\S+)\s+Nseqs2=(\S+)\s+Neff2=(\S+)"),
           "scores": re.compile("Smith-Waterman score = (\S+)\s+Evalue = (\S+)"),
           "start": re.compile("(\d+)"),
           "align": re.compile("^.{15}(\S+)"),
           "positive_alignment": re.compile("^.{15}(.+)"),
          }


def __read_names(record, line):
    # Ali1: 60456.blo.gz.aln  Ali2: allscop//14984.blo.gz.aln
    #       ------query-----        -------hit-------------
    if "Ali1:" not in line:
        raise ValueError("Line does not contain 'Ali1:':\n%s" % line)
    m = __regex["names"].search(line)
    record.query = m.group(1)
    record.hit = m.group(2)


def __read_threshold(record, line):
    if not line.startswith("Threshold"):
        raise ValueError("Line does not start with 'Threshold':\n%s" % line)
    m = __regex["threshold"].search(line)
    record.gap_threshold = float(m.group(1))


def __read_lengths(record, line):
    if not line.startswith("length1="):
        raise ValueError("Line does not start with 'length1=':\n%s" % line)
    m = __regex["lengths"].search(line)
    record.query_length = int(m.group(1))
    record.query_filtered_length = float(m.group(2))
    record.hit_length = int(m.group(3))
    record.hit_filtered_length = float(m.group(4))


def __read_profilewidth(record, line):
    if "Nseqs1" not in line:
        raise ValueError("Line does not contain 'Nseqs1':\n%s" % line)
    m = __regex["profilewidth"].search(line)
    record.query_nseqs = int(m.group(1))
    record.query_neffseqs = float(m.group(2))
    record.hit_nseqs = int(m.group(3))
    record.hit_neffseqs = float(m.group(4))


def __read_scores(record, line):
    if not line.startswith("Smith-Waterman"):
        raise ValueError("Line does not start with 'Smith-Waterman':\n%s" % line)
    m = __regex["scores"].search(line)
    if m:
        record.sw_score = int(m.group(1))
        record.evalue = float(m.group(2))
    else:
        record.sw_score = 0
        record.evalue = -1.0


def __read_query_alignment(record, line):
    m = __regex["start"].search(line)
    if m:
        record.query_start = int(m.group(1))
    m = __regex["align"].match(line)
    assert m is not None, "invalid match"
    record.query_aln += m.group(1)


def __read_positive_alignment(record, line):
    m = __regex["positive_alignment"].match(line)
    assert m is not None, "invalid match"
    record.positives += m.group(1)


def __read_hit_alignment(record, line):
    m = __regex["start"].search(line)
    if m:
        record.hit_start = int(m.group(1))
    m = __regex["align"].match(line)
    assert m is not None, "invalid match"
    record.hit_aln += m.group(1)
