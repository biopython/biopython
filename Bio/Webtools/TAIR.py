# Copyright 2012 by Kevin Murray.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Functions to get Arabidopsis sequences by their AGI identifier.
"""
from Bio.Webtools.multiparthandler import multiparthandler
import urllib2
from StringIO import StringIO
from Bio import SeqIO
import re


def _sanitise_agis(agis):
    """Takes a list of agis, and returns a list of only those which are valid.
    """
    clean_agis = []
    agi_re = re.compile(r"AT[12345CM]G\d{5}(\.\d)?")
    for agi in agis:
        agi_match = agi_re.match(agi)
        if agi_match is not None:
            clean_agis.append(agi)
    return clean_agis


# functions to get sequences directly from arabidopsis.org
TAIR_DATASETS = {
    "transcript": "At_transcripts",
    "cds": "ATH1_cds",
    "gene": "ATH1_seq",
    "genomic": "ATH1_seq",
    "protein": "ATH1_pep",
    "peptide": "ATH1_pep",
    "coding_sequence": "ATH1_cds",
    "genomic_locus_sequence": "ATH1_seq",
    "upstream_500": "At_upstream_500",
    "upstream_1000": "At_upstream_1000",
    "upstream_3000": "At_upstream_3000",
    "downstream_500": "At_downstream_500",
    "downstream_1000": "At_downstream_1000",
    "downstream_3000": "At_downstream_3000",
    "intergenic": "At_intergenic",
    "intron": "At_intron",
    "3prime_utr": "ATH1_3_UTR",
    "5prime_utr": "ATH1_5_UTR"
    }

TAIR_TARGETS = {
    "representative": "rep_gene",
    "all": "both",
    "specified": "genemodel"
    }


def get(agis, dataset, target):
    """Get TAIR sequence(s) from AGI from the arabidopsis.org server directly.
    """
    bad_agi_exception = ValueError(
        "Must specify AGIs as a iterable, list or tuple"
        )
    agis = _sanitise_agis(agis)
    if agis is None:
        raise bad_agi_exception

    # Check dataset
    if dataset in TAIR_DATASETS:
        dataset = TAIR_DATASETS[dataset]
    elif dataset in TAIR_DATASETS.values():
        pass  # dataset is already equal to required value
    else:
        raise ValueError("%s is an invalid TAIR dataset" % dataset)

    # Check dataset
    if target in TAIR_TARGETS:
        target = TAIR_TARGETS[target]
    elif target in TAIR_TARGETS.values():
        pass  # target is already equal to required value
    else:
        raise ValueError("%s is an invalid TAIR target" % target)

    # Prepare http request params
    url = 'http://www.arabidopsis.org/cgi-bin/bulk/sequences/getseq.pl'
    params = {
        'loci': "\r\n".join(agis),
        'file': "",
        "dataset": dataset,
        "search_against": target,
        "outputformat": "outputformat"
        }

    # Prepare and get request
    opener = urllib2.build_opener(multiparthandler)
    raw_seq_string = opener.open(url, params).read()

    # Trim returned fasta file, responce has leading non-fasta text
    in_fasta_text = False
    seq_string = ""
    for line in raw_seq_string.splitlines():
        if len(line) < 1:
            continue
        if line[0] == ">":
            in_fasta_text = True
        if in_fasta_text:
            if line[0] == "-":
                break
            seq_string += line
            seq_string += "\n"
    return SeqIO.parse(StringIO(seq_string), "fasta")
