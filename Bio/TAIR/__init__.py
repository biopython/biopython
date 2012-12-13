# Copyright 2012 by Kevin Murray.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from MultipartPostHandler import MultipartPostHandler
import urllib2
from StringIO import StringIO
from Bio import SeqIO


class TAIRDirect:
    """
    """
    def __init__(self):
        self.datasets = {
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
        self.targets = {
            "representative": "rep_gene",
            "all": "both",
            "specified": "genemodel"
            }

    def _getSeqFastaText(self, agis, dataset, target):
        if agis is None:
            raise ValueError(
                "Must specify AGIs, specify agis='Demo' for a demo run"
                )

        # Check dataset
        if dataset in self.datasets.keys():
            dataset = self.datasets[dataset]
        elif dataset in self.datasets.values():
            pass  # dataset is already equal to required value
        else:
            raise ValueError("%s is an invalid TAIR dataset" % dataset)

        # Check dataset
        if target in self.targets.keys():
            target = self.targets[target]
        elif target in self.targets.values():
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
        opener = urllib2.build_opener(MultipartPostHandler)
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
                seq_string += line
                seq_string += "\n"
        return seq_string

    def get(self, agis, dataset="gene", target="rep_gene"):
        return SeqIO.parse(
            StringIO(self._getSeqFastaText(agis, dataset, target)),
            "fasta"
            )


def get(agis, dataset="gene", target="rep_gene"):
    tair = TAIRDirect()
    return tair.get(agis, dataset, target)
