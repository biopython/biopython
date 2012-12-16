# Copyright 2012 by Kevin Murray.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Functions to get sequences by their Arabidopsis Genome Initiative (AGI)
identifier, which are commonly used by biologists working with Arabidopsis.
Sequences can be retreived either directly from The Arabidopsis Information
Resource (TAIR) at arabidopsis.org, or via GenBank/RefSeq, using
Bio.Entrez.efetch().
"""

from MultipartPostHandler import MultipartPostHandler
import urllib2
from StringIO import StringIO
from Bio import SeqIO, Entrez
from Bio.TAIR._ncbi import ncbi_prot, ncbi_rna
import re

NCBI_RNA = 1
NCBI_PROTEIN = 2


def _sanitise_agis(agis):
    clean_agis = []
    agi_re = re.compile(r"AT[12345CM]G\d{5}(\.\d)?")
    for agi in agis:
        agi_match = agi_re.match(agi)
        if agi_match is not None:
            clean_agis.append(agi)
    return clean_agis


class TAIRDirect(object):
    """
    Class of functions to get sequences directly from arabidopsis.org, by their
    AGI identifier.
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

    def _get_fasta_text(self, agis, dataset, target):
        bad_agi_exception = ValueError(
            "Must specify AGIs as a iterable, list or tuple"
            )
        agis = _sanitise_agis(agis)
        if agis is None:
            raise bad_agi_exception

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
                if line[0] == "-":
                    break
                seq_string += line
                seq_string += "\n"
        return seq_string

    def get(self, agis, dataset="gene", target="rep_gene"):
        return SeqIO.parse(
            StringIO(self._get_fasta_text(agis, dataset, target)),
            "fasta"
            )


class TAIRNCBI(object):
    """
    Class of functions to get sequences from NCBI using Efetch, by their AGI
    identifier.
    """
    def _agi_to_rna(self, agis):
        rna_ids = []
        for agi in agis:
            try:
                rna_ids.append(ncbi_rna[agi])
            except LookupError:
                next
        return rna_ids

    def _agi_to_protein(self, agis):
        protein_ids = []
        for agi in agis:
            try:
                protein_ids.append(ncbi_prot[agi])
            except LookupError:
                next
        return protein_ids

    def get_rna_from_ncbi(self, agis):
        Entrez.email = ""
        entrez_handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(self._agi_to_rna(agis)),
                rettype="gb",
                retmode="text"
                )
        return SeqIO.parse(entrez_handle, "gb")

    def get_protein_from_ncbi(self, agis):
        Entrez.email = ""
        entrez_handle = Entrez.efetch(
                db="protein",
                id=",".join(self._agi_to_protein(agis)),
                rettype="gb",
                retmode="text"
                )
        return SeqIO.parse(entrez_handle, "gb")


def get(agis, dataset="gene", target="rep_gene"):
    tair = TAIRDirect()
    return tair.get(_sanitise_agis(agis), dataset, target)


def get_from_ncbi(agis, mode):
    tair = TAIRNCBI()
    if mode == NCBI_RNA:
        return tair.get_rna_from_ncbi(agis)
    elif mode == NCBI_PROTEIN:
        return tair.get_protein_from_ncbi(agis)
