# Copyright 2006-2013,2020 by Peter Cock.
# Revisions copyright 2008-2009 by Michiel de Hoon.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "swiss" (aka SwissProt/UniProt) file format.

You are expected to use this module via the Bio.SeqIO functions.
See also the Bio.SwissProt module which offers more than just accessing
the sequences as SeqRecord objects.

See also Bio.SeqIO.UniprotIO.py which supports the "uniprot-xml" format.
"""
from Bio import SeqFeature
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _make_position(location_string, offset=0):
    """Turn a Swiss location position into a SeqFeature position object (PRIVATE).

    An offset of -1 is used with a start location to make it pythonic.
    """
    if location_string == "?":
        return SeqFeature.UnknownPosition()
    # Hack so that feature from 0 to 0 becomes 0 to 0, not -1 to 0.
    try:
        return SeqFeature.ExactPosition(max(0, offset + int(location_string)))
    except ValueError:
        pass
    if location_string.startswith("<"):
        try:
            return SeqFeature.BeforePosition(max(0, offset + int(location_string[1:])))
        except ValueError:
            pass
    elif location_string.startswith(">"):  # e.g. ">13"
        try:
            return SeqFeature.AfterPosition(max(0, offset + int(location_string[1:])))
        except ValueError:
            pass
    elif location_string.startswith("?"):  # e.g. "?22"
        try:
            return SeqFeature.UncertainPosition(
                max(0, offset + int(location_string[1:]))
            )
        except ValueError:
            pass
    raise NotImplementedError("Cannot parse location '%s'" % location_string)


def SwissIterator(source):
    """Break up a Swiss-Prot/UniProt file into SeqRecord objects.

    Argument source is a file-like object or a path to a file.

    Every section from the ID line to the terminating // becomes
    a single SeqRecord with associated annotation and features.

    This parser is for the flat file "swiss" format as used by:
     - Swiss-Prot aka SwissProt
     - TrEMBL
     - UniProtKB aka UniProt Knowledgebase

    For consistency with BioPerl and EMBOSS we call this the "swiss"
    format. See also the SeqIO support for "uniprot-xml" format.

    Rather than calling it directly, you are expected to use this
    parser via Bio.SeqIO.parse(..., format="swiss") instead.
    """
    swiss_records = SwissProt.parse(source)

    for swiss_record in swiss_records:
        # Convert the SwissProt record to a SeqRecord
        record = SeqRecord(
            Seq(swiss_record.sequence),
            id=swiss_record.accessions[0],
            name=swiss_record.entry_name,
            description=swiss_record.description,
            features=swiss_record.features,
        )
        for cross_reference in swiss_record.cross_references:
            if len(cross_reference) < 2:
                continue
            database, accession = cross_reference[:2]
            dbxref = "%s:%s" % (database, accession)
            if dbxref not in record.dbxrefs:
                record.dbxrefs.append(dbxref)
        annotations = record.annotations
        annotations["molecule_type"] = "protein"
        annotations["accessions"] = swiss_record.accessions
        if swiss_record.protein_existence:
            annotations["protein_existence"] = swiss_record.protein_existence
        if swiss_record.created:
            date, version = swiss_record.created
            annotations["date"] = date
            annotations["sequence_version"] = version
        if swiss_record.sequence_update:
            date, version = swiss_record.sequence_update
            annotations["date_last_sequence_update"] = date
            annotations["sequence_version"] = version
        if swiss_record.annotation_update:
            date, version = swiss_record.annotation_update
            annotations["date_last_annotation_update"] = date
            annotations["entry_version"] = version
        if swiss_record.gene_name:
            annotations["gene_name"] = swiss_record.gene_name
        annotations["organism"] = swiss_record.organism.rstrip(".")
        annotations["taxonomy"] = swiss_record.organism_classification
        annotations["ncbi_taxid"] = swiss_record.taxonomy_id
        if swiss_record.host_organism:
            annotations["organism_host"] = swiss_record.host_organism
        if swiss_record.host_taxonomy_id:
            annotations["host_ncbi_taxid"] = swiss_record.host_taxonomy_id
        if swiss_record.comments:
            annotations["comment"] = "\n".join(swiss_record.comments)
        if swiss_record.references:
            annotations["references"] = []
            for reference in swiss_record.references:
                feature = SeqFeature.Reference()
                feature.comment = " ".join("%s=%s;" % k_v for k_v in reference.comments)
                for key, value in reference.references:
                    if key == "PubMed":
                        feature.pubmed_id = value
                    elif key == "MEDLINE":
                        feature.medline_id = value
                    elif key == "DOI":
                        pass
                    elif key == "AGRICOLA":
                        pass
                    else:
                        raise ValueError("Unknown key %s found in references" % key)
                feature.authors = reference.authors
                feature.title = reference.title
                feature.journal = reference.location
                annotations["references"].append(feature)
        if swiss_record.keywords:
            record.annotations["keywords"] = swiss_record.keywords
        yield record
