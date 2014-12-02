# Copyright 2006-2013 by Peter Cock.
# Revisions copyright 2008-2009 by Michiel de Hoon.
# All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "swiss" (aka SwissProt/UniProt) file format.

You are expected to use this module via the Bio.SeqIO functions.
See also the Bio.SwissProt module which offers more than just accessing
the sequences as SeqRecord objects.

See also Bio.SeqIO.UniprotIO.py which supports the "uniprot-xml" format.
"""

from __future__ import print_function

from Bio import Seq
from Bio import SeqRecord
from Bio import Alphabet
from Bio import SeqFeature
from Bio import SwissProt

__docformat__ = "restructuredtext en"


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
            return SeqFeature.UncertainPosition(max(0, offset + int(location_string[1:])))
        except ValueError:
            pass
    raise NotImplementedError("Cannot parse location '%s'" % location_string)


def _make_seqfeature(name, from_res, to_res, description, ft_id):
    """Construct SeqFeature from feature data from parser (PRIVATE)."""
    loc = SeqFeature.FeatureLocation(_make_position(from_res, -1),
                                     _make_position(to_res, 0))
    if not ft_id:
        ft_id = "<unknown id>"  # The default in SeqFeature object
    return SeqFeature.SeqFeature(loc, type=name, id=ft_id,
                                 qualifiers={"description": description})


def SwissIterator(handle):
    """Breaks up a Swiss-Prot/UniProt file into SeqRecord objects.

    Every section from the ID line to the terminating // becomes
    a single SeqRecord with associated annotation and features.

    This parser is for the flat file "swiss" format as used by:
     - Swiss-Prot aka SwissProt
     - TrEMBL
     - UniProtKB aka UniProt Knowledgebase

    For consistency with BioPerl and EMBOSS we call this the "swiss"
    format. See also the SeqIO support for "uniprot-xml" format.
    """
    swiss_records = SwissProt.parse(handle)
    for swiss_record in swiss_records:
        # Convert the SwissProt record to a SeqRecord
        seq = Seq.Seq(swiss_record.sequence, Alphabet.generic_protein)
        record = SeqRecord.SeqRecord(seq,
                                     id=swiss_record.accessions[0],
                                     name=swiss_record.entry_name,
                                     description=swiss_record.description,
                                     features=[_make_seqfeature(*f) for f
                                               in swiss_record.features],
                                     )
        record.description = swiss_record.description
        for cross_reference in swiss_record.cross_references:
            if len(cross_reference) < 2:
                continue
            database, accession = cross_reference[:2]
            dbxref = "%s:%s" % (database, accession)
            if dbxref not in record.dbxrefs:
                record.dbxrefs.append(dbxref)
        annotations = record.annotations
        annotations['accessions'] = swiss_record.accessions
        if swiss_record.created:
            annotations['date'] = swiss_record.created[0]
        if swiss_record.sequence_update:
            annotations[
                'date_last_sequence_update'] = swiss_record.sequence_update[0]
        if swiss_record.annotation_update:
            annotations['date_last_annotation_update'] = swiss_record.annotation_update[0]
        if swiss_record.gene_name:
            annotations['gene_name'] = swiss_record.gene_name
        annotations['organism'] = swiss_record.organism.rstrip(".")
        annotations['taxonomy'] = swiss_record.organism_classification
        annotations['ncbi_taxid'] = swiss_record.taxonomy_id
        if swiss_record.host_organism:
            annotations['organism_host'] = swiss_record.host_organism
        if swiss_record.host_taxonomy_id:
            annotations['host_ncbi_taxid'] = swiss_record.host_taxonomy_id
        if swiss_record.comments:
            annotations['comment'] = "\n".join(swiss_record.comments)
        if swiss_record.references:
            annotations['references'] = []
            for reference in swiss_record.references:
                feature = SeqFeature.Reference()
                feature.comment = " ".join("%s=%s;" % k_v for k_v in reference.comments)
                for key, value in reference.references:
                    if key == 'PubMed':
                        feature.pubmed_id = value
                    elif key == 'MEDLINE':
                        feature.medline_id = value
                    elif key == 'DOI':
                        pass
                    elif key == 'AGRICOLA':
                        pass
                    else:
                        raise ValueError(
                            "Unknown key %s found in references" % key)
                feature.authors = reference.authors
                feature.title = reference.title
                feature.journal = reference.location
                annotations['references'].append(feature)
        if swiss_record.keywords:
            record.annotations['keywords'] = swiss_record.keywords
        yield record

if __name__ == "__main__":
    print("Quick self test...")

    example_filename = "../../Tests/SwissProt/sp008"

    import os
    if not os.path.isfile(example_filename):
        print("Missing test file %s" % example_filename)
    else:
        # Try parsing it!
        with open(example_filename) as handle:
            records = SwissIterator(handle)
            for record in records:
                print(record.name)
                print(record.id)
                print(record.annotations['keywords'])
                print(repr(record.annotations['organism']))
                print(str(record.seq)[:20] + "...")
                for f in record.features:
                    print(f)
