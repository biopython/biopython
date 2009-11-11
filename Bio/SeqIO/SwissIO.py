# Copyright 2006-2009 by Peter Cock and Michiel de Hoon.
# All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "swiss" (aka SwissProt/UniProt) file format.

You are expected to use this module via the Bio.SeqIO functions.
See also the Bio.SwissProt module which offers more than just accessing
the sequences as SeqRecord objects."""

from Bio import Seq
from Bio import SeqRecord
from Bio import Alphabet
from Bio import SeqFeature
from Bio import SwissProt
    
#This is a generator function!
def SwissIterator(handle):
    """Breaks up a Swiss-Prot/UniProt file into SeqRecord objects.

    Every section from the ID line to the terminating // becomes
    a single SeqRecord with associated annotation and features.

    This parser is for the flat file "swiss" format as used by:
     * Swiss-Prot aka SwissProt
     * TrEMBL
     * UniProtKB aka UniProt Knowledgebase

    It does NOT read their new XML file format.
    http://www.expasy.org/sprot/

    For consistency with BioPerl and EMBOSS we call this the "swiss"
    format.
    """
    swiss_records = SwissProt.parse(handle)
    for swiss_record in swiss_records:
        # Convert the SwissProt record to a SeqRecord
        seq = Seq.Seq(swiss_record.sequence, Alphabet.generic_protein)
        record = SeqRecord.SeqRecord(seq,
                                     id=swiss_record.accessions[0],
                                     name=swiss_record.entry_name,
                                     description=swiss_record.description,
                                    )
        record.description = swiss_record.description
        for cross_reference in swiss_record.cross_references:
            if len(cross_reference) < 2: continue
            database, accession = cross_reference[:2]
            dbxref = "%s:%s" % (database, accession)
            if not dbxref in record.dbxrefs:
                record.dbxrefs.append(dbxref)
        annotations = record.annotations
        annotations['accessions'] = swiss_record.accessions
        annotations['date'] = swiss_record.created[0]
        annotations['date_last_sequence_update'] = swiss_record.sequence_update[0]
        if swiss_record.annotation_update:
            annotations['date_last_annotation_update'] = swiss_record.annotation_update[0]
        if swiss_record.gene_name:
            annotations['gene_name'] = swiss_record.gene_name
        annotations['organism'] = swiss_record.organism.rstrip(".")
        annotations['taxonomy'] = swiss_record.organism_classification
        annotations['ncbi_taxid'] = swiss_record.taxonomy_id
        if swiss_record.host_organism:
            annotations['organism_host'] = [word.rstrip(".") for word in swiss_record.host_organism]
        if swiss_record.comments:
            annotations['comment'] = "\n".join(swiss_record.comments)
        if swiss_record.references:
            annotations['references'] = []
            for reference in swiss_record.references:
                feature = SeqFeature.Reference()
                feature.comment = " ".join(["%s=%s;" % (key, value) for key, value in reference.comments])
                for key, value in reference.references:
                    if key=='PubMed':
                        feature.pubmed_id = value
                    elif key=='MEDLINE':
                        feature.medline_id = value
                    elif key=='DOI':
                        pass
                    elif key=='AGRICOLA':
                        pass
                    else:
                        raise ValueError, "Unknown key %s found in references" % key
                feature.authors = reference.authors
                feature.title = reference.title
                feature.journal = reference.location
                annotations['references'].append(feature)
        if swiss_record.keywords:
            record.annotations['keywords'] = swiss_record.keywords
        yield record

if __name__ == "__main__":
    print "Quick self test..."

    example_filename = "../../Tests/SwissProt/sp008"

    import os
    if not os.path.isfile(example_filename):
        print "Missing test file %s" % example_filename
    else:
        #Try parsing it!
        handle = open(example_filename)
        records = SwissIterator(handle)
        for record in records:
            print record.name
            print record.id
            print record.annotations['keywords']
            print repr(record.annotations['organism'])
            print record.seq.tostring()[:20] + "..."
        handle.close()
