# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for fetching NCBI qblast.

Uses Bio.Blast.NCBIWWW.qblast() to run some online blast queries, get XML
blast results back, and then checks Bio.Blast.NCBIXML.parse() can read them.

Goals:
    Make sure that all retrieval is working as expected.
    Make sure we can parse the latest XML format being used by the NCBI.
"""
import sys
import requires_internet
requires_internet.check()
from Bio import MissingExternalDependencyError 
from urllib2 import HTTPError

#We want to test these:
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


#####################################################################

#List of qblast requests stored as a tuple of parameters:
# - program
# - database
# - query identifier or sequence
# - expectation value threshold
# - Entrez filter string (or None)
# - list of hit identifiers expected to be found (or None if expect 0)
tests = [ \
    #Simple protein blast filtered for rat only, using protein GI:160837788
    #the actin related protein 2/3 complex, subunit 1B [Mus musculus]
    ("blastp", "nr", "160837788", 0.001,
     "rat [ORGN]", ['9506405','13592137','37589612','149064087','56912225']),
    #This next example finds PCR primer matches in Chimpanzees, e.g. BRCA1:
    ("blastn", "nr", "GTACCTTGATTTCGTATTC"+("N"*30)+"GACTCTACTACCTTTACCC",
     10, "pan [ORGN]", ["37953274","51104367","51104367","51104367"]),
    #Try an orchid EST (nucleotide) sequence against NR using BLASTX
    ("blastx", "nr", """>gi|116660609|gb|EG558220.1|EG558220 CR02019H04 Leaf CR02 cDNA library Catharanthus roseus cDNA clone CR02019H04 5', mRNA sequence
CTCCATTCCCTCTCTATTTTCAGTCTAATCAAATTAGAGCTTAAAAGAATGAGATTTTTAACAAATAAAA
AAACATAGGGGAGATTTCATAAAAGTTATATTAGTGATTTGAAGAATATTTTAGTCTATTTTTTTTTTTT
TCTTTTTTTGATGAAGAAAGGGTATATAAAATCAAGAATCTGGGGTGTTTGTGTTGACTTGGGTCGGGTG
TGTATAATTCTTGATTTTTTCAGGTAGTTGAAAAGGTAGGGAGAAAAGTGGAGAAGCCTAAGCTGATATT
GAAATTCATATGGATGGAAAAGAACATTGGTTTAGGATTGGATCAAAAAATAGGTGGACATGGAACTGTA
CCACTACGTCCTTACTATTTTTGGCCGAGGAAAGATGCTTGGGAAGAACTTAAAACAGTTTTAGAAAGCA
AGCCATGGATTTCTCAGAAGAAAATGATTATACTTCTTAATCAGGCAACTGATATTATCAATTTATGGCA
GCAGAGTGGTGGCTCCTTGTCCCAGCAGCAGTAATTACTTTTTTTTCTCTTTTTGTTTCCAAATTAAGAA
ACATTAGTATCATATGGCTATTTGCTCAATTGCAGATTTCTTTCTTTTGTGAATG""",
     0.0000001, None, ["21554275","18409071","296087288"]),
]

print "Checking Bio.Blast.NCBIWWW.qblast() with various queries"
for program,database,query,e_value,entrez_filter,expected_hits in tests:
    print "qblast('%s', '%s', %s, ...)" % (program, database, repr(query))
    try:
        if program=="blastn":
            #Check the megablast parameter is accepted
            handle = NCBIWWW.qblast(program, database, query, \
                                    alignments=10, descriptions=10, \
                                    hitlist_size=10, \
                                    entrez_query=entrez_filter,
                                    expect=e_value, megablast="FALSE")
        else:
            handle = NCBIWWW.qblast(program, database, query, \
                                    alignments=10, descriptions=10, \
                                    hitlist_size=10, \
                                    entrez_query=entrez_filter,
                                    expect=e_value)
    except HTTPError:
        #e.g. a proxy error
        raise MissingExternalDependencyError("internet connection failed")
    record = NCBIXML.read(handle)

    if record.query == "No definition line":
        #We used a sequence as the query
        assert len(query) == record.query_letters
    elif query.startswith(">"):
        #We used a FASTA record as the query
        assert query[1:].split("\n",1)[0] == (record.query)
    else:
        #We used an identifier as the query
        assert query in record.query_id.split("|")

    #Check the recorded input parameters agree with those requested
    assert float(record.expect) == e_value
    assert record.application.lower() == program
    assert len(record.alignments) <= 10
    assert len(record.descriptions) <= 10

    #Check the expected result(s) are found in the alignments
    if expected_hits is None:
        assert len(record.alignments)==0, "Expected no alignments!"
    else:
        assert len(record.alignments) > 0, "Expected some alignments!"
        found_result = False
        for expected_hit in expected_hits:
            for alignment in record.alignments:
                if expected_hit in alignment.hit_id.split("|"):
                    found_result = True
                    break
        if len(expected_hits)==1:
            print "Update this test to have some redundancy..."
            for alignment in record.alignments:
                print alignment.hit_id
        assert found_result, "Missing all of %s in alignments" \
               % ", ".join(expected_hits)

    #Check the expected result(s) are found in the descriptions
    if expected_hits is None:
        assert len(record.descriptions)==0, "Expected no descriptions!"
    else:
        assert len(record.descriptions) > 0, "Expected some descriptions!"
        found_result = False
        for expected_hit in expected_hits:
            for descr in record.descriptions:
                if expected_hit == descr.accession \
                or expected_hit in descr.title.split(None,1)[0].split("|"):
                    found_result = True
                    break
        assert found_result, "Missing all of %s in descriptions" % expected_hit

print "Done"
