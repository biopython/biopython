# Copyright 2005 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from Bio.Blast import NCBIXML


all_tests = [
    'xbt001.xml', # BLASTP 2.2.12, gi|49176427|ref|NP_418280.3|
    'xbt002.xml', # BLASTN 2.2.12, gi|1348916|gb|G26684.1|G26684
    'xbt003.xml', # BLASTX 2.2.12, gi|1347369|gb|G25137.1|G25137
    'xbt004.xml', # TBLASTN 2.2.12, gi|729325|sp|P39483|DHG2_BACME
    'xbt005.xml', # TBLASTX 2.2.12, gi|1348853|gb|G26621.1|G26621, BLOSUM80
    ]

detailed_tests = [
    'xbt001.xml', # BLASTP 2.2.12, gi|49176427|ref|NP_418280.3|
    'xbt002.xml', # BLASTN 2.2.12, gi|1348916|gb|G26684.1|G26684
    'xbt003.xml', # BLASTX 2.2.12, gi|1347369|gb|G25137.1|G25137
    ]

### NCBIXML.BlastParser

print "Running tests on NCBIXML.BlastParser"

parser = NCBIXML.BlastParser()
for test in all_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    input = open(datafile)

    result = parser.parse(input)

    # In 2006, in order to cope with the new style XML files
    # for multiple queries introduced in Blast 2.2.14 the
    # parser was changed to return a list or records
    # (instead of one single record).
    if isinstance(result, list) :
        assert len(result)==1
        result = result[0]

    if not test in detailed_tests: continue

    alignments = result.alignments
    print 'Blast of', result.query
    print 'Found %s alignments with a total of %s HSPs' % (len(alignments),
              reduce(lambda a,b: a+b, [len(a.hsps) for a in alignments]))
    for alignment in alignments:
        print alignment.title[:50], alignment.length, 'bp', len(alignment.hsps), 'HSPs'

    # Cookbook example
    E_VALUE_THRESH = 0.04
    for alignment in alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print '*****'
                print 'sequence', alignment.title
                print 'length', alignment.length
                # The following is a quick and dirty hack to get the
                # number of digits in the exponent to match that on record
                # as the expected output.
                f = str(hsp.expect)
                if f.find("e-") <> -1 :
                    matissa, exponent = str(f).split("e-")
                    print 'e value %se-%02i' % (matissa, int(exponent))
                else :
                    print 'e value', hsp.expect
                print hsp.query[:75] + '...'
                print hsp.match[:75] + '...'
                print hsp.sbjct[:75] + '...'
