# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Consumers for AlignACE and CompareACE parsers.
"""

class AlignAceScanner:
    """Scannner for AlignACE output

    Methods:
    feed     Feed data into the scanner.

    The scanner generates (and calls the consumer) the following types of events:

    noevent - blank line

    version - AlignACE version number
    command_line - AlignACE command line string
    parameters - the begining of the parameters
    parameter - the line containing a parameter
    sequences - the begining of the sequences list
    sequence - line containing the name of the input sequence (and a respective number)
    motif - the begining of the motif (contains the number)
    motif_hit - one hit for a motif
    motif_mask - mask of the motif (space - gap, asterisk - significant position)
    motif_score - MAP score of the motif - approx. N * log R, where R == (num. of actual occur.) / (num. of occur. expected by random.)
    
    """
    def feed(self, handle, consumer):
        """S.feed(handle, consumer)

        Feed in a AlignACE report for scanning.  handle is a file-like
        object that contains the AlignACE report.  consumer is a Consumer
        object that will receive events as the report is scanned.
        """
        consumer.version(handle.readline())
        consumer.command_line(handle.readline())
        for line in handle:
            if line.strip() == "":
                consumer.noevent(line)
            elif line[:4]=="Para":
                consumer.parameters(line)
            elif line[0]=="#":
                consumer.sequence(line)
            elif "=" in line:
                consumer.parameter(line)
            elif line[:5]=="Input":
                consumer.sequences(line)
            elif line[:5]=="Motif":
                consumer.motif(line)
            elif line[:3]=="MAP":
                consumer.motif_score(line)
            elif len(line.split("\t"))==4:
                consumer.motif_hit(line)
            elif "*" in line:
                consumer.motif_mask(line)
            else:
                raise ValueError(line)

class CompareAceScanner:
    """Scannner for CompareACE output

    Methods:
    feed     Feed data into the scanner.

    The scanner generates (and calls the consumer) the following types of events:

    motif_score - CompareACE score of motifs

    ###### TO DO #############3
    extend the scanner to include other, more complex outputs.
    """
    def feed(self, handle, consumer):
        """S.feed(handle, consumer)

        Feed in a CompareACE report for scanning.  handle is a file-like
        object that contains the CompareACE report.  consumer is a Consumer
        object that will receive events as the report is scanned.
        """
        consumer.motif_score(handle.readline())
