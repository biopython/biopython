# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Record.py

Classes:
Heavyweight
Lightweight

PSIBlast
MasterSlave

DescriptionHit
AlignmentHit
HSP

"""


class Heavyweight:
    """Saves the results from a blast search.

    Members:
    application
    version
    date
    reference

    query
    query_letters
    
    database
    database_sequences
    database_letters

    descriptions  A list of DescriptionHit objects.
    alignment     A list of AlignmentHit objects.
    
    """
    def __init__(self):
        self.application = None
        self.version = None
        self.reference = None

        self.database = None
        self.db_seqs = None
        self.db_letters = None

        self.descriptions = []
        self.alignments = []

class DescriptionHit:
    """Stores information about one hit in the descriptions section.

    Members:
    title
    score       High Score.
    p           P value.
    
    """
    def __init__(self):
        self.name = None
        self.score = None
        self.p = None

class AlignmentHit:
    """Stores information about one hit in the alignments section.

    Members:
    title      Name.
    length     Length.
    hsps       A list of HSP objects.

    """
    def __init__(self):
        self.name = None
        self.length = None
        self.hsps = []

class HSP:
    """Stores information about one hsp in an alignment hit.

    Members:
    score        BLAST score of hit
    bits         Number of bits for that score
    expect       Expect value
    identities   number of identities
    positives    number of positives
    gaps         number of gaps (gapped BLAST only)

    query
    query_start
    match
    sbjct
    sbjct_start
    
    Not all flavors of BLAST return values for every attribute:
              score     expect     identities   positives    strand  frame
    BLASTP     X          X            X            X
    BLASTN     X          X            X            X          X
    BLASTX     X          X            X            X                  X
    TBLASTN    X          X            X            X                  X
    TBLASTX    X          X            X            X                 X/X

    Note: for BLASTX, the query sequence is shown as a protein sequence,
    but the numbering is based on the nucleotides.  Thus, the numbering
    is 3x larger than the number of amino acid residues.  A similar effect
    can be seen for the sbjct sequence in TBLASTN, and for both sequences
    in TBLASTX.

    Also, for negative frames, the sequence numbering starts from
    query_start and counts down.

    """
    def __init__(self):
        self.score = None
        self.bits = None
        self.expect = None
        self.identities = None
        self.positives = None
        self.strand = None
        self.frame = None
        self.gaps = None
        
        self.query = None
        self.query_start = None
        self.match = None
        self.sbjct = None
        self.sbjct_start = None
