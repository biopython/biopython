# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Record.py

Classes:
Comprehensive   Holds all the information from a BLAST record.
Lightweight     Not implemented.
PSIBlast        Not implemented.
MasterSlave     Not implemented.

Description     Holds information about a hit description.
Alignment       Holds information about an alignment hit.
HSP             Holds information about 1 HSP.

"""

class Comprehensive:
    """Saves the results from a blast search.

    Members:
    application         The name of the BLAST flavor that generated this data.
    version             Version of blast used.
    date                Date this data was generated.
    reference           Reference for blast.

    query               Name of query sequence.
    query_letters       Number of letters in the query sequence.  (int)
    
    database            Name of the database.
    database_sequences  Number of sequences in the database.  (int)
    database_letters    Number of letters in the database.  (int)

    descriptions        A list of Description objects.
    alignments          A list of Alignment objects.

    posted_date         The date the database was posted.
    ka_params           A tuple of (lambda, k, h) values.  (floats)
    gapped         # XXX this isn't set right!
    ka_params_gap       A tuple of (lambda, k, h) values.  (floats)
    
    matrix              Name of the matrix.
    gap_penalties       Tuple of (open, extend) penalties.  (floats)
    num_hits            Number of hits to the database.  (int)
    num_sequences       Number of sequences.  (int)
    num_good_extends    Number of extensions.  (int)
    num_seqs_better_e   Number of sequences better than e-value.  (int)
    hsps_no_gap         Number of HSP's better, without gapping.  (int)
    hsps_prelim_gapped  Number of HSP's gapped in prelim test.  (int)
    hsps_prelim_gapped_attemped  Number of HSP's attempted in prelim.  (int)
    hsps_gapped         Total number of HSP's gapped.  (int)
    query_length        Length of the query.  (int)
    database_length     Number of letters in the database.  (int)
    effective_hsp_length         Effective HSP length.  (int)
    effective_query_length       Effective length of query.  (int)
    effective_database_length    Effective length of database.  (int)
    effective_search_space       Effective search space.  (int)
    effective_search_space_used  Effective search space used.  (int)
    frameshift          Frameshift window.  Tuple of (int, float)
    threshold           Threshold.  (int)
    window_size         Window size.  (int)
    dropoff_1st_pass    Tuple of (score, bits).  (int, float)
    gap_x_dropoff       Tuple of (score, bits).  (int, float)
    gap_x_dropoff_final Tuple of (score, bits).  (int, float)
    gap_trigger         Tuple of (score, bits).  (int, float)
    blast_cutoff        Tuple of (score, bits).  (int, float)
    
    """
    def __init__(self):
        self.application = ''
        self.version = ''
        self.date = ''
        self.reference = ''

        self.query = ''
        self.query_letters = None

        self.database = ''
        self.database_sequences = None
        self.database_letters = None

        self.descriptions = []
        self.alignments = []

        self.posted_date = ''
        self.ka_params = (None, None, None)
        self.gapped = None # XXX ???
        self.ka_params_gap = (None, None, None)

        self.matrix = ''
        self.gap_penalties = (None, None)
        self.num_hits = None
        self.num_sequences = None
        self.num_good_extends = None
        self.num_seqs_better_e = None
        self.hsps_no_gap = None
        self.hsps_prelim_gapped = None
        self.hsps_prelim_gapped_attemped = None
        self.hsps_gapped = None
        self.query_length = None
        self.database_length = None
        self.effective_hsp_length = None
        self.effective_query_length = None
        self.effective_database_length = None
        self.effective_search_space = None
        self.effective_search_space_used = None
        self.frameshift = (None, None)
        self.threshold = None
        self.window_size = None
        self.dropoff_1st_pass = (None, None)
        self.gap_x_dropoff = (None, None)
        self.gap_x_dropoff_final = (None, None)
        self.gap_trigger = (None, None)
        self.blast_cutoff = (None, None)

class Description:
    """Stores information about one hit in the descriptions section.

    Members:
    title       Title of the hit.
    score       Number of bits.  (int)
    p           P value.  (float)
    
    """
    def __init__(self):
        self.title = ''
        self.score = None
        self.p = None

class Alignment:
    """Stores information about one hit in the alignments section.

    Members:
    title      Name.
    length     Length.  (int)
    hsps       A list of HSP objects.

    """
    def __init__(self):
        self.title = ''
        self.length = None
        self.hsps = []

class HSP:
    """Stores information about one hsp in an alignment hit.

    Members:
    score        BLAST score of hit.  (float)
    bits         Number of bits for that score.  (float)
    expect       Expect value.  (float)
    identities   Number of identities.  (int)
    positives    Number of positives.  (int)
    gaps         Number of gaps (gapped BLAST only).  XXX ???
    strand       Tuple of (query, target) strand.
    frame        Tuple of 1 or 2 frame shifts, depending on the flavor.

    query        List of query sequences.
    query_start  List of where the query sequences start.
    match        List of the match sequences.
    sbjct        List of the sbjct sequences.
    sbjct_start  List of where the sbjct sequences start.
    
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
        self.gaps = None
        self.strand = ()
        self.frame = ()
        
        self.query = []
        self.query_start = []
        self.match = []
        self.sbjct = []
        self.sbjct_start = []
