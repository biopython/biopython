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
    application
    version
    date
    reference

    query
    query_letters
    
    database
    database_sequences
    database_letters

    descriptions  A list of Description objects.
    alignments    A list of Alignment objects.

    posted_date
    ka_params
    gapped        # XXX this isn't set right!
    ka_params_gap
    
    matrix
    gap_penalties
    num_hits
    num_sequences
    num_good_extends
    num_seqs_better_e
    hsps_no_gap
    hsps_prelim_gapped
    hsps_prelim_gapped_attemped
    hsps_gapped
    query_length
    database_length
    effective_hsp_length
    effective_query_length
    effective_database_length
    effective_search_space
    effective_search_space_used
    frameshift
    threshold
    window_size
    dropoff_1st_pass
    gap_x_dropoff
    gap_x_dropoff_final
    gap_trigger
    blast_cutoff
    
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

        self.posted_date = ''
        self.ka_params = (None, None, None)
        self.gapped = None # XXX ???
        self.ka_params_gap = (None, None, None)

        self._matrix = ''
        self._gap_penalties = ()
        self._num_hits = None
        self._num_sequences = None
        self._num_good_extends = None
        self._num_seqs_better_e = None
        self._hsps_no_gap = None
        self._hsps_prelim_gapped = None
        self._hsps_prelim_gapped_attemped = None
        self._hsps_gapped = None
        self._query_length = None
        self._database_length = None
        self._effective_hsp_length = None
        self._effective_query_length = None
        self._effective_database_length = None
        self._effective_search_space = None
        self._effective_search_space_used = None
        self._frameshift = ()
        self._threshold = None
        self._window_size = None
        self._dropoff_1st_pass = ()
        self._gap_x_dropoff = ()
        self._gap_x_dropoff_final = ()
        self._gap_trigger = ()
        self._blast_cutoff = ()

class Description:
    """Stores information about one hit in the descriptions section.

    Members:
    title       Title of the hit.
    score       Number of bits.  (integer)
    p           P value.  (float)
    
    """
    def __init__(self):
        self.title = None
        self.score = None
        self.p = None

class Alignment:
    """Stores information about one hit in the alignments section.

    Members:
    title      Name.
    length     Length.
    hsps       A list of HSP objects.

    """
    def __init__(self):
        self.title = None
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
    strand       tuple of (query, target) strand
    frame        tuple of 1 or 2 frame shifts, depending on the flavor

    # XXX change these to strings!
    query        list of query sequences
    query_start  list of where the query sequences start
    match        list of the match sequences
    sbjct        list of the sbjct sequences
    sbjct_start  list of where the sbjct sequences start
    
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
        self.score = ''
        self.bits = ''
        self.expect = ''
        self.identities = ''
        self.positives = ''
        self.strand = None
        self.frame = None
        self.gaps = ''
        
        self.query = []
        self.query_start = []
        self.match = []
        self.sbjct = []
        self.sbjct_start = []
