# Copyright 2006 by Sean Davis, National Cancer Institute, NIH.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Parse Unigene flat file format files such as the Hs.data file.

Here is an overview of the flat file format that this parser deals with:

   Line types/qualifiers::

       ID           UniGene cluster ID
       TITLE        Title for the cluster
       GENE         Gene symbol
       CYTOBAND     Cytological band
       EXPRESS      Tissues of origin for ESTs in cluster
       RESTR_EXPR   Single tissue or development stage contributes
                    more than half the total EST frequency for this gene.
       GNM_TERMINUS genomic confirmation of presence of a 3' terminus;
                    T if a non-templated polyA tail is found among
                    a cluster's sequences; else
                    I if templated As are found in genomic sequence or
                    S if a canonical polyA signal is found on
                      the genomic sequence
       GENE_ID      Entrez gene identifier associated with at least one
                    sequence in this cluster;
                    to be used instead of LocusLink.
       LOCUSLINK    LocusLink identifier associated with at least one
                    sequence in this cluster;
                    deprecated in favor of GENE_ID
       HOMOL        Homology;
       CHROMOSOME   Chromosome.  For plants, CHROMOSOME refers to mapping
                    on the arabidopsis genome.
       STS          STS
            ACC=         GenBank/EMBL/DDBJ accession number of STS
                         [optional field]
            UNISTS=      identifier in NCBI's UNISTS database
       TXMAP        Transcript map interval
            MARKER=      Marker found on at least one sequence in this
                         cluster
            RHPANEL=     Radiation Hybrid panel used to place marker
       PROTSIM      Protein Similarity data for the sequence with
                    highest-scoring protein similarity in this cluster
            ORG=         Organism
            PROTGI=      Sequence GI of protein
            PROTID=      Sequence ID of protein
            PCT=         Percent alignment
            ALN=         length of aligned region (aa)
       SCOUNT       Number of sequences in the cluster
       SEQUENCE     Sequence
            ACC=         GenBank/EMBL/DDBJ accession number of sequence
            NID=         Unique nucleotide sequence identifier (gi)
            PID=         Unique protein sequence identifier (used for
                         non-ESTs)
            CLONE=       Clone identifier (used for ESTs only)
            END=         End (5'/3') of clone insert read (used for
                         ESTs only)
            LID=         Library ID; see Hs.lib.info for library name
                         and tissue
            MGC=         5' CDS-completeness indicator; if present, the
                         clone associated with this sequence is believed
                         CDS-complete. A value greater than 511 is the gi
                         of the CDS-complete mRNA matched by the EST,
                         otherwise the value is an indicator of the
                         reliability of the test indicating CDS
                         completeness; higher values indicate more
                         reliable CDS-completeness predictions.
           SEQTYPE=      Description of the nucleotide sequence.
                         Possible values are mRNA, EST and HTC.
           TRACE=        The Trace ID of the EST sequence, as provided by
                         NCBI Trace Archive

"""


class SequenceLine:
    """Store the information for one SEQUENCE line from a Unigene file.

    Initialize with the text part of the SEQUENCE line, or nothing.

    Attributes and descriptions (access as LOWER CASE):
     - ACC=         GenBank/EMBL/DDBJ accession number of sequence
     - NID=         Unique nucleotide sequence identifier (gi)
     - PID=         Unique protein sequence identifier (used for non-ESTs)
     - CLONE=       Clone identifier (used for ESTs only)
     - END=         End (5'/3') of clone insert read (used for ESTs only)
     - LID=         Library ID; see Hs.lib.info for library name and tissue
     - MGC=         5' CDS-completeness indicator; if present,
       the clone associated with this sequence
       is believed CDS-complete. A value greater than 511
       is the gi of the CDS-complete mRNA matched by the EST,
       otherwise the value is an indicator of the reliability
       of the test indicating CDS completeness;
       higher values indicate more reliable CDS-completeness
       predictions.
     - SEQTYPE=     Description of the nucleotide sequence. Possible values
       are mRNA, EST and HTC.
     - TRACE=       The Trace ID of the EST sequence, as provided by NCBI
       Trace Archive

    """

    def __init__(self, text=None):
        """Initialize the class."""
        self.acc = ""
        self.nid = ""
        self.lid = ""
        self.pid = ""
        self.clone = ""
        self.image = ""
        self.is_image = False
        self.end = ""
        self.mgc = ""
        self.seqtype = ""
        self.trace = ""
        if text is not None:
            self.text = text
            self._init_from_text(text)

    def _init_from_text(self, text):
        parts = text.split("; ")
        for part in parts:
            key, val = part.split("=")
            if key == "CLONE":
                if val[:5] == "IMAGE":
                    self.is_image = True
                    self.image = val[6:]
            setattr(self, key.lower(), val)

    def __repr__(self):
        """Return UniGene SequenceLine object as a string."""
        return self.text


class ProtsimLine:
    """Store the information for one PROTSIM line from a Unigene file.

    Initialize with the text part of the PROTSIM line, or nothing.

    Attributes and descriptions (access as LOWER CASE)
    ORG=         Organism
    PROTGI=      Sequence GI of protein
    PROTID=      Sequence ID of protein
    PCT=         Percent alignment
    ALN=         length of aligned region (aa)
    """

    def __init__(self, text=None):
        """Initialize the class."""
        self.org = ""
        self.protgi = ""
        self.protid = ""
        self.pct = ""
        self.aln = ""
        if text is not None:
            self.text = text
            self._init_from_text(text)

    def _init_from_text(self, text):
        parts = text.split("; ")

        for part in parts:
            key, val = part.split("=")
            setattr(self, key.lower(), val)

    def __repr__(self):
        """Return UniGene ProtsimLine object as a string."""
        return self.text


class STSLine:
    """Store the information for one STS line from a Unigene file.

    Initialize with the text part of the STS line, or nothing.

    Attributes and descriptions (access as LOWER CASE)

    ACC=         GenBank/EMBL/DDBJ accession number of STS [optional field]
    UNISTS=      identifier in NCBI's UNISTS database
    """

    def __init__(self, text=None):
        """Initialize the class."""
        self.acc = ""
        self.unists = ""
        if text is not None:
            self.text = text
            self._init_from_text(text)

    def _init_from_text(self, text):
        parts = text.split(" ")

        for part in parts:
            key, val = part.split("=")
            setattr(self, key.lower(), val)

    def __repr__(self):
        """Return UniGene STSLine object as a string."""
        return self.text


class Record:
    """Store a Unigene record.

    Here is what is stored::

        self.ID           = ''  # ID line
        self.species      = ''  # Hs, Bt, etc.
        self.title        = ''  # TITLE line
        self.symbol       = ''  # GENE line
        self.cytoband     = ''  # CYTOBAND line
        self.express      = []  # EXPRESS line, parsed on ';'
                                # Will be an array of strings
        self.restr_expr   = ''  # RESTR_EXPR line
        self.gnm_terminus = ''  # GNM_TERMINUS line
        self.gene_id      = ''  # GENE_ID line
        self.locuslink    = ''  # LOCUSLINK line
        self.homol        = ''  # HOMOL line
        self.chromosome   = ''  # CHROMOSOME line
        self.protsim      = []  # PROTSIM entries, array of Protsims
                                # Type ProtsimLine
        self.sequence     = []  # SEQUENCE entries, array of Sequence entries
                                # Type SequenceLine
        self.sts          = []  # STS entries, array of STS entries
                                # Type STSLine
        self.txmap        = []  # TXMAP entries, array of TXMap entries

    """

    def __init__(self):
        """Initialize the class."""
        self.ID = ""  # ID line
        self.species = ""  # Hs, Bt, etc.
        self.title = ""  # TITLE line
        self.symbol = ""  # GENE line
        self.cytoband = ""  # CYTOBAND line
        self.express = []  # EXPRESS line, parsed on ';'
        self.restr_expr = ""  # RESTR_EXPR line
        self.gnm_terminus = ""  # GNM_TERMINUS line
        self.gene_id = ""  # GENE_ID line
        self.locuslink = ""  # LOCUSLINK line
        self.homol = ""  # HOMOL line
        self.chromosome = ""  # CHROMOSOME line
        self.protsim = []  # PROTSIM entries, array of Protsims
        self.sequence = []  # SEQUENCE entries, array of Sequence entries
        self.sts = []  # STS entries, array of STS entries
        self.txmap = []  # TXMAP entries, array of TXMap entries

    def __repr__(self):
        """Represent the UniGene Record object as a string for debugging."""
        return f"<{self.__class__.__name__}> {self.ID} {self.symbol} {self.title}"


def parse(handle):
    """Read and load a UniGene records, for files containing multiple records."""
    while True:
        record = _read(handle)
        if not record:
            return
        yield record


def read(handle):
    """Read and load a UniGene record, one record per file."""
    record = _read(handle)
    if not record:
        raise ValueError("No SwissProt record found")
    # We should have reached the end of the record by now
    remainder = handle.read()
    if remainder:
        raise ValueError("More than one SwissProt record found")
    return record


# Everything below is private


def _read(handle):
    UG_INDENT = 12
    record = None
    for line in handle:
        tag, value = line[:UG_INDENT].rstrip(), line[UG_INDENT:].rstrip()
        line = line.rstrip()
        if tag == "ID":
            record = Record()
            record.ID = value
            record.species = record.ID.split(".")[0]
        elif tag == "TITLE":
            record.title = value
        elif tag == "GENE":
            record.symbol = value
        elif tag == "GENE_ID":
            record.gene_id = value
        elif tag == "LOCUSLINK":
            record.locuslink = value
        elif tag == "HOMOL":
            if value == "YES":
                record.homol = True
            elif value == "NO":
                record.homol = True
            else:
                raise ValueError(f"Cannot parse HOMOL line {line}")
        elif tag == "EXPRESS":
            record.express = [word.strip() for word in value.split("|")]
        elif tag == "RESTR_EXPR":
            record.restr_expr = [word.strip() for word in value.split("|")]
        elif tag == "CHROMOSOME":
            record.chromosome = value
        elif tag == "CYTOBAND":
            record.cytoband = value
        elif tag == "PROTSIM":
            protsim = ProtsimLine(value)
            record.protsim.append(protsim)
        elif tag == "SCOUNT":
            scount = int(value)
        elif tag == "SEQUENCE":
            sequence = SequenceLine(value)
            record.sequence.append(sequence)
        elif tag == "STS":
            sts = STSLine(value)
            record.sts.append(sts)
        elif tag == "//":
            if len(record.sequence) != scount:
                raise ValueError(
                    "The number of sequences specified in the record "
                    "(%d) does not agree with the number of sequences found (%d)"
                    % (scount, len(record.sequence))
                )
            return record
        else:
            raise ValueError(f"Unknown tag {tag}")
    if record:
        raise ValueError("Unexpected end of stream.")
