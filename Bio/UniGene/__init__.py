# Copyright 2006 by Sean Davis.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# $Id: __init__.py,v 1.12 2009-04-24 12:03:45 mdehoon Exp $
# Sean Davis <sdavis2 at mail dot nih dot gov>
# National Cancer Institute
# National Institutes of Health
# Bethesda, MD, USA
#

"""Parse Unigene flat file format files such as the Hs.data file.

Here is an overview of the flat file format that this parser deals with:
   Line types/qualifiers:

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
    """Store the information for one SEQUENCE line from a Unigene file

    Initialize with the text part of the SEQUENCE line, or nothing.

    Attributes and descriptions (access as LOWER CASE)
    ACC=         GenBank/EMBL/DDBJ accession number of sequence
    NID=         Unique nucleotide sequence identifier (gi)
    PID=         Unique protein sequence identifier (used for non-ESTs)
    CLONE=       Clone identifier (used for ESTs only)
    END=         End (5'/3') of clone insert read (used for ESTs only) 
    LID=         Library ID; see Hs.lib.info for library name and tissue
    MGC=         5' CDS-completeness indicator; if present, 
                 the clone associated with this sequence  
                 is believed CDS-complete. A value greater than 511
                 is the gi of the CDS-complete mRNA matched by the EST,
                 otherwise the value is an indicator of the reliability
                 of the test indicating CDS completeness;
                 higher values indicate more reliable CDS-completeness
                 predictions. 
    SEQTYPE=     Description of the nucleotide sequence. Possible values
                 are mRNA, EST and HTC.
    TRACE=       The Trace ID of the EST sequence, as provided by NCBI
                 Trace Archive
    """
    
    def __init__(self,text=None):
        self.acc = ''
        self.nid = ''
        self.lid = ''
        self.pid = ''
        self.clone = ''
        self.image = ''
        self.is_image = False
        self.end = ''
        self.mgc = ''
        self.seqtype = ''
        self.trace = ''
        if not text==None:
            self.text=text
            self._init_from_text(text)

    def _init_from_text(self,text):
        parts = text.split('; ');
        for part in parts:
            key, val = part.split("=")
            if key=='CLONE':
                if val[:5]=='IMAGE':
                    self.is_image=True
                    self.image = val[6:]
            setattr(self,key.lower(),val)

    def __repr__(self):
        return self.text
        

class ProtsimLine:
    """Store the information for one PROTSIM line from a Unigene file

    Initialize with the text part of the PROTSIM line, or nothing.

    Attributes and descriptions (access as LOWER CASE)
    ORG=         Organism
    PROTGI=      Sequence GI of protein
    PROTID=      Sequence ID of protein
    PCT=         Percent alignment
    ALN=         length of aligned region (aa)
    """

    def __init__(self,text=None):
        self.org = ''
        self.protgi = ''
        self.protid = ''
        self.pct = ''
        self.aln = ''
        if not text==None:
            self.text=text
            self._init_from_text(text)

    def _init_from_text(self,text):
        parts = text.split('; ');
        
        for part in parts:
            key, val = part.split("=")
            setattr(self,key.lower(),val)

    def __repr__(self):
        return self.text
        

class STSLine:
    """Store the information for one STS line from a Unigene file

    Initialize with the text part of the STS line, or nothing.

    Attributes and descriptions (access as LOWER CASE)

    ACC=         GenBank/EMBL/DDBJ accession number of STS [optional field]
    UNISTS=      identifier in NCBI's UNISTS database
    """

    def __init__(self,text=None):
        self.acc = ''
        self.unists = ''
        if not text==None:
            self.text=text
            self._init_from_text(text)

    def _init_from_text(self,text):
        parts = text.split(' ');
        
        for part in parts:
            key, val = part.split("=")
            setattr(self,key.lower(),val)

    def __repr__(self):
        return self.text
        

class Record:
    """Store a Unigene record

    Here is what is stored:
    
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
        self.ID           = ''  # ID line
        self.species      = ''  # Hs, Bt, etc.
        self.title        = ''  # TITLE line
        self.symbol       = ''  # GENE line
        self.cytoband     = ''  # CYTOBAND line
        self.express      = []  # EXPRESS line, parsed on ';'
        self.restr_expr   = ''  # RESTR_EXPR line
        self.gnm_terminus = ''  # GNM_TERMINUS line
        self.gene_id      = ''  # GENE_ID line
        self.locuslink    = ''  # LOCUSLINK line
        self.homol        = ''  # HOMOL line
        self.chromosome   = ''  # CHROMOSOME line
        self.protsim      = []  # PROTSIM entries, array of Protsims
        self.sequence     = []  # SEQUENCE entries, array of Sequence entries
        self.sts          = []  # STS entries, array of STS entries
        self.txmap        = []  # TXMAP entries, array of TXMap entries

    def __repr__(self):
        return "<%s> %s %s\n%s" % (self.__class__.__name__,
                          self.ID, self.symbol, self.title)

def parse(handle):
    while True:
        record = _read(handle)
        if not record:
            return
        yield record


def read(handle):
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
        if tag=="ID":
            record = Record()
            record.ID = value
            record.species = record.ID.split('.')[0]
        elif tag=="TITLE":
            record.title = value
        elif tag=="GENE":
            record.symbol = value
        elif tag=="GENE_ID":
            record.gene_id = value
        elif tag=="LOCUSLINK":
            record.locuslink = value
        elif tag=="HOMOL":
            if value=="YES":
                record.homol = True
            elif value=="NO":
                record.homol = True
            else:
                raise ValueError, "Cannot parse HOMOL line %s" % line
        elif tag=="EXPRESS":
            record.express = [word.strip() for word in value.split("|")]
        elif tag=="RESTR_EXPR":
            record.restr_expr = [word.strip() for word in value.split("|")]
        elif tag=="CHROMOSOME":
            record.chromosome = value
        elif tag=="CYTOBAND":
            record.cytoband = value
        elif tag=="PROTSIM":
            protsim = ProtsimLine(value)
            record.protsim.append(protsim)
        elif tag=="SCOUNT":
            scount = int(value)
        elif tag=="SEQUENCE":
            sequence = SequenceLine(value)
            record.sequence.append(sequence)
        elif tag=="STS":
            sts = STSLine(value)
            record.sts.append(sts)
        elif tag=='//':
            if len(record.sequence)!=scount:
                raise ValueError, "The number of sequences specified in the record (%d) does not agree with the number of sequences found (%d)" % (scount, len(record.sequence))
            return record
        else:
            raise ValueError, "Unknown tag %s" % tag
    if record:
        raise ValueError("Unexpected end of stream.")


# Everything below is deprecated


from Bio.ParserSupport import *
import re
import Bio

#
# CONSTANTS
#
UG_INDENT=12

class UnigeneSequenceRecord:
    """Store the information for one SEQUENCE line from a Unigene file
    (DEPRECATED).

    Initialize with the text part of the SEQUENCE line, or nothing.

    Attributes and descriptions (access as LOWER CASE)
    ACC=         GenBank/EMBL/DDBJ accession number of sequence
    NID=         Unique nucleotide sequence identifier (gi)
    PID=         Unique protein sequence identifier (used for non-ESTs)
    CLONE=       Clone identifier (used for ESTs only)
    END=         End (5'/3') of clone insert read (used for ESTs only) 
    LID=         Library ID; see Hs.lib.info for library name and tissue
    MGC=         5' CDS-completeness indicator; if present, 
                 the clone associated with this sequence  
                 is believed CDS-complete. A value greater than 511
                 is the gi of the CDS-complete mRNA matched by the EST,
                 otherwise the value is an indicator of the reliability
                 of the test indicating CDS comleteness;
                 higher values indicate more reliable CDS-completeness predictions. 
    SEQTYPE=     Description of the nucleotide sequence. Possible values are
                 mRNA, EST and HTC.
    TRACE=       The Trace ID of the EST sequence, as provided by NCBI Trace Archive
    PERIPHERAL=  Indicator that the sequence is a suboptimal 
                 representative of the gene represented by this cluster.
                 Peripheral sequences are those that are in a cluster
                 which represents a spliced gene without sharing a
                 splice junction with any other sequence.  In many
                 cases, they are unspliced transcripts originating
                 from the gene.

    This class is DEPRECATED; please use the read() function in this module
    instead.
    """
    
    def __init__(self,text=None):
        import warnings
        warnings.warn("Bio.UniGene.UnigeneSequenceRecord is deprecated; please use the read() function in this module instead", Bio.BiopythonDeprecationWarning)
        self.acc = ''
        self.nid = ''
        self.lid = ''
        self.pid = ''
        self.clone = ''
        self.image = ''
        self.is_image = False
        self.end = ''
        self.mgc = ''
        self.seqtype = ''
        self.Trace = ''
        self.peripheral = ''
        if not text==None:
            self.text=text
            return self._init_from_text(text)

    def _init_from_text(self,text):
        parts = text.split('; ');
        for part in parts:
            key,val = re.match('(\w+)=(\S+)',part).groups()
            if key=='CLONE':
                if val[:5]=='IMAGE':
                    self.is_image=True
                    self.image = val[6:]
            setattr(self,key.lower(),val)

    def __repr__(self):
        return self.text
        

class UnigeneProtsimRecord:
    """Store the information for one PROTSIM line from a Unigene file
    (DEPRECATED).

    Initialize with the text part of the PROTSIM line, or nothing.

    Attributes and descriptions (access as LOWER CASE)
    ORG=         Organism
    PROTGI=      Sequence GI of protein
    PROTID=      Sequence ID of protein
    PCT=         Percent alignment
    ALN=         length of aligned region (aa)

    This class is DEPRECATED; please use the read() function in this module
    instead.
    """

    def __init__(self,text=None):
        import warnings
        warnings.warn("Bio.UniGene.UnigeneProtsimRecord is deprecated; please use the read() function in this module instead", Bio.BiopythonDeprecationWarning)
        self.org = ''
        self.protgi = ''
        self.protid = ''
        self.pct = ''
        self.aln = ''
        if not text==None:
            self.text=text
            return self._init_from_text(text)

    def _init_from_text(self,text):
        parts = text.split('; ');
        
        for part in parts:
            key,val = re.match('(\w+)=(\S+)',part).groups()
            setattr(self,key.lower(),val)

    def __repr__(self):
        return self.text
        

class UnigeneSTSRecord:
    """Store the information for one STS line from a Unigene file
    (DEPRECATED).

    Initialize with the text part of the STS line, or nothing.

    Attributes and descriptions (access as LOWER CASE)

    NAME=        Name of STS
    ACC=         GenBank/EMBL/DDBJ accession number of STS [optional field]
    DSEG=        GDB Dsegment number [optional field]
    UNISTS=      identifier in NCBI's UNISTS database

    This class is DEPRECATED; please use the read() function in this module
    instead.
    """

    def __init__(self,text=None):
        import warnings
        warnings.warn("Bio.UniGene.UnigeneSTSRecord is deprecated; please use the read() function in this module instead", Bio.BiopythonDeprecationWarning)
        self.name = ''
        self.acc = ''
        self.dseg = ''
        self.unists = ''
        if not text==None:
            self.text=text
            return self._init_from_text(text)

    def _init_from_text(self,text):
        parts = text.split(' ');
        
        for part in parts:
            key,val = re.match('(\w+)=(\S+)',part).groups()
            setattr(self,key.lower(),val)

    def __repr__(self):
        return self.text
        

class UnigeneRecord:
    """Store a Unigene record (DEPRECATED).

    Here is what is stored:
    
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
        self.chromosome   = ''  # CHROMOSOME
        self.protsim      = []  # PROTSIM entries, array of Protsims
                                # Type UnigeneProtsimRecord
        self.sequence     = []  # SEQUENCE entries, array of Sequence entries
                                # Type UnigeneSequenceRecord
        self.sts          = []  # STS entries, array of STS entries
                                # Type UnigeneSTSRecord
        self.txmap        = []  # TXMAP entries, array of TXMap entries

    This class is DEPRECATED; please use the read() function in this module
    instead.
    """

    def __init__(self):
        import warnings
        warnings.warn("Bio.UniGene.UnigeneRecord is deprecated; please use the read() function in this module instead", Bio.BiopythonDeprecationWarning)
        self.ID           = ''  # ID line
        self.species      = ''  # Hs, Bt, etc.
        self.title        = ''  # TITLE line
        self.symbol       = ''  # GENE line
        self.cytoband     = ''  # CYTOBAND line
        self.express      = []  # EXPRESS line, parsed on ';'
        self.restr_expr   = ''  # RESTR_EXPR line
        self.gnm_terminus = ''  # GNM_TERMINUS line
        self.gene_id      = ''  # GENE_ID line
        self.chromosome   = ''  # CHROMOSOME
        self.protsim      = []  # PROTSIM entries, array of Protsims
        self.sequence     = []  # SEQUENCE entries, array of Sequence entries
        self.sts          = []  # STS entries, array of STS entries
        self.txmap        = []  # TXMAP entries, array of TXMap entries

    def __repr__(self):
        return "<%s> %s %s\n%s" % (self.__class__.__name__,
                          self.ID, self.symbol, self.title)


class _RecordConsumer(AbstractConsumer):
    """This class is DEPRECATED; please use the read() function in this module
    instead."""

    def __init__(self):
        import warnings
        warnings.warn("Bio.UniGene._RecordConsumer is deprecated; please use the read() function in this module instead", Bio.BiopythonDeprecationWarning)
        self.unigene_record = UnigeneRecord()
    def ID(self,line):
        self.unigene_record.ID = self._get_single_entry(line)
        self.unigene_record.species = self.unigene_record.ID.split('.')[0]
    def TITLE(self,line):
        self.unigene_record.title = self._get_single_entry(line)
    def GENE(self,line):
        self.unigene_record.symbol = self._get_single_entry(line)
    def EXPRESS(self,line):
        self.unigene_record.express = self._get_array_entry(line,split_on='; ')
    def RESTR_EXPR(self,line):
        self.unigene_record.restr_expr = self._get_single_entry(line)
    def GENE_ID(self,line):
        self.unigene_record.gene_id = self._get_single_entry(line)
    def CHROMOSOME(self,line):
        self.unigene_record.chromosome = self._get_single_entry(line)
    def GENE_ID(self,line):
        self.unigene_record.gene_id = self._get_single_entry(line)
    def SEQUENCE(self,line):
        ug_seqrecord = UnigeneSequenceRecord(self._get_single_entry(line))
        self.unigene_record.sequence.append(ug_seqrecord)
    def PROTSIM(self,line):
        ug_protsimrecord = UnigeneProtsimRecord(self._get_single_entry(line))
        self.unigene_record.protsim.append(ug_protsimrecord)
    def STS(self,line):
        ug_stsrecord = UnigeneSTSRecord(self._get_single_entry(line))
        self.unigene_record.sts.append(ug_stsrecord)
    

    def _get_single_entry(self,line):
        """Consume a single-value line
        """
        return line[UG_INDENT:]

    def _get_array_entry(self,line,split_on):
        """Consume a multi-value line by splitting on split_on
        """
        return line[UG_INDENT:].split(split_on)
    

class _Scanner:
    """Scans a Unigene Flat File Format file (DEPRECATED).

    This class is DEPRECATED; please use the read() function in this module
    instead.
    """

    def __init__(self):
        import warnings
        warnings.warn("Bio.UniGene._Scanner is deprecated; please use the read() function in this module instead", Bio.BiopythonDeprecationWarning)

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed events from parsing a Unigene file to a consumer.
        handle is a file-like object, and consumer is a consumer object
        that will receive events as the file is scanned

        """
        consumer.start_record()
        for line in handle:
            tag = line.split(' ')[0]
            line = line.rstrip()
            if line=='//':
                consumer.end_record()
                break
            try:
                f = getattr(consumer, tag)
            except AttributeError:
                print 'no method called', tag
            else:
                if callable(f):
                    f(line)

        
class RecordParser(AbstractParser):
    """This class is DEPRECATED; please use the read() function in this module
    instead."""

    def __init__(self):
        import warnings
        warnings.warn("Bio.UniGene._RecordParser is deprecated; please use the read() function in this module instead", Bio.BiopythonDeprecationWarning)
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
            self._scanner.feed(uhandle, self._consumer)
        return self._consumer.unigene_record

class Iterator:
    """This class is DEPRECATED; please use the parse() function in this module
    instead."""

    def __init__(self, handle, parser=None):
        import warnings
        warnings.warn("Bio.UniGene.Iterator is deprecated; please use the parse() function in this module instead", Bio.BiopythonDeprecationWarning)
        self._uhandle = File.UndoHandle(handle)

    def next(self):
        self._parser = RecordParser()
        lines = []
        while True:
            line = self._uhandle.readline()
            if not line: break
            if line[:2] == '//':
                break
            lines.append(line)
        if not lines:
            return None
        lines.append('//')
        data = ''.join(lines)
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

    def __iter__(self):
        return iter(self.next, None)
