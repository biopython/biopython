# Copyright 2006 by Sean Davis.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# $Id: __init__.py,v 1.11 2008-07-19 12:43:40 peterc Exp $
# Sean Davis <sdavis2 at mail dot nih dot gov>
# National Cancer Institute
# National Institutes of Health
# Bethesda, MD, USA
#

"""
Parse Unigene flat file format files such as the Hs.data file.

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
       GENE_ID      Entrez gene identifier associated with at least one sequence in this cluster; 
	            to be used instead of LocusLink.  
       LOCUSLINK    LocusLink identifier associated with at least one sequence in this cluster;  
		    deprecated in favor of GENE_ID
       CHROMOSOME   Chromosome.  For plants, CHROMOSOME refers to mapping on the arabidopsis genome.
       STS          STS
            NAME=        Name of STS
            ACC=         GenBank/EMBL/DDBJ accession number of STS [optional field]
            DSEG=        GDB Dsegment number [optional field]
            UNISTS=      identifier in NCBI's UNISTS database
       TXMAP        Transcript map interval
            MARKER=      Marker found on at least one sequence in this cluster
            RHPANEL=     Radiation Hybrid panel used to place marker
       PROTSIM      Protein Similarity data for the sequence with highest-scoring protein similarity in this cluster
            ORG=         Organism
            PROTGI=      Sequence GI of protein
            PROTID=      Sequence ID of protein
            PCT=         Percent alignment
            ALN=         length of aligned region (aa)
       SCOUNT       Number of sequences in the cluster
       SEQUENCE     Sequence
            ACC=         GenBank/EMBL/DDBJ accession number of sequence
            NID=         Unique nucleotide sequence identifier (gi)
            PID=         Unique protein sequence identifier (used for non-ESTs)
            CLONE=       Clone identifier (used for ESTs only)
            END=         End (5'/3') of clone insert read (used for ESTs only) 
            LID=         Library ID; see Hs.lib.info for library name and tissue  	
            MGC=	 5' CDS-completeness indicator; if present, 
			 the clone associated with this sequence  
			 is believed CDS-complete. A value greater than 511
			 is the gi of the CDS-complete mRNA matched by the EST,
 	 		 otherwise the value is an indicator of the reliability
                         of the test indicating CDS comleteness;
 			 higher values indicate more reliable CDS-completeness predictions. 
           SEQTYPE=	 Description of the nucleotide sequence. Possible values are
			 mRNA, EST and HTC.
           TRACE=	 The Trace ID of the EST sequence, as provided by NCBI Trace Archive
           PERIPHERAL=   Indicator that the sequence is a suboptimal 
	                 representative of the gene represented by this cluster.
                         Peripheral sequences are those that are in a cluster
                         which represents a spliced gene without sharing a
                         splice junction with any other sequence.  In many
                         cases, they are unspliced transcripts originating
                         from the gene.

       //           End of record
"""
from Bio.ParserSupport import *
import re

#
# CONSTANTS
#
UG_INDENT=12

class UnigeneSequenceRecord:
    """Store the information for one SEQUENCE line from a Unigene file

    Initialize with the text part of the SEQUENCE line, or nothing.

    Attributes and descriptions (access as LOWER CASE)
    ACC=         GenBank/EMBL/DDBJ accession number of sequence
    NID=         Unique nucleotide sequence identifier (gi)
    PID=         Unique protein sequence identifier (used for non-ESTs)
    CLONE=       Clone identifier (used for ESTs only)
    END=         End (5'/3') of clone insert read (used for ESTs only) 
    LID=         Library ID; see Hs.lib.info for library name and tissue  	
    MGC=	 5' CDS-completeness indicator; if present, 
                 the clone associated with this sequence  
                 is believed CDS-complete. A value greater than 511
                 is the gi of the CDS-complete mRNA matched by the EST,
                 otherwise the value is an indicator of the reliability
                 of the test indicating CDS comleteness;
                 higher values indicate more reliable CDS-completeness predictions. 
    SEQTYPE=	 Description of the nucleotide sequence. Possible values are
                 mRNA, EST and HTC.
    TRACE=	 The Trace ID of the EST sequence, as provided by NCBI Trace Archive
    PERIPHERAL=   Indicator that the sequence is a suboptimal 
                  representative of the gene represented by this cluster.
                  Peripheral sequences are those that are in a cluster
                  which represents a spliced gene without sharing a
                  splice junction with any other sequence.  In many
                  cases, they are unspliced transcripts originating
                  from the gene.
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

    Initialize with the text part of the STS line, or nothing.

    Attributes and descriptions (access as LOWER CASE)

    NAME=        Name of STS
    ACC=         GenBank/EMBL/DDBJ accession number of STS [optional field]
    DSEG=        GDB Dsegment number [optional field]
    UNISTS=      identifier in NCBI's UNISTS database
    """

    def __init__(self,text=None):
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
        self.chromosome   = ''  # CHROMOSOME
        self.protsim      = []  # PROTSIM entries, array of Protsims
                                # Type UnigeneProtsimRecord
        self.sequence     = []  # SEQUENCE entries, array of Sequence entries
                                # Type UnigeneSequenceRecord
        self.sts          = []  # STS entries, array of STS entries
                                # Type UnigeneSTSRecord
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
        self.chromosome   = ''  # CHROMOSOME
        self.protsim      = []  # PROTSIM entries, array of Protsims
        self.sequence     = []  # SEQUENCE entries, array of Sequence entries
        self.sts          = []  # STS entries, array of STS entries
        self.txmap        = []  # TXMAP entries, array of TXMap entries

    def __repr__(self):
        return "<%s> %s %s\n%s" % (self.__class__.__name__,
                          self.ID, self.symbol, self.title)


class _RecordConsumer(AbstractConsumer):

    def __init__(self):
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
    """Scans a Unigene Flat File Format file
    """

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
	def __init__(self):
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
	def __init__(self, handle, parser=None):
		self._uhandle = File.UndoHandle(handle)

	def next(self):
		self._parser = RecordParser()
		lines = []
		while 1:
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
