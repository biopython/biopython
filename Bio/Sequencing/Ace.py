"""
Parser for ACE files output by PHRAP.

version 1.0
Written by Frank Kauff (fkauff@duke.edu) and
Cymon J. Cox (cymon@duke.edu)

Tags RT, CT and WA are *not* supported
Uses the Biopython Parser interface: ParserSupport.py

Usage:

There are two ways of reading an ace file: The ACEParser() reads
the whole file at once and the RecordParser() reads contig after
contig.
    
1) Parse whole ace file at once:
        aceparser=Ace.ACEParser()
        acefilerecord=aceparser.parse(open('my_ace_file.ace','r'))

gives you 
        acefilerecord.contigs (the number of contigs in the ace file)
        acefilerecord.reads (the number of reads in the ace file)
        acefilerecord.records[] (one record for each contig)

The latter one includes the data for each contig. All the info in the
        co and bq tags is unique for each contig: contig_name, bases, reads
        segments, coru, sequence, quality, like in 
        
        acefilerecord.records[0].sequences (the contig sequence)

        af, bs, qa, ds appear multiple times, like in 

        acefilerecord.records[0].af[1].padded_start
        acefilerecord.records[1].rd[17].read_tags

        see _RecordConsumer for details.


2) Or you can iterate over the contigs of an ace file one by one in the ususal way:        
        recordparser=Ace.RecordParser()
        iterator=Ace.Iterator(open('my_ace_file.ace','r'),recordparser)
        while 1:
            contig=iterator.next()
            if contig is None:
                break
            ...

"""

import os
from types import *

from Bio import File
from Bio import Index
#from Bio import Seq
#from Bio import SeqRecord
from Bio.ParserSupport import *
from Bio.Alphabet import IUPAC


class _rd:
    def __init__(self):
        self.name=''
        self.padded_bases=None
        self.info_items=None
        self.read_tags=None
        self.sequence=''

class _qa:
    def __init__(self):
        self.qual_clipping_start=None
        self.qual_clipping_end=None
        self.align_clipping_start=None
        self.align_clipping_end=None

class _ds:
    def __init__(self):
        self.chromat_file=''
        self.phd_file=''
        self.time=''
        self.chem=''
        self.dye=''
        self.template=''
        self.direction=''
    
class _af:
    def __init__(self):
        self.name=''
        self.coru=None
        self.padded_start=None

class _bs:
    def __init__(self):
        self.name=''
        self.padded_start=None
        self.padded_end=None
        
class Record:
    """Hold information from a ACE record

    """
    def __init__(self):
        self.contig_name = ''
        self.bases=None
        self.reads=None
        self.segments=None
        self.uorc=None
        self.sequence=None
        self.quality=None
        self.af=[]
        self.bs=[]
        self.rd=[]
        self.qa=[]
        self.ds=[]
        
class Iterator:
    """Iterates over a ACE-file with multiple contigs
    
    Methods: 
    next    Return the next record from the stream, or None.
    """

    def __init__(self, handle, parser=None):
        """__init__(self, handle, parser=None)
        
        Create a new iterator.  handle is a file-like object.  parser
        is an optional Parser object to change the results into another form.
        If set to None, then the raw contents of the file will be returned.
        """

        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise ValueError, "I expected a file handle or file-like object"
        self._uhandle = File.UndoHandle(handle)
        self._parser = parser

    def next(self):
        """next(self) -> object

        Return the next contig record from the file. If no more records
        return None.
        """

        lines = []
        while 1: 
            # if at beginning, skip the AS and look for first CO command
            line=self._uhandle.readline()
            if not line:                    # empty or corrupt file
                return None
            if line[:2]=='CO':
                lines.append(line)
                break
        while 1:
            line = self._uhandle.readline()
            if not line:
                break
            # If a new record, then put the line back and stop.
            if lines and line[:2] == 'CO':
                self._uhandle.saveline(line)
                break
            lines.append(line)

        if not lines:
            return None

        data = ''.join(lines)
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

class RecordParser(AbstractParser):
    """Parses ACE file data into a Record object

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        self._scanner.feed(uhandle, self._consumer)
        return self._consumer.data

class ACEFileRecord:
    """Holds data of ACE file.
    """
    def __init__(self):
        self.records=[]
        self.contigs=None
        self.reads=None

class ACEParser(AbstractParser):
    """Parses full ACE file in list of records.
    """

    def __init__(self):
        self.data=ACEFileRecord()
        
    def parse(self,handle):
        firstline=handle.readline()
        if firstline[:2]!='AS':
            raise SyntaxError, "File does not start with 'AS'."
        self.data.contigs=eval(firstline.split()[1])
        self.data.reads=eval(firstline.split()[2])
        recparser=RecordParser()
        iter=Iterator(handle,recparser)
        while 1:
            rec=iter.next()
            if not rec:
                break
            self.data.records.append(rec)
        return self.data
        
class _Scanner:
    """Scans a ACE-formatted file
    
    Methods:
    feed - Feed one ACE record.
    """
    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in ACE data for scanning.  handle is a file-like object
        containing ACE data.  consumer is a Consumer object that will
        receive events as the ACE data is scanned.
        """
        assert isinstance(handle, File.UndoHandle), \
            "handle must be an UndoHandle"
        if handle.peekline():
            self._scan_record(handle, consumer)

    def _scan_record(self, uhandle, consumer):
        consumer.begin_contig()
        read_and_call(uhandle,consumer.co_header,start='CO ')
        self._scan_co_data(uhandle, consumer)
        read_and_call_while(uhandle,consumer.noevent,blank=1)
        read_and_call(uhandle,consumer.bq_header,start='BQ')
        self._scan_bq_data(uhandle, consumer)
        read_and_call_while(uhandle,consumer.noevent,blank=1)
        read_and_call_while(uhandle,consumer.af,start='AF ')
        read_and_call_while(uhandle,consumer.noevent,blank=1)
        read_and_call_while(uhandle,consumer.bs,start='BS ')
        while 1:
            read_and_call_until(uhandle,consumer.noevent,start='RD ')
            read_and_call(uhandle,consumer.rd_header,start='RD ')
            self._scan_rd_data(uhandle, consumer)
            read_and_call_while(uhandle,consumer.noevent,blank=1)
            read_and_call(uhandle,consumer.qa,start='QA ')
            read_and_call_while(uhandle,consumer.noevent,blank=1)
            read_and_call(uhandle,consumer.ds,start='DS ')
            while 1:
                try:
                    line=safe_peekline(uhandle)
                except SyntaxError:
                    break
                else:
                    if line[:2]=='RD':
                        break
                    line=safe_readline(uhandle)
            if not line[:2]=='RD':
                break
        consumer.end_contig()
    
    def _scan_co_data(self, uhandle, consumer):
        consumer.co_data(self._scan_sequence_data(uhandle))
            
    def _scan_rd_data(self, uhandle, consumer):
        consumer.rd_data(self._scan_sequence_data(uhandle))
    
    def _scan_bq_data(self, uhandle, consumer):
        qual=''
        while 1:
            line=uhandle.readline()
            if is_blank_line(line):
                uhandle.saveline(line)
                break
            qual+=' '+line
        consumer.bq_data(qual)
    
    def _scan_sequence_data(self,uhandle):
        seq=''
        while 1:
            line=uhandle.readline()
            if is_blank_line(line):
                uhandle.saveline(line)
                break
            seq+=line.strip()
        return seq
        
class _RecordConsumer(AbstractConsumer):
    """Reads the ace tags into data records."""
    
    def __init__(self):
        self.data = None

    def begin_contig(self):
        self.data = Record()

    def end_contig(self):
        pass

    def bq_header(self,line):
        pass
    
    def bq_data(self,qual):
        self.data.quality=map(eval,qual.split())

    def co_header(self,line):
        header=line.split()
        self.data.contig_name=header[1]
        self.data.bases=eval(header[2])
        self.data.reads=eval(header[3])
        self.data.segments=eval(header[4])
        self.data.uorc=header[5]

    def co_data(self,seq):
        self.data.sequence=seq
    
    def af(self,line):
        header=line.split()
        afdata=_af()
        afdata.name=header[1]
        afdata.coru=header[2]
        afdata.padded_start=eval(header[3])
        self.data.af.append(afdata)
        
    def bs(self,line):
        header=line.split()
        bsdata=_bs()
        bsdata.padded_start=eval(header[1])
        bsdata.padded_end=eval(header[2])
        bsdata.name=header[3]
        self.data.bs.append(bsdata)
        
    def rd_header(self,line):
        header=line.split()
        rddata=_rd()
        rddata.name=header[1]
        rddata.padded_bases=eval(header[2])
        rddata.info_items=eval(header[3])
        rddata.read_tags=eval(header[4])
        self.data.rd.append(rddata)

    def rd_data(self,seq):
        self.data.rd[-1].sequence=seq
        
    def qa(self,line):
        header=map(eval,line.split()[1:])
        qadata=_qa()
        qadata.qual_clipping_start=header[0]
        qadata.qual_clipping_end=header[1]
        qadata.align_clipping_start=header[2]
        qadata.align_clipping_end=header[3]
        self.data.qa.append(qadata)

    def ds(self,line):
        dsdata=_ds()
        tags=['CHROMAT_FILE','PHD_FILE','TIME','CHEM','DYE','TEMPLATE','DIRECTION']
        poss=map(line.find,tags)
        tagpos=dict(zip(poss,tags))
        if tagpos.has_key(-1):
            del tagpos[-1]
        ps=tagpos.keys()
        ps.sort()
        for (p1,p2) in zip(ps,ps[1:]+[len(line)+1]):
            setattr(dsdata,tagpos[p1].lower(),line[p1+len(tagpos[p1])+1:p2])   
        self.data.ds.append(dsdata)

