"""
Parser for (new) ACE files output by PHRAP.

version 1.3, 05/06/2004
Written by Frank Kauff (fkauff@duke.edu) and
Cymon J. Cox (cymon@duke.edu)

Uses the Biopython Parser interface: ParserSupport.py

Usage:

There are two ways of reading an ace file: The ACEParser() reads
the whole file at once and the RecordParser() reads contig after
contig.
    
1) Parse whole ace file at once:
        aceparser=Ace.ACEParser()
        acefilerecord=aceparser.parse(open('my_ace_file.ace','r'))

gives you 
        acefilerecord.ncontigs (the number of contigs in the ace file)
        acefilerecord.nreads (the number of reads in the ace file)
        acefilerecord.contigs[] (one instance of the Contig class for each contig)
        The Contig class holds the info of the CO tag, CT and WA tags, and all the reads used
        for this contig in a list of instances of the Read class, e.g.:
        contig3=acefilerecord.contigs[2]
        read4=contig3.reads[3]
        RD_of_read4=read4.rd
        DS_of_read4=read4.ds

        CT, WA, RT tags from the end of the file can appear anywhere are automatically
        sorted into the right place.

        see _RecordConsumer for details.

2) *** DEPRECATED: not entirely suitable for ACE files! 
        Or you can iterate over the contigs of an ace file one by one in the ususal way:        
        recordparser=Ace.RecordParser()
        iterator=Ace.Iterator(open('my_ace_file.ace','r'),recordparser)
        while 1:
            contig=iterator.next()
            if contig is None:
                break
            ...

        if WA, CT, RT, WR tags are at the end and the iterator is used, they will be returned
        with the last contig record. This is is necessarily the case when using an interator.
        Thus an ace file does not entirerly suit the concept of iterating. If WA, CT, RT, WR tags
        are needed, the ACEParser instead of the RecordParser might be appropriate.
        
"""
import os
from types import *

from Bio import File
from Bio import Index
#from Bio import Seq
#from Bio import SeqRecord
from Bio.ParserSupport import *
from Bio.Alphabet import IUPAC


class rd:
    def __init__(self):
        self.name=''
        self.padded_bases=None
        self.info_items=None
        self.read_tags=None
        self.sequence=''

class qa:
    def __init__(self):
        self.qual_clipping_start=None
        self.qual_clipping_end=None
        self.align_clipping_start=None
        self.align_clipping_end=None

class ds:
    def __init__(self):
        self.chromat_file=''
        self.phd_file=''
        self.time=''
        self.chem=''
        self.dye=''
        self.template=''
        self.direction=''
    
class af:
    def __init__(self):
        self.name=''
        self.coru=None
        self.padded_start=None

class bs:
    def __init__(self):
        self.name=''
        self.padded_start=None
        self.padded_end=None

class rt:
    def __init__(self):
        self.name=''
        self.tag_type=''
        self.program=''
        self.padded_start=None
        self.padded_end=None
        self.date=''

class ct:
    def __init__(self):
        self.name=''
        self.tag_type=''
        self.program=''
        self.padded_start=None
        self.padded_end=None
        self.date=''
        self.notrans=''
        self.info=[]

class wa:
    def __init__(self):
        self.tag_type=''
        self.program=''
        self.date=''
        self.info=[]

class wr:
    def __init__(self):
        self.name=''
        self.aligned=''
        self.program=''
        self.date=[]

class Reads:
    def __init__(self):
        self.rd=None    # one per read
        self.qa=None    # one per read
        self.ds=None    # none or one per read
        self.rt=None    # none or many per read
        self.wr=None    # none or many per read
        
class Contig:
    """Hold information from a ACE record """
    def __init__(self):
        self.name = ''
        self.nbases=None
        self.nreads=None
        self.nsegments=None
        self.uorc=None
        self.sequence=None
        self.quality=None
        self.af=[]
        self.bs=[]
        self.reads=[]
        self.ct=None    # none or many
        self.wa=None    # none or many
        
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
        self.ncontigs=None
        self.nreads=None
        self.contigs=[]
        self.wa=None    # none or many

    def sort(self):
        """Sorts wr, rt and ct tags into the appropriate contig / read instance, if possible.  """
       
        ct=[]
        rt=[]
        wr=[]
        # search for tags that aren't in the right position
        for i in range(len(self.contigs)):
            c = self.contigs[i]
            if c.wa:
                if not self.wa:
                    self.wa=[]
                self.wa.extend(c.wa)
            if c.ct:
                newcts=[ct_tag for ct_tag in c.ct if ct_tag.name!=c.name]
                map(self.contigs[i].ct.remove,newcts)
                ct.extend(newcts)
            for j in range(len(c.reads)):
                r = c.reads[j]
                if r.rt:
                    newrts=[rt_tag for rt_tag in r.rt if rt_tag.name!=r.rd.name]
                    map(self.contigs[i].reads[j].rt.remove,newrts)
                    rt.extend(newrts)
                if r.wr:
                    newwrs=[wr_tag for wr_tag in r.wr if wr_tag.name!=r.rd.name]
                    map(self.contigs[i].reads[j].wr.remove,newwrs)
                    wr.extend(newwrs)
        # now sort them into their proper place
        for i in range(len(self.contigs)):
            c = self.contigs[i]
            for ct_tag in ct:
                if ct_tag.name==c.name:
                    if self.contigs[i].ct is None:
                        self.contigs[i].ct=[]
                    self.contigs[i].ct.append(ct_tag)
            if rt or wr:
                for j in range(len(c.reads)):
                    r = c.reads[j]
                    for rt_tag in rt:
                        if rt_tag.name==r.rd.name:
                            if self.contigs[i].reads[j].rt is None:
                                self.contigs[i].reads[j].rt=[]
                            self.contigs[i].reads[j].rt.append(rt_tag)
                    for wr_tag in wr:
                        if wr_tag.name==r.rd.name:
                            if self.contigs[i].reads[j].wr is None:
                                self.contigs[i].reads[j].wr=[]
                            self.contigs[i].reads[j].wr.append(wr_tag)
       
class ACEParser(AbstractParser):
    """Parses full ACE file in list of contigs.

    """

    def __init__(self):
        self.data=ACEFileRecord()
        
    def parse(self,handle):
        firstline=handle.readline()
        # check if the file starts correctly
        if firstline[:2]!='AS':
            raise ValueError, "File does not start with 'AS'."
        self.data.ncontigs=eval(firstline.split()[1])
        self.data.nreads=eval(firstline.split()[2])
        # now read all the records
        recparser=RecordParser()
        iter=Iterator(handle,recparser)
        while 1:
            rec=iter.next()
            if not rec:
                break
            self.data.contigs.append(rec)
        # wa, ct, rt rags are usually at the end of the file, but not necessarily (correct?).
        # If the iterator is used, the tags are returned with the contig or the read after which they appear,
        # if all tags are at the end, they are read with the last contig. The concept of an
        # iterator leaves no other choice. But if the user uses the ACEParser, we can check
        # them and put them into the appropriate contig/read instance.
        # Conclusion: An ACE file is not a filetype for which iteration is 100% suitable...
        self.data.sort()
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
        consumer.co_data(self._scan_sequence_data(uhandle))
        read_and_call_while(uhandle,consumer.noevent,blank=1)
        read_and_call(uhandle,consumer.bq_header,start='BQ')
        consumer.bq_data(self._scan_bq_data(uhandle, consumer))
        read_and_call_while(uhandle,consumer.noevent,blank=1)
        read_and_call_while(uhandle,consumer.af,start='AF ')
        read_and_call_while(uhandle,consumer.noevent,blank=1)
        read_and_call_while(uhandle,consumer.bs,start='BS ')
        # now read all the read data
        # it starts with a 'RD', and then a mandatory QA
        # then follows an optional DS
        # CT,RT,WA,WR may or may not be there in unlimited quantity. They might refer to the actial read or contig,
        # or, if encountered at the end of file, to any previous read or contig. the sort() method deals
        # with that later.
        while 1:
            # each read must have a rd and qa
            read_and_call_until(uhandle,consumer.noevent,start='RD ')
            read_and_call(uhandle,consumer.rd_header,start='RD ')
            consumer.rd_data(self._scan_sequence_data(uhandle))
            read_and_call_while(uhandle,consumer.noevent,blank=1)
            read_and_call(uhandle,consumer.qa,start='QA ')
            # now one ds can follow
            try:
                read_and_call_while(uhandle,consumer.noevent,blank=1)
                attempt_read_and_call(uhandle,consumer.ds,start='DS ')
            except ValueError:
                # file ends
                consumer.end_contig()
                return
            # the file could just end, or there's some more stuff. In ace files, everything can happen.
            # the following tags are interspersed between reads and can ap[pear multiple times. 
            while 1:
                # something left 
                try:
                    read_and_call_while(uhandle,consumer.noevent,blank=1)
                except ValueError:
                    # file ends here
                    consumer.end_contig()
                    return
                else:
                    if attempt_read_and_call(uhandle,consumer.rt_start,start='RT'):
                        consumer.rt_data(self._scan_bracket_tags(uhandle))
                    elif attempt_read_and_call(uhandle,consumer.wr_start,start='WR'):
                        consumer.wr_data(self._scan_bracket_tags(uhandle))
                    elif attempt_read_and_call(uhandle,consumer.wa_start,start='WA'):
                        consumer.wa_data(self._scan_bracket_tags(uhandle))
                    elif attempt_read_and_call(uhandle,consumer.ct_start,start='CT'):
                        consumer.ct_data(self._scan_bracket_tags(uhandle))
                    else:
                        line=safe_peekline(uhandle)
                        break
            if not line.startswith('RD'): # another read?
                consumer.end_contig()
                break    
    
    def _scan_bq_data(self, uhandle, consumer):
        """Scans multiple lines of quality data and concatenates them."""
        
        qual=''
        while 1:
            line=uhandle.readline()
            if is_blank_line(line):
                uhandle.saveline(line)
                break
            qual+=' '+line
        return qual
   
    def _scan_sequence_data(self,uhandle):
        """Scans multiple lines of sequence data and concatenates them."""
        
        seq=''
        while 1:
            line=uhandle.readline()
            if is_blank_line(line):
                uhandle.saveline(line)
                break
            seq+=line.strip()
        return seq
     
    def _scan_bracket_tags(self,uhandle):
        """Reads the data lines of a {} tag."""
        
        fulltag=[]
        while 1:
            line=uhandle.readline().strip()
            fulltag.append(line)
            if line.endswith('}'):
                fulltag[-1]=fulltag[-1][:-1]    # delete the ending }
                if fulltag[-1]=='':
                    fulltag=fulltag[:-1]        # delete empty line
                break
        return fulltag
            
class _RecordConsumer(AbstractConsumer):
    """Reads the ace tags into data records."""
    
    def __init__(self):
        self.data = None

    def begin_contig(self):
        self.data = Contig()

    def end_contig(self):
        pass

    def co_header(self,line):
        header=line.split()
        self.data.name=header[1]
        self.data.nbases=eval(header[2])
        self.data.nreads=eval(header[3])
        self.data.nsegments=eval(header[4])
        self.data.uorc=header[5]

    def co_data(self,seq):
        self.data.sequence=seq
    
    def bq_header(self,line):
        pass
    
    def bq_data(self,qual):
        self.data.quality=map(eval,qual.split())

    def af(self,line):
        header=line.split()
        afdata=af()
        afdata.name=header[1]
        afdata.coru=header[2]
        afdata.padded_start=eval(header[3])
        self.data.af.append(afdata)
        
    def bs(self,line):
        header=line.split()
        bsdata=bs()
        bsdata.padded_start=eval(header[1])
        bsdata.padded_end=eval(header[2])
        bsdata.name=header[3]
        self.data.bs.append(bsdata)
        
    def rd_header(self,line):
        header=line.split()
        # Reads start with rd, so we create a new read record here
        self.data.reads.append(Reads())
        rddata=rd()
        rddata.name=header[1]
        rddata.padded_bases=eval(header[2])
        rddata.info_items=eval(header[3])
        rddata.read_tags=eval(header[4])
        self.data.reads[-1].rd=rddata

    def rd_data(self,seq):
        self.data.reads[-1].rd.sequence=seq
        
    def qa(self,line):
        header=map(eval,line.split()[1:])
        qadata=qa()
        qadata.qual_clipping_start=header[0]
        qadata.qual_clipping_end=header[1]
        qadata.align_clipping_start=header[2]
        qadata.align_clipping_end=header[3]
        self.data.reads[-1].qa=qadata

    def ds(self,line):
        dsdata=ds()
        tags=['CHROMAT_FILE','PHD_FILE','TIME','CHEM','DYE','TEMPLATE','DIRECTION']
        poss=map(line.find,tags)
        tagpos=dict(zip(poss,tags))
        if tagpos.has_key(-1):
            del tagpos[-1]
        ps=tagpos.keys()
        ps.sort()
        for (p1,p2) in zip(ps,ps[1:]+[len(line)+1]):
            setattr(dsdata,tagpos[p1].lower(),line[p1+len(tagpos[p1])+1:p2].strip())   
        self.data.reads[-1].ds=dsdata

    def ct_start(self,line):
        if not line.strip().endswith('{'):
            print line
            raise ValueError, 'CT tag does not start with CT{'
        ctdata=ct()
        if self.data.ct is None:
            self.data.ct=[]
        self.data.ct.append(ctdata)   
    
    def ct_data(self,taglines):
        if len(taglines)<1:
            raise ValueError, 'Missing header line in CT tag'
        header=taglines[0].split()
        self.data.ct[-1].name=header[0]
        self.data.ct[-1].tag_type=header[1]
        self.data.ct[-1].program=header[2]
        self.data.ct[-1].padded_start=eval(header[3])
        self.data.ct[-1].padded_end=eval(header[4])
        self.data.ct[-1].date=header[5]
        if len(header)==7:
            self.data.ct[-1].notrans=header[6]
        self.data.ct[-1].info=taglines[1:] 

    def rt_start(self,line):
        if not line.strip().endswith('{'):
            raise ValueError, 'RT tag does not start with RT{'
        rtdata=rt()
        # now if we're at the end of the file, this rt could belong to a previous read, not the actual one
        # we store it here were it appears, the user can sort later. 
        if self.data.reads[-1].rt is None:
            self.data.reads[-1].rt=[]
        self.data.reads[-1].rt.append(rtdata)
   
    def rt_data(self,taglines):
        if len(taglines)<1:
            raise ValueError, 'Missing header line in RT tag'
        header=taglines[0].split()
        self.data.reads[-1].rt[-1].name=header[0]
        self.data.reads[-1].rt[-1].tag_type=header[1]
        self.data.reads[-1].rt[-1].program=header[2]
        self.data.reads[-1].rt[-1].padded_start=eval(header[3])
        self.data.reads[-1].rt[-1].padded_end=eval(header[4])
        self.data.reads[-1].rt[-1].date=header[5]
    
    def wa_start(self,line):
        if not line.strip().endswith('{'):
            raise ValueError, 'WA tag does not start with WA{'
        wadata=wa()
        if self.data.wa is None:
            self.data.wa=[]
        self.data.wa.append(wadata)
    
    def wa_data(self,taglines):
        if len(taglines)<1:
            raise ValueError, 'Missing header line in WA tag'
        header=taglines[0].split()
        self.data.wa[-1].tag_type=header[0]
        self.data.wa[-1].program=header[1]
        self.data.wa[-1].date=header[2]
        self.data.wa[-1].info=taglines[1:] 
    
    def wr_start(self,line):
        if not line.strip().endswith('{'):
            raise ValueError, 'WR tag does not start with WR{'
        wrdata=wr()
        if self.data.reads[-1].wr is None:
            self.data.reads[-1].wr=[]
        self.data.reads[-1].wr.append(wrdata)
    
    def wr_data(self,taglines):
        if len(taglines)<1:
            raise ValueError, 'Missing header line in WR tag'
        header=taglines[0].split()
        self.data.reads[-1].wr[-1].name=header[0]
        self.data.reads[-1].wr[-1].aligned=header[1]
        self.data.reads[-1].wr[-1].program=header[2]
        self.data.reads[-1].wr[-1].date=header[3]
    
    
