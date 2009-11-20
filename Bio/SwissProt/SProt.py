# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with the sprotXX.dat file from
Utilities for working with FASTA-formatted sequences (DEPRECATED).
http://www.expasy.ch/sprot/sprot-top.html

Please see Bio.SwissProt for alternatives for the functionality in this module.

Tested with:
Release 37, Release 38, Release 39

Limited testing with:
Release 51, 54


Classes:
Record             Holds SwissProt data.
Reference          Holds reference data from a SwissProt entry.
Iterator           Iterates over entries in a SwissProt file.
Dictionary         Accesses a SwissProt file using a dictionary interface.
RecordParser       Parses a SwissProt record into a Record object.
SequenceParser     Parses a SwissProt record into a SeqRecord object.

_Scanner           Scans SwissProt-formatted data.
_RecordConsumer    Consumes SwissProt data to a SProt.Record object.
_SequenceConsumer  Consumes SwissProt data to a SeqRecord object.


Functions:
index_file         Index a SwissProt file for a Dictionary.

"""
import warnings
warnings.warn("Bio.SwissProt.SProt is deprecated. Please use the functions Bio.SwissProt.parse or Bio.SwissProt.read if you want to get a SwissProt.Record, or Bio.SeqIO.parse or Bio.SeqIO.read if you want to get a SeqRecord. If these solutions do not work for you, please get in contact with the Biopython developers (biopython-dev@biopython.org).",
              DeprecationWarning)

from types import *
import os
from Bio import File
from Bio import Index
from Bio import Alphabet
from Bio import Seq
from Bio import SeqRecord
from Bio.ParserSupport import *

# The parse(), read() functions can probably be simplified if we don't
# use the "parser = RecordParser(); parser.parse(handle)" approach.
def parse(handle):
    from SProt import RecordParser
    import cStringIO
    parser = RecordParser()
    text = ""
    for line in handle:
        text += line
        if line[:2]=='//':
            handle = cStringIO.StringIO(text)
            record = parser.parse(handle)
            text = ""
            yield record

def read(handle):
    from SProt import RecordParser
    parser = RecordParser()
    try:
        record = parser.parse(handle)
    except ValueError, error:
        if error.message.startswith("Line does not start with 'ID':"): 
            raise ValueError("No SwissProt record found")
        else:
            raise error
    # We should have reached the end of the record by now
    remainder = handle.read()
    if remainder:
        raise ValueError("More than one SwissProt record found")
    return record


_CHOMP = " \n\r\t.,;" #whitespace and trailing punctuation

class Record:
    """Holds information from a SwissProt record.

    Members:
    entry_name        Name of this entry, e.g. RL1_ECOLI.
    data_class        Either 'STANDARD' or 'PRELIMINARY'.
    molecule_type     Type of molecule, 'PRT',
    sequence_length   Number of residues.

    accessions        List of the accession numbers, e.g. ['P00321']
    created           A tuple of (date, release).
    sequence_update   A tuple of (date, release).
    annotation_update A tuple of (date, release).

    description       Free-format description.
    gene_name         Gene name.  See userman.txt for description.
    organism          The source of the sequence.
    organelle         The origin of the sequence.
    organism_classification  The taxonomy classification.  List of strings.
                             (http://www.ncbi.nlm.nih.gov/Taxonomy/)
    taxonomy_id       A list of NCBI taxonomy id's.
    host_organism     A list of NCBI taxonomy id's of the hosts of a virus,
                      if any.
    references        List of Reference objects.
    comments          List of strings.
    cross_references  List of tuples (db, id1[, id2][, id3]).  See the docs.
    keywords          List of the keywords.
    features          List of tuples (key name, from, to, description).
                      from and to can be either integers for the residue
                      numbers, '<', '>', or '?'

    seqinfo           tuple of (length, molecular weight, CRC32 value)
    sequence          The sequence.
    
    """
    def __init__(self):
        self.entry_name = None
        self.data_class = None
        self.molecule_type = None
        self.sequence_length = None
        
        self.accessions = []
        self.created = None
        self.sequence_update = None
        self.annotation_update = None
        
        self.description = ''
        self.gene_name = ''
        self.organism = ''
        self.organelle = ''
        self.organism_classification = []
        self.taxonomy_id = []
        self.host_organism = []
        self.references = []
        self.comments = []
        self.cross_references = []
        self.keywords = []
        self.features = []
        
        self.seqinfo = None
        self.sequence = ''

class Reference:
    """Holds information from 1 references in a SwissProt entry.

    Members:
    number      Number of reference in an entry.
    positions   Describes extent of work.  list of strings.
    comments    Comments.  List of (token, text).
    references  References.  List of (dbname, identifier)
    authors     The authors of the work.
    title       Title of the work.
    location    A citation for the work.
    
    """
    def __init__(self):
        self.number = None
        self.positions = []
        self.comments = []
        self.references = []
        self.authors = ''
        self.title = ''
        self.location = ''


class Dictionary:
    """Accesses a SwissProt file using a dictionary interface.

    """
    __filename_key = '__filename'
    
    def __init__(self, indexname, parser=None):
        """__init__(self, indexname, parser=None)

        Open a SwissProt Dictionary.  indexname is the name of the
        index for the dictionary.  The index should have been created
        using the index_file function.  parser is an optional Parser
        object to change the results into another form.  If set to None,
        then the raw contents of the file will be returned.

        """
        self._index = Index.Index(indexname)
        self._handle = open(self._index[self.__filename_key])
        self._parser = parser

    def __len__(self):
        return len(self._index)

    def __getitem__(self, key):
        start, len = self._index[key]
        self._handle.seek(start)
        data = self._handle.read(len)
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

    def __getattr__(self, name):
        return getattr(self._index, name)

    def keys(self):
        # I only want to expose the keys for SwissProt.
        k = self._index.keys()
        k.remove(self.__filename_key)
        return k


class RecordParser(AbstractParser):
    """Parses SwissProt data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class SequenceParser(AbstractParser):
    """Parses SwissProt data into a standard SeqRecord object.
    """
    def __init__(self, alphabet = Alphabet.generic_protein):
        """Initialize a SequenceParser.

        Arguments:
        o alphabet - The alphabet to use for the generated Seq objects. If
        not supplied this will default to the generic protein alphabet.
        """
        self._scanner = _Scanner()
        self._consumer = _SequenceConsumer(alphabet)

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans SwissProt-formatted data.

    Tested with:
    Release 37
    Release 38
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in SwissProt data for scanning.  handle is a file-like
        object that contains swissprot data.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        self._scan_record(uhandle, consumer)
        
    def _skip_starstar(self, uhandle):
        """Ignores any lines starting **"""
        #See Bug 2353, some files from the EBI have extra lines
        #starting "**" (two asterisks/stars), usually between the
        #features and sequence but not all the time.  They appear
        #to be unofficial automated annotations. e.g.
        #**
        #**   #################    INTERNAL SECTION    ##################
        #**HA SAM; Annotated by PicoHamap 1.88; MF_01138.1; 09-NOV-2003.
        while "**" == uhandle.peekline()[:2]:
            skip = uhandle.readline()
            #print "Skipping line: %s" % skip.rstrip()

    def _scan_record(self, uhandle, consumer):
        consumer.start_record()
        for fn in self._scan_fns:
            self._skip_starstar(uhandle)
            fn(self, uhandle, consumer)

            # In Release 38, ID N33_HUMAN has a DR buried within comments.
            # Check for this and do more comments, if necessary.
            # XXX handle this better
            if fn is self._scan_dr.im_func:
                self._scan_cc(uhandle, consumer)
                self._scan_dr(uhandle, consumer)
        consumer.end_record()

    def _scan_line(self, line_type, uhandle, event_fn,
                   exactly_one=None, one_or_more=None, any_number=None,
                   up_to_one=None):
        # Callers must set exactly one of exactly_one, one_or_more, or
        # any_number to a true value.  I do not explicitly check to
        # make sure this function is called correctly.
        
        # This does not guarantee any parameter safety, but I
        # like the readability.  The other strategy I tried was have
        # parameters min_lines, max_lines.
        
        if exactly_one or one_or_more:
            read_and_call(uhandle, event_fn, start=line_type)
        if one_or_more or any_number:
            while 1:
                if not attempt_read_and_call(uhandle, event_fn,
                                             start=line_type):
                    break
        if up_to_one:
            attempt_read_and_call(uhandle, event_fn, start=line_type)

    def _scan_id(self, uhandle, consumer):
        self._scan_line('ID', uhandle, consumer.identification, exactly_one=1)

    def _scan_ac(self, uhandle, consumer):
        # Until release 38, this used to match exactly_one.
        # However, in release 39, 1A02_HUMAN has 2 AC lines, and the
        # definition needed to be expanded.
        self._scan_line('AC', uhandle, consumer.accession, any_number=1)
    
    def _scan_dt(self, uhandle, consumer):
        self._scan_line('DT', uhandle, consumer.date, exactly_one=1)
        self._scan_line('DT', uhandle, consumer.date, exactly_one=1)
        # IPI doesn't necessarily contain the third line about annotations
        self._scan_line('DT', uhandle, consumer.date, up_to_one=1)

    def _scan_de(self, uhandle, consumer):
        # IPI can be missing a DE line
        self._scan_line('DE', uhandle, consumer.description, any_number=1)
    
    def _scan_gn(self, uhandle, consumer):
        self._scan_line('GN', uhandle, consumer.gene_name, any_number=1)
    
    def _scan_os(self, uhandle, consumer):
        self._scan_line('OS', uhandle, consumer.organism_species,
                        one_or_more=1)
    
    def _scan_og(self, uhandle, consumer):
        self._scan_line('OG', uhandle, consumer.organelle, any_number=1)
    
    def _scan_oc(self, uhandle, consumer):
        self._scan_line('OC', uhandle, consumer.organism_classification,
                        one_or_more=1)

    def _scan_ox(self, uhandle, consumer):
        self._scan_line('OX', uhandle, consumer.taxonomy_id,
                        any_number=1)

    def _scan_oh(self, uhandle, consumer):
        # viral host organism. introduced after SwissProt 39.
        self._scan_line('OH', uhandle, consumer.organism_host, any_number=1)

    def _scan_reference(self, uhandle, consumer):
        while True:
            if safe_peekline(uhandle)[:2] != 'RN':
                break
            self._scan_rn(uhandle, consumer)
            self._scan_rp(uhandle, consumer)
            self._scan_rc(uhandle, consumer)
            self._scan_rx(uhandle, consumer)
            # ws:2001-12-05 added, for record with RL before RA
            self._scan_rl(uhandle, consumer)
            self._scan_ra(uhandle, consumer)
            #EBI copy of P72010 is missing the RT line, and has one
            #of their ** lines in its place noting "**   /NO TITLE."
            #See also bug 2353
            self._skip_starstar(uhandle) 
            self._scan_rt(uhandle, consumer)
            self._scan_rl(uhandle, consumer)
    
    def _scan_rn(self, uhandle, consumer):
        self._scan_line('RN', uhandle, consumer.reference_number,
                        exactly_one=1)
    
    def _scan_rp(self, uhandle, consumer):
        self._scan_line('RP', uhandle, consumer.reference_position,
                        one_or_more=1)
    
    def _scan_rc(self, uhandle, consumer):
        self._scan_line('RC', uhandle, consumer.reference_comment,
                        any_number=1)
    
    def _scan_rx(self, uhandle, consumer):
        self._scan_line('RX', uhandle, consumer.reference_cross_reference,
                        any_number=1)
    
    def _scan_ra(self, uhandle, consumer):
        # In UniProt release 1.12 of 6/21/04, there is a new RG
        # (Reference Group) line, which references a group instead of
        # an author.  Each block must have at least 1 RA or RG line.
        self._scan_line('RA', uhandle, consumer.reference_author,
                        any_number=1)
        self._scan_line('RG', uhandle, consumer.reference_author,
                        any_number=1)
        # PRKN_HUMAN has RG lines, then RA lines.  The best solution
        # is to write code that accepts either of the line types.
        # This is the quick solution...
        self._scan_line('RA', uhandle, consumer.reference_author,
                        any_number=1)
    
    def _scan_rt(self, uhandle, consumer):
        self._scan_line('RT', uhandle, consumer.reference_title,
                        any_number=1)
    
    def _scan_rl(self, uhandle, consumer):
        # This was one_or_more, but P82909 in TrEMBL 16.0 does not
        # have one.
        self._scan_line('RL', uhandle, consumer.reference_location,
                        any_number=1)
    
    def _scan_cc(self, uhandle, consumer):
        self._scan_line('CC', uhandle, consumer.comment, any_number=1)
    
    def _scan_dr(self, uhandle, consumer):
        self._scan_line('DR', uhandle, consumer.database_cross_reference,
                        any_number=1)
    
    def _scan_kw(self, uhandle, consumer):
        self._scan_line('KW', uhandle, consumer.keyword, any_number=1)
    
    def _scan_ft(self, uhandle, consumer):
        self._scan_line('FT', uhandle, consumer.feature_table, any_number=1)

    def _scan_pe(self, uhandle, consumer):
        self._scan_line('PE', uhandle, consumer.protein_existence, any_number=1)

    def _scan_sq(self, uhandle, consumer):
        self._scan_line('SQ', uhandle, consumer.sequence_header, exactly_one=1)
    
    def _scan_sequence_data(self, uhandle, consumer):
        self._scan_line('  ', uhandle, consumer.sequence_data, one_or_more=1)
    
    def _scan_terminator(self, uhandle, consumer):
        self._scan_line('//', uhandle, consumer.terminator, exactly_one=1)

    _scan_fns = [
        _scan_id,
        _scan_ac,
        _scan_dt,
        _scan_de,
        _scan_gn,
        _scan_os,
        _scan_og,
        _scan_oc,
        _scan_ox,
        _scan_oh,
        _scan_reference,
        _scan_cc,
        _scan_dr,
        _scan_pe,
        _scan_kw,
        _scan_ft,
        _scan_sq,
        _scan_sequence_data,
        _scan_terminator
        ]

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a SwissProt record to a Record object.

    Members:
    data    Record with SwissProt data.

    """
    def __init__(self):
        self.data = None
        
    def __repr__(self):
        return "Bio.SwissProt.SProt._RecordConsumer()"
        
    def start_record(self):
        self.data = Record()
        self._sequence_lines = []
        
    def end_record(self):
        self._clean_record(self.data)
        self.data.sequence = "".join(self._sequence_lines)

    def identification(self, line):
        cols = line.split()
        #Prior to release 51, included with MoleculeType:
        #ID   EntryName DataClass; MoleculeType; SequenceLength.
        #
        #Newer files lack the MoleculeType:
        #ID   EntryName DataClass; SequenceLength.
        #
        #Note that cols is split on white space, so the length
        #should become two fields (number and units)
        if len(cols) == 6:
            self.data.entry_name = cols[1]
            self.data.data_class = cols[2].rstrip(_CHOMP)    # don't want ';'
            self.data.molecule_type = cols[3].rstrip(_CHOMP) # don't want ';'
            self.data.sequence_length = int(cols[4])
        elif len(cols) == 5:
            self.data.entry_name = cols[1]
            self.data.data_class = cols[2].rstrip(_CHOMP)    # don't want ';'
            self.data.molecule_type = None
            self.data.sequence_length = int(cols[3])
        else:
            #Should we print a warning an continue?
            raise ValueError("ID line has unrecognised format:\n"+line)
        
        # data class can be 'STANDARD' or 'PRELIMINARY'
        # ws:2001-12-05 added IPI
        # pjc:2006-11-02 added 'Reviewed' and 'Unreviewed'
        if self.data.data_class not in ['STANDARD', 'PRELIMINARY', 'IPI',
                                        'Reviewed', 'Unreviewed']: 
            raise ValueError("Unrecognized data class %s in line\n%s" % \
                  (self.data.data_class, line))
        # molecule_type should be 'PRT' for PRoTein
        # Note that has been removed in recent releases (set to None)
        if self.data.molecule_type is not None \
        and self.data.molecule_type != 'PRT':
            raise ValueError("Unrecognized molecule type %s in line\n%s" % \
                  (self.data.molecule_type, line))
    
    def accession(self, line):
        cols = line[5:].rstrip(_CHOMP).strip().split(';')
        for ac in cols:
            if ac.strip():
                #remove any leading or trailing white space
                self.data.accessions.append(ac.strip())
    
    def date(self, line):
        uprline = line.upper()
        cols = line.rstrip().split()
                
        if uprline.find('CREATED') >= 0 \
        or uprline.find('LAST SEQUENCE UPDATE') >= 0 \
        or uprline.find('LAST ANNOTATION UPDATE') >= 0:
            # Old style DT line
            # =================
            # e.g.
            # DT   01-FEB-1995 (Rel. 31, Created)
            # DT   01-FEB-1995 (Rel. 31, Last sequence update)
            # DT   01-OCT-2000 (Rel. 40, Last annotation update)
            #
            # or:
            # DT   08-JAN-2002 (IPI Human rel. 2.3, Created)
            # ... 
            
            # find where the version information will be located
            # This is needed for when you have cases like IPI where
            # the release verison is in a different spot:
            # DT   08-JAN-2002 (IPI Human rel. 2.3, Created)
            uprcols = uprline.split()
            rel_index = -1
            for index in range(len(uprcols)):
                if uprcols[index].find("REL.") >= 0:
                    rel_index = index
            assert rel_index >= 0, \
                    "Could not find Rel. in DT line: %s" % (line)
            version_index = rel_index + 1
            # get the version information
            str_version = cols[version_index].rstrip(_CHOMP)
            # no version number
            if str_version == '':
                version = 0
            # dot versioned
            elif str_version.find(".") >= 0:
                version = str_version
            # integer versioned
            else:
                version = int(str_version)

            if uprline.find('CREATED') >= 0:
                self.data.created = cols[1], version
            elif uprline.find('LAST SEQUENCE UPDATE') >= 0:
                self.data.sequence_update = cols[1], version
            elif uprline.find( 'LAST ANNOTATION UPDATE') >= 0:
                self.data.annotation_update = cols[1], version
            else:
                assert False, "Shouldn't reach this line!"
        elif uprline.find('INTEGRATED INTO') >= 0 \
        or uprline.find('SEQUENCE VERSION') >= 0 \
        or uprline.find('ENTRY VERSION') >= 0:
            # New style DT line
            # =================
            # As of UniProt Knowledgebase release 7.0 (including
            # Swiss-Prot release 49.0 and TrEMBL release 32.0) the
            # format of the DT lines and the version information
            # in them was changed - the release number was dropped.
            #
            # For more information see bug 1948 and
            # http://ca.expasy.org/sprot/relnotes/sp_news.html#rel7.0
            #
            # e.g.
            # DT   01-JAN-1998, integrated into UniProtKB/Swiss-Prot.
            # DT   15-OCT-2001, sequence version 3.
            # DT   01-APR-2004, entry version 14.
            #
            #This is a new style DT line...

            # The date should be in string cols[1]
            # Get the version number if there is one.
            # For the three DT lines above: 0, 3, 14
            try:
                version = int(cols[-1])
            except ValueError:
                version = 0

            # Re-use the historical property names, even though
            # the meaning has changed slighty:
            if uprline.find("INTEGRATED") >= 0:
                self.data.created = cols[1], version
            elif uprline.find('SEQUENCE VERSION') >= 0:
                self.data.sequence_update = cols[1], version
            elif uprline.find( 'ENTRY VERSION') >= 0:
                self.data.annotation_update = cols[1], version
            else:
                assert False, "Shouldn't reach this line!"
        else:
            raise ValueError("I don't understand the date line %s" % line)
    
    def description(self, line):
        self.data.description += line[5:].strip() + " "
    
    def gene_name(self, line):
        self.data.gene_name += line[5:]
    
    def organism_species(self, line):
        self.data.organism += line[5:]
    
    def organelle(self, line):
        self.data.organelle += line[5:]
    
    def organism_classification(self, line):
        line = line[5:].rstrip(_CHOMP)
        cols = line.split(';')
        for col in cols:
            self.data.organism_classification.append(col.lstrip())

    def taxonomy_id(self, line):
        # The OX line is in the format:
        # OX   DESCRIPTION=ID[, ID]...;
        # If there are too many id's to fit onto a line, then the ID's
        # continue directly onto the next line, e.g.
        # OX   DESCRIPTION=ID[, ID]...
        # OX   ID[, ID]...;
        # Currently, the description is always "NCBI_TaxID".

        # To parse this, I need to check to see whether I'm at the
        # first line.  If I am, grab the description and make sure
        # it's an NCBI ID.  Then, grab all the id's.
        line = line[5:].rstrip(_CHOMP)
        index = line.find('=')
        if index >= 0:
            descr = line[:index]
            assert descr == "NCBI_TaxID", "Unexpected taxonomy type %s" % descr
            ids = line[index+1:].split(',')
        else:
            ids = line.split(',')
        self.data.taxonomy_id.extend([id.strip() for id in ids])

    def organism_host(self, line):
        # Line type OH (Organism Host) for viral hosts
        # same code as in taxonomy_id()
        line = line[5:].rstrip(_CHOMP)
        index = line.find('=')
        if index >= 0:
            descr = line[:index]
            assert descr == "NCBI_TaxID", "Unexpected taxonomy type %s" % descr
            ids = line[index+1:].split(',')
        else:
            ids = line.split(',')
        self.data.host_organism.extend([id.strip() for id in ids])
        
    def reference_number(self, line):
        rn = line[5:].rstrip()
        assert rn[0] == '[' and rn[-1] == ']', "Missing brackets %s" % rn
        ref = Reference()
        ref.number = int(rn[1:-1])
        self.data.references.append(ref)
    
    def reference_position(self, line):
        assert self.data.references, "RP: missing RN"
        self.data.references[-1].positions.append(line[5:].rstrip())
    
    def reference_comment(self, line):
        assert self.data.references, "RC: missing RN"
        cols = line[5:].rstrip().split( ';')
        ref = self.data.references[-1]
        for col in cols:
            if not col:  # last column will be the empty string
                continue
            # The token is everything before the first '=' character.
            index = col.find('=')
            token, text = col[:index], col[index+1:]
            # According to the spec, there should only be 1 '='
            # character.  However, there are too many exceptions to
            # handle, so we'll ease up and allow anything after the
            # first '='.
            #if col == ' STRAIN=TISSUE=BRAIN':
            #    # from CSP_MOUSE, release 38
            #    token, text = "TISSUE", "BRAIN"
            #elif col == ' STRAIN=NCIB 9816-4, AND STRAIN=G7 / ATCC 17485':
            #    # from NDOA_PSEPU, release 38
            #    token, text = "STRAIN", "NCIB 9816-4 AND G7 / ATCC 17485"
            #elif col == ' STRAIN=ISOLATE=NO 27, ANNO 1987' or \
            #     col == ' STRAIN=ISOLATE=NO 27 / ANNO 1987':
            #    # from NU3M_BALPH, release 38, release 39
            #    token, text = "STRAIN", "ISOLATE NO 27, ANNO 1987"
            #else:
            #    token, text = string.split(col, '=')
            ref.comments.append((token.lstrip(), text))
    
    def reference_cross_reference(self, line):
        assert self.data.references, "RX: missing RN"
        # The basic (older?) RX line is of the form:
        # RX   MEDLINE; 85132727.
        # but there are variants of this that need to be dealt with (see below)
        
        # CLD1_HUMAN in Release 39 and DADR_DIDMA in Release 33
        # have extraneous information in the RX line.  Check for
        # this and chop it out of the line.
        # (noticed by katel@worldpath.net)
        ind = line.find('[NCBI, ExPASy, Israel, Japan]')
        if ind >= 0:
            line = line[:ind]

        # RX lines can also be used of the form
        # RX   PubMed=9603189;
        # reported by edvard@farmasi.uit.no
        # and these can be more complicated like:
        # RX   MEDLINE=95385798; PubMed=7656980;
        # RX   PubMed=15060122; DOI=10.1136/jmg 2003.012781;
        # We look for these cases first and deal with them
        if line.find("=") != -1:
            cols = line[2:].split("; ")
            cols = [x.strip() for x in cols]
            cols = [x for x in cols if x]
            for col in cols:
                x = col.split("=")
                assert len(x) == 2, "I don't understand RX line %s" % line
                key, value = x[0].rstrip(_CHOMP), x[1].rstrip(_CHOMP)
                ref = self.data.references[-1].references
                ref.append((key, value))
        # otherwise we assume we have the type 'RX   MEDLINE; 85132727.'
        else:
            cols = line.split()
            # normally we split into the three parts
            assert len(cols) == 3, "I don't understand RX line %s" % line
            self.data.references[-1].references.append(
                (cols[1].rstrip(_CHOMP), cols[2].rstrip(_CHOMP)))
    
    def reference_author(self, line):
        assert self.data.references, "RA: missing RN"
        ref = self.data.references[-1]
        ref.authors += line[5:]
    
    def reference_title(self, line):
        assert self.data.references, "RT: missing RN"
        ref = self.data.references[-1]
        ref.title += line[5:]
    
    def reference_location(self, line):
        assert self.data.references, "RL: missing RN"
        ref = self.data.references[-1]
        ref.location += line[5:]
    
    def comment(self, line):
        if line[5:8] == '-!-':   # Make a new comment
            self.data.comments.append(line[9:])
        elif line[5:8] == '   ': # add to the previous comment
            if not self.data.comments:
                # TCMO_STRGA in Release 37 has comment with no topic
                self.data.comments.append(line[9:])
            else:
                self.data.comments[-1] += line[9:]
        elif line[5:8] == '---':
            # If there are no comments, and it's not the closing line,
            # make a new comment.
            if not self.data.comments or self.data.comments[-1][:3] != '---':
                self.data.comments.append(line[5:])
            else:
                self.data.comments[-1] += line[5:]
        else:  # copyright notice
            self.data.comments[-1] += line[5:]
    
    def database_cross_reference(self, line):
        # From CLD1_HUMAN, Release 39:
        # DR   EMBL; [snip]; -. [EMBL / GenBank / DDBJ] [CoDingSequence]
        # DR   PRODOM [Domain structure / List of seq. sharing at least 1 domai
        # DR   SWISS-2DPAGE; GET REGION ON 2D PAGE.
        line = line[5:]
        # Remove the comments at the end of the line
        i = line.find('[')
        if i >= 0:
            line = line[:i]
        cols = line.rstrip(_CHOMP).split(';')
        cols = [col.lstrip() for col in cols]
        self.data.cross_references.append(tuple(cols))
    
    def keyword(self, line):
        cols = line[5:].rstrip(_CHOMP).split(';')
        self.data.keywords.extend([c.lstrip() for c in cols])

    def feature_table(self, line):
        line = line[5:]    # get rid of junk in front
        name = line[0:8].rstrip()
        try:
            from_res = int(line[9:15])
        except ValueError:
            from_res = line[9:15].lstrip()
        try:
            to_res = int(line[16:22])
        except ValueError:
            to_res = line[16:22].lstrip()
        description = line[29:70].rstrip()
        #if there is a feature_id (FTId), store it away
        if line[29:35]==r"/FTId=":
            ft_id = line[35:70].rstrip()[:-1]
        else:
            ft_id =""
        if not name:  # is continuation of last one
            assert not from_res and not to_res
            name, from_res, to_res, old_description,old_ft_id = self.data.features[-1]
            del self.data.features[-1]
            description = "%s %s" % (old_description, description)
            
            # special case -- VARSPLIC, reported by edvard@farmasi.uit.no
            if name == "VARSPLIC":
                description = self._fix_varsplic_sequences(description)
        self.data.features.append((name, from_res, to_res, description,ft_id))

    def _fix_varsplic_sequences(self, description):
        """Remove unwanted spaces in sequences.

        During line carryover, the sequences in VARSPLIC can get mangled
        with unwanted spaces like:
        'DISSTKLQALPSHGLESIQT -> PCRATGWSPFRRSSPC LPTH'
        We want to check for this case and correct it as it happens.
        """
        descr_cols = description.split(" -> ")
        if len(descr_cols) == 2:
            first_seq = descr_cols[0]
            second_seq = descr_cols[1]
            extra_info = ''
            # we might have more information at the end of the
            # second sequence, which should be in parenthesis
            extra_info_pos = second_seq.find(" (")
            if extra_info_pos != -1:
                extra_info = second_seq[extra_info_pos:]
                second_seq = second_seq[:extra_info_pos]

            # now clean spaces out of the first and second string
            first_seq = first_seq.replace(" ", "")
            second_seq = second_seq.replace(" ", "")

            # reassemble the description
            description = first_seq + " -> " + second_seq + extra_info

        return description
    
    def protein_existence(self, line):
        #TODO - Record this information?
        pass
    
    def sequence_header(self, line):
        cols = line.split()
        assert len(cols) == 8, "I don't understand SQ line %s" % line
        # Do more checking here?
        self.data.seqinfo = int(cols[2]), int(cols[4]), cols[6]
    
    def sequence_data(self, line):
        #It should be faster to make a list of strings, and join them at the end.
        self._sequence_lines.append(line.replace(" ", "").rstrip())
    
    def terminator(self, line):
        pass

    #def _clean(self, line, rstrip=1):
    #    if rstrip:
    #        return string.rstrip(line[5:])
    #    return line[5:]

    def _clean_record(self, rec):
        # Remove trailing newlines
        members = ['description', 'gene_name', 'organism', 'organelle']
        for m in members:
            attr = getattr(rec, m)
            setattr(rec, m, attr.rstrip())
        for ref in rec.references:
            self._clean_references(ref)

    def _clean_references(self, ref):
        # Remove trailing newlines
        members = ['authors', 'title', 'location']
        for m in members:
            attr = getattr(ref, m)
            setattr(ref, m, attr.rstrip())

class _SequenceConsumer(AbstractConsumer):
    """Consumer that converts a SwissProt record to a SeqRecord object.

    Members:
    data      Record with SwissProt data.
    alphabet  The alphabet the generated Seq objects will have.
    """
    #TODO - Cope with references as done for GenBank
    def __init__(self, alphabet = Alphabet.generic_protein):
        """Initialize a Sequence Consumer

        Arguments:
        o alphabet - The alphabet to use for the generated Seq objects. If
        not supplied this will default to the generic protein alphabet.
        """
        self.data = None
        self.alphabet = alphabet
        
    def start_record(self):
        seq = Seq.Seq("", self.alphabet)
        self.data = SeqRecord.SeqRecord(seq)
        self.data.description = ""
        self.data.name = ""
        self._current_ref = None
        self._sequence_lines = []
        
    def end_record(self):
        if self._current_ref is not None:
            self.data.annotations['references'].append(self._current_ref)
            self._current_ref = None
        self.data.description = self.data.description.rstrip()
        self.data.seq = Seq.Seq("".join(self._sequence_lines), self.alphabet)
        self.data.annotations['organism'] = self.data.annotations['organism'].rstrip(_CHOMP)

    def identification(self, line):
        cols = line.split()
        self.data.name = cols[1]

    def accession(self, line):
        #Note that files can and often do contain multiple AC lines.
        ids = line[5:].strip().split(';')
        #Remove any white space
        ids = [x.strip() for x in ids if x.strip()]
        
        #Use the first as the ID, but record them ALL in the annotations
        try:
            self.data.annotations['accessions'].extend(ids)
        except KeyError:
            self.data.annotations['accessions'] = ids
            
        #Use the FIRST accession as the ID, not the first on this line!
        self.data.id = self.data.annotations['accessions'][0]
        #self.data.id = ids[0]

    def description(self, line):
        self.data.description += line[5:].strip() + " "
        
    def sequence_data(self, line):
        #It should be faster to make a list of strings, and join them at the end.
        self._sequence_lines.append(line.replace(" ", "").rstrip())

    def gene_name(self, line):
        #We already store the identification/accession as the records name/id
        try:
            self.data.annotations['gene_name'] += " " + line[5:].rstrip()
        except KeyError:
            self.data.annotations['gene_name'] =  line[5:].rstrip()

    def comment(self, line):
        #Try and agree with SeqRecord convention from the GenBank parser,
        #which stores the comments as a long string with newlines
        #with key 'comment'
        #TODO - Follow SwissProt conventions more closely?
        prefix = line[5:8]
        text = line[9:].rstrip()
        if prefix == '-!-':   # Make a new comment
            try:
                self.data.annotations['comment'] += "\n" + text
            except KeyError:
                self.data.annotations['comment'] =  text
        elif prefix == '   ':
            try:
                # add to the previous comment
                self.data.annotations['comment'] += " " + text
            except KeyError:
                # TCMO_STRGA in Release 37 has comment with no topic
                self.data.annotations['comment'] =  text

    def database_cross_reference(self, line):
        #Format of the line is described in the manual dated 04-Dec-2007 as:
        #DR   DATABASE; PRIMARY; SECONDARY[; TERTIARY][; QUATERNARY].
        #However, some older files only seem to have a single identifier:
        #DR   DATABASE; PRIMARY.
        #
        #Also must cope with things like this from Tests/SwissProt/sp007,
        #DR   PRODOM [Domain structure / List of seq. sharing at least 1 domain]
        #
        #Store these in the dbxref list, but for consistency with
        #the GenBank parser and with what BioSQL can cope with,
        #store only DATABASE_IDENTIFIER:PRIMARY_IDENTIFIER
        parts = [x.strip() for x in line[5:].strip(_CHOMP).split(";")]
        if len(parts) > 1:
            value = "%s:%s" % (parts[0], parts[1])
            #Avoid duplicate entries
            if value not in self.data.dbxrefs:
                self.data.dbxrefs.append(value)
        #else:
            #print "Bad DR line:\n%s" % line
            
                    
    def date(self, line):
        date_str = line.split()[1]
        uprline = line.upper()
        if uprline.find('CREATED') >= 0:
            #Try and agree with SeqRecord convention from the GenBank parser,
            #which stores the submitted date as 'date'
            self.data.annotations['date'] = date_str
        elif uprline.find('LAST SEQUENCE UPDATE') >= 0:
            #There is no existing convention from the GenBank SeqRecord parser
            self.data.annotations['date_last_sequence_update'] = date_str
        elif uprline.find('LAST ANNOTATION UPDATE') >= 0:
            #There is no existing convention from the GenBank SeqRecord parser
            self.data.annotations['date_last_annotation_update'] = date_str
        elif uprline.find('INTEGRATED INTO') >= 0:
            self.data.annotations['date'] = date_str.rstrip(",")
        elif uprline.find('SEQUENCE VERSION') >= 0:
            self.data.annotations['date_last_sequence_update'] = date_str.rstrip(",")
        elif uprline.find('ENTRY VERSION') >= 0:
            self.data.annotations['date_last_annotation_update'] = date_str.rstrip(",")

    def keyword(self, line):
        #Try and agree with SeqRecord convention from the GenBank parser,
        #which stores a list as 'keywords'
        cols = line[5:].rstrip(_CHOMP).split(';')
        cols = [c.strip() for c in cols]
        cols = filter(None, cols)
        try:
            #Extend any existing list of keywords
            self.data.annotations['keywords'].extend(cols)
        except KeyError:
            #Create the list of keywords
            self.data.annotations['keywords'] = cols

    def organism_species(self, line):
        #Try and agree with SeqRecord convention from the GenBank parser,
        #which stores the organism as a string with key 'organism'
        data = line[5:].rstrip()
        try:
            #Append to any existing data split over multiple lines
            self.data.annotations['organism'] += " " + data
        except KeyError:
            self.data.annotations['organism'] = data

    def organism_host(self, line):
        #There is no SeqRecord convention from the GenBank parser,
        data = line[5:].rstrip(_CHOMP)
        index = data.find('=')
        if index >= 0:
            descr = data[:index]
            assert descr == "NCBI_TaxID", "Unexpected taxonomy type %s" % descr
            ids = data[index+1:].split(',')
        else:
            ids = data.split(',')

        try:
            #Append to any existing data
            self.data.annotations['organism_host'].extend(ids)
        except KeyError:
            self.data.annotations['organism_host'] = ids

    def organism_classification(self, line):
        #Try and agree with SeqRecord convention from the GenBank parser,
        #which stores this taxonomy lineage ese as a list of strings with
        #key 'taxonomy'.
        #Note that 'ncbi_taxid' is used for the taxonomy ID (line OX)
        line = line[5:].rstrip(_CHOMP)
        cols = [col.strip() for col in line.split(';')]
        try:
            #Append to any existing data
            self.data.annotations['taxonomy'].extend(cols)
        except KeyError:
            self.data.annotations['taxonomy'] = cols

    def taxonomy_id(self, line):
        #Try and agree with SeqRecord convention expected in BioSQL
        #the NCBI taxon id with key 'ncbi_taxid'.
        #Note that 'taxonomy' is used for the taxonomy lineage
        #(held as a list of strings, line type OC)

        line = line[5:].rstrip(_CHOMP)
        index = line.find('=')
        if index >= 0:
            descr = line[:index]
            assert descr == "NCBI_TaxID", "Unexpected taxonomy type %s" % descr
            ids = line[index+1:].split(',')
        else:
            ids = line.split(',')

        try:
            #Append to any existing data
            self.data.annotations['ncbi_taxid'].extend(ids)
        except KeyError:
            self.data.annotations['ncbi_taxid'] = ids

    def reference_number(self, line):
        """RN line, reference number (start of new reference)."""
        from Bio.SeqFeature import Reference
        # if we have a current reference that hasn't been added to
        # the list of references, add it.
        if self._current_ref is not None:
            self.data.annotations['references'].append(self._current_ref)
        else:
            self.data.annotations['references'] = []

        self._current_ref = Reference()

    def reference_position(self, line):
        """RP line, reference position."""
        assert self._current_ref is not None, "RP: missing RN"
        #Should try and store this in self._current_ref.location
        #but the SwissProt locations don't match easily to the
        #format used in GenBank...
        pass

    def reference_cross_reference(self, line):
        """RX line, reference cross-references."""
        assert self._current_ref is not None, "RX: missing RN"
        # The basic (older?) RX line is of the form:
        # RX   MEDLINE; 85132727.
        # or more recently:
        # RX   MEDLINE=95385798; PubMed=7656980;
        # RX   PubMed=15060122; DOI=10.1136/jmg 2003.012781;
        # We look for these cases first and deal with them
        if line.find("=") != -1:
            cols = line[2:].split("; ")
            cols = [x.strip() for x in cols]
            cols = [x for x in cols if x]
            for col in cols:
                x = col.split("=")
                assert len(x) == 2, "I don't understand RX line %s" % line
                key, value = x[0].rstrip(_CHOMP), x[1].rstrip(_CHOMP)
                if key == "MEDLINE":
                    self._current_ref.medline_id = value
                elif key == "PubMed":
                    self._current_ref.pubmed_id = value
                else:
                    #Sadly the SeqFeature.Reference object doesn't
                    #support anything else (yet)
                    pass
        # otherwise we assume we have the type 'RX   MEDLINE; 85132727.'
        else:
            # CLD1_HUMAN in Release 39 and DADR_DIDMA in Release 33
            # have extraneous information in the RX line.  Check for
            # this and chop it out of the line.
            # (noticed by katel@worldpath.net)
            ind = line.find('[NCBI, ExPASy, Israel, Japan]')
            if ind >= 0:
                line = line[:ind]
            cols = line.split()
            # normally we split into the three parts
            assert len(cols) == 3, "I don't understand RX line %s" % line
            key = cols[1].rstrip(_CHOMP)
            value = cols[2].rstrip(_CHOMP)
            if key == "MEDLINE":
                self._current_ref.medline_id = value
            elif key == "PubMed":
                self._current_ref.pubmed_id = value
            else:
                #Sadly the SeqFeature.Reference object doesn't
                #support anything else (yet)
                pass

    def reference_author(self, line):
        """RA line, reference author(s)."""
        assert self._current_ref is not None, "RA: missing RN"
        if self._current_ref.authors:
            self._current_ref.authors += " "
        self._current_ref.authors += line[5:].rstrip("\n")

    def reference_title(self, line):
        """RT line, reference title."""
        assert self._current_ref is not None, "RT: missing RN"
        if self._current_ref.title:
            self._current_ref.title += " "
        self._current_ref.title += line[5:].rstrip("\n")
    
    def reference_location(self, line):
        """RL line, reference 'location' - journal, volume, pages, year."""
        assert self._current_ref is not None, "RL: missing RN"
        if self._current_ref.journal:
            self._current_ref.journal += " "
        self._current_ref.journal += line[5:].rstrip("\n")

    def reference_comment(self, line):
        """RC line, reference comment."""
        assert self._current_ref is not None, "RC: missing RN"
        #This has a key=value; structure...
        #Can we do a better job with the current Reference class?
        if self._current_ref.comment:
            self._current_ref.comment += " "
        self._current_ref.comment += line[5:].rstrip("\n")

def index_file(filename, indexname, rec2key=None):
    """index_file(filename, indexname, rec2key=None)

    Index a SwissProt file.  filename is the name of the file.
    indexname is the name of the dictionary.  rec2key is an
    optional callback that takes a Record and generates a unique key
    (e.g. the accession number) for the record.  If not specified,
    the entry name will be used.

    """
    from Bio.SwissProt import parse
    if not os.path.exists(filename):
        raise ValueError("%s does not exist" % filename)
    
    index = Index.Index(indexname, truncate=1)
    index[Dictionary._Dictionary__filename_key] = filename
    
    handle = open(filename)
    records = parse(handle)
    end = 0L
    for record in records:
        start = end
        end = long(handle.tell())
        length = end - start
        
        if rec2key is not None:
            key = rec2key(record)
        else:
            key = record.entry_name

        if not key:
            raise KeyError("empty sequence key was produced")
        elif key in index:
            raise KeyError("duplicate key %s found" % key)

        index[key] = start, length
