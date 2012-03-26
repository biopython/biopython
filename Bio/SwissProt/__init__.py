# Copyright 2007 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
This module provides code to work with the sprotXX.dat file from
SwissProt.
http://www.expasy.ch/sprot/sprot-top.html

Tested with:
Release 56.9, 03-March-2009.


Classes:
Record             Holds SwissProt data.
Reference          Holds reference data from a SwissProt record.

Functions:
read               Read one SwissProt record
parse              Read multiple SwissProt records

"""

from Bio._py3k import _as_string

class Record(object):
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
    host_organism     A list of names of the hosts of a virus, if any.
    host_taxonomy_id  A list of NCBI taxonomy id's of the hosts, if any.
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
        
        self.description = []
        self.gene_name = ''
        self.organism = []
        self.organelle = ''
        self.organism_classification = []
        self.taxonomy_id = []
        self.host_organism = []
        self.host_taxonomy_id = []
        self.references = []
        self.comments = []
        self.cross_references = []
        self.keywords = []
        self.features = []
        
        self.seqinfo = None
        self.sequence = ''


class Reference(object):
    """Holds information from one reference in a SwissProt entry.

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
        self.authors = []
        self.title = []
        self.location = []


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

 
# Everything below is considered private


def _read(handle):
    record = None
    unread = ""
    for line in handle:
        #This is for Python 3 to cope with a binary handle (byte strings),
        #or a text handle (unicode strings):
        line = _as_string(line)
        key, value = line[:2], line[5:].rstrip()
        if unread:
            value = unread + " " + value
            unread = ""
        if key=='**':
            #See Bug 2353, some files from the EBI have extra lines
            #starting "**" (two asterisks/stars).  They appear
            #to be unofficial automated annotations. e.g.
            #**
            #**   #################    INTERNAL SECTION    ##################
            #**HA SAM; Annotated by PicoHamap 1.88; MF_01138.1; 09-NOV-2003.
            pass
        elif key=='ID':
            record = Record()
            _read_id(record, line)
            _sequence_lines = []
        elif key=='AC':
            accessions = [word for word in value.rstrip(";").split("; ")]
            record.accessions.extend(accessions)
        elif key=='DT':
            _read_dt(record, line)
        elif key=='DE':
            record.description.append(value.strip())
        elif key=='GN':
            if record.gene_name:
                record.gene_name += " "
            record.gene_name += value
        elif key=='OS':
            record.organism.append(value)
        elif key=='OG':
            record.organelle += line[5:]
        elif key=='OC':
            cols = [col for col in value.rstrip(";.").split("; ")]
            record.organism_classification.extend(cols)
        elif key=='OX':
            _read_ox(record, line)
        elif key=='OH':
            _read_oh(record, line)
        elif key=='RN':
            reference = Reference()
            _read_rn(reference, value)
            record.references.append(reference)
        elif key=='RP':
            assert record.references, "RP: missing RN"
            record.references[-1].positions.append(value)
        elif key=='RC':
            assert record.references, "RC: missing RN"
            reference = record.references[-1]
            unread = _read_rc(reference, value)
        elif key=='RX':
            assert record.references, "RX: missing RN"
            reference = record.references[-1]
            _read_rx(reference, value)
        elif key=='RL':
            assert record.references, "RL: missing RN"
            reference = record.references[-1]
            reference.location.append(value)
        # In UniProt release 1.12 of 6/21/04, there is a new RG
        # (Reference Group) line, which references a group instead of
        # an author.  Each block must have at least 1 RA or RG line.
        elif key=='RA':
            assert record.references, "RA: missing RN"
            reference = record.references[-1]
            reference.authors.append(value)
        elif key=='RG':
            assert record.references, "RG: missing RN"
            reference = record.references[-1]
            reference.authors.append(value)
        elif key=="RT":
            assert record.references, "RT: missing RN"
            reference = record.references[-1]
            reference.title.append(value)
        elif key=='CC':
            _read_cc(record, line)
        elif key=='DR':
            _read_dr(record, value)
        elif key=='PE':
            #TODO - Record this information?
            pass
        elif key=='KW':
            cols = value.rstrip(";.").split('; ')
            record.keywords.extend(cols)
        elif key=='FT':
            _read_ft(record, line)
        elif key=='SQ':
            cols = value.split()
            assert len(cols) == 7, "I don't understand SQ line %s" % line
            # Do more checking here?
            record.seqinfo = int(cols[1]), int(cols[3]), cols[5]
        elif key=='  ':
            _sequence_lines.append(value.replace(" ", "").rstrip())
        elif key=='//':
            # Join multiline data into one string
            record.description = " ".join(record.description)
            record.organism = " ".join(record.organism)
            record.organelle   = record.organelle.rstrip()
            for reference in record.references:
                reference.authors = " ".join(reference.authors).rstrip(";")
                reference.title = " ".join(reference.title).rstrip(";")
                if reference.title.startswith('"') and reference.title.endswith('"'):
                    reference.title = reference.title[1:-1] #remove quotes
                reference.location = " ".join(reference.location)
            record.sequence = "".join(_sequence_lines)
            return record
        else:
            raise ValueError("Unknown keyword '%s' found" % key)
    if record:
        raise ValueError("Unexpected end of stream.")


def _read_id(record, line):
    cols = line[5:].split()
    #Prior to release 51, included with MoleculeType:
    #ID   EntryName DataClass; MoleculeType; SequenceLength AA.
    #
    #Newer files lack the MoleculeType:
    #ID   EntryName DataClass; SequenceLength AA.
    if len(cols) == 5:
        record.entry_name = cols[0]
        record.data_class = cols[1].rstrip(";")
        record.molecule_type = cols[2].rstrip(";")
        record.sequence_length = int(cols[3])
    elif len(cols) == 4:
        record.entry_name = cols[0]
        record.data_class = cols[1].rstrip(";")
        record.molecule_type = None
        record.sequence_length = int(cols[2])
    else:
        raise ValueError("ID line has unrecognised format:\n"+line)
    # check if the data class is one of the allowed values
    allowed = ('STANDARD', 'PRELIMINARY', 'IPI', 'Reviewed', 'Unreviewed')
    if record.data_class not in allowed:
        raise ValueError("Unrecognized data class %s in line\n%s" % \
              (record.data_class, line))
    # molecule_type should be 'PRT' for PRoTein
    # Note that has been removed in recent releases (set to None)
    if record.molecule_type not in (None, 'PRT'):
        raise ValueError("Unrecognized molecule type %s in line\n%s" % \
              (record.molecule_type, line))


def _read_dt(record, line):
    value = line[5:]
    uprline = value.upper()
    cols = value.rstrip().split()
    if 'CREATED' in uprline \
    or 'LAST SEQUENCE UPDATE' in uprline \
    or 'LAST ANNOTATION UPDATE' in uprline:
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
                "Could not find Rel. in DT line: %s" % line
        version_index = rel_index + 1
        # get the version information
        str_version = cols[version_index].rstrip(",")
        # no version number
        if str_version == '':
            version = 0
        # dot versioned
        elif str_version.find(".") >= 0:
            version = str_version
        # integer versioned
        else:
            version = int(str_version)
        date = cols[0]

        if 'CREATED' in uprline:
            record.created = date, version
        elif 'LAST SEQUENCE UPDATE' in uprline:
            record.sequence_update = date, version
        elif 'LAST ANNOTATION UPDATE' in uprline:
            record.annotation_update = date, version
        else:
            assert False, "Shouldn't reach this line!"
    elif 'INTEGRATED INTO' in uprline \
    or 'SEQUENCE VERSION' in uprline \
    or 'ENTRY VERSION' in uprline:
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
        date = cols[0].rstrip(",")

        # Re-use the historical property names, even though
        # the meaning has changed slighty:
        if "INTEGRATED"  in uprline:
            record.created = date, version
        elif 'SEQUENCE VERSION' in uprline:
            record.sequence_update = date, version
        elif 'ENTRY VERSION' in uprline:
            record.annotation_update = date, version
        else:
            assert False, "Shouldn't reach this line!"
    else:
        raise ValueError("I don't understand the date line %s" % line)


def _read_ox(record, line):
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
    if record.taxonomy_id:
        ids = line[5:].rstrip().rstrip(";")
    else:
        descr, ids = line[5:].rstrip().rstrip(";").split("=")
        assert descr == "NCBI_TaxID", "Unexpected taxonomy type %s" % descr
    record.taxonomy_id.extend(ids.split(', '))


def _read_oh(record, line):
    # Line type OH (Organism Host) for viral hosts
    assert line[5:].startswith("NCBI_TaxID="), "Unexpected %s" % line
    line = line[16:].rstrip()
    assert line[-1]=="." and line.count(";")==1, line
    taxid, name = line[:-1].split(";")
    record.host_taxonomy_id.append(taxid.strip())
    record.host_organism.append(name.strip())


def _read_rn(reference, rn):
    assert rn[0] == '[' and rn[-1] == ']', "Missing brackets %s" % rn
    reference.number = int(rn[1:-1])


def _read_rc(reference, value):
    cols = value.split(';')
    if value[-1]==';':
        unread = ""
    else:
        cols, unread = cols[:-1], cols[-1]
    for col in cols:
        if not col:  # last column will be the empty string
            return
        # The token is everything before the first '=' character.
        i = col.find("=")
        if i>=0:
            token, text = col[:i], col[i+1:]
            comment = token.lstrip(), text
            reference.comments.append(comment)
        else:
            comment = reference.comments[-1]
            comment = "%s %s" % (comment, col)
            reference.comments[-1] = comment
    return unread


def _read_rx(reference, value):
    # The basic (older?) RX line is of the form:
    # RX   MEDLINE; 85132727.
    # but there are variants of this that need to be dealt with (see below)

    # CLD1_HUMAN in Release 39 and DADR_DIDMA in Release 33
    # have extraneous information in the RX line.  Check for
    # this and chop it out of the line.
    # (noticed by katel@worldpath.net)
    value = value.replace(' [NCBI, ExPASy, Israel, Japan]','')

    # RX lines can also be used of the form
    # RX   PubMed=9603189;
    # reported by edvard@farmasi.uit.no
    # and these can be more complicated like:
    # RX   MEDLINE=95385798; PubMed=7656980;
    # RX   PubMed=15060122; DOI=10.1136/jmg 2003.012781;
    # We look for these cases first and deal with them
    warn = False
    if "=" in value:
        cols = value.split("; ")
        cols = [x.strip() for x in cols]
        cols = [x for x in cols if x]
        for col in cols:
            x = col.split("=")
            if len(x) != 2 or x == ("DOI", "DOI"):
                warn = True
                break
            assert len(x) == 2, "I don't understand RX line %s" % value
            reference.references.append((x[0], x[1].rstrip(";")))
    # otherwise we assume we have the type 'RX   MEDLINE; 85132727.'
    else:
        cols = value.split("; ")
        # normally we split into the three parts
        if len(cols) != 2:
            warn = True
        else:
            reference.references.append((cols[0].rstrip(";"), cols[1].rstrip(".")))
    if warn:
        import warnings
        from Bio import BiopythonParserWarning
        warnings.warn("Possibly corrupt RX line %r" % value,
                      BiopythonParserWarning)

def _read_cc(record, line):
    key, value = line[5:8], line[9:].rstrip()
    if key=='-!-':   # Make a new comment
        record.comments.append(value)
    elif key=='   ': # add to the previous comment
        if not record.comments:
            # TCMO_STRGA in Release 37 has comment with no topic
            record.comments.append(value)
        else:
            record.comments[-1] += " " + value


def _read_dr(record, value):
    # Remove the comments at the end of the line
    i = value.find(' [')
    if i >= 0:
        value = value[:i]
    cols = value.rstrip(".").split('; ')
    record.cross_references.append(tuple(cols))


def _read_ft(record, line):
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
    #if there is a feature_id (FTId), store it away
    if line[29:35]==r"/FTId=":
        ft_id = line[35:70].rstrip()[:-1]
        description = ""
    else:
        ft_id =""
        description = line[29:70].rstrip()
    if not name:  # is continuation of last one
        assert not from_res and not to_res
        name, from_res, to_res, old_description,old_ft_id = record.features[-1]
        del record.features[-1]
        description = ("%s %s" % (old_description, description)).strip()

        # special case -- VARSPLIC, reported by edvard@farmasi.uit.no
        if name == "VARSPLIC":
            # Remove unwanted spaces in sequences.
            # During line carryover, the sequences in VARSPLIC can get mangled
            # with unwanted spaces like:
            # 'DISSTKLQALPSHGLESIQT -> PCRATGWSPFRRSSPC LPTH'
            # We want to check for this case and correct it as it happens.
            descr_cols = description.split(" -> ")
            if len(descr_cols) == 2:
                first_seq, second_seq = descr_cols
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
    record.features.append((name, from_res, to_res, description,ft_id))


if __name__ == "__main__":
    print "Quick self test..."

    example_filename = "../../Tests/SwissProt/sp008"

    import os
    if not os.path.isfile(example_filename):
        print "Missing test file %s" % example_filename
    else:
        #Try parsing it!
        
        handle = open(example_filename)
        records = parse(handle)
        for record in records:
            print record.entry_name
            print ",".join(record.accessions)
            print record.keywords
            print repr(record.organism)
            print record.sequence[:20] + "..."
        handle.close()
