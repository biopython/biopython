"""Hold GenBank data in a straightforward format.

classes:
o Record - All of the information in a GenBank record.
o Reference - hold reference data for a record.
o Feature - Hold the information in a Feature Table.
o Qualifier - Qualifiers on a Feature.
"""
# standard modules
import string

def _wrapped_genbank(information, indent, wrap_space = 1):
    """Write a line of GenBank info that can wrap over multiple lines.

    This takes a line of information which can potentially wrap over
    multiple lines, and breaks it up with carriage returns and
    indentation so it fits properly.

    Arguments:

    o information - The string holding the information we want
    wrapped in GenBank method.

    o indent - The indentation on the lines we are writing.

    o wrap_space - Whether or not to wrap only on spaces in the
    information. 
    """
    info_length = Record.GB_LINE_LENGTH - indent

    line_start = 0
    output_info = ""
    while 1:
        if wrap_space:
            if line_start + info_length >= len(information):
                line_end = len(information)
            else:
                line_end = information.rfind(" ", line_start,
                                             line_start + info_length)
            # if we can't find a space before the end, we are forced
            # to wrap with what we've got
            if line_end == -1:
                line_end = line_start + info_length
        else:
            line_end = line_start + info_length

        cur_info_line = information[line_start:line_end]

        if not(cur_info_line):
            break

        # only add indent spaces if we are not at the start of the info
        if line_start != 0:
            output_info += " " * indent
        # add the actual information
        output_info += "%s\n" % string.lstrip(cur_info_line)

        # update where we are at in the information
        line_start = line_end

    return output_info

class Record:
    """Hold GenBank information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at GenBank data.

    Attributes:
    o locus - The name specified after the LOCUS keyword in the GenBank
    record. This may be the accession number, or a clone id or something else.
    o size - The size of the record.
    o residue_type - The type of residues making up the sequence in this
    record. Normally something like RNA, DNA or PROTEIN, but may be as
    esoteric as 'ss-RNA circular'.
    o data_file_division - The division this record is stored under in
    GenBank (ie. PLN -> plants; PRI -> humans, primates; BCT -> bacteria...)
    o date - The date of submission of the record, in a form like '28-JUL-1998'
    o accession - list of all accession numbers for the sequence.
    o nid - Nucleotide identifier number.
    o version - The accession number + version (ie. AB01234.2)
    o gi - The NCBI gi identifier for the record.
    o keywords - A list of keywords related to the record.
    o segment - If the record is one of a series, this is info about which
    segment this record is (something like '1 of 6').
    o source - The source of material where the sequence came from.
    o organism - The genus and species of the organism (ie. 'Homo sapiens')
    o taxonomy - A listing of the taxonomic classification of the organism,
    starting general and getting more specific.
    o references - A list of Reference objects.
    o comment - Text with any kind of comment about the record.
    o features - A listing of Features making up the feature table.
    o base_counts - A string with the counts of bases for the sequence.
    o origin - A string specifying info about the origin of the sequence.
    o sequence - A string with the sequence itself.
    """
    # constants for outputting GenBank information
    GB_LINE_LENGTH = 79
    GB_BASE_INDENT = 12
    GB_FEATURE_INDENT = 21
    GB_INTERNAL_INDENT = 2
    GB_FEATURE_INTERNAL_INDENT = 5
    GB_SEQUENCE_INDENT = 9

    BASE_FORMAT = "%-" + str(GB_BASE_INDENT) + "s"
    INTERNAL_FORMAT = " " * GB_INTERNAL_INDENT + "%-" + \
                      str(GB_BASE_INDENT - GB_INTERNAL_INDENT) + "s"

    BASE_FEATURE_FORMAT = "%-" + str(GB_FEATURE_INDENT) + "s"
    INTERNAL_FEATURE_FORMAT = " " * GB_FEATURE_INTERNAL_INDENT + "%-" + \
                              str(GB_FEATURE_INDENT -
                                  GB_FEATURE_INTERNAL_INDENT) + "s"
    SEQUENCE_FORMAT = "%" + str(GB_SEQUENCE_INDENT) + "s"
    
    def __init__(self):
        self.locus = ''
        self.size = ''
        self.residue_type = ''
        self.data_file_division = ''
        self.date = ''
        self.definition = ''
        self.accession = []
        self.nid = ''
        self.version = ''
        self.gi = ''
        self.keywords = []
        self.segment = ''
        self.source = ''
        self.organism = ''
        self.taxonomy = []
        self.references = []
        self.comment = ''
        self.features = []
        self.base_counts = ''
        self.origin = ''
        self.sequence = ''

    def __str__(self):
        """Provide a GenBank formatted output option for a Record.

        The objective of this is to provide an easy way to read in a GenBank
        record, modify it somehow, and then output it in 'GenBank format.'
        We are striving to make this work so that a parsed Record that is
        output using this function will look exactly like the original
        record.

        Much of the output is based on format description info at:

        ftp://ncbi.nlm.nih.gov/genbank/gbrel.txt
        """
        output = self._locus_line()
        output += self._definition_line()
        output += self._accession_line()
        output += self._version_line()
        output += self._nid_line()
        output += self._keywords_line()
        output += self._segment_line()
        output += self._source_line()
        output += self._organism_line()
        for reference in self.references:
            output += str(reference)
        output += self._comment_line()
        output += self._features_line()
        for feature in self.features:
            output += str(feature)
        output += self._base_count_line()
        output += self._origin_line()
        output += self._sequence_line()

        output += "//"
        return output
            
    def _locus_line(self):
        """Provide the output string for the LOCUS line.
        """
        output = "LOCUS"
        output += " " * 7 # 6-12 spaces
        output += "%9s" % self.locus
        output += " " # 22 space
        output += "%7s" % self.size
        output += " bp"
        # treat circular types differently, since they'll have long residue
        # types
        if self.residue_type.find("circular") >= 0:
            output += "%18s" % self.residue_type
        else:
            output += "%9s" % self.residue_type
            output += "         " # spaces for circular

        output += " " # space at 52
        output += "%3s" % self.data_file_division
        output += " " * 8 # spaces for 56-63
        output += "%11s" % self.date
        output += "\n"
        return output

    def _definition_line(self):
        """Provide output for the DEFINITION line.
        """
        output = Record.BASE_FORMAT % "DEFINITION"
        output += _wrapped_genbank(self.definition, Record.GB_BASE_INDENT)
        return output

    def _accession_line(self):
        """Output for the ACCESSION line.
        """
        output = Record.BASE_FORMAT % "ACCESSION"

        acc_info = ""
        for accession in self.accession:
            acc_info += "%s " % accession
        output += _wrapped_genbank(acc_info, Record.GB_BASE_INDENT)
        
        return output

    def _version_line(self):
        """Output for the VERSION line.
        """
        output = Record.BASE_FORMAT % "VERSION"
        output += self.version
        output += " GI:"
        output += "%s\n" % self.gi
        return output

    def _nid_line(self):
        """Output for the NID line. Use of NID is obsolete in GenBank files.
        """
        if self.nid:
            output = Record.BASE_FORMAT % "NID"
            output += "%s\n" % self.nid
        else:
            output = ""
        return output

    def _keywords_line(self):
        """Output for the KEYWORDS line.
        """
        output = ""
        if len(self.keywords) >= 0:
            output +=  Record.BASE_FORMAT % "KEYWORDS"
            keyword_info = ""
            for keyword in self.keywords:
                keyword_info += "%s; " % keyword
            # replace the ; at the end with a period
            keyword_info = keyword_info[:-2]
            keyword_info += "."
            
            output += _wrapped_genbank(keyword_info,
                                       Record.GB_BASE_INDENT)

        return output

    def _segment_line(self):
        """Output for the SEGMENT line.
        """
        output = ""
        if self.segment:
            output += Record.BASE_FORMAT % "SEGMENT"
        return output

    def _source_line(self):
        """Output for SOURCE line on where the sample came from.
        """
        output = Record.BASE_FORMAT % "SOURCE"
        output += _wrapped_genbank(self.source, Record.GB_BASE_INDENT)
        return output
    
    def _organism_line(self):
        """Output for ORGANISM line with taxonomy info.
        """
        output = Record.INTERNAL_FORMAT % "ORGANISM"
        output += "%s\n" % self.organism
        output += " " * Record.GB_BASE_INDENT
        taxonomy_info = ""
        for tax in self.taxonomy:
            taxonomy_info += "%s; " % tax
        # replace the ; at the end with a period
        taxonomy_info = taxonomy_info[:-2]
        taxonomy_info += "."
        output += _wrapped_genbank(taxonomy_info, Record.GB_BASE_INDENT)

        return output
            
    def _comment_line(self):
        """Output for the COMMENT lines.
        """
        output = ""
        if self.comment:
            output += Record.BASE_FORMAT % "COMMENT"
            output += _wrapped_genbank(self.comment,
                                       Record.GB_BASE_INDENT)
        return output

    def _features_line(self):
        """Output for the FEATURES line.
        """
        output = ""
        if len(self.features) > 0:
            output += Record.BASE_FEATURE_FORMAT % "FEATURES"
            output += "Location/Qualifiers\n"
        return output

    def _base_count_line(self):
        """Output for the BASE COUNT line with base information.
        """
        output = Record.BASE_FORMAT % "BASE_COUNT"
        output += "%s\n" % self.base_counts 
        return output

    def _origin_line(self):
        """Output for the ORIGIN line
        """
        output = Record.BASE_FORMAT % "ORIGIN"
        if self.origin:
            output += _wrapped_genbank(self.origin,
                                       Record.GB_BASE_INDENT)
        else:
            output += "\n"
        return output

    def _sequence_line(self):
        """Output for all of the sequence.
        """
        output = ""
        cur_seq_pos = 0
        while cur_seq_pos < len(self.sequence):
            output += Record.SEQUENCE_FORMAT % str(cur_seq_pos + 1)

            for section in range(6):
                start_pos = cur_seq_pos + section * 10
                end_pos = start_pos + 10
                seq_section = self.sequence[start_pos:end_pos]
                output += " %s" % seq_section.lower()

                # stop looping if we are out of sequence
                if end_pos > len(self.sequence):
                    break
                
            output += "\n"
            cur_seq_pos += 60
        return output
        
class Reference:
    """Hold information from a GenBank reference.

    Attributes:
    o number - The number of the reference in the listing of references.
    o bases - The bases in the sequence the reference refers to.
    o authors - String with all of the authors.
    o title - The title of the reference.
    o journal - Information about the journal where the reference appeared.
    o medline_id - The medline id for the reference.
    o pubmed_id - The pubmed_id for the reference.
    o remark - Free-form remarks about the reference.
    """
    def __init__(self):
        self.number = ''
        self.bases = ''
        self.authors = ''
        self.title = ''
        self.journal = ''
        self.medline_id = ''
        self.pubmed_id = ''
        self.remark = ''

    def __str__(self):
        output = self._reference_line()
        output += self._authors_line()
        output += self._title_line()
        output += self._journal_line()
        output += self._medline_line()
        output += self._pubmed_line()
        output += self._remark_line()
        
        return output

    def _reference_line(self):
        """Output for REFERENCE lines.
        """
        output = Record.BASE_FORMAT % "REFERENCE"
        if self.number:
            output += "%s" % self.number
        if self.number and self.bases:
            output += " "
        if self.bases:
            output += "%s" % self.bases
        output += "\n"
        return output

    def _authors_line(self):
        """Output for AUTHORS information.
        """
        output = ""
        if self.authors:
            output += Record.INTERNAL_FORMAT % "AUTHORS"
            output += _wrapped_genbank(self.authors, Record.GB_BASE_INDENT)
        return output

    def _title_line(self):
        """Output for TITLE information.
        """
        output = ""
        if self.title:
            output += Record.INTERNAL_FORMAT % "TITLE"
            output += _wrapped_genbank(self.title, Record.GB_BASE_INDENT)
        return output

    def _journal_line(self):
        """Output for JOURNAL information.
        """
        output = ""
        if self.journal:
            output += Record.INTERNAL_FORMAT % "JOURNAL"
            output += _wrapped_genbank(self.journal, Record.GB_BASE_INDENT)
        return output

    def _medline_line(self):
        """Output for MEDLINE information.
        """
        output = ""
        if self.medline_id:
            output += Record.INTERNAL_FORMAT % "MEDLINE"
            output += self.medline_id
        return output
    
    def _pubmed_line(self):
        """Output for PUBMED information.
        """
        output = ""
        if self.pubmed_id:
            output += Record.INTERNAL_FORMAT % "PUBMED"
            output += self.pubmed_id
        return output
    
    def _remark_line(self):
        """Output for REMARK information.
        """
        output = ""
        if self.remark:
            output += Record.INTERNAL_FORMAT % "RECORD"
            output += _wrapped_genbank(self.remark, Record.GB_BASE_INDENT)
        return output
    
class Feature:
    """Hold information about a Feature in the Feature Table of GenBank record.

    Attributes:
    o key - The key name of the featue (ie. source)
    o location - The string specifying the location of the feature.
    o qualfiers - A listing Qualifier objects in the feature.
    """
    def __init__(self):
        self.key = ''
        self.location = ''
        self.qualifiers = []

    def __str__(self):
        output = Record.INTERNAL_FEATURE_FORMAT % self.key
        output += _wrapped_genbank(self.location, Record.GB_FEATURE_INDENT)
        for qualifier in self.qualifiers:
            output += " " * Record.GB_FEATURE_INDENT
            output += qualifier.key
            output += _wrapped_genbank(qualifier.value,
                                       Record.GB_FEATURE_INDENT)
        return output

class Qualifier:
    """Hold information about a qualifier in a GenBank feature.

    Attributes:
    o key - The key name of the qualifier (ie. /organim=)
    o value - The value of the qualifier ("Dictyostelium discoideum").
    """
    def __init__(self):
        self.key = ''
        self.value = ''
