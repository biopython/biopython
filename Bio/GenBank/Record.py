# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Hold GenBank data in a straightforward format.

Classes:
 - Record - All of the information in a GenBank record.
 - Reference - hold reference data for a record.
 - Feature - Hold the information in a Feature Table.
 - Qualifier - Qualifiers on a Feature.

"""

import Bio.GenBank


def _wrapped_genbank(information, indent, wrap_space=1, split_char=" "):
    """Write a line of GenBank info that can wrap over multiple lines (PRIVATE).

    This takes a line of information which can potentially wrap over
    multiple lines, and breaks it up with carriage returns and
    indentation so it fits properly into a GenBank record.

    Arguments:
     - information - The string holding the information we want
       wrapped in GenBank method.
     - indent - The indentation on the lines we are writing.
     - wrap_space - Whether or not to wrap only on spaces in the
       information.
     - split_char - A specific character to split the lines on. By default
       spaces are used.

    """
    info_length = Record.GB_LINE_LENGTH - indent

    if not information:
        # GenBank files use "." for missing data
        return ".\n"

    if wrap_space:
        info_parts = information.split(split_char)
    else:
        cur_pos = 0
        info_parts = []
        while cur_pos < len(information):
            info_parts.append(information[cur_pos : cur_pos + info_length])
            cur_pos += info_length

    # first get the information string split up by line
    output_parts = []
    cur_part = ""
    for info_part in info_parts:
        if len(cur_part) + 1 + len(info_part) > info_length:
            if cur_part:
                if split_char != " ":
                    cur_part += split_char
                output_parts.append(cur_part)
            cur_part = info_part
        else:
            if cur_part == "":
                cur_part = info_part
            else:
                cur_part += split_char + info_part

    # add the last bit of information to the output
    if cur_part:
        output_parts.append(cur_part)

    # now format the information string for return
    output_info = output_parts[0] + "\n"
    for output_part in output_parts[1:]:
        output_info += " " * indent + output_part + "\n"

    return output_info


def _indent_genbank(information, indent):
    """Write out information with the specified indent (PRIVATE).

    Unlike _wrapped_genbank, this function makes no attempt to wrap
    lines -- it assumes that the information already has newlines in the
    appropriate places, and will add the specified indent to the start of
    each line.
    """
    # split the info into lines based on line breaks
    info_parts = information.split("\n")

    # the first line will have no indent
    output_info = info_parts[0] + "\n"
    for info_part in info_parts[1:]:
        output_info += " " * indent + info_part + "\n"

    return output_info


class Record:
    """Hold GenBank information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at GenBank data.

    Attributes:
     - locus - The name specified after the LOCUS keyword in the GenBank
       record. This may be the accession number, or a clone id or something else.
     - size - The size of the record.
     - residue_type - The type of residues making up the sequence in this
       record. Normally something like RNA, DNA or PROTEIN, but may be as
       esoteric as 'ss-RNA circular'.
     - data_file_division - The division this record is stored under in
       GenBank (ie. PLN -> plants; PRI -> humans, primates; BCT -> bacteria...)
     - date - The date of submission of the record, in a form like '28-JUL-1998'
     - accession - list of all accession numbers for the sequence.
     - nid - Nucleotide identifier number.
     - pid - Proteint identifier number
     - version - The accession number + version (ie. AB01234.2)
     - db_source - Information about the database the record came from
     - gi - The NCBI gi identifier for the record.
     - keywords - A list of keywords related to the record.
     - segment - If the record is one of a series, this is info about which
       segment this record is (something like '1 of 6').
     - source - The source of material where the sequence came from.
     - organism - The genus and species of the organism (ie. 'Homo sapiens')
     - taxonomy - A listing of the taxonomic classification of the organism,
       starting general and getting more specific.
     - references - A list of Reference objects.
     - comment - Text with any kind of comment about the record.
     - features - A listing of Features making up the feature table.
     - base_counts - A string with the counts of bases for the sequence.
     - origin - A string specifying info about the origin of the sequence.
     - sequence - A string with the sequence itself.
     - contig - A string of location information for a CONTIG in a RefSeq file
     - project - The genome sequencing project numbers
       (will be replaced by the dblink cross-references in 2009).
     - dblinks - The genome sequencing project number(s) and other links.
       (will replace the project information in 2009).

    """

    # constants for outputting GenBank information
    GB_LINE_LENGTH = 79
    GB_BASE_INDENT = 12
    GB_FEATURE_INDENT = 21
    GB_INTERNAL_INDENT = 2
    GB_OTHER_INTERNAL_INDENT = 3
    GB_FEATURE_INTERNAL_INDENT = 5
    GB_SEQUENCE_INDENT = 9

    BASE_FORMAT = "%-" + str(GB_BASE_INDENT) + "s"
    INTERNAL_FORMAT = (
        " " * GB_INTERNAL_INDENT + "%-" + str(GB_BASE_INDENT - GB_INTERNAL_INDENT) + "s"
    )
    OTHER_INTERNAL_FORMAT = (
        " " * GB_OTHER_INTERNAL_INDENT
        + "%-"
        + str(GB_BASE_INDENT - GB_OTHER_INTERNAL_INDENT)
        + "s"
    )

    BASE_FEATURE_FORMAT = "%-" + str(GB_FEATURE_INDENT) + "s"
    INTERNAL_FEATURE_FORMAT = (
        " " * GB_FEATURE_INTERNAL_INDENT
        + "%-"
        + str(GB_FEATURE_INDENT - GB_FEATURE_INTERNAL_INDENT)
        + "s"
    )
    SEQUENCE_FORMAT = "%" + str(GB_SEQUENCE_INDENT) + "s"

    def __init__(self):
        """Initialize the class."""
        self.accession = []
        self.base_counts = ""
        self.comment = ""
        self.contig = ""
        self.data_file_division = ""
        self.date = ""
        self.db_source = ""
        self.dblinks = []
        self.definition = ""
        self.features = []
        self.gi = ""
        self.keywords = []
        self.locus = ""
        self.molecule_type = ""
        self.nid = ""
        self.organism = ""
        self.origin = ""
        self.pid = ""
        self.primary = []
        self.projects = []
        self.references = []
        self.residue_type = ""
        self.segment = ""
        self.sequence = ""
        self.size = ""
        self.source = ""
        self.taxonomy = []
        self.topology = ""
        self.version = ""
        self.wgs = ""
        self.wgs_scafld = []

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
        output += self._project_line()
        output += self._dblink_line()
        output += self._nid_line()
        output += self._pid_line()
        output += self._keywords_line()
        output += self._db_source_line()
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
        output += self._wgs_line()
        output += self._wgs_scafld_line()
        output += self._contig_line()
        output += "//"
        return output

    def _locus_line(self):
        """Provide the output string for the LOCUS line (PRIVATE)."""
        output = "LOCUS"
        output += " " * 7  # 6-12 spaces
        output += "%-9s" % self.locus
        output += " "  # 22 space
        output += "%7s" % self.size
        if "PROTEIN" in self.residue_type:
            output += " aa"
        else:
            output += " bp "

        # treat circular types differently, since they'll have long residue
        # types
        if "circular" in self.residue_type:
            output += "%17s" % self.residue_type
        # second case: ss-DNA types of records
        elif "-" in self.residue_type:
            output += "%7s" % self.residue_type
            output += " " * 10  # spaces for circular
        else:
            output += " " * 3  # spaces for stuff like ss-
            output += "%-4s" % self.residue_type
            output += " " * 10  # spaces for circular

        output += " " * 2
        output += "%3s" % self.data_file_division
        output += " " * 7  # spaces for 56-63
        output += "%11s" % self.date
        output += "\n"
        return output

    def _definition_line(self):
        """Provide output for the DEFINITION line (PRIVATE)."""
        output = Record.BASE_FORMAT % "DEFINITION"
        output += _wrapped_genbank(self.definition + ".", Record.GB_BASE_INDENT)
        return output

    def _accession_line(self):
        """Output for the ACCESSION line (PRIVATE)."""
        if self.accession:
            output = Record.BASE_FORMAT % "ACCESSION"

            acc_info = ""
            for accession in self.accession:
                acc_info += f"{accession} "
            # strip off an extra space at the end
            acc_info = acc_info.rstrip()
            output += _wrapped_genbank(acc_info, Record.GB_BASE_INDENT)
        else:
            output = ""

        return output

    def _version_line(self):
        """Output for the VERSION line (PRIVATE)."""
        if self.version:
            output = Record.BASE_FORMAT % "VERSION"
            output += self.version
            output += "  GI:"
            output += f"{self.gi}\n"
        else:
            output = ""
        return output

    def _project_line(self):
        output = ""
        if len(self.projects) > 0:
            output = Record.BASE_FORMAT % "PROJECT"
            output += f"{'  '.join(self.projects)}\n"
        return output

    def _dblink_line(self):
        output = ""
        if len(self.dblinks) > 0:
            output = Record.BASE_FORMAT % "DBLINK"
            dblink_info = "\n".join(self.dblinks)
            output += _wrapped_genbank(dblink_info, Record.GB_BASE_INDENT)
        return output

    def _nid_line(self):
        """Output for the NID line. Use of NID is obsolete in GenBank files (PRIVATE)."""
        if self.nid:
            output = Record.BASE_FORMAT % "NID"
            output += f"{self.nid}\n"
        else:
            output = ""
        return output

    def _pid_line(self):
        """Output for PID line. Presumedly, PID usage is also obsolete (PRIVATE)."""
        if self.pid:
            output = Record.BASE_FORMAT % "PID"
            output += f"{self.pid}\n"
        else:
            output = ""
        return output

    def _keywords_line(self):
        """Output for the KEYWORDS line (PRIVATE)."""
        output = ""
        if self.keywords:
            output += Record.BASE_FORMAT % "KEYWORDS"
            keyword_info = ""
            for keyword in self.keywords:
                keyword_info += f"{keyword}; "
            # replace the ; at the end with a period
            keyword_info = keyword_info[:-2]
            keyword_info += "."

            output += _wrapped_genbank(keyword_info, Record.GB_BASE_INDENT)

        return output

    def _db_source_line(self):
        """Output for DBSOURCE line (PRIVATE)."""
        if self.db_source:
            output = Record.BASE_FORMAT % "DBSOURCE"
            output += f"{self.db_source}\n"
        else:
            output = ""
        return output

    def _segment_line(self):
        """Output for the SEGMENT line (PRIVATE)."""
        output = ""
        if self.segment:
            output += Record.BASE_FORMAT % "SEGMENT"
            output += _wrapped_genbank(self.segment, Record.GB_BASE_INDENT)
        return output

    def _source_line(self):
        """Output for SOURCE line on where the sample came from (PRIVATE)."""
        output = Record.BASE_FORMAT % "SOURCE"
        output += _wrapped_genbank(self.source, Record.GB_BASE_INDENT)
        return output

    def _organism_line(self):
        """Output for ORGANISM line with taxonomy info (PRIVATE)."""
        output = Record.INTERNAL_FORMAT % "ORGANISM"
        # Now that species names can be too long, this line can wrap (Bug 2591)
        output += _wrapped_genbank(self.organism, Record.GB_BASE_INDENT)
        output += " " * Record.GB_BASE_INDENT
        taxonomy_info = ""
        for tax in self.taxonomy:
            taxonomy_info += f"{tax}; "
        # replace the ; at the end with a period
        taxonomy_info = taxonomy_info[:-2]
        taxonomy_info += "."
        output += _wrapped_genbank(taxonomy_info, Record.GB_BASE_INDENT)

        return output

    def _comment_line(self):
        """Output for the COMMENT lines (PRIVATE)."""
        output = ""
        if self.comment:
            output += Record.BASE_FORMAT % "COMMENT"
            output += _indent_genbank(self.comment, Record.GB_BASE_INDENT)
        return output

    def _features_line(self):
        """Output for the FEATURES line (PRIVATE)."""
        output = ""
        if len(self.features) > 0:
            output += Record.BASE_FEATURE_FORMAT % "FEATURES"
            output += "Location/Qualifiers\n"
        return output

    def _base_count_line(self):
        """Output for the BASE COUNT line with base information (PRIVATE)."""
        output = ""
        if self.base_counts:
            output += Record.BASE_FORMAT % "BASE COUNT  "
            # split up the base counts into their individual parts
            count_parts = self.base_counts.split(" ")
            while "" in count_parts:
                count_parts.remove("")
            # deal with the standard case, with a normal origin line
            # like: 474 a    356 c    428 g    364 t
            if len(count_parts) % 2 == 0:
                while len(count_parts) > 0:
                    count_info = count_parts.pop(0)
                    count_type = count_parts.pop(0)

                    output += f"{count_info:>7} {count_type}"
            # deal with ugly ORIGIN lines like:
            # 1311257 a2224835 c2190093 g1309889 t
            # by just outputting the raw information
            else:
                output += self.base_counts
            output += "\n"
        return output

    def _origin_line(self):
        """Output for the ORIGIN line (PRIVATE)."""
        output = ""
        # only output the ORIGIN line if we have a sequence
        if self.sequence:
            output += Record.BASE_FORMAT % "ORIGIN"
            if self.origin:
                output += _wrapped_genbank(self.origin, Record.GB_BASE_INDENT)
            else:
                output += "\n"
        return output

    def _sequence_line(self):
        """Output for all of the sequence (PRIVATE)."""
        output = ""
        if self.sequence:
            cur_seq_pos = 0
            while cur_seq_pos < len(self.sequence):
                output += Record.SEQUENCE_FORMAT % str(cur_seq_pos + 1)

                for section in range(6):
                    start_pos = cur_seq_pos + section * 10
                    end_pos = start_pos + 10
                    seq_section = self.sequence[start_pos:end_pos]
                    output += f" {seq_section.lower()}"

                    # stop looping if we are out of sequence
                    if end_pos > len(self.sequence):
                        break

                output += "\n"
                cur_seq_pos += 60
        return output

    def _wgs_line(self):
        output = ""
        if self.wgs:
            output += Record.BASE_FORMAT % "WGS"
            output += self.wgs
        return output

    def _wgs_scafld_line(self):
        output = ""
        if self.wgs_scafld:
            output += Record.BASE_FORMAT % "WGS_SCAFLD"
            output += self.wgs_scafld
        return output

    def _contig_line(self):
        """Output for CONTIG location information from RefSeq (PRIVATE)."""
        output = ""
        if self.contig:
            output += Record.BASE_FORMAT % "CONTIG"
            output += _wrapped_genbank(
                self.contig, Record.GB_BASE_INDENT, split_char=","
            )
        return output


class Reference:
    """Hold information from a GenBank reference.

    Attributes:
     - number - The number of the reference in the listing of references.
     - bases - The bases in the sequence the reference refers to.
     - authors - String with all of the authors.
     - consrtm - Consortium the authors belong to.
     - title - The title of the reference.
     - journal - Information about the journal where the reference appeared.
     - medline_id - The medline id for the reference.
     - pubmed_id - The pubmed_id for the reference.
     - remark - Free-form remarks about the reference.

    """

    def __init__(self):
        """Initialize the class."""
        self.number = ""
        self.bases = ""
        self.authors = ""
        self.consrtm = ""
        self.title = ""
        self.journal = ""
        self.medline_id = ""
        self.pubmed_id = ""
        self.remark = ""

    def __str__(self):
        """Convert the reference to a GenBank format string."""
        output = self._reference_line()
        output += self._authors_line()
        output += self._consrtm_line()
        output += self._title_line()
        output += self._journal_line()
        output += self._medline_line()
        output += self._pubmed_line()
        output += self._remark_line()

        return output

    def _reference_line(self):
        """Output for REFERENCE lines (PRIVATE)."""
        output = Record.BASE_FORMAT % "REFERENCE"
        if self.number:
            if self.bases:
                output += "%-3s" % self.number
                output += f"{self.bases}"
            else:
                output += f"{self.number}"

        output += "\n"
        return output

    def _authors_line(self):
        """Output for AUTHORS information (PRIVATE)."""
        output = ""
        if self.authors:
            output += Record.INTERNAL_FORMAT % "AUTHORS"
            output += _wrapped_genbank(self.authors, Record.GB_BASE_INDENT)
        return output

    def _consrtm_line(self):
        """Output for CONSRTM information (PRIVATE)."""
        output = ""
        if self.consrtm:
            output += Record.INTERNAL_FORMAT % "CONSRTM"
            output += _wrapped_genbank(self.consrtm, Record.GB_BASE_INDENT)
        return output

    def _title_line(self):
        """Output for TITLE information (PRIVATE)."""
        output = ""
        if self.title:
            output += Record.INTERNAL_FORMAT % "TITLE"
            output += _wrapped_genbank(self.title, Record.GB_BASE_INDENT)
        return output

    def _journal_line(self):
        """Output for JOURNAL information (PRIVATE)."""
        output = ""
        if self.journal:
            output += Record.INTERNAL_FORMAT % "JOURNAL"
            output += _wrapped_genbank(self.journal, Record.GB_BASE_INDENT)
        return output

    def _medline_line(self):
        """Output for MEDLINE information (PRIVATE)."""
        output = ""
        if self.medline_id:
            output += Record.INTERNAL_FORMAT % "MEDLINE"
            output += self.medline_id + "\n"
        return output

    def _pubmed_line(self):
        """Output for PUBMED information (PRIVATE)."""
        output = ""
        if self.pubmed_id:
            output += Record.OTHER_INTERNAL_FORMAT % "PUBMED"
            output += self.pubmed_id + "\n"
        return output

    def _remark_line(self):
        """Output for REMARK information (PRIVATE)."""
        output = ""
        if self.remark:
            output += Record.INTERNAL_FORMAT % "REMARK"
            output += _wrapped_genbank(self.remark, Record.GB_BASE_INDENT)
        return output


class Feature:
    """Hold information about a Feature in the Feature Table of GenBank record.

    Attributes:
     - key - The key name of the feature (ie. source)
     - location - The string specifying the location of the feature.
     - qualifiers - A list of Qualifier objects in the feature.

    """

    def __init__(self, key="", location=""):
        """Initialize the class."""
        self.key = key
        self.location = location
        self.qualifiers = []

    def __repr__(self):
        """Representation of the object for debugging or logging."""
        return f"Feature(key={self.key!r}, location={self.location!r})"

    def __str__(self):
        """Return feature as a GenBank format string."""
        output = Record.INTERNAL_FEATURE_FORMAT % self.key
        output += _wrapped_genbank(
            self.location, Record.GB_FEATURE_INDENT, split_char=","
        )
        for qualifier in self.qualifiers:
            output += str(qualifier)
        return output


class Qualifier:
    """Hold information about a qualifier in a GenBank feature.

    Attributes:
     - key - The key name of the qualifier (ie. /organism=)
     - value - The value of the qualifier ("Dictyostelium discoideum").

    """

    def __init__(self, key="", value=""):
        """Initialize the class."""
        self.key = key
        self.value = value

    def __repr__(self):
        """Representation of the object for debugging or logging."""
        return f"Qualifier(key={self.key!r}, value={self.value!r})"

    def __str__(self):
        """Return feature qualifier as a GenBank format string."""
        output = " " * Record.GB_FEATURE_INDENT
        # determine whether we can wrap on spaces
        space_wrap = 1
        for no_space_key in Bio.GenBank._BaseGenBankConsumer.remove_space_keys:
            if no_space_key in self.key:
                space_wrap = 0
        # return double quotes as-is, leave it to the user to escape them
        return output + _wrapped_genbank(
            self.key + self.value, Record.GB_FEATURE_INDENT, space_wrap
        )
