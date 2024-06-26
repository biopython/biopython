# Copyright 2000 by Jeffrey Chang, Brad Chapman.  All rights reserved.
# Copyright 2006-2017 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to work with GenBank formatted files.

Rather than using Bio.GenBank, you are now encouraged to use Bio.SeqIO with
the "genbank" or "embl" format names to parse GenBank or EMBL files into
SeqRecord and SeqFeature objects (see the Biopython tutorial for details).

Using Bio.GenBank directly to parse GenBank files is only useful if you want
to obtain GenBank-specific Record objects, which is a much closer
representation to the raw file contents than the SeqRecord alternative from
the FeatureParser (used in Bio.SeqIO).

To use the Bio.GenBank parser, there are two helper functions:

    - read                  Parse a handle containing a single GenBank record
      as Bio.GenBank specific Record objects.
    - parse                 Iterate over a handle containing multiple GenBank
      records as Bio.GenBank specific Record objects.

The following internal classes are not intended for direct use and may
be deprecated in a future release.

Classes:
 - Iterator              Iterate through a file of GenBank entries
 - FeatureParser         Parse GenBank data in SeqRecord and SeqFeature objects.
 - RecordParser          Parse GenBank data into a Record object.

Exceptions:
 - ParserFailureError    Exception indicating a failure in the parser (ie.
   scanner or consumer)

"""

import re
import warnings

from Bio import BiopythonParserWarning
from Bio.Seq import Seq
from Bio.SeqFeature import Location
from Bio.SeqFeature import Reference
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import LocationParserError

# other Bio.GenBank stuff
from .utils import FeatureValueCleaner
from .Scanner import GenBankScanner


# Constants used to parse GenBank header lines
GENBANK_INDENT = 12
GENBANK_SPACER = " " * GENBANK_INDENT

# Constants for parsing GenBank feature lines
FEATURE_KEY_INDENT = 5
FEATURE_QUALIFIER_INDENT = 21
FEATURE_KEY_SPACER = " " * FEATURE_KEY_INDENT
FEATURE_QUALIFIER_SPACER = " " * FEATURE_QUALIFIER_INDENT


class Iterator:
    """Iterator interface to move over a file of GenBank entries one at a time (OBSOLETE).

    This class is likely to be deprecated in a future release of Biopython.
    Please use Bio.SeqIO.parse(..., format="gb") or Bio.GenBank.parse(...)
    for SeqRecord and GenBank specific Record objects respectively instead.
    """

    def __init__(self, handle, parser=None):
        """Initialize the iterator.

        Arguments:
         - handle - A handle with GenBank entries to iterate through.
         - parser - An optional parser to pass the entries through before
           returning them. If None, then the raw entry will be returned.

        """
        self.handle = handle
        self._parser = parser

    def __next__(self):
        """Return the next GenBank record from the handle.

        Will return None if we ran out of records.
        """
        if self._parser is None:
            lines = []
            while True:
                line = self.handle.readline()
                if not line:
                    return None  # Premature end of file?
                lines.append(line)
                if line.rstrip() == "//":
                    break
            return "".join(lines)
        try:
            return self._parser.parse(self.handle)
        except StopIteration:
            return None

    def __iter__(self):
        """Iterate over the records."""
        return iter(self.__next__, None)


class ParserFailureError(ValueError):
    """Failure caused by some kind of problem in the parser."""


_cleaner = FeatureValueCleaner()


class FeatureParser:
    """Parse GenBank files into Seq + Feature objects (OBSOLETE).

    Direct use of this class is discouraged, and may be deprecated in
    a future release of Biopython.

    Please use Bio.SeqIO.parse(...) or Bio.SeqIO.read(...) instead.
    """

    def __init__(self, debug_level=0, use_fuzziness=1, feature_cleaner=None):
        """Initialize a GenBank parser and Feature consumer.

        Arguments:
         - debug_level - An optional argument that species the amount of
           debugging information the parser should spit out. By default we have
           no debugging info (the fastest way to do things), but if you want
           you can set this as high as two and see exactly where a parse fails.
         - use_fuzziness - Specify whether or not to use fuzzy representations.
           The default is 1 (use fuzziness).
         - feature_cleaner - A class which will be used to clean out the
           values of features. This class must implement the function
           clean_value. GenBank.utils has a "standard" cleaner class, which
           is used by default.

        """
        self._scanner = GenBankScanner(debug_level)
        self.use_fuzziness = use_fuzziness
        if feature_cleaner:
            self._cleaner = feature_cleaner
        else:
            self._cleaner = _cleaner  # default

    def parse(self, handle):
        """Parse the specified handle."""
        _consumer = _FeatureConsumer(self.use_fuzziness, self._cleaner)
        self._scanner.feed(handle, _consumer)
        return _consumer.data


class RecordParser:
    """Parse GenBank files into Record objects (OBSOLETE).

    Direct use of this class is discouraged, and may be deprecated in
    a future release of Biopython.

    Please use the Bio.GenBank.parse(...) or Bio.GenBank.read(...) functions
    instead.
    """

    def __init__(self, debug_level=0):
        """Initialize the parser.

        Arguments:
         - debug_level - An optional argument that species the amount of
           debugging information the parser should spit out. By default we have
           no debugging info (the fastest way to do things), but if you want
           you can set this as high as two and see exactly where a parse fails.

        """
        self._scanner = GenBankScanner(debug_level)

    def parse(self, handle):
        """Parse the specified handle into a GenBank record."""
        _consumer = _RecordConsumer()

        self._scanner.feed(handle, _consumer)
        return _consumer.data


class _BaseGenBankConsumer:
    """Abstract GenBank consumer providing useful general functions (PRIVATE).

    This just helps to eliminate some duplication in things that most
    GenBank consumers want to do.
    """

    # Special keys in GenBank records that we should remove spaces from
    # For instance, \translation keys have values which are proteins and
    # should have spaces and newlines removed from them. This class
    # attribute gives us more control over specific formatting problems.
    remove_space_keys = ["translation"]

    def __init__(self):
        pass

    @staticmethod
    def _split_keywords(keyword_string):
        """Split a string of keywords into a nice clean list (PRIVATE)."""
        # process the keywords into a python list
        if keyword_string == "" or keyword_string == ".":
            keywords = ""
        elif keyword_string[-1] == ".":
            keywords = keyword_string[:-1]
        else:
            keywords = keyword_string
        keyword_list = keywords.split(";")
        return [x.strip() for x in keyword_list]

    @staticmethod
    def _split_accessions(accession_string):
        """Split a string of accession numbers into a list (PRIVATE)."""
        # first replace all line feeds with spaces
        # Also, EMBL style accessions are split with ';'
        accession = accession_string.replace("\n", " ").replace(";", " ")

        return [x.strip() for x in accession.split() if x.strip()]

    @staticmethod
    def _split_taxonomy(taxonomy_string):
        """Split a string with taxonomy info into a list (PRIVATE)."""
        if not taxonomy_string or taxonomy_string == ".":
            # Missing data, no taxonomy
            return []

        if taxonomy_string[-1] == ".":
            tax_info = taxonomy_string[:-1]
        else:
            tax_info = taxonomy_string
        tax_list = tax_info.split(";")
        new_tax_list = []
        for tax_item in tax_list:
            new_items = tax_item.split("\n")
            new_tax_list.extend(new_items)
        while "" in new_tax_list:
            new_tax_list.remove("")
        return [x.strip() for x in new_tax_list]

    @staticmethod
    def _clean_location(location_string):
        """Clean whitespace out of a location string (PRIVATE).

        The location parser isn't a fan of whitespace, so we clean it out
        before feeding it into the parser.
        """
        # Originally this imported string.whitespace and did a replace
        # via a loop.  It's simpler to just split on whitespace and rejoin
        # the string - and this avoids importing string too.  See Bug 2684.
        return "".join(location_string.split())

    @staticmethod
    def _remove_newlines(text):
        """Remove any newlines in the passed text, returning the new string (PRIVATE)."""
        # get rid of newlines in the qualifier value
        newlines = ["\n", "\r"]
        for ws in newlines:
            text = text.replace(ws, "")

        return text

    @staticmethod
    def _normalize_spaces(text):
        """Replace multiple spaces in the passed text with single spaces (PRIVATE)."""
        # get rid of excessive spaces
        return " ".join(x for x in text.split(" ") if x)

    @staticmethod
    def _remove_spaces(text):
        """Remove all spaces from the passed text (PRIVATE)."""
        return text.replace(" ", "")

    @staticmethod
    def _convert_to_python_numbers(start, end):
        """Convert a start and end range to python notation (PRIVATE).

        In GenBank, starts and ends are defined in "biological" coordinates,
        where 1 is the first base and [i, j] means to include both i and j.

        In python, 0 is the first base and [i, j] means to include i, but
        not j.

        So, to convert "biological" to python coordinates, we need to
        subtract 1 from the start, and leave the end and things should
        be converted happily.
        """
        new_start = start - 1
        new_end = end

        return new_start, new_end


class _FeatureConsumer(_BaseGenBankConsumer):
    """Create a SeqRecord object with Features to return (PRIVATE).

    Attributes:
     - use_fuzziness - specify whether or not to parse with fuzziness in
       feature locations.
     - feature_cleaner - a class that will be used to provide specialized
       cleaning-up of feature values.

    """

    def __init__(self, use_fuzziness, feature_cleaner=None):
        from Bio.SeqRecord import SeqRecord

        _BaseGenBankConsumer.__init__(self)
        self.data = SeqRecord(None, id=None)
        self.data.id = None
        self.data.description = ""

        self._use_fuzziness = use_fuzziness
        self._feature_cleaner = feature_cleaner

        self._seq_type = ""
        self._seq_data = []
        self._cur_reference = None
        self._cur_feature = None
        self._expected_size = None

    def locus(self, locus_name):
        """Set the locus name is set as the name of the Sequence."""
        self.data.name = locus_name

    def size(self, content):
        """Record the sequence length."""
        self._expected_size = int(content)

    def residue_type(self, type):
        """Record the sequence type (SEMI-OBSOLETE).

        This reflects the fact that the topology (linear/circular) and
        molecule type (e.g. DNA vs RNA) were a single field in early
        files. Current GenBank/EMBL files have two fields.
        """
        self._seq_type = type.strip()

    def topology(self, topology):
        """Validate and record sequence topology.

        The topology argument should be "linear" or "circular" (string).
        """
        if topology:
            if topology not in ["linear", "circular"]:
                raise ParserFailureError(
                    f"Unexpected topology {topology!r} should be linear or circular"
                )
            self.data.annotations["topology"] = topology

    def molecule_type(self, mol_type):
        """Validate and record the molecule type (for round-trip etc)."""
        if mol_type:
            if "circular" in mol_type or "linear" in mol_type:
                raise ParserFailureError(
                    f"Molecule type {mol_type!r} should not include topology"
                )

            # Writing out records will fail if we have a lower case DNA
            # or RNA string in here, so upper case it.
            # This is a bit ugly, but we don't want to upper case e.g.
            # the m in mRNA, but thanks to the strip we lost the spaces
            # so we need to index from the back
            if mol_type[-3:].upper() in ("DNA", "RNA") and not mol_type[-3:].isupper():
                warnings.warn(
                    f"Non-upper case molecule type in LOCUS line: {mol_type}",
                    BiopythonParserWarning,
                )

            self.data.annotations["molecule_type"] = mol_type

    def data_file_division(self, division):
        self.data.annotations["data_file_division"] = division

    def date(self, submit_date):
        self.data.annotations["date"] = submit_date

    def definition(self, definition):
        """Set the definition as the description of the sequence."""
        if self.data.description:
            # Append to any existing description
            # e.g. EMBL files with two DE lines.
            self.data.description += " " + definition
        else:
            self.data.description = definition

    def accession(self, acc_num):
        """Set the accession number as the id of the sequence.

        If we have multiple accession numbers, the first one passed is
        used.
        """
        new_acc_nums = self._split_accessions(acc_num)

        # Also record them ALL in the annotations
        try:
            # On the off chance there was more than one accession line:
            for acc in new_acc_nums:
                # Prevent repeat entries
                if acc not in self.data.annotations["accessions"]:
                    self.data.annotations["accessions"].append(acc)
        except KeyError:
            self.data.annotations["accessions"] = new_acc_nums

        # if we haven't set the id information yet, add the first acc num
        if not self.data.id:
            if len(new_acc_nums) > 0:
                # self.data.id = new_acc_nums[0]
                # Use the FIRST accession as the ID, not the first on this line!
                self.data.id = self.data.annotations["accessions"][0]

    def tls(self, content):
        self.data.annotations["tls"] = content.split("-")

    def tsa(self, content):
        self.data.annotations["tsa"] = content.split("-")

    def wgs(self, content):
        self.data.annotations["wgs"] = content.split("-")

    def add_wgs_scafld(self, content):
        self.data.annotations.setdefault("wgs_scafld", []).append(content.split("-"))

    def nid(self, content):
        self.data.annotations["nid"] = content

    def pid(self, content):
        self.data.annotations["pid"] = content

    def version(self, version_id):
        # Want to use the versioned accession as the record.id
        # This comes from the VERSION line in GenBank files, or the
        # obsolete SV line in EMBL.  For the new EMBL files we need
        # both the version suffix from the ID line and the accession
        # from the AC line.
        if version_id.count(".") == 1 and version_id.split(".")[1].isdigit():
            self.accession(version_id.split(".")[0])
            self.version_suffix(version_id.split(".")[1])
        elif version_id:
            # For backwards compatibility...
            self.data.id = version_id

    def project(self, content):
        """Handle the information from the PROJECT line as a list of projects.

        e.g.::

            PROJECT     GenomeProject:28471

        or::

            PROJECT     GenomeProject:13543  GenomeProject:99999

        This is stored as dbxrefs in the SeqRecord to be consistent with the
        projected switch of this line to DBLINK in future GenBank versions.
        Note the NCBI plan to replace "GenomeProject:28471" with the shorter
        "Project:28471" as part of this transition.
        """
        content = content.replace("GenomeProject:", "Project:")
        self.data.dbxrefs.extend(p for p in content.split() if p)

    def dblink(self, content):
        """Store DBLINK cross references as dbxrefs in our record object.

        This line type is expected to replace the PROJECT line in 2009. e.g.

        During transition::

            PROJECT     GenomeProject:28471
            DBLINK      Project:28471
                        Trace Assembly Archive:123456

        Once the project line is dropped::

            DBLINK      Project:28471
                        Trace Assembly Archive:123456

        Note GenomeProject -> Project.

        We'll have to see some real examples to be sure, but based on the
        above example we can expect one reference per line.

        Note that at some point the NCBI have included an extra space, e.g.::

            DBLINK      Project: 28471

        """
        # During the transition period with both PROJECT and DBLINK lines,
        # we don't want to add the same cross reference twice.
        while ": " in content:
            content = content.replace(": ", ":")
        if content.strip() not in self.data.dbxrefs:
            self.data.dbxrefs.append(content.strip())

    def version_suffix(self, version):
        """Set the version to overwrite the id.

        Since the version provides the same information as the accession
        number, plus some extra info, we set this as the id if we have
        a version.
        """
        # e.g. GenBank line:
        # VERSION     U49845.1  GI:1293613
        # or the obsolete EMBL line:
        # SV   U49845.1
        # Scanner calls consumer.version("U49845.1")
        # which then calls consumer.version_suffix(1)
        #
        # e.g. EMBL new line:
        # ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
        # Scanner calls consumer.version_suffix(1)
        assert version.isdigit()
        self.data.annotations["sequence_version"] = int(version)

    def db_source(self, content):
        self.data.annotations["db_source"] = content.rstrip()

    def gi(self, content):
        self.data.annotations["gi"] = content

    def keywords(self, content):
        if "keywords" in self.data.annotations:
            # Multi-line keywords, append to list
            # Note EMBL states "A keyword is never split between lines."
            self.data.annotations["keywords"].extend(self._split_keywords(content))
        else:
            self.data.annotations["keywords"] = self._split_keywords(content)

    def segment(self, content):
        self.data.annotations["segment"] = content

    def source(self, content):
        # Note that some software (e.g. VectorNTI) may produce an empty
        # source (rather than using a dot/period as might be expected).
        if content == "":
            source_info = ""
        elif content[-1] == ".":
            source_info = content[:-1]
        else:
            source_info = content
        self.data.annotations["source"] = source_info

    def organism(self, content):
        self.data.annotations["organism"] = content

    def taxonomy(self, content):
        """Record (another line of) the taxonomy lineage."""
        lineage = self._split_taxonomy(content)
        try:
            self.data.annotations["taxonomy"].extend(lineage)
        except KeyError:
            self.data.annotations["taxonomy"] = lineage

    def reference_num(self, content):
        """Signal the beginning of a new reference object."""
        # if we have a current reference that hasn't been added to
        # the list of references, add it.
        if self._cur_reference is not None:
            self.data.annotations["references"].append(self._cur_reference)
        else:
            self.data.annotations["references"] = []

        self._cur_reference = Reference()

    def reference_bases(self, content):
        """Attempt to determine the sequence region the reference entails.

        Possible types of information we may have to deal with:

        (bases 1 to 86436)
        (sites)
        (bases 1 to 105654; 110423 to 111122)
        1  (residues 1 to 182)
        """
        # first remove the parentheses
        assert content.endswith(")"), content
        ref_base_info = content[1:-1]

        all_locations = []
        # parse if we've got 'bases' and 'to'
        if "bases" in ref_base_info and "to" in ref_base_info:
            # get rid of the beginning 'bases'
            ref_base_info = ref_base_info[5:]
            locations = self._split_reference_locations(ref_base_info)
            all_locations.extend(locations)
        elif "residues" in ref_base_info and "to" in ref_base_info:
            residues_start = ref_base_info.find("residues")
            # get only the information after "residues"
            ref_base_info = ref_base_info[(residues_start + len("residues ")) :]
            locations = self._split_reference_locations(ref_base_info)
            all_locations.extend(locations)

        # make sure if we are not finding information then we have
        # the string 'sites' or the string 'bases'
        elif ref_base_info == "sites" or ref_base_info.strip() == "bases":
            pass
        # otherwise raise an error
        else:
            raise ValueError(
                f"Could not parse base info {ref_base_info} in record {self.data.id}"
            )

        self._cur_reference.location = all_locations

    def _split_reference_locations(self, location_string):
        """Get reference locations out of a string of reference information (PRIVATE).

        The passed string should be of the form::

            1 to 20; 20 to 100

        This splits the information out and returns a list of location objects
        based on the reference locations.
        """
        # split possibly multiple locations using the ';'
        all_base_info = location_string.split(";")

        new_locations = []
        for base_info in all_base_info:
            start, end = base_info.split("to")
            new_start, new_end = self._convert_to_python_numbers(
                int(start.strip()), int(end.strip())
            )
            this_location = SimpleLocation(new_start, new_end)
            new_locations.append(this_location)
        return new_locations

    def authors(self, content):
        if self._cur_reference.authors:
            self._cur_reference.authors += " " + content
        else:
            self._cur_reference.authors = content

    def consrtm(self, content):
        if self._cur_reference.consrtm:
            self._cur_reference.consrtm += " " + content
        else:
            self._cur_reference.consrtm = content

    def title(self, content):
        if self._cur_reference is None:
            warnings.warn(
                "GenBank TITLE line without REFERENCE line.", BiopythonParserWarning
            )
        elif self._cur_reference.title:
            self._cur_reference.title += " " + content
        else:
            self._cur_reference.title = content

    def journal(self, content):
        if self._cur_reference.journal:
            self._cur_reference.journal += " " + content
        else:
            self._cur_reference.journal = content

    def medline_id(self, content):
        self._cur_reference.medline_id = content

    def pubmed_id(self, content):
        self._cur_reference.pubmed_id = content

    def remark(self, content):
        """Deal with a reference comment."""
        if self._cur_reference.comment:
            self._cur_reference.comment += " " + content
        else:
            self._cur_reference.comment = content

    def comment(self, content):
        try:
            self.data.annotations["comment"] += "\n" + "\n".join(content)
        except KeyError:
            self.data.annotations["comment"] = "\n".join(content)

    def structured_comment(self, content):
        self.data.annotations["structured_comment"] = content

    def features_line(self, content):
        """Get ready for the feature table when we reach the FEATURE line."""
        self.start_feature_table()

    def start_feature_table(self):
        """Indicate we've got to the start of the feature table."""
        # make sure we've added on our last reference object
        if self._cur_reference is not None:
            self.data.annotations["references"].append(self._cur_reference)
            self._cur_reference = None

    def feature_key(self, content):
        # start a new feature
        self._cur_feature = SeqFeature()
        self._cur_feature.type = content
        self.data.features.append(self._cur_feature)

    def location(self, content):
        """Parse out location information from the location string.

        This uses simple Python code with some regular expressions to do the
        parsing, and then translates the results into appropriate objects.
        """
        # clean up newlines and other whitespace inside the location before
        # parsing - locations should have no whitespace whatsoever
        location_line = self._clean_location(content)

        # Older records have junk like replace(266,"c") in the
        # location line. Newer records just replace this with
        # the number 266 and have the information in a more reasonable
        # place. So we'll just grab out the number and feed this to the
        # parser. We shouldn't really be losing any info this way.
        if "replace" in location_line:
            comma_pos = location_line.find(",")
            location_line = location_line[8:comma_pos]

        length = self._expected_size
        # Check if the sequence is circular for features that span the origin
        is_circular = "circular" in self.data.annotations.get("topology", "").lower()
        stranded = "PROTEIN" not in self._seq_type.upper()

        try:
            location = Location.fromstring(location_line, length, is_circular, stranded)
        except LocationParserError as e:
            warnings.warn(
                f"{e}; setting feature location to None.", BiopythonParserWarning
            )
            location = None
        self._cur_feature.location = location

    def feature_qualifier(self, key, value):
        """When we get a qualifier key and its value.

        Can receive None, since you can have valueless keys such as /pseudo
        """
        # Hack to try to preserve historical behaviour of /pseudo etc
        if value is None:
            # if the key doesn't exist yet, add an empty string
            if key not in self._cur_feature.qualifiers:
                self._cur_feature.qualifiers[key] = [""]
                return
            # otherwise just skip this key
            return

        # Remove enclosing quotation marks
        if len(value) > 1 and value[0] == '"' and value[-1] == '"':
            value = value[1:-1]

        # Handle NCBI escaping
        # Warn if escaping is not according to standard
        if re.search(r'[^"]"[^"]|^"[^"]|[^"]"$', value):
            warnings.warn(
                'The NCBI states double-quote characters like " should be escaped as "" '
                "(two double - quotes), but here it was not: %r" % value,
                BiopythonParserWarning,
            )
        # Undo escaping, repeated double quotes -> one double quote
        value = value.replace('""', '"')

        if self._feature_cleaner is not None:
            value = self._feature_cleaner.clean_value(key, value)

        # if the qualifier name exists, append the value
        if key in self._cur_feature.qualifiers:
            self._cur_feature.qualifiers[key].append(value)
        # otherwise start a new list of the key with its values
        else:
            self._cur_feature.qualifiers[key] = [value]

    def feature_qualifier_name(self, content_list):
        """Use feature_qualifier instead (OBSOLETE)."""
        raise NotImplementedError("Use the feature_qualifier method instead.")

    def feature_qualifier_description(self, content):
        """Use feature_qualifier instead (OBSOLETE)."""
        raise NotImplementedError("Use the feature_qualifier method instead.")

    def contig_location(self, content):
        """Deal with CONTIG information."""
        # Historically this was stored as a SeqFeature object, but it was
        # stored under record.annotations["contig"] and not under
        # record.features with the other SeqFeature objects.
        #
        # The CONTIG location line can include additional tokens like
        # Gap(), Gap(100) or Gap(unk100) which are not used in the feature
        # location lines, so storing it using SeqFeature based location
        # objects is difficult.
        #
        # We now store this a string, which means for BioSQL we are now in
        # much better agreement with how BioPerl records the CONTIG line
        # in the database.
        #
        # NOTE - This code assumes the scanner will return all the CONTIG
        # lines already combined into one long string!
        self.data.annotations["contig"] = content

    def origin_name(self, content):
        pass

    def base_count(self, content):
        pass

    def base_number(self, content):
        pass

    def sequence(self, content):
        """Add up sequence information as we get it.

        To try and make things speedier, this puts all of the strings
        into a list of strings, and then uses string.join later to put
        them together. Supposedly, this is a big time savings
        """
        assert " " not in content
        self._seq_data.append(content.upper())

    def record_end(self, content):
        """Clean up when we've finished the record."""
        # Try and append the version number to the accession for the full id
        if not self.data.id:
            if "accessions" in self.data.annotations:
                raise ValueError(
                    "Problem adding version number to accession: "
                    + str(self.data.annotations["accessions"])
                )
            self.data.id = self.data.name  # Good fall back?
        elif self.data.id.count(".") == 0:
            try:
                self.data.id += ".%i" % self.data.annotations["sequence_version"]
            except KeyError:
                pass

        # add the sequence information

        sequence = "".join(self._seq_data)

        if (
            self._expected_size is not None
            and len(sequence) != 0
            and self._expected_size != len(sequence)
        ):
            warnings.warn(
                "Expected sequence length %i, found %i (%s)."
                % (self._expected_size, len(sequence), self.data.id),
                BiopythonParserWarning,
            )

        molecule_type = None
        if self._seq_type:
            # mRNA is really also DNA, since it is actually cDNA
            if "DNA" in self._seq_type.upper() or "MRNA" in self._seq_type.upper():
                molecule_type = "DNA"
            # are there ever really RNA sequences in GenBank?
            elif "RNA" in self._seq_type.upper():
                # Even for data which was from RNA, the sequence string
                # is usually given as DNA (T not U).  Bug 3010
                molecule_type = "RNA"
            elif (
                "PROTEIN" in self._seq_type.upper() or self._seq_type == "PRT"
            ):  # PRT is used in EMBL-bank for patents
                molecule_type = "protein"
            # work around ugly GenBank records which have circular or
            # linear but no indication of sequence type
            elif self._seq_type in ["circular", "linear", "unspecified"]:
                pass
            # we have a bug if we get here
            else:
                raise ValueError(
                    f"Could not determine molecule_type for seq_type {self._seq_type}"
                )
        # Don't overwrite molecule_type
        if molecule_type is not None:
            self.data.annotations["molecule_type"] = self.data.annotations.get(
                "molecule_type", molecule_type
            )
        if not sequence and self._expected_size:
            self.data.seq = Seq(None, length=self._expected_size)
        else:
            self.data.seq = Seq(sequence)


class _RecordConsumer(_BaseGenBankConsumer):
    """Create a GenBank Record object from scanner generated information (PRIVATE)."""

    def __init__(self):
        _BaseGenBankConsumer.__init__(self)
        from . import Record

        self.data = Record.Record()

        self._seq_data = []
        self._cur_reference = None
        self._cur_feature = None
        self._cur_qualifier = None

    def tls(self, content):
        self.data.tls = content.split("-")

    def tsa(self, content):
        self.data.tsa = content.split("-")

    def wgs(self, content):
        self.data.wgs = content.split("-")

    def add_wgs_scafld(self, content):
        self.data.wgs_scafld.append(content.split("-"))

    def locus(self, content):
        self.data.locus = content

    def size(self, content):
        self.data.size = content

    def residue_type(self, content):
        # Be lenient about parsing, but technically lowercase residue types are malformed.
        if "dna" in content or "rna" in content:
            warnings.warn(
                f"Invalid seq_type ({content}): DNA/RNA should be uppercase.",
                BiopythonParserWarning,
            )
        self.data.residue_type = content

    def data_file_division(self, content):
        self.data.data_file_division = content

    def date(self, content):
        self.data.date = content

    def definition(self, content):
        self.data.definition = content

    def accession(self, content):
        for acc in self._split_accessions(content):
            if acc not in self.data.accession:
                self.data.accession.append(acc)

    def molecule_type(self, mol_type):
        """Validate and record the molecule type (for round-trip etc)."""
        if mol_type:
            if "circular" in mol_type or "linear" in mol_type:
                raise ParserFailureError(
                    f"Molecule type {mol_type!r} should not include topology"
                )

            # Writing out records will fail if we have a lower case DNA
            # or RNA string in here, so upper case it.
            # This is a bit ugly, but we don't want to upper case e.g.
            # the m in mRNA, but thanks to the strip we lost the spaces
            # so we need to index from the back
            if mol_type[-3:].upper() in ("DNA", "RNA") and not mol_type[-3:].isupper():
                warnings.warn(
                    f"Non-upper case molecule type in LOCUS line: {mol_type}",
                    BiopythonParserWarning,
                )

            self.data.molecule_type = mol_type

    def topology(self, topology):
        """Validate and record sequence topology.

        The topology argument should be "linear" or "circular" (string).
        """
        if topology:
            if topology not in ["linear", "circular"]:
                raise ParserFailureError(
                    f"Unexpected topology {topology!r} should be linear or circular"
                )
            self.data.topology = topology

    def nid(self, content):
        self.data.nid = content

    def pid(self, content):
        self.data.pid = content

    def version(self, content):
        self.data.version = content

    def db_source(self, content):
        self.data.db_source = content.rstrip()

    def gi(self, content):
        self.data.gi = content

    def keywords(self, content):
        self.data.keywords = self._split_keywords(content)

    def project(self, content):
        self.data.projects.extend(p for p in content.split() if p)

    def dblink(self, content):
        self.data.dblinks.append(content)

    def segment(self, content):
        self.data.segment = content

    def source(self, content):
        self.data.source = content

    def organism(self, content):
        self.data.organism = content

    def taxonomy(self, content):
        self.data.taxonomy = self._split_taxonomy(content)

    def reference_num(self, content):
        """Grab the reference number and signal the start of a new reference."""
        # check if we have a reference to add
        if self._cur_reference is not None:
            self.data.references.append(self._cur_reference)

        from . import Record

        self._cur_reference = Record.Reference()
        self._cur_reference.number = content

    def reference_bases(self, content):
        self._cur_reference.bases = content

    def authors(self, content):
        self._cur_reference.authors = content

    def consrtm(self, content):
        self._cur_reference.consrtm = content

    def title(self, content):
        if self._cur_reference is None:
            warnings.warn(
                "GenBank TITLE line without REFERENCE line.", BiopythonParserWarning
            )
            return
        self._cur_reference.title = content

    def journal(self, content):
        self._cur_reference.journal = content

    def medline_id(self, content):
        self._cur_reference.medline_id = content

    def pubmed_id(self, content):
        self._cur_reference.pubmed_id = content

    def remark(self, content):
        self._cur_reference.remark = content

    def comment(self, content):
        self.data.comment += "\n".join(content)

    def structured_comment(self, content):
        self.data.structured_comment = content

    def primary_ref_line(self, content):
        """Save reference data for the PRIMARY line."""
        self.data.primary.append(content)

    def primary(self, content):
        pass

    def features_line(self, content):
        """Get ready for the feature table when we reach the FEATURE line."""
        self.start_feature_table()

    def start_feature_table(self):
        """Signal the start of the feature table."""
        # we need to add on the last reference
        if self._cur_reference is not None:
            self.data.references.append(self._cur_reference)

    def feature_key(self, content):
        """Grab the key of the feature and signal the start of a new feature."""
        # first add on feature information if we've got any
        self._add_feature()

        from . import Record

        self._cur_feature = Record.Feature()
        self._cur_feature.key = content

    def _add_feature(self):
        """Add a feature to the record, with relevant checks (PRIVATE).

        This does all of the appropriate checking to make sure we haven't
        left any info behind, and that we are only adding info if it
        exists.
        """
        if self._cur_feature is not None:
            # if we have a left over qualifier, add it to the qualifiers
            # on the current feature
            if self._cur_qualifier is not None:
                self._cur_feature.qualifiers.append(self._cur_qualifier)

            self._cur_qualifier = None
            self.data.features.append(self._cur_feature)

    def location(self, content):
        self._cur_feature.location = self._clean_location(content)

    def feature_qualifier(self, key, value):
        self.feature_qualifier_name([key])
        if value is not None:
            self.feature_qualifier_description(value)

    def feature_qualifier_name(self, content_list):
        """Deal with qualifier names.

        We receive a list of keys, since you can have valueless keys such as
        /pseudo which would be passed in with the next key (since no other
        tags separate them in the file)
        """
        from . import Record

        for content in content_list:
            # the record parser keeps the /s -- add them if we don't have 'em
            if not content.startswith("/"):
                content = f"/{content}"
            # add on a qualifier if we've got one
            if self._cur_qualifier is not None:
                self._cur_feature.qualifiers.append(self._cur_qualifier)

            self._cur_qualifier = Record.Qualifier()
            self._cur_qualifier.key = content

    def feature_qualifier_description(self, content):
        # if we have info then the qualifier key should have a ='s
        if "=" not in self._cur_qualifier.key:
            self._cur_qualifier.key = f"{self._cur_qualifier.key}="
        cur_content = self._remove_newlines(content)
        # remove all spaces from the value if it is a type where spaces
        # are not important
        for remove_space_key in self.__class__.remove_space_keys:
            if remove_space_key in self._cur_qualifier.key:
                cur_content = self._remove_spaces(cur_content)
        self._cur_qualifier.value = self._normalize_spaces(cur_content)

    def base_count(self, content):
        self.data.base_counts = content

    def origin_name(self, content):
        self.data.origin = content

    def contig_location(self, content):
        """Signal that we have contig information to add to the record."""
        self.data.contig = self._clean_location(content)

    def sequence(self, content):
        """Add sequence information to a list of sequence strings.

        This removes spaces in the data and uppercases the sequence, and
        then adds it to a list of sequences. Later on we'll join this
        list together to make the final sequence. This is faster than
        adding on the new string every time.
        """
        assert " " not in content
        self._seq_data.append(content.upper())

    def record_end(self, content):
        """Signal the end of the record and do any necessary clean-up."""
        # add together all of the sequence parts to create the
        # final sequence string
        self.data.sequence = "".join(self._seq_data)
        # add on the last feature
        self._add_feature()


def parse(handle):
    """Iterate over GenBank formatted entries as Record objects.

    >>> from Bio import GenBank
    >>> with open("GenBank/NC_000932.gb") as handle:
    ...     for record in GenBank.parse(handle):
    ...         print(record.accession)
    ['NC_000932']

    To get SeqRecord objects use Bio.SeqIO.parse(..., format="gb")
    instead.
    """
    return iter(Iterator(handle, RecordParser()))


def read(handle):
    """Read a handle containing a single GenBank entry as a Record object.

    >>> from Bio import GenBank
    >>> with open("GenBank/NC_000932.gb") as handle:
    ...     record = GenBank.read(handle)
    ...     print(record.accession)
    ['NC_000932']

    To get a SeqRecord object use Bio.SeqIO.read(..., format="gb")
    instead.
    """
    iterator = parse(handle)
    try:
        record = next(iterator)
    except StopIteration:
        raise ValueError("No records found in handle") from None
    try:
        next(iterator)
        raise ValueError("More than one record found in handle")
    except StopIteration:
        pass
    return record


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
