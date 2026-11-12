"""Bio.SeqIO support for GenBank's 5-column tab-delimited feature table

This format is used when preparing GenBank submissions. It contains features and
no actual sequences.

Documentation:
    - https://www.ncbi.nlm.nih.gov/genbank/table2asn/
    - https://www.ncbi.nlm.nih.gov/genbank/feature_table/
"""

import re
import warnings

from Bio import BiopythonWarning, BiopythonParserWarning
from Bio.File import as_handle
from Bio.Seq import Seq
from Bio.SeqFeature import Reference, Position, SimpleLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator, SequenceWriter


def _split_by_tabs(line):
    """Split a line by tab characters or at least 2 consecutive spaces (PRIVATE)."""
    if "\t" in line:
        return line.split("\t")

    # Some files use spaces for identation
    parts = [""]
    space_ctr = 0
    for ch in line:
        if ch == " ":
            parts[-1] += ch
            space_ctr += 1
        else:
            if space_ctr >= 2:
                parts[-1] = parts[-1][:-space_ctr]
                parts.append(ch)
            else:
                parts[-1] += ch
            space_ctr = 0
    return parts


def _make_simple_location(start_loc_str, end_loc_str, offset):
    """Constructs a SimpleLocation from start location and end location strings (PRIVATE)."""
    start_loc = Position.fromstring(start_loc_str)
    end_loc = Position.fromstring(end_loc_str)
    if start_loc < end_loc:
        return SimpleLocation(start_loc - 1 + offset, end_loc + offset, strand=1)
    else:
        return SimpleLocation(end_loc - 1 + offset, start_loc + offset, strand=-1)


def _fix_references(rec):
    """Convert REFERENCE features from SeqFeature to annotations (PRIVATE)."""
    ref_features = [feature for feature in rec.features if feature.type == "REFERENCE"]
    for feature in ref_features:
        rec.features.remove(feature)

    refs = []
    for ref_feature in ref_features:
        ref = Reference()
        ref.location = ref_feature.location.parts
        if "PubMed" in ref_feature.qualifiers:
            if len(ref_feature.qualifiers["PubMed"]) > 1:
                raise ValueError(f"Reference has multiple PubMed IDs: {ref_feature}")
            ref.pubmed_id = ref_feature.qualifiers["PubMed"][0]
        if "Medline" in ref_feature.qualifiers:
            warnings.warn(
                f"Reference has Medline ID, but only PubMed IDs are officially supported in .tbl files: {ref_feature}",
                BiopythonParserWarning,
            )
            if len(ref_feature.qualifiers["Medline"]) > 1:
                raise ValueError(f"Reference has multiple Medline IDs: {ref_feature}")
            ref.medline_id = ref_feature.qualifiers["Medline"][0]
        if ref.pubmed_id == "" and ref.medline_id == "":
            warnings.warn(
                f"Ignoring reference without PubMed or Medline ID: {ref_feature}",
                BiopythonParserWarning,
            )
            continue
        refs.append(ref)
    if len(refs) > 0:
        rec.annotations["references"] = refs


def _set_seq_len(rec):
    """If there is only one reference, set the sequence length to the length of the reference (PRIVATE)."""
    if "references" in rec.annotations and len(rec.annotations["references"]) == 1:
        ref = rec.annotations["references"][0]
        rec.seq = Seq(None, length=len(sum(ref.location)))
    pass


def _add_codon_start(rec):
    """Add codon_start=1 to CDS features without a codon_start qualifier (PRIVATE)."""
    for feature in rec.features:
        if feature.type == "CDS" and "codon_start" not in feature.qualifiers:
            feature.qualifiers["codon_start"] = ["1"]


def _find_overlap(rec):
    """Add gene qualifiers to features without one but overlap gene features (PRIVATE)."""
    genes = [
        feature
        for feature in rec.features
        if feature.type == "gene" and "gene" in feature.qualifiers
    ]
    for feature in rec.features:
        loc = feature.location
        if feature.type in ["mRNA", "tRNA", "CDS"] and "gene" not in feature.qualifiers:
            for gene in genes:
                if loc.start in gene.location or loc.end in gene.location:
                    feature.qualifiers["gene"] = gene.qualifiers["gene"]


def _record_fixups(rec):
    """Modify the record as needed before it is yielded (PRIVATE)."""
    _fix_references(rec)
    _set_seq_len(rec)
    _add_codon_start(rec)
    _find_overlap(rec)


class FeatureTableIterator(SequenceIterator):
    """Parser for GenBank's 5-column tab-delimited feature table"""

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
    ) -> None:
        """Iterate over a .tbl file as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
        """
        super().__init__(source, fmt="5-column tab-delimited feature table")
        self.rec = None
        self.offset = 0

    def __next__(self):
        for line in self.stream:
            if line == "\n":
                continue

            if line.startswith(">Feature"):
                # First line of a new table
                parts = line.strip().split(" ")
                if len(parts) < 2:
                    raise ValueError(
                        f"First line of table has no sequence identifier: {line}"
                    )
                new_rec = SeqRecord(None, id=parts[1], name=parts[1])

                if self.rec is not None:
                    self.offset = 0
                    _record_fixups(self.rec)
                    old_rec = self.rec
                    self.rec = new_rec
                    return old_rec
                else:
                    self.rec = new_rec
            elif line.startswith("[offset="):
                matches = re.findall(r"\d+", line)
                if len(matches) != 1:
                    raise ValueError(f"Malformed offset: {line}")
                self.offset = int(matches[0])
            elif not (line.startswith("\t") or line.startswith(" ")):
                # New feature in the current table
                if self.rec is None:
                    raise ValueError(f"Feature found outside of a table: {line}")

                parts = _split_by_tabs(line.strip())
                if len(parts) == 3:
                    start_loc, end_loc, feature_key = parts
                    loc = _make_simple_location(start_loc, end_loc, self.offset)
                    self.rec.features.append(SeqFeature(loc, type=feature_key))
                elif len(parts) == 2:
                    # Add new location to current feature
                    start_loc, end_loc = parts
                    self.rec.features[-1].location += _make_simple_location(
                        start_loc, end_loc, self.offset
                    )
                else:
                    raise ValueError(f"Malformed first line of feature: {line}")
            else:
                # New qualifier in the current feature
                parts = _split_by_tabs(line.strip())
                if len(parts) != 2:
                    raise ValueError(f"Expected qualifier key-value pair: {line}")
                key, value = parts

                feature = self.rec.features[-1]
                if (
                    "product" in feature.qualifiers
                    and len(feature.qualifiers["product"]) == 1
                    and key == "product"
                ):
                    # Documentation says that a second product must be
                    # converted into a note
                    key = "note"

                if key not in feature.qualifiers:
                    feature.qualifiers[key] = [value]
                else:
                    feature.qualifiers[key].append(value)

        if self.rec is not None:
            _record_fixups(self.rec)
            rec = self.rec
            self.rec = None
            return rec
        raise StopIteration


class FeatureTableWriter(SequenceWriter):
    """Writer for GenBank's 5-column tab-delimited feature table"""

    modes = "t"

    def _simplelocation_to_string(self, loc):
        """Convert a SimpleLocation into a tab-delimited string (PRIVATE)."""
        if loc.strand == 1:
            return f"{loc.start + 1}\t{loc.end}"
        elif loc.strand == -1:
            return f"{loc.end}\t{loc.start + 1}"
        else:
            warnings.warn(
                "Feature or reference location is on mixed strands. Writing as if stand=1.",
                BiopythonParserWarning,
            )
            return f"{loc.start + 1}\t{loc.end}"

    def _write_feature_header(self, loc, name):
        """Write the header of a feature/references (PRIVATE).

        If name is "REFERENCE", then loc may be a list of SimpleLocations
        instead of one CompoundLocation, but the list must be non-empty.
        """
        if isinstance(loc, list):
            # Reference's locations is a list of locations, not one
            # CompoundLocation, so make it a CompoundLocation
            assert name == "REFERENCE"
            assert len(loc) > 0
            loc = sum(loc)

        loc_str = self._simplelocation_to_string(loc.parts[0])
        self.handle.write(f"{loc_str}\t{name}\n")
        for parent_loc in loc.parts[1:]:
            loc_str = self._simplelocation_to_string(parent_loc)
            self.handle.write(f"{loc_str}\n")

    def write_record(self, rec):
        """Write a single record to the output file."""

        seq_id = self.clean(rec.id)
        self.handle.write(f">Feature {seq_id}\n")

        if "references" in rec.annotations:
            for ref in rec.annotations["references"]:
                if ref.location == []:
                    warnings.warn(
                        f"Ignoring reference with no location: {ref}", BiopythonWarning
                    )
                    continue
                if ref.pubmed_id == "" and ref.medline_id == "":
                    warnings.warn(
                        f"Ignoring reference with no pubmed_id or medline_id: {ref}",
                        BiopythonWarning,
                    )
                    continue

                self._write_feature_header(ref.location, "REFERENCE")
                if ref.pubmed_id != "":
                    self.handle.write(f"\t\t\tPubMed\t{self.clean(ref.pubmed_id)}\n")
                elif ref.medline_id != "":
                    warnings.warn(
                        "Reference does not have a PubMed ID. Using the Medline"
                        " ID, but 5-column tab-delimited feature tables only"
                        " officially support PubMed IDs: {ref}",
                        BiopythonParserWarning,
                    )
                    self.handle.write(f"\t\t\tMedline\t{self.clean(ref.medline_id)}\n")

        for feature in rec.features:
            if feature.type == "":
                raise TypeError(f"Feature has no type: {feature}")
            if feature.location is None:
                raise TypeError(f"Feature has no location: {feature}")

            if feature.type == "source":
                # Features tables don't care about source
                continue
            self._write_feature_header(feature.location, feature.type)
            for qualifier in feature.qualifiers:
                if qualifier == "translation":
                    # Features tables don't care about translations
                    continue
                elif self.clean(qualifier).strip() == "":
                    raise TypeError(f"Qualifier is empty when cleaned: '{qualifier}'")
                for value in feature.qualifiers[qualifier]:
                    if self.clean(value).strip() == "":
                        raise TypeError(f"Value is empty when cleaned: '{value}'")
                    self.handle.write(
                        f"\t\t\t{self.clean(qualifier)}\t{self.clean(value)}\n"
                    )
