"""Bio.SeqIO support for GenBank's 5-column tab-delimited feature table

This format is used when preparing GenBank submissions. It contains features and
no actual sequences.

Documentation: https://www.ncbi.nlm.nih.gov/genbank/feature_table/
"""

import re
import warnings

from Bio import BiopythonWarning
from Bio.File import as_handle
from Bio.SeqFeature import Reference, Position, SimpleLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from .Interfaces import SequenceWriter


def _make_simple_location(start_loc_str, end_loc_str, offset):
    """Constructs a SimpleLocation from start location and end location strings (PRIVATE)."""
    start_loc = Position.fromstring(start_loc_str, -1 + offset)
    end_loc = Position.fromstring(end_loc_str, offset)
    if start_loc < end_loc:
        return SimpleLocation(start_loc, end_loc, strand=1)
    else:
        return SimpleLocation(end_loc - 1, start_loc + 1, strand=-1)


def FeatureTableIterator(source):
    """Parser for GenBank's 5-column tab-delimited feature table"""
    with as_handle(source) as handle:
        rec = None
        offset = 0

        for line in handle:
            if line == "\n":
                continue

            if line.startswith(">Feature"):
                # First line of a new table
                if rec is not None:
                    offset = 0
                    yield rec

                parts = line.strip().split(" ")
                if len(parts) < 2:
                    raise ValueError(
                        f"First line of table has no sequence identifier: {line}"
                    )
                rec = SeqRecord(None, id=parts[1], name=parts[1])
            elif line.startswith("[offset="):
                try:
                    offset = re.findall(r"\d+", "[offset=2000]")[0]
                except IndexError:
                    raise ValueError(f"Malformed offset: {line}")
            elif not line.startswith("\t"):
                # New feature in the current table
                if rec is None:
                    raise ValueError(f"Feature found outside of a table: {line}")

                parts = line.strip().split("\t")
                if len(parts) == 3:
                    start_loc, end_loc, feature_key = parts
                    loc = _make_simple_location(start_loc, end_loc, offset)

                    if feature_key == "REFERENCE":
                        # REFERENCEs are annotations, not features, so they need
                        # special handling
                        try:
                            ref_line = next(handle)
                        except StopIteration:
                            raise ValueError(
                                f"Expected qualifier key-value pair line to follow REFERENCE: {line}"
                            )
                        parts = ref_line.strip().split("\t")
                        if len(parts) != 2:
                            raise ValueError(
                                f"Expected qualifier key-value pair: {ref_line}"
                            )
                        if parts[0] != "PubMed":
                            warnings.warn(
                                f"Reference qualifier key is not PubMed: {ref_line}",
                                BiopythonWarning,
                            )

                        ref = Reference()
                        ref.location = loc
                        if parts[0] == "PubMed":
                            ref.pubmed_id = parts[1]
                        else:
                            ref.title = f"{parts[0]} {parts[1]}"
                        if "references" not in rec.annotations:
                            rec.annotations["references"] = [ref]
                        else:
                            rec.annotations["references"].append(ref)
                    else:
                        rec.features.append(SeqFeature(loc, type=feature_key))
                    pass
                elif len(parts) == 2:
                    # add new location to current feature
                    start_loc, end_loc = parts
                    rec.features[-1].location += _make_simple_location(
                        start_loc, end_loc, offset
                    )
                else:
                    raise ValueError(f"Malformed first line of feature: {line}")
            else:
                # New qualifier in the current feature
                parts = line.strip().split("\t")
                if len(parts) != 2:
                    raise ValueError(f"Expected qualifier key-value pair: {line}")
                key, value = parts

                if key not in rec.features[-1].qualifiers:
                    rec.features[-1].qualifiers[key] = [value]
                else:
                    rec.features[-1].qualifiers[key].append(value)

        if rec is not None:
            yield rec


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
                BiopythonWarning,
            )
            return f"{loc.start + 1}\t{loc.end}"

    def _write_feature_header(self, loc, name):
        """Write the header of a feature/references (PRIVATE)."""
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
                self._write_feature_header(ref.location, "REFERENCE")
                if ref.pubmed_id != "":
                    self.handle.write(f"\t\t\tPubMed\t{self.clean(ref.pubmed_id)}\n")
                elif ref.medline_id != "":
                    warnings.warn(
                        "Reference does not have a PubMed ID. Using the Medline"
                        " ID, but 5-column tab-delimited feature tables only"
                        " officially support PubMed IDs.",
                        BiopythonWarning,
                    )
                    self.handle.write(f"\t\t\tMedline\t{self.clean(ref.medline_id)}\n")
                else:
                    raise TypeError("Reference has no medline_id or pubmed_id.")

        for feature in rec.features:
            self._write_feature_header(feature.location, feature.type)
            for qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    self.handle.write(f"\t\t\t{self.clean(qualifier)}\t{str(value)}\n")
