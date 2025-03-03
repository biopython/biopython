# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Bio.SearchIO base classes for Infernal-related code."""


from Bio.SearchIO._model import Hit


class _BaseInfernalParser:
    """Abstract class for Infernal parser objects."""

    def _add_hit_to_dict(self, hit_attrs, hsp, hit_dict):
        """Add a hit information to the hit container (PRIVATE)."""
        # add the hit to the container dict, creating a new
        # entry if it does not exist yet
        hid = hit_attrs["id"]
        if hid not in hit_dict:
            hit_dict[hid] = {"attrs": hit_attrs, "hsps": []}
        else:
            assert hit_dict[hid]["attrs"]["query_id"] == hit_attrs["query_id"]
            assert hit_dict[hid]["attrs"]["description"] == hit_attrs["description"]

        hit_dict[hid]["hsps"].append(hsp)

    def _hit_to_list(self, hit_dict):
        """Create a Hit list from the hit container (PRIVATE)."""

        hit_list = []
        for hit_info in hit_dict.values():
            hit = Hit(hit_info["hsps"])
            for attr, value in hit_info["attrs"].items():
                setattr(hit, attr, value)
            hit_list.append(hit)
        return hit_list
