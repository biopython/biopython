# Copyright 2014 Joe Cora.
# Revisions copyright 2017 Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Objects to represent NEXUS standard data type matrix coding."""


class NexusError(Exception):
    """Provision for the management of Nexus exceptions."""


class StandardData:
    """Create a StandardData iterable object.

    Each coding specifies t [type] => (std [standard], multi [multistate] or
    uncer [uncertain]) and d [data]
    """

    def __init__(self, data):
        """Initialize the class."""
        self._data = []
        self._current_pos = 0

        # Enforce string data requirement
        if not isinstance(data, str):
            raise NexusError(
                "The coding data given to a StandardData object should be a string"
            )

        # Transfer each coding to a position within a sequence
        multi_coding = False
        uncertain_coding = False
        coding_list = {"t": "std", "d": []}

        for pos, coding in enumerate(data):
            # Check if in a multiple coded or uncertain character
            if multi_coding:
                # End multicoding if close parenthesis
                if coding == ")":
                    multi_coding = False
                else:
                    # Add current coding to list and advance to next coding
                    coding_list["d"].append(coding)
                    continue
            elif uncertain_coding:
                # End multicoding if close parenthesis
                if coding == "}":
                    uncertain_coding = False
                else:
                    # Add current coding to list and advance to next coding
                    coding_list["d"].append(coding)
                    continue
            else:
                # Check if a multiple coded or uncertain character is starting
                if coding == "(":
                    multi_coding = True
                    coding_list["t"] = "multi"
                    continue
                elif coding == "{":
                    uncertain_coding = True
                    coding_list["t"] = "uncer"
                    continue
                elif coding in [")", "}"]:
                    raise NexusError(
                        "Improper character %s at position %i of a coding sequence."
                        % (coding, pos)
                    )
                else:
                    coding_list["d"].append(coding)

            # Add character coding to data
            self._data.append(coding_list.copy())
            coding_list = {"t": "std", "d": []}

    def __len__(self):
        """Return the length of the coding, use len(my_coding)."""
        return len(self._data)

    def __getitem__(self, arg):
        """Pull out child by index."""
        return self._data[arg]

    def __iter__(self):
        """Iterate over the items."""
        return self

    def __next__(self):
        """Return next item."""
        try:
            return_coding = self._data[self._current_pos]
        except IndexError:
            self._current_pos = 0
            raise StopIteration from None
        else:
            self._current_pos += 1
            return return_coding

    def raw(self):
        """Return the full coding as a python list."""
        return self._data

    def __str__(self):
        """Return the full coding as a python string, use str(my_coding)."""
        str_return = ""
        for coding in self._data:
            if coding["t"] == "multi":
                str_return += "(" + "".join(coding["d"]) + ")"
            elif coding["t"] == "uncer":
                str_return += "{" + "".join(coding["d"]) + "}"
            else:
                str_return += coding["d"][0]
        return str_return
