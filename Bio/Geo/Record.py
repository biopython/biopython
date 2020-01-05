# Copyright 2001 by Katharine Lindner.  All rights reserved.
# Copyright 2006 by PeterC.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Hold GEO data in a straightforward format.

classes:
o Record - All of the information in an GEO record.

See http://www.ncbi.nlm.nih.gov/geo/
"""


class Record:
    """Hold GEO information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at GEO data.

    Attributes:
    entity_type
    entity_id
    entity_attributes
    col_defs
    table_rows

    """

    def __init__(self):
        """Initialize the class."""
        self.entity_type = ""
        self.entity_id = ""
        self.entity_attributes = {}
        self.col_defs = {}
        self.table_rows = []

    def __str__(self):
        """Return the GEO record as a string."""
        output = ""
        output += "GEO Type: %s\n" % self.entity_type
        output += "GEO Id: %s\n" % self.entity_id
        att_keys = sorted(self.entity_attributes)
        for key in att_keys:
            contents = self.entity_attributes[key]
            if isinstance(contents, list):
                for item in contents:
                    try:
                        output += "%s: %s\n" % (key, item[:40])
                        output += out_block(item[40:])
                    except Exception:  # TODO: IndexError?
                        pass
            elif isinstance(contents, str):
                output += "%s: %s\n" % (key, contents[:40])
                output += out_block(contents[40:])
            else:
                print(contents)
                output += "%s: %s\n" % (key, contents[:40])
                output += out_block(contents[40:])
        col_keys = sorted(self.col_defs)
        output += "Column Header Definitions\n"
        for key in col_keys:
            val = self.col_defs[key]
            output += "    %s: %s\n" % (key, val[:40])
            output += out_block(val[40:], "    ")
        # May have to display VERY large tables,
        # so only show the first 20 lines of data
        MAX_ROWS = 20 + 1  # include header in count
        for row in self.table_rows[0:MAX_ROWS]:
            output += "%s: " % self.table_rows.index(row)
            for col in row:
                output += "%s\t" % col
            output += "\n"
        if len(self.table_rows) > MAX_ROWS:
            output += "...\n"
            row = self.table_rows[-1]
            output += "%s: " % self.table_rows.index(row)
            for col in row:
                output += "%s\t" % col
            output += "\n"

        return output


def out_block(text, prefix=""):
    """Format text in blocks of 80 chars with an additional optional prefix."""
    output = ""
    for j in range(0, len(text), 80):
        output += "%s%s\n" % (prefix, text[j : j + 80])
    output += "\n"
    return output
