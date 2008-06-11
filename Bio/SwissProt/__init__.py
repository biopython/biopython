# Copyright 2007 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parsers for file formats from the SwissProt database.
"""

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
            raise ValueError, "No SwissProt record found"
        else:
            raise error
    # We should have reached the end of the record by now
    remainder = handle.read()
    if remainder:
        raise ValueError, "More than one SwissProt record found"
    return record


if __name__ == "__main__" :
    print "Quick self test..."

    example_filename = "../../Tests/SwissProt/sp008"

    import os
    if not os.path.isfile(example_filename):
        print "Missing test file %s" % example_filename
    else :
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
