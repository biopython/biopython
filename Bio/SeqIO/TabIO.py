# Copyright 2008-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "tab" (simple tab separated) file format.

You are expected to use this module via the Bio.SeqIO functions.

The "tab" format is an ad-hoc plain text file format where each sequence is
on one (long) line.  Each line contains the identifier/description, followed
by a tab, followed by the sequence.  For example, consider the following
short FASTA format file:

>ID123456 possible binding site?
CATCNAGATGACACTACGACTACGACTCAGACTAC
>ID123457 random sequence
ACACTACGACTACGACTCAGACTACAAN

Apart from the descriptions, this can be represented in the simple two column
tab separated format as follows:

ID123456(tab)CATCNAGATGACACTACGACTACGACTCAGACTAC
ID123457(tab)ACACTACGACTACGACTCAGACTACAAN

When reading this file, "ID123456" or "ID123457" will be taken as the record's
.id and .name property.  There is no other information to record.

Similarly, when writing to this format, Biopython will ONLY record the record's
.id and .seq (and not the description or any other information) as in the
example above.
"""

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequentialSequenceWriter

#This is a generator function!
def TabIterator(handle, alphabet = single_letter_alphabet):
    """Iterates over tab separated lines (as SeqRecord objects).

    Each line of the file should contain one tab only, dividing the line
    into an identifier and the full sequence.

    handle - input file
    alphabet - optional alphabet

    The first field is taken as the record's .id and .name (regardless of
    any spaces within the text) and the second field is the sequence.

    Any blank lines are ignored.
    """
    for line in handle:
        try:
            title, seq = line.split("\t") #will fail if more than one tab!
        except:
            if line.strip() == "":
                #It's a blank line, ignore it
                continue
            raise ValueError("Each line should have one tab separating the" + \
                             " title and sequence, this line has %i tabs: %s" \
                             % (line.count("\t"), repr(line)))
        title = title.strip()
        seq = seq.strip() #removes the trailing new line
        yield SeqRecord(Seq(seq, alphabet),
                        id=title, name=title,
                        description="")

class TabWriter(SequentialSequenceWriter):
    """Class to write simple tab separated format files.

    Each line consists of "id(tab)sequence" only.

    Any description, name or other annotation is not recorded.
    """
    def write_record(self, record):
        """Write a single tab line to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True
        
        title = self.clean(record.id)
        seq = self._get_seq_string(record) #Catches sequence being None
        assert "\t" not in title
        assert "\n" not in title
        assert "\r" not in title
        assert "\t" not in seq
        assert "\n" not in seq
        assert "\r" not in seq
        self.handle.write("%s\t%s\n" % (title, seq))


if __name__ == "__main__":
    print "Running quick self test"
    from StringIO import StringIO

    #This example has a trailing blank line which should be ignored
    handle = StringIO("Alpha\tAAAAAAA\nBeta\tCCCCCCC\n\n")
    records = list(TabIterator(handle))
    assert len(records) == 2

    handle = StringIO("Alpha\tAAAAAAA\tExtra\nBeta\tCCCCCCC\n")
    try:
        records = list(TabIterator(handle))
        assert False, "Should have reject this invalid example!"
    except ValueError:
        #Good!
        pass

    print "Done"    
