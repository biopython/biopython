# Copyright 2006-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading and writing FASTA format files as SeqRecord
# objects.  The code is partly inspired  by earlier Biopython modules,
# Bio.Fasta.* and the now deprecated Bio.SeqIO.FASTA

"""Bio.SeqIO support for the "fasta" (aka FastA or Pearson) file format.

You are expected to use this module via the Bio.SeqIO functions."""

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequentialSequenceWriter

#This is a generator function!
def FastaIterator(handle, alphabet = single_letter_alphabet, title2ids = None):
    """Generator function to iterate over Fasta records (as SeqRecord objects).

    handle - input file
    alphabet - optional alphabet
    title2ids - A function that, when given the title of the FASTA
    file (without the beginning >), will return the id, name and
    description (in that order) for the record as a tuple of strings.

    If this is not given, then the entire title line will be used
    as the description, and the first word as the id and name.

    Note that use of title2ids matches that of Bio.Fasta.SequenceParser
    but the defaults are slightly different.
    """
    #Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0]!=">":
            raise ValueError("Records in Fasta files should start with '>' character")
        if title2ids:
            id, name, descr = title2ids(line[1:].rstrip())
        else:
            descr = line[1:].rstrip()
            try:
                id = descr.split()[0]
            except IndexError:
                assert not descr, repr(line)
                #Should we use SeqRecord default for no ID?
                id = ""
            name = id

        lines = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">": break
            lines.append(line.rstrip())
            line = handle.readline()

        #Remove trailing whitespace, and any internal spaces
        #(and any embedded \r which are possible in mangled files
        #when not opened in universal read lines mode)
        result = "".join(lines).replace(" ", "").replace("\r", "")

        #Return the record and then continue...
        yield SeqRecord(Seq(result, alphabet),
                         id = id, name = name, description = descr)

        if not line : return #StopIteration
    assert False, "Should not reach this line"

class FastaWriter(SequentialSequenceWriter):
    """Class to write Fasta format files."""
    def __init__(self, handle, wrap=60, record2title=None):
        """Create a Fasta writer.

        handle - Handle to an output file, e.g. as returned
                 by open(filename, "w")
        wrap -   Optional line length used to wrap sequence lines.
                 Defaults to wrapping the sequence at 60 characters
                 Use zero (or None) for no wrapping, giving a single
                 long line for the sequence.
        record2title - Optional function to return the text to be
                 used for the title line of each record.  By default the
                 a combination of the record.id and record.description
                 is used.  If the record.description starts with the
                 record.id, then just the record.description is used.

        You can either use:

        myWriter = FastaWriter(open(filename,"w"))
        writer.write_file(myRecords)

        Or, follow the sequential file writer system, for example:

        myWriter = FastaWriter(open(filename,"w"))
        writer.write_header() # does nothing for Fasta files
        ...
        Multiple calls to writer.write_record() and/or writer.write_records()
        ...
        writer.write_footer() # does nothing for Fasta files
        writer.close()
        """
        SequentialSequenceWriter.__init__(self, handle)
        #self.handle = handle
        self.wrap = None
        if wrap:
            if wrap < 1:
                raise ValueError
        self.wrap = wrap
        self.record2title = record2title

    def write_record(self, record):
        """Write a single Fasta record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        if self.record2title:
            title=self.clean(self.record2title(record))
        else:
            id = self.clean(record.id)
            description = self.clean(record.description)

            #if description[:len(id)]==id:
            if description and description.split(None,1)[0]==id:
                #The description includes the id at the start
                title = description
            elif description:
                title = "%s %s" % (id, description)
            else:
                title = id

        assert "\n" not in title
        assert "\r" not in title
        self.handle.write(">%s\n" % title)

        data = self._get_seq_string(record) #Catches sequence being None

        assert "\n" not in data
        assert "\r" not in data

        if self.wrap:
            for i in range(0, len(data), self.wrap):
                self.handle.write(data[i:i+self.wrap] + "\n")
        else:
            self.handle.write(data + "\n")

if __name__ == "__main__":
    print "Running quick self test"

    import os
    from Bio.Alphabet import generic_protein, generic_nucleotide

    #Download the files from here:
    #ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Nanoarchaeum_equitans
    fna_filename = "NC_005213.fna"
    faa_filename = "NC_005213.faa"

    def genbank_name_function(text):
        text, descr = text.split(None,1)
        id = text.split("|")[3]
        name = id.split(".",1)[0]
        return id, name, descr

    def print_record(record):
        #See also bug 2057
        #http://bugzilla.open-bio.org/show_bug.cgi?id=2057
        print "ID:" + record.id
        print "Name:" + record.name
        print "Descr:" + record.description
        print record.seq
        for feature in record.annotations:
            print '/%s=%s' % (feature, record.annotations[feature])
        if record.dbxrefs:
            print "Database cross references:"
            for x in record.dbxrefs : print " - %s" % x

    if os.path.isfile(fna_filename):
        print "--------"
        print "FastaIterator (single sequence)"
        iterator = FastaIterator(open(fna_filename, "r"), alphabet=generic_nucleotide, title2ids=genbank_name_function)
        count=0
        for record in iterator:
            count=count+1
            print_record(record)
        assert count == 1
        print str(record.__class__)

    if os.path.isfile(faa_filename):
        print "--------"
        print "FastaIterator (multiple sequences)"
        iterator = FastaIterator(open(faa_filename, "r"), alphabet=generic_protein, title2ids=genbank_name_function)
        count=0
        for record in iterator:
            count=count+1
            print_record(record)
            break
        assert count>0
        print str(record.__class__)

    from cStringIO import StringIO
    print "--------"
    print "FastaIterator (empty input file)"
    #Just to make sure no errors happen
    iterator = FastaIterator(StringIO(""))
    count = 0
    for record in iterator:
        count = count+1
    assert count==0

    print "Done"
