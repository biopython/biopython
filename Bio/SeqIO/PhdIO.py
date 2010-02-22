# Copyright 2008-2010 by Peter Cock.  All rights reserved.
# Revisions copyright 2009 by Cymon J. Cox.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "phd" file format.

PHD files are output by PHRED and used by PHRAP and CONSED.

You are expected to use this module via the Bio.SeqIO functions, under the
format name "phd". See also the underlying Bio.Sequencing.Phd module.

For example, using Bio.SeqIO we can read in one of the example PHRED files
from the Biopython unit tests:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse(open("Phd/phd1"), "phd"):
    ...     print record.id
    ...     print record.seq[:10], "..."
    ...     print record.letter_annotations["phred_quality"][:10], "..."
    34_222_(80-A03-19).b.ab1
    ctccgtcgga ...
    [9, 9, 10, 19, 22, 37, 28, 28, 24, 22] ...
    425_103_(81-A03-19).g.ab1
    cgggatccca ...
    [14, 17, 22, 10, 10, 10, 15, 8, 8, 9] ...
    425_7_(71-A03-19).b.ab1
    acataaatca ...
    [10, 10, 10, 10, 8, 8, 6, 6, 6, 6] ...

Since PHRED files contain quality scores, you can save them as FASTQ or as
QUAL files, for example using Bio.SeqIO.write(...), or simply with the format
method of the SeqRecord object:

    >>> print record[:50].format("fastq")
    @425_7_(71-A03-19).b.ab1
    acataaatcaaattactnaccaacacacaaaccngtctcgcgtagtggag
    +
    ++++))'''')(''')$!$''')''''(+.''$!$))))+)))'''''''
    <BLANKLINE>

Or,

    >>> print record[:50].format("qual")
    >425_7_(71-A03-19).b.ab1
    10 10 10 10 8 8 6 6 6 6 8 7 6 6 6 8 3 0 3 6 6 6 8 6 6 6 6 7
    10 13 6 6 3 0 3 8 8 8 8 10 8 8 8 6 6 6 6 6 6 6
    <BLANKLINE>

Note these examples only show the first 50 bases to keep the output short.
"""

from Bio.SeqRecord import SeqRecord
from Bio.Sequencing import Phd
from Bio.SeqIO.Interfaces import SequentialSequenceWriter
from Bio.SeqIO import QualityIO
    
#This is a generator function!
def PhdIterator(handle):
    """Returns SeqRecord objects from a PHD file.

    This uses the Bio.Sequencing.Phd module to do the hard work.
    """
    phd_records = Phd.parse(handle)
    for phd_record in phd_records:
        #Convert the PHY record into a SeqRecord...
        #The "filename" can contain spaces, e.g. 'HWI-EAS94_4_1_1_602_99 1'
        #from unit test example file phd_solexa.
        #This will cause problems if used as the record identifier
        #(e.g. output for FASTQ format).
        name = phd_record.file_name.split(None, 1)[0]
        seq_record = SeqRecord(phd_record.seq,
                               id = name, name = name,
                               description= phd_record.file_name)
        #Just re-use the comments dictionary as the SeqRecord's annotations
        seq_record.annotations = phd_record.comments
        #And store the qualities and peak locations as per-letter-annotation
        seq_record.letter_annotations["phred_quality"] = \
                [int(site[1]) for site in phd_record.sites]
        try:
            seq_record.letter_annotations["peak_location"] = \
                    [int(site[2]) for site in phd_record.sites]
        except IndexError:
            # peak locations are not always there according to
            # David Gordon (the Consed author)
            pass
        yield seq_record 
    #All done

class PhdWriter(SequentialSequenceWriter):
    """Class to write Phd format files"""

    def __init__(self, handle):
        SequentialSequenceWriter.__init__(self, handle)

    def write_record(self, record):
        """Write a single Phd record to the file."""
        assert record.seq, "No sequence present in SeqRecord"
        # This method returns the 'phred_quality' scores or converted
        # 'solexa_quality' scores if present, else raises a value error
        phred_qualities = QualityIO._get_phred_quality(record)
        peak_locations = record.letter_annotations.get("peak_location", None)
        assert len(record.seq) == len(phred_qualities), "Number of " + \
                "phd quality scores does not match length of sequence"
        if peak_locations:
            assert len(record.seq) == len(peak_locations), "Number " + \
                    "of peak location scores does not match length of sequence"
        if None in phred_qualities:
            raise ValueError("A quality value of None was found")
        if record.description.startswith("%s " % record.id):
            title = record.description
        else:
            title = "%s %s" % (record.id, record.description)
        self.handle.write("BEGIN_SEQUENCE %s\nBEGIN_COMMENT\n" \
                          % self.clean(title))
        for annot in [k.lower() for k in Phd.CKEYWORDS]:
            value = None
            if annot == "trim":
                if record.annotations.get("trim", None):
                    value = "%s %s %.4f" % record.annotations["trim"]
            elif annot == "trace_peak_area_ratio":
                if record.annotations.get("trace_peak_area_ratio", None):
                    value = "%.4f" % record.annotations["trace_peak_area_ratio"]
            else:
                value = record.annotations.get(annot, None)
            if value or value == 0:
                self.handle.write("%s: %s\n" % (annot.upper(), value))

        self.handle.write("END_COMMENT\nBEGIN_DNA\n")
        for i, site in enumerate(record.seq):
            if peak_locations:
                self.handle.write("%s %i %i\n" % (
                        site,
                        round(phred_qualities[i]),
                        peak_locations[i])
                        )
            else:
                self.handle.write("%s %i\n" % (
                        site,
                        round(phred_qualities[i]))
                        )

        self.handle.write("END_DNA\nEND_SEQUENCE\n")

def _test():
    """Run the Bio.SeqIO.PhdIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..", "..", "Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "..", "Tests"))
        assert os.path.isfile("Phd/phd1")
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
        
if __name__ == "__main__":
    _test()

        
    
