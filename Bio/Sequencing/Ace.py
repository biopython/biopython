# Copyright 2004 by Frank Kauff and Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for ACE files output by PHRAP.

Written by Frank Kauff (fkauff@duke.edu) and
Cymon J. Cox (cymon@duke.edu)

Usage:

There are two ways of reading an ace file:

1. The function 'read' reads the whole file at once;
2. The function 'parse' reads the file contig after contig.

First option, parse whole ace file at once::

        from Bio.Sequencing import Ace
        acefilerecord = Ace.read(open('my_ace_file.ace'))

This gives you:
 - acefilerecord.ncontigs (the number of contigs in the ace file)
 - acefilerecord.nreads (the number of reads in the ace file)
 - acefilerecord.contigs[] (one instance of the Contig class for each contig)

The Contig class holds the info of the CO tag, CT and WA tags, and all the reads used
for this contig in a list of instances of the Read class, e.g.::

        contig3 = acefilerecord.contigs[2]
        read4 = contig3.reads[3]
        RD_of_read4 = read4.rd
        DS_of_read4 = read4.ds

CT, WA, RT tags from the end of the file can appear anywhere are automatically
sorted into the right place.

see _RecordConsumer for details.

The second option is to  iterate over the contigs of an ace file one by one
in the ususal way::

    from Bio.Sequencing import Ace
    contigs = Ace.parse(open('my_ace_file.ace'))
    for contig in contigs:
        print(contig.name)
        ...

Please note that for memory efficiency, when using the iterator approach, only one
contig is kept in memory at once.  However, there can be a footer to the ACE file
containing WA, CT, RT or WR tags which contain additional meta-data on the contigs.
Because the parser doesn't see this data until the final record, it cannot be added to
the appropriate records.  Instead these tags will be returned with the last contig record.
Thus an ace file does not entirerly suit the concept of iterating. If WA, CT, RT, WR tags
are needed, the 'read' function rather than the 'parse' function might be more appropriate.
"""


class rd:
    """RD (reads), store a read with its name, sequence etc.

    The location and strand each read is mapped to is held in the AF lines.
    """

    def __init__(self):
        """Initialize the class."""
        self.name = ""
        self.padded_bases = None
        self.info_items = None
        self.read_tags = None
        self.sequence = ""


class qa:
    """QA (read quality), including which part if any was used as the consensus."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.qual_clipping_start = None
        self.qual_clipping_end = None
        self.align_clipping_start = None
        self.align_clipping_end = None
        if line:
            header = line.split()
            self.qual_clipping_start = int(header[1])
            self.qual_clipping_end = int(header[2])
            self.align_clipping_start = int(header[3])
            self.align_clipping_end = int(header[4])


class ds:
    """DS lines, include file name of a read's chromatogram file."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.chromat_file = ""
        self.phd_file = ""
        self.time = ""
        self.chem = ""
        self.dye = ""
        self.template = ""
        self.direction = ""
        if line:
            tags = [
                "CHROMAT_FILE",
                "PHD_FILE",
                "TIME",
                "CHEM",
                "DYE",
                "TEMPLATE",
                "DIRECTION",
            ]
            poss = [line.find(x) for x in tags]
            tagpos = dict(zip(poss, tags))
            if -1 in tagpos:
                del tagpos[-1]
            ps = sorted(tagpos)  # the keys
            for (p1, p2) in zip(ps, ps[1:] + [len(line) + 1]):
                setattr(
                    self,
                    tagpos[p1].lower(),
                    line[p1 + len(tagpos[p1]) + 1 : p2].strip(),
                )


class af:
    """AF lines, define the location of the read within the contig.

    Note attribute coru is short for complemented (C) or uncomplemented (U),
    since the strand information is stored in an ACE file using either the
    C or U character.
    """

    def __init__(self, line=None):
        """Initialize the class."""
        self.name = ""
        self.coru = None
        self.padded_start = None
        if line:
            header = line.split()
            self.name = header[1]
            self.coru = header[2]
            self.padded_start = int(header[3])


class bs:
    """BS (base segment), which read was chosen as the consensus at each position."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.name = ""
        self.padded_start = None
        self.padded_end = None
        if line:
            header = line.split()
            self.padded_start = int(header[1])
            self.padded_end = int(header[2])
            self.name = header[3]


class rt:
    """RT (transient read tags), generated by crossmatch and phrap."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.name = ""
        self.tag_type = ""
        self.program = ""
        self.padded_start = None
        self.padded_end = None
        self.date = ""
        self.comment = []
        if line:
            header = line.split()
            self.name = header[0]
            self.tag_type = header[1]
            self.program = header[2]
            self.padded_start = int(header[3])
            self.padded_end = int(header[4])
            self.date = header[5]


class ct:
    """CT (consensus tags)."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.name = ""
        self.tag_type = ""
        self.program = ""
        self.padded_start = None
        self.padded_end = None
        self.date = ""
        self.notrans = ""
        self.info = []
        self.comment = []
        if line:
            header = line.split()
            self.name = header[0]
            self.tag_type = header[1]
            self.program = header[2]
            self.padded_start = int(header[3])
            self.padded_end = int(header[4])
            self.date = header[5]
            if len(header) == 7:
                self.notrans = header[6]


class wa:
    """WA (whole assembly tag), holds the assembly program name, version, etc."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.tag_type = ""
        self.program = ""
        self.date = ""
        self.info = []
        if line:
            header = line.split()
            self.tag_type = header[0]
            self.program = header[1]
            self.date = header[2]


class wr:
    """WR lines."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.name = ""
        self.aligned = ""
        self.program = ""
        self.date = []
        if line:
            header = line.split()
            self.name = header[0]
            self.aligned = header[1]
            self.program = header[2]
            self.date = header[3]


class Reads:
    """Holds information about a read supporting an ACE contig."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.rd = None  # one per read
        self.qa = None  # one per read
        self.ds = None  # none or one per read
        self.rt = None  # none or many per read
        self.wr = None  # none or many per read
        if line:
            self.rd = rd()
            header = line.split()
            self.rd.name = header[1]
            self.rd.padded_bases = int(header[2])
            self.rd.info_items = int(header[3])
            self.rd.read_tags = int(header[4])


class Contig:
    """Holds information about a contig from an ACE record."""

    def __init__(self, line=None):
        """Initialize the class."""
        self.name = ""
        self.nbases = None
        self.nreads = None
        self.nsegments = None
        self.uorc = None
        self.sequence = ""
        self.quality = []
        self.af = []
        self.bs = []
        self.reads = []
        self.ct = None  # none or many
        self.wa = None  # none or many
        if line:
            header = line.split()
            self.name = header[1]
            self.nbases = int(header[2])
            self.nreads = int(header[3])
            self.nsegments = int(header[4])
            self.uorc = header[5]


def parse(source):
    """Iterate of ACE file contig by contig.

    Argument source is a file-like object or a path to a file.

    This function returns an iterator that allows you to iterate
    over the ACE file record by record::

        records = parse(source)
        for record in records:
            # do something with the record

    where each record is a Contig object.
    """
    try:
        handle = open(source)
    except TypeError:
        handle = source
        if handle.read(0) != "":
            raise ValueError("Ace files must be opened in text mode.") from None

    try:
        line = ""
        while True:
            # at beginning, skip the AS and look for first CO command
            try:
                while True:
                    if line.startswith("CO"):
                        break
                    line = next(handle)
            except StopIteration:
                return

            record = Contig(line)

            for line in handle:
                line = line.strip()
                if not line:
                    break
                record.sequence += line

            for line in handle:
                if line.strip():
                    break
            if not line.startswith("BQ"):
                raise ValueError("Failed to find BQ line")

            for line in handle:
                if not line.strip():
                    break
                record.quality.extend(int(x) for x in line.split())

            for line in handle:
                if line.strip():
                    break

            while True:
                if not line.startswith("AF "):
                    break
                record.af.append(af(line))
                try:
                    line = next(handle)
                except StopIteration:
                    raise ValueError("Unexpected end of AF block") from None

            while True:
                if line.strip():
                    break
                try:
                    line = next(handle)
                except StopIteration:
                    raise ValueError("Unexpected end of file") from None

            while True:
                if not line.startswith("BS "):
                    break
                record.bs.append(bs(line))
                try:
                    line = next(handle)
                except StopIteration:
                    raise ValueError("Failed to find end of BS block") from None

            # now read all the read data
            # it starts with a 'RD', and then a mandatory QA
            # then follows an optional DS
            # CT,RT,WA,WR may or may not be there in unlimited quantity.
            # They might refer to the actual read or contig, or, if
            # encountered at the end of file, to any previous read or contig.
            # The sort() method deals with that later.
            while True:

                # each read must have a rd and qa
                try:
                    while True:
                        # If I've met the condition, then stop reading the line.
                        if line.startswith("RD "):
                            break
                        line = next(handle)
                except StopIteration:
                    raise ValueError("Failed to find RD line") from None

                record.reads.append(Reads(line))

                for line in handle:
                    line = line.strip()
                    if not line:
                        break
                    record.reads[-1].rd.sequence += line

                for line in handle:
                    if line.strip():
                        break
                if not line.startswith("QA "):
                    raise ValueError("Failed to find QA line")
                record.reads[-1].qa = qa(line)

                # now one ds can follow
                for line in handle:
                    if line.strip():
                        break
                else:
                    break

                if line.startswith("DS "):
                    record.reads[-1].ds = ds(line)
                    line = ""
                # the file could just end, or there's some more stuff.
                # In ace files, anything can happen.
                # the following tags are interspersed between reads and can appear multiple times.
                while True:
                    # something left
                    try:
                        while True:
                            if line.strip():
                                break
                            line = next(handle)
                    except StopIteration:
                        # file ends here
                        break
                    if line.startswith("RT{"):
                        # now if we're at the end of the file, this rt could
                        # belong to a previous read, not the actual one.
                        # we store it here were it appears, the user can sort later.
                        if record.reads[-1].rt is None:
                            record.reads[-1].rt = []
                        for line in handle:
                            line = line.strip()
                            # if line=="COMMENT{":
                            if line.startswith("COMMENT{"):
                                if line[8:].strip():
                                    # MIRA 3.0.5 would miss the new line out :(
                                    record.reads[-1].rt[-1].comment.append(line[8:])
                                for line in handle:
                                    line = line.strip()
                                    if line.endswith("C}"):
                                        break
                                    record.reads[-1].rt[-1].comment.append(line)
                            elif line == "}":
                                break
                            else:
                                record.reads[-1].rt.append(rt(line))
                        line = ""
                    elif line.startswith("WR{"):
                        if record.reads[-1].wr is None:
                            record.reads[-1].wr = []
                        for line in handle:
                            line = line.strip()
                            if line == "}":
                                break
                            record.reads[-1].wr.append(wr(line))
                        line = ""
                    elif line.startswith("WA{"):
                        if record.wa is None:
                            record.wa = []
                        try:
                            line = next(handle)
                        except StopIteration:
                            raise ValueError("Failed to read WA block") from None
                        record.wa.append(wa(line))
                        for line in handle:
                            line = line.strip()
                            if line == "}":
                                break
                            record.wa[-1].info.append(line)
                        line = ""
                    elif line.startswith("CT{"):
                        if record.ct is None:
                            record.ct = []
                        try:
                            line = next(handle)
                        except StopIteration:
                            raise ValueError("Failed to read CT block") from None
                        record.ct.append(ct(line))
                        for line in handle:
                            line = line.strip()
                            if line == "COMMENT{":
                                for line in handle:
                                    line = line.strip()
                                    if line.endswith("C}"):
                                        break
                                    record.ct[-1].comment.append(line)
                            elif line == "}":
                                break
                            else:
                                record.ct[-1].info.append(line)
                        line = ""
                    else:
                        break

                if not line.startswith("RD"):  # another read?
                    break

            yield record

    finally:
        if handle is not source:
            handle.close()


class ACEFileRecord:
    """Holds data of an ACE file."""

    def __init__(self):
        """Initialize the class."""
        self.ncontigs = None
        self.nreads = None
        self.contigs = []
        self.wa = None  # none or many

    def sort(self):
        """Sorts wr, rt and ct tags into the appropriate contig / read instance, if possible."""
        ct = []
        rt = []
        wr = []
        # search for tags that aren't in the right position
        for i, c in enumerate(self.contigs):
            if c.wa:
                if not self.wa:
                    self.wa = []
                self.wa.extend(c.wa)
            if c.ct:
                newcts = [ct_tag for ct_tag in c.ct if ct_tag.name != c.name]
                for x in newcts:
                    self.contigs[i].ct.remove(x)
                ct.extend(newcts)
            for j, r in enumerate(c.reads):
                if r.rt:
                    newrts = [rt_tag for rt_tag in r.rt if rt_tag.name != r.rd.name]
                    for x in newrts:
                        self.contigs[i].reads[j].rt.remove(x)
                    rt.extend(newrts)
                if r.wr:
                    newwrs = [wr_tag for wr_tag in r.wr if wr_tag.name != r.rd.name]
                    for x in newwrs:
                        self.contigs[i].reads[j].wr.remove(x)
                    wr.extend(newwrs)
        # now sort them into their proper place
        for i, c in enumerate(self.contigs):
            for ct_tag in ct:
                if ct_tag.name == c.name:
                    if self.contigs[i].ct is None:
                        self.contigs[i].ct = []
                    self.contigs[i].ct.append(ct_tag)
            if rt or wr:
                for j, r in enumerate(c.reads):
                    for rt_tag in rt:
                        if rt_tag.name == r.rd.name:
                            if self.contigs[i].reads[j].rt is None:
                                self.contigs[i].reads[j].rt = []
                            self.contigs[i].reads[j].rt.append(rt_tag)
                    for wr_tag in wr:
                        if wr_tag.name == r.rd.name:
                            if self.contigs[i].reads[j].wr is None:
                                self.contigs[i].reads[j].wr = []
                            self.contigs[i].reads[j].wr.append(wr_tag)


def read(handle):
    """Parse a full ACE file into a list of contigs."""
    handle = iter(handle)

    record = ACEFileRecord()

    try:
        line = next(handle)
    except StopIteration:
        raise ValueError("Premature end of file") from None

    # check if the file starts correctly
    if not line.startswith("AS"):
        raise ValueError("File does not start with 'AS'.")

    words = line.split()
    record.ncontigs = int(words[1])
    record.nreads = int(words[2])

    # now read all the records
    record.contigs = list(parse(handle))
    # wa, ct, rt rags are usually at the end of the file, but not necessarily (correct?).
    # If the iterator is used, the tags are returned with the contig or the read after which they appear,
    # if all tags are at the end, they are read with the last contig. The concept of an
    # iterator leaves no other choice. But if the user uses the ACEParser, we can check
    # them and put them into the appropriate contig/read instance.
    # Conclusion: An ACE file is not a filetype for which iteration is 100% suitable...
    record.sort()
    return record
