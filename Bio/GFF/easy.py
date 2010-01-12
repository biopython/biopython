#!/usr/bin/env python
#
# Copyright 2002 by Michael Hoffman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.GFF.easy: some functions to ease the use of Biopython (DEPRECATED)

This is part of the "old" Bio.GFF module by Michael Hoffman, which offered
access to a MySQL database holding GFF data loaded by BioPerl. This code has
now been deprecated, and will probably be removed in order to free the Bio.GFF
namespace for a new GFF parser in Biopython (including GFF3 support).

Some of the more useful ideas of Bio.GFF.easy may be reworked for Bio.GenBank,
using the standard SeqFeature objects used elsewhere in Biopython.
"""

import copy
import re
import sys

import Bio
from Bio import GenBank
from Bio.Data import IUPACData
from Bio.Seq import Seq

from Bio import SeqIO
from Bio import SeqUtils

import GenericTools

class FeatureDict(dict):
    """ JH:  accessing feature.qualifiers as a list is stupid.  Here's a dict that does it"""
    def __init__(self, feature_list, default=None):
        dict.__init__(self)
        self.default = default
        key_re = re.compile(r'/(\S+)=')

        for i in feature_list:
            key = key_re.match(i.key).group(1)
            val = i.value.replace('"','')
            self[key] = val
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return self.default

class Location(GenericTools.VerboseList):
    """
    this is really best interfaced through LocationFromString
    fuzzy: < or >
    join: {0 = no join, 1 = join, 2 = order}
    
    >>> location = Location([Location([339]), Location([564])]) # zero-based
    >>> location
    Location(Location(339), Location(564))
    >>> print location # one-based
    340..565
    >>> print location.five_prime()
    340
    >>> location_rev = Location([Location([339]), Location([564])], 1)
    >>> print location_rev
    complement(340..565)
    >>> print location_rev.five_prime()
    565
    """
    def __init__(self, the_list, complement=0, seqname=None):
        self.complement = complement
        self.join = 0
        self.fuzzy = None
        self.seqname = seqname
        list.__init__(self, the_list)

    def _joinstr(self):
        if self.join == 1:
            label = 'join'
        elif self.join == 2:
            label = 'order'
        return "%s(%s)" % (label, ",".join(map(str, self)))
    
    def __str__(self):
        if self.seqname:
            format = "%s:%%s" % self.seqname
        else:
            format = "%s"

        if self.complement:
            format = format % "complement(%s)"

        if self.join:
            return format % self._joinstr()
            
        elif isinstance(self[0], list):
            return format % "%s..%s" % (str(self[0]), str(self[1]))
        else:
            if self.fuzzy:
                format = format % self.fuzzy + "%s"
            return format % str(self[0] + 1)

    def __repr__(self):
        return "Location(%s)" % ", ".join(map(repr, self))

    direction2index = {1: 0, -1: -1}
    def direction_and_index(self, direction):
        """
        1: 5'
        -1: 3'

        >>> loc1 = LocationFromString("join(1..3,complement(5..6))")
        >>> loc1.direction_and_index(1)
        (1, 0)
        >>> loc1.direction_and_index(-1)
        (-1, -1)
        >>> loc1.reverse()
        >>> print loc1
        complement(join(1..3,complement(5..6)))
        >>> loc1.direction_and_index(1)
        (-1, -1)
        """
        if self.complement:
            direction = direction * -1
        index = self.direction2index[direction]
        return direction, index
    
    def findside(self, direction):
        """
        >>> loc = LocationFromString("complement(join(1..5,complement(6..10)))")
        >>> loc.findside(1)
        Location(5)
        >>> loc.findside(-1)
        Location(0)
        """
        direction, index = self.direction_and_index(direction)
        if self.join or isinstance(self[0], list):
            return self[index].findside(direction)
        else:
            return self

    def findseqname_3prime(self):
        """
        >>> loc = LocationFromString("complement(join(MOOCOW:1..5,SEQ:complement(6..10)))")
        >>> loc.findseqname_3prime()
        'MOOCOW'
        """
        return self.findseqname(-1)

    def findseqname(self, direction=1): # find 5' seqname
        """
        >>> loc = LocationFromString("complement(join(MOOCOW:1..5,SEQ:complement(6..10)))")
        >>> loc.findseqname()
        'SEQ'
        >>> loc.findseqname(-1)
        'MOOCOW'
        """
        direction, index = self.direction_and_index(direction)
        if self.seqname:
            return self.seqname
        elif self.join:
            return self[index].findseqname(direction)
        else:
            raise AttributeError('no sequence name')

    def five_prime(self):
        return self.findside(1)
    def three_prime(self):
        return self.findside(-1)

    def length(self):
        """
        WARNING: doesn't deal with joins!!!!
        """
        return self.end()-self.start()

    def intersection(self, other):
        """
        WARNING: doesn't deal with joins!!!!

        >>> location1 = LocationFromString("1..50")
        >>> location2 = LocationFromString("25..200")
        >>> print location1.intersection(location2)
        25..50
        >>> print location1.intersection(location2)
        25..50
        """
        if self.start() >= other.start():
            start = self.start()
        else:
            start = other.start()
        if self.end() <= other.end():
            end = self.end()
        else:
            end = other.end()
        return Location([Location([start]), Location([end])])

    def start(self):
        # zero-based
        if self.complement:
            return self.three_prime()[0]
        else:
            return self.five_prime()[0]

    def end(self):
        # zero-based
        if self.complement:
            return self.five_prime()[0]
        else:
            return self.three_prime()[0]

    def three_prime_range(self, window):
        three_prime_loc = self.three_prime()
        if self.complement:
            return Location([three_prime_loc-window, three_prime_loc], complement=1)
        else:
            return Location([three_prime_loc, three_prime_loc+window])

    def sublocation(self, sub_location):
        """
        >>> fwd_location = LocationFromString('X:5830132..5831528')
        >>> print fwd_location.sublocation(LocationFromString('1..101'))
        X:5830132..5830232
        >>> print fwd_location.sublocation(LocationFromString('1267..1286'))
        X:5831398..5831417
        >>> rev_location = LocationFromString('I:complement(8415686..8416216)')
        >>> print rev_location.sublocation(LocationFromString('1..101'))
        I:complement(8416116..8416216)
        >>> print rev_location.sublocation(LocationFromString('100..200'))
        I:complement(8416017..8416117)
        """
        
        absolute_location = copy.deepcopy(self)
        for i in xrange(2):
            absolute_location[i] = self.five_prime().add(sub_location[i], self.complement)
        if absolute_location.complement:
            list.reverse(absolute_location)
        return absolute_location

    def __add__(self, addend):
        return self.add(addend)

    def add(self, addend, complement=0):
        self_copy = copy.deepcopy(self)
        if isinstance(addend, Location):
            addend = addend[0]
        if complement:
            addend *= -1
        self_copy[0] += addend
        return self_copy

    def __sub__(self, subtrahend):
        return self + -subtrahend

    def reverse(self):
        self.complement = [1, 0][self.complement]

    def reorient(self):
        """
        >>> loc1 = LocationFromString("join(I:complement(1..9000),I:complement(9001..10000))")
        >>> loc1.reorient()
        >>> print loc1
        complement(join(I:1..9000,I:9001..10000))
        >>> loc2 = LocationFromString("join(I:complement(1..9000),I:9001..10000)")
        >>> loc2.reorient()
        >>> print loc2
        join(I:complement(1..9000),I:9001..10000)
        """
        if self.join:
            if len([x for x in self if x.complement]) == len(self):
                self.reverse()
                for segment in self:
                    segment.reverse()

    def bounding(self):
        """
        works for single level non-complex joins

        >>> LOC = LocationFromString
        >>> l1 = LOC("join(alpha:1..30,alpha:50..70)")
        >>> print l1.bounding()
        join(alpha:1..70)
        >>> l2 = LOC("join(alpha:1..30,alpha:complement(50..70))")
        >>> print l2.bounding()
        join(alpha:1..30,alpha:complement(50..70))
        >>> l3 = LOC("join(alpha:1..30,alpha:complement(50..70),beta:6..20,alpha:25..45)")
        >>> print l3.bounding()
        join(alpha:1..45,alpha:complement(50..70),beta:6..20)

        """
        if not self.join:
            return

        seqdict = {}
        seqkeys = []
        for subloc in self:
            assert subloc.seqname
            assert not subloc.join
            try:
                seqdict[_hashname(subloc)].append(subloc)
            except KeyError:
                key = _hashname(subloc)
                seqdict[key] = [subloc]
                seqkeys.append(key)

        res = LocationJoin()
        for key in seqkeys:
            locations = seqdict[key]
            coords = []
            for subloc in locations:
                coords.append(subloc.start())
                coords.append(subloc.end())
            res.append(LocationFromCoords(min(coords), max(coords), locations[0].complement, locations[0].seqname))
        return res

def _hashname(location):
    return str(location.complement) + location.seqname

class LocationJoin(Location):
    """
    >>> join = LocationJoin([LocationFromCoords(339, 564, 1), LocationFromString("complement(100..339)")])
    >>> appendloc = LocationFromString("order(complement(66..99),complement(5..55))")
    >>> join.append(appendloc)
    >>> print join
    join(complement(340..565),complement(100..339),order(complement(66..99),complement(5..55)))
    >>> join2 = LocationJoin()
    >>> join2.append(LocationFromString("complement(66..99)"))
    >>> join2.append(LocationFromString("complement(5..55)"))
    >>> print join2
    join(complement(66..99),complement(5..55))
    """
    def __init__(self, the_list = [], complement=0, seqname=None):
        self.complement = complement
        self.join = 1
        self.fuzzy = None
        self.seqname = seqname
        list.__init__(self, the_list)

class LocationFromCoords(Location):
    """
    >>> print LocationFromCoords(339, 564)
    340..565
    >>> print LocationFromCoords(339, 564, seqname="I")
    I:340..565
    >>> print LocationFromCoords(999, 3234, "-", seqname="NC_343434")
    NC_343434:complement(1000..3235)
    """
    def __init__(self, start, end, strand=0, seqname=None):
        if strand == "+":
            strand = 0
        elif strand == "-":
            strand = 1
        Location.__init__(self, [Location([start]), Location([end])], strand, seqname)

# see http://www.ncbi.nlm.nih.gov/collab/FT/index.html#backus-naur
# for how this should actually be implemented
re_complement = re.compile(r"^complement\((.*)\)$")
re_seqname = re.compile(r"^(?!join|order|complement)([^\:]+?):(.*)$") # not every character is allowed by spec
re_join = re.compile(r"^(join|order)\((.*)\)$")
re_dotdot = re.compile(r"^([><]*\d+)\.\.([><]*\d+)$")
re_fuzzy = re.compile(r"^([><])(\d+)")
class LocationFromString(Location):
    """
    >>> # here are some tests from http://www.ncbi.nlm.nih.gov/collab/FT/index.html#location
    >>> print LocationFromString("467")
    467
    >>> print LocationFromString("340..565")
    340..565
    >>> print LocationFromString("<345..500")
    <345..500
    >>> print LocationFromString("<1..888")
    <1..888
    >>> # (102.110) and 123^124 syntax unimplemented
    >>> print LocationFromString("join(12..78,134..202)")
    join(12..78,134..202)
    >>> print LocationFromString("complement(join(2691..4571,4918..5163))")
    complement(join(2691..4571,4918..5163))
    >>> print LocationFromString("join(complement(4918..5163),complement(2691..4571))")
    join(complement(4918..5163),complement(2691..4571))
    >>> print LocationFromString("order(complement(4918..5163),complement(2691..4571))")
    order(complement(4918..5163),complement(2691..4571))
    >>> print LocationFromString("NC_001802x.fna:73..78")
    NC_001802x.fna:73..78
    >>> print LocationFromString("J00194:100..202")
    J00194:100..202

    >>> print LocationFromString("join(117505..118584,1..609)")
    join(117505..118584,1..609)
    >>> print LocationFromString("join(test3:complement(4..6),test3:complement(1..3))")
    join(test3:complement(4..6),test3:complement(1..3))
    >>> print LocationFromString("test3:join(test1:complement(1..3),4..6)")
    test3:join(test1:complement(1..3),4..6)
    """
    def __init__(self, location_str):
        match_seqname = re_seqname.match(location_str)
        if match_seqname:
            self.seqname = match_seqname.group(1)
            location_str = match_seqname.group(2)
        else:
            self.seqname = None
        match_complement = re_complement.match(location_str)
        if match_complement:
            self.complement = 1
            location_str = match_complement.group(1)
        else:
            self.complement = 0
        match_join = re_join.match(location_str)
        if match_join:
            self.join = {'join':1, 'order':2}[match_join.group(1)]
            list.__init__(self, map(lambda x: LocationFromString(x), match_join.group(2).split(",")))
        else:
            self.join = 0
            match_dotdot = re_dotdot.match(location_str)
            if match_dotdot:
                list.__init__(self, map(lambda x: LocationFromString(match_dotdot.group(x)), (1, 2)))
            else:
                match_fuzzy = re_fuzzy.match(location_str)
                if match_fuzzy:
                    self.fuzzy = match_fuzzy.group(1)
                    location_str = match_fuzzy.group(2)
                else:
                    self.fuzzy = None
                    
                list.__init__(self, [int(location_str)-1]) # zero based, nip it in the bud

def open_file(filename):
    if filename:
        return open(filename)
    else:
        return sys.stdin

def fasta_single(filename=None, string=None):
    """
    >>> record = fasta_single(string=\"""
    ... >gi|9629360|ref|NP_057850.1| Gag [Human immunodeficiency virus type 1]
    ... MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQT
    ... GSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQG
    ... QMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAA
    ... EWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPT
    ... SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTAC
    ... QGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEG
    ... HQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLR
    ... SLFGNDPSSQ
    ... \""")
    >>> record.id
    'gi|9629360|ref|NP_057850.1|'
    >>> record.description
    'gi|9629360|ref|NP_057850.1| Gag [Human immunodeficiency virus type 1]'
    >>> record.seq[0:5]
    Seq('MGARA', SingleLetterAlphabet())
    """
    #Returns the first record in a fasta file as a SeqRecord,
    #or None if there are no records in the file.
    if string:
        import cStringIO
        handle = cStringIO.StringIO(string)
    else:
        handle = open_file(filename)
    try:
        record = SeqIO.parse(handle, format="fasta").next()
    except StopIteration:
        record = None
    return record

def fasta_multi(filename=None):
    #Simple way is just:
    #return SeqIO.parse(open_file(filename), format="fasta")
    #However, for backwards compatibility make sure we raise
    #the StopIteration exception rather than returning None.
    reader = SeqIO.parse(open_file(filename), format="fasta")
    while True:
        record = reader.next()
        if record is None:
            raise StopIteration
        else:
            yield record

def fasta_readrecords(filename=None):
    """
    >>> records = fasta_readrecords('GFF/multi.fna')
    >>> records[0].id
    'test1'
    >>> records[2].seq
    Seq('AAACACAC', SingleLetterAlphabet())
    """
    return list(SeqIO.parse(open_file(filename), format="fasta"))

def fasta_write(filename, records):
    handle = open(filename, "w")
    SeqIO.write(records, handle, format="fasta")
    handle.close()

def genbank_single(filename):
    """
    >>> record = genbank_single("GFF/NC_001422.gbk")
    >>> record.taxonomy
    ['Viruses', 'ssDNA viruses', 'Microviridae', 'Microvirus']
    >>> cds = record.features[-4]
    >>> cds.key
    'CDS'
    >>> location = LocationFromString(cds.location)
    >>> print location
    2931..3917
    >>> subseq = record_subseq(record, location)
    >>> subseq[0:20]
    Seq('ATGTTTGGTGCTATTGCTGG', Alphabet())
    """
    return GenBank.RecordParser().parse(open(filename))

def record_subseq(record, location, *args, **keywds):
    """
    >>> from Bio.SeqRecord import SeqRecord    
    >>> record = SeqRecord(Seq("gagttttatcgcttccatga"),
    ...                    "ref|NC_001422",
    ...                    "Coliphage phiX174, complete genome",
    ...                    "bases 1-11")
    >>> record_subseq(record, LocationFromString("1..4")) # one-based
    Seq('GAGT', Alphabet())
    >>> record_subseq(record, LocationFromString("complement(1..4)")) # one-based
    Seq('ACTC', Alphabet())
    >>> record_subseq(record, LocationFromString("join(complement(1..4),1..4)")) # what an idea!
    Seq('ACTCGAGT', Alphabet())
    >>> loc = LocationFromString("complement(join(complement(5..7),1..4))")
    >>> print loc
    complement(join(complement(5..7),1..4))
    >>> record_subseq(record, loc)
    Seq('ACTCTTT', Alphabet())
    >>> print loc
    complement(join(complement(5..7),1..4))
    >>> loc.reverse()
    >>> record_subseq(record, loc)
    Seq('AAAGAGT', Alphabet())
    >>> record_subseq(record, loc, upper=1)
    Seq('AAAGAGT', Alphabet())
    """
    if location.join:
        subsequence_list = []
        if location.complement:
            location_copy = copy.copy(location)
            list.reverse(location_copy)
        else:
            location_copy = location
        for sublocation in location_copy:
            if location.complement:
                sublocation_copy = copy.copy(sublocation)
                sublocation_copy.reverse()
            else:
                sublocation_copy = sublocation
            subsequence_list.append(record_subseq(record, sublocation_copy, *args, **keywds).tostring())
        return Seq(''.join(subsequence_list), record_sequence(record).alphabet)
    else:
        return record_coords(record, location.start(), location.end()+1, location.complement, *args, **keywds)

def record_sequence(record):
    """
    returns the sequence of a record

    can be Bio.SeqRecord.SeqRecord or Bio.GenBank.Record.Record
    """
    if isinstance(record, Bio.SeqRecord.SeqRecord):
        return record.seq
    elif isinstance(record, Bio.GenBank.Record.Record):
        return Seq(record.sequence)
    else:
        raise TypeError('not Bio.SeqRecord.SeqRecord or Bio.GenBank.Record.Record')

def record_coords(record, start, end, strand=0, upper=0):
    """
    >>> from Bio.SeqRecord import SeqRecord
    >>> record = SeqRecord(Seq("gagttttatcgcttccatga"),
    ...                    "ref|NC_001422",
    ...                    "Coliphage phiX174, complete genome",
    ...                    "bases 1-11")
    >>> record_coords(record, 0, 4) # zero-based
    Seq('GAGT', Alphabet())
    >>> record_coords(record, 0, 4, "-") # zero-based
    Seq('ACTC', Alphabet())
    >>> record_coords(record, 0, 4, "-", upper=1) # zero-based
    Seq('ACTC', Alphabet())
    """

    subseq = record_sequence(record)[start:end]
    subseq_str = subseq.tostring()
    subseq_str = subseq_str.upper()
    subseq = Seq(subseq_str, subseq.alphabet)
    if strand == '-' or strand == 1:
        return subseq.reverse_complement()
    else:
        return subseq

def _test():
    """Run the Bio.GFF.easy module's doctests (PRIVATE).

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests")) :
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    if __debug__:
        _test()
