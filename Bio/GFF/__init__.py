#!/usr/bin/env python

"""
GFF.py: Access to General Feature Format databases created with Bio::DB:GFF

based on documentation for Lincoln Stein's Perl Bio::DB::GFF
"""

__version__ = "$Revision: 1.1 $"
# $Source: /home/bartek/cvs2bzr/biopython_fastimport/cvs_repo/biopython/Bio/GFF/__init__.py,v $

import exceptions
import operator
import os.path
import sys
import types

import MySQLdb

from Bio.Alphabet import IUPAC
from Bio import DocSQL
import Bio.easy
from Bio import GenericTools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Translate

import binning

DEFAULT_ALPHABET = IUPAC.unambiguous_dna
standard_translator = Translate.unambiguous_dna_by_id[1]

class Segment(object):
    """
    this will only work for the simplest of Bio.easy.Location objects
    """
    def __init__(self, gff, location=None):
        self.gff = gff
        if location is None:
            self.seqname = None
            self.coords = None
            self.complement = None
        else:
            self.seqname = location.seqname
            self.coords = [location.start(), location.end()]
            self.complement = location.complement

    def features(self, *args, **keywds):
        return FeatureQuery(self.seqname, self.coords, self.complement, connection=self.gff.db, *args, **keywds)

class Connection(Segment):
    """
    Connection to GFF database

    also functions as whole segment
    """
    
    def __init__(self, *args, **keywds):
        try:
            RetrieveSeqname._dir = keywds['fastadir']
            del keywds['fastadir']
        except KeyError:
            RetrieveSeqname._dir = '.'
        self.db = MySQLdb.connect(*args, **keywds)
        Segment.__init__(self, self)

    def segment(self, *args, **keywds):
        return Segment(self, *args, **keywds)

class RetrieveSeqname(GenericTools.Surrogate, SeqRecord):
    """
    Singleton: contain records of loaded FASTA files

    >>> RetrieveSeqname._dir = '/home/grouse/lib/python/t/Profiler/BLAST'
    >>> RetrieveSeqname._diagnostics = 1
    >>> sys.stderr = sys.stdout # dump diagnostics to stdout so doctest can see them
    >>> record1 = RetrieveSeqname('NC_001802.fna')
    Loading /home/grouse/lib/python/t/Profiler/BLAST/NC_001802.fna... Done.
    >>> record1.id
    'gi|9629357|ref|NC_001802.1|'
    >>> record2 = RetrieveSeqname('NC_001802.fna')
    >>> record1.seq is record2.seq # should be same space in memory
    1
    """
    __records = {}
    _diagnostics = 0

    def __init__(self, seqname):
        try:
            GenericTools.Surrogate.__init__(self, self.__records[seqname])
        except KeyError:
            filename = os.path.join(self._dir, seqname)
            if self._diagnostics:
                sys.stderr.write("Loading %s..." % filename)
            record = Bio.easy.fasta_single(filename)
            self.__records[seqname] = record
            GenericTools.Surrogate.__init__(self, self.__records[seqname])
            if self._diagnostics:
                print >>sys.stderr, " Done."

class Feature(object):
    """
    strand may be:
    +/0 = Watson
    -/1 = Crick

    I propose that we start calling these the Rosalind and Franklin strands
    
    >>> RetrieveSeqname._dir = '/home/grouse/lib/python/grouse-misc-python/tests/'
    >>> feature = Feature("NC_001802x.fna", 73, 78) # not having the x will interfere with the RetrieveSequence test
    >>> feature.seq()
    Seq('AATAAA', Alphabet())
    >>> print feature.location()
    NC_001802x.fna:73..78
    >>> from Bio.SeqIO.FASTA import FastaWriter
    >>> writer = FastaWriter(sys.stdout)
    >>> writer.write(feature.record())
    > NC_001802x.fna:73..78
    AATAAA
    >>> feature2 = Feature(location=Bio.easy.LocationFromString("NC_001802x.fna:73..78"))
    >>> writer.write(feature2.record())
    > NC_001802x.fna:73..78
    AATAAA
    >>> location3 = Bio.easy.LocationFromString("NC_001802x.fna:complement(73..78)")
    >>> feature3 = Feature(location=location3)
    >>> writer.write(feature3.record())
    > NC_001802x.fna:complement(73..78)
    TTTATT
    >>> location4 = Bio.easy.LocationFromString("NC_001802x.fna:336..1631")
    >>> feature4 = Feature(location=location4, frame=0)
    >>> feature4.frame
    0
    >>> feature4.translate()[:7]
    Seq('MGARASV', HasStopCodon(IUPACProtein(), '*'))
    >>> feature4.frame = 6 # can't happen, but a useful demonstration
    >>> feature4.translate()[:5]
    Seq('ARASV', HasStopCodon(IUPACProtein(), '*'))
    >>> feature4.frame = 1
    >>> feature4.translate()[:5]
    Seq('WVRER', HasStopCodon(IUPACProtein(), '*'))
    >>> location5 = Bio.easy.LocationFromString("NC_001802lc.fna:336..1631") # lowercase data
    >>> feature5 = Feature(location=location5, frame=0)
    >>> feature5.translate()[:7]
    Seq('MGARASV', HasStopCodon(IUPACProtein(), '*'))
    >>> location6 = Bio.easy.LocationFromString("NC_001802lc.fna:335..351")
    >>> feature6 = Feature(location=location6, frame=1)
    >>> feature6.translate()
    Seq('MGARA', HasStopCodon(IUPACProtein(), '*'))
    """
    def __init__(self,
                 seqname=None,
                 start=None,
                 end=None,
                 strand="+",
                 frame=None,
                 location=None,
                 alphabet=DEFAULT_ALPHABET):
        if not isinstance(seqname, (types.StringType, types.NoneType)):
            raise exceptions.TypeError, "seqname needs to be string"
        self.frame = frame
        self.alphabet = alphabet
        if location is None:
            self.seqname = seqname
            self.start = start
            self.end = end
            self.strand = strand
            return
        else:
            self.seqname = location.seqname
            self.start = location.start() + 1
            self.end = location.end() + 1
            self.strand = location.complement

    def seq(self):
        rec = RetrieveSeqname(self.seqname)
        return Bio.easy.record_subseq(rec, self.location(), upper=1)

    def translate(self):
        seq = self.seq()
        try:
            seq = Seq(seq.tostring()[self.frame:], self.alphabet)
        except TypeError:
            seq.alphabet = self.alphabet
        try:
            return standard_translator.translate(seq)
        except AssertionError:
            # if the feature was pickled then we have problems
            import cPickle
            if cPickle.dumps(seq.alphabet) == cPickle.dumps(DEFAULT_ALPHABET):
                seq.alphabet = DEFAULT_ALPHABET
                return standard_translator.translate(seq)
            else:
                raise

    def location(self):
        return Bio.easy.LocationFromCoords(self.start-1, self.end-1, self.strand, seqname=self.seqname)

    def target_location(self):
        if self.target_start <= self.target_end:
            return Bio.easy.LocationFromCoords(self.target_start-1, self.target_end-1, 0, seqname=self.gname)
        else:
            # switch start and end and make it complement:
            return Bio.easy.LocationFromCoords(self.target_end-1, self.target_start-1, 1, seqname=self.gname)

    def id(self):
        try:
            return "%s:%s" % (self.gclass, self.gname)
        except AttributeError:
            return ""

    def record(self):
        return SeqRecord(self.seq(), self.id(), "", self.location())

    def __str__(self):
        return "Feature(%s)" % self.location()

class FeatureQueryRow(DocSQL.QueryRow, Feature):
    """
    row of FeatureQuery results
    works like a Feature
    """
    def __init__(self, *args, **keywds):
        DocSQL.QueryRow.__init__(self, *args, **keywds)
        try:
            self.frame = int(self.frame)
        except ValueError:
            self.frame = ''
        except TypeError:
            self.frame = None
        self.alphabet = DEFAULT_ALPHABET

class FeatureQuery(DocSQL.Query):
    """
    SELECT fdata.fref AS seqname,
        ftype.fsource AS source,
        ftype.fmethod AS feature,
        fdata.fstart AS start,
        fdata.fstop AS end,
        fdata.fscore AS score,
        fdata.fstrand AS strand,
        fdata.fphase AS frame,
        fdata.ftarget_start AS target_start,
        fdata.ftarget_stop AS target_end,
        fdata.gid,
        fgroup.gclass,
        fgroup.gname
    FROM fdata, ftype, fgroup
    WHERE fdata.ftypeid = ftype.ftypeid
        AND fdata.gid = fgroup.gid
    """
    
    row_class = FeatureQueryRow
    def __init__(self,
                 seqname=None,
                 coords=None,
                 complement=0,
                 method=None,
                 source=None,
                 object=None, # "class:name"
                 *args,
                 **keywds):
        
        DocSQL.Query.__init__(self, *args, **keywds)

        if seqname is not None:
            self.statement += ' AND fref="%s"\n' % seqname
        if coords is not None:
            self.statement += " AND (%s)\n" % binning.query(coords[0], coords[1])
            if coords[0] is not None:
                self.statement += (" AND fstop >= %s\n" % coords[0])
            if coords[1] is not None:
                self.statement += (" AND fstart <= %s\n" % coords[1])
        if method is not None:
            self.statement += ' AND fmethod LIKE "%s"\n' % method # LIKE allows SQL "%" wildcards
        if source is not None:
            self.statement += ' AND fsource LIKE "%s"\n' % source
        if object is not None:
            self.statement += ' AND %s\n' % object2fgroup_sql(object)
            
        if complement:
            self.statement += " ORDER BY 0-fstart, 0-fstop"
        else:
            self.statement += " ORDER BY fstart, fstop"

    def aggregate(self):
        return FeatureAggregate(self)

def object2fgroup_sql(object):
    return 'gclass LIKE "%s" AND gname LIKE "%s"' % object_split(object)

class FeatureAggregate(list, Feature):
    """
    >>> feature1_1 = Feature(location=Bio.easy.LocationFromString("NC_001802x.fna:336..1631"), frame=0) # gag-pol
    >>> feature1_2 = Feature(location=Bio.easy.LocationFromString("NC_001802x.fna:1631..4642"), frame=0) # slippage
    >>> aggregate = FeatureAggregate([feature1_1, feature1_2])
    >>> print aggregate.location()
    join(NC_001802x.fna:336..1631,NC_001802x.fna:1631..4642)
    >>> xlate_str = aggregate.translate().tostring()
    >>> xlate_str[:5], xlate_str[-5:]
    ('MGARA', 'RQDED')
    
    >>> location1 = Bio.easy.LocationFromString("NC_001802x.fna:complement(1..6)")
    >>> location2 = Bio.easy.LocationFromString("NC_001802x.fna:complement(7..12)")
    >>> feature2_1 = Feature(location=location1, frame=0)
    >>> feature2_2 = Feature(location=location2, frame=0)
    >>> aggregate2 = FeatureAggregate([feature2_1, feature2_2])
    >>> print aggregate2.location()
    complement(join(NC_001802x.fna:1..6,NC_001802x.fna:7..12))
    >>> print aggregate2.translate()
    Seq('TRET', HasStopCodon(IUPACProtein(), '*'))
    >>> location1.reverse()
    >>> location2.reverse()
    >>> aggregate3 = FeatureAggregate([Feature(location=x, frame=0) for x in [location1, location2]])
    >>> print aggregate3.location()
    join(NC_001802x.fna:1..6,NC_001802x.fna:7..12)
    >>> print aggregate3.translate()
    Seq('GLSG', HasStopCodon(IUPACProtein(), '*'))
    >>> aggregate3[0].frame = 3
    >>> print aggregate3.translate()
    Seq('LSG', HasStopCodon(IUPACProtein(), '*'))
    
    >>> aggregate4 = FeatureAggregate()
    >>> aggregate4.append(Feature(location=Bio.easy.LocationFromString("NC_001802x.fna:1..5"), frame=0))
    >>> aggregate4.append(Feature(location=Bio.easy.LocationFromString("NC_001802x.fna:6..12"), frame=2))
    >>> aggregate4.seq()
    Seq('GGTCTCTCTGGT', Alphabet())
    >>> aggregate4.translate()
    Seq('GLSG', HasStopCodon(IUPACProtein(), '*'))
    """
    def __init__(self, feature_query=None):
        if feature_query is None:
            list.__init__(self, [])
        else:
            list.__init__(self, map(lambda x: x, feature_query))

    def location(self):
        loc = Bio.easy.LocationJoin(map(lambda x: x.location(), self))
        loc.reorient()
        return loc

    def map(self, func):
        mapped = map(func, self)
        if self.location().complement:
            mapped.reverse()
        return reduce(operator.add, mapped)

    def seq(self):
        return self.map(lambda x: x.seq())

    def translate(self):
        attributes = ['alphabet', 'frame']
        index = [0, -1][self.location().complement]
        for attribute in attributes:
            self.__dict__[attribute] = self[index].__dict__[attribute]
        translation = Feature.translate(self)
        try:
            assert translation.tostring().index("*") == len(translation) - 1
            return translation[:translation.tostring().index("*")]
        except ValueError:
            return translation

def object_split(object):
    """
    >>> object_split("Sequence:F02E9.2a")
    ('Sequence', 'F02E9.2a')
    """
    return tuple(object.split(":"))

def _test(*args, **keywds):
    import doctest, sys
    doctest.testmod(sys.modules[__name__], *args, **keywds)

if __name__ == "__main__":
    if __debug__:
        _test()
