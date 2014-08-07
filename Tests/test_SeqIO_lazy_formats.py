import unittest

import sys
import imp
import os
import tempfile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqIO import InsdcIO, UniprotIO, FastaIO, _lazy

test_files = [ \
    ("GenBank", ["brca_FJ940752.gb"], 'genbank'),
    ("GenBank", ["pri1.gb"], 'genbank'),
    ("GenBank", ["cor6_6.gb"], 'genbank'),
    ("GenBank", ["cor6_6.gb", "pri1.gb", "brca_FJ940752.gb"], 'genbank'),
    ("GenBank", ["1MRR_A.gp"], 'genbank'),
    ("EMBL", ["AAA03323.embl"], "embl"),
    ("EMBL", ["SC10H5.embl"], "embl"),
    ("EMBL", ["AE017046.embl"], "embl"),
    ("EMBL", ["patents.embl"], "embl"),
    #this is important because the default cache size if 5
    ("Fasta", ["f001", "f002", "elderberry.nu", \
               "loveliesbleeding.pro", "rose.pro", "rosemary.pro"], "fasta"),
    ("Fasta", ["f001"], "fasta"),
    ("Fasta", ["f002"], "fasta"),
    ("Fasta", ["elderberry.nu"], "fasta"),
    ("SwissProt", ["uni001"], "uniprot-xml")]

test_only_lazy = [ \
    ("SwissProt", ["uni001", "uni002", "uni003"], "uniprot-xml")]

class TestMultipleFormats(unittest.TestCase):

    _returnclasses = {"genbank": InsdcIO.GenbankSeqRecProxy,
                      "uniprot-xml": UniprotIO.UniprotXMLSeqRecProxy,
                      "fasta": FastaIO.FastaSeqRecProxy,
                      "embl": InsdcIO.EmblSeqRecProxy }

    def _multifile_iter_old(self, files, format):
        for filename in files:
            recorditer = SeqIO.parse(filename, format)
            for record in recorditer:
                yield record

    def features_to_strings(self, featurelist):
        returnlist = [(f.location.nofuzzy_start,
                       f.location.nofuzzy_end,
                       f.type) for f in featurelist]
        #returnlist.sort(key = lambda t: t[0]+(t[1]/100.0) + \
        #                  sum(ord(c) for c in t[2])/100000)

        returnlist.sort(key = lambda t: t[0])
        returnlist.sort(key = lambda t: t[1])
        returnlist.sort(key = lambda t: str(t[2]))
        return returnlist

    def _test_iter(self, folder, files, format):
        files = [os.path.join(folder, f) for f in files]
        memoryparser = list(self._multifile_iter_old(files, format))
        returnclass = self._returnclasses[format]
        lazyparser = iter(_lazy.LazyIterator(files, returnclass))
        for lazy, default in zip(lazyparser, memoryparser):
            self.compare_lazy_and_default(lazy, default)

    def compare_lazy_and_default(self, lazy, default):
        #check name, id
        self.assertEqual(default.id, lazy.id)
        self.assertEqual(default.name, lazy.name)
        #check self reported length
        self.assertEqual(len(default), len(lazy))
        #check seq slicing and reading
        if len(default) > 100:
            self.assertEqual(str(default[50:100].seq),
                             str(lazy[50:100].seq))
            self.assertEqual(str(default[:65].seq),
                             str(lazy[:65].seq))
            self.assertEqual(str(default[50:].seq),
                             str(lazy[50:].seq))
            self.assertEqual(str(default.seq), str(lazy.seq))
            self.assertEqual(str(default[50:100].seq),
                             str(lazy[50:100].seq))
            self.assertEqual(str(default[:65].seq),
                             str(lazy[:65].seq))
            self.assertEqual(str(default[50:].seq),
                             str(lazy[50:].seq))
            #check alphabet
            self.assertTrue(default.seq.alphabet is not None)
            self.assertTrue(lazy.seq.alphabet is not None)
            self.assertEqual(repr(default.seq.alphabet),
                             repr(lazy.seq.alphabet))
        else:
            self.assertEqual(str(default.seq), str(lazy.seq))
        #now check features
        fttostrings = self.features_to_strings
        if len(default) > 100:
            self.features_to_strings
            self.assertEqual(fttostrings(lazy[:].features),
                             fttostrings(default[:].features))
            self.assertEqual(fttostrings(lazy[99:].features),
                             fttostrings(default[99:].features))
            self.assertEqual(fttostrings(lazy[:99].features),
                             fttostrings(default[:99].features))
            self.assertEqual(fttostrings(lazy[:].features),
                             fttostrings(default[:].features))
            self.assertEqual(fttostrings(lazy[:].features),
                             fttostrings(default[:].features))
            self.assertEqual(fttostrings(lazy[99:].features),
                             fttostrings(default[99:].features))
            self.assertEqual(fttostrings(lazy[:99].features),
                             fttostrings(default[:99].features))
        else:
            self.assertEqual(fttostrings(lazy[:].features),
                             fttostrings(default[:].features))
        #check annotations
        try:
            self.assertEqual(len(default.annotations), len(lazy.annotations))
            oldkeys = [repr(v) for v in default.annotations.keys()]
            lzykeys = [repr(v) for v in lazy.annotations.keys()]
            oldkeys.sort()
            lzykeys.sort()
            self.assertEqual(oldkeys, lzykeys)
            oldvals = [repr(v) for v in default.annotations.values()]
            lzyvals = [repr(v) for v in lazy.annotations.values()]
            oldvals.sort()
            lzyvals.sort
        except TypeError:
            #in some cases None is a valid value for annotations,
            # The type error will force comparison for simple cases
            self.assertEqual(default.annotations, lazy.annotations)

    def _test_lazy_db_only(self, folder, files, format):
        oshandlelazy, db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        returnclass = self._returnclasses[format]
        files = [os.path.join(folder, f) for f in files]
        lazyparser = _lazy.LazyIterator(files, returnclass, index=db_name)

        keys = lazyparser.keys()
        keys.sort()
        for k in keys:
            rec = lazyparser[k]
            self.assertTrue(isinstance(rec, returnclass))
            self.assertTrue(isinstance(rec[0:99].seq, Seq))
            self.assertTrue(isinstance(rec.seq, Seq))
            self.assertTrue(isinstance(rec[0:99].seq, Seq))
            for f in rec[0:99].features:
                self.assertTrue(isinstance(f, SeqFeature))
        #cleanup
        del lazyparser
        os.remove(db_name)



    def _test_index_db(self, folder, files, format):
        oshandlelazy, lazy_db_name = tempfile.mkstemp()
        oshandleold, old_db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        os.close(oshandleold)
        os.remove(old_db_name)

        files = [os.path.join(folder, f) for f in files]
        old_db = SeqIO.index_db(old_db_name, files, format)
        returnclass = self._returnclasses[format]
        lazyparser = _lazy.LazyIterator(files, returnclass, index=lazy_db_name)
        newkeylist = [k for k in lazyparser.keys()]
        oldkeylist = [k for k in old_db.keys()]
        newkeylist.sort()
        oldkeylist.sort()
        self.assertEqual(newkeylist, oldkeylist)
        for key in newkeylist:
            lazy = lazyparser[key]
            default = old_db[key]
            self.assertTrue(isinstance(default, SeqRecord))
            self.assertTrue(isinstance(lazy, returnclass))
            self.compare_lazy_and_default(lazy, default)
        del(lazy)
        del(default)
        del(lazyparser)
        del(old_db)
        os.remove(lazy_db_name)
        os.remove(old_db_name)

    def _test_premade_db_use(self, folder, files, format):
        pass
        """    def setUp(self):
        #db setup: make tempfile and close os handle to help with cleanup
        oshandle, self.dbfilename = tempfile.mkstemp()
        recordfile = os.path.join("GenBank", self.recordfile)
        os.close(oshandle)
        #iter setup lazy
        tempopen = open(recordfile, 'rb')
        #make index
        returncls = SeqIO.InsdcIO.GenbankSeqRecProxy
        SeqIO._lazy._make_index_db(handle = tempopen,
                             return_class = returncls,
                             indexdb = self.dbfilename,
                             format = 'genbank' )
        tempopen.close()
        self.file = open(recordfile, 'rb')
        self.handle = self.file
        self.parser = lambda handle: iter(SeqIO._lazy.LazyIterator(handle, \
                                        returncls, index=self.dbfilename))
        self.lazy_iter = self.parser(self.file)
        self.standard_iter = SeqIO.parse(recordfile, 'genbank')
        self.oldrec = next(self.standard_iter)"""

#data sets for comparitive testing
for (folder, files, format) in test_files:
    name = "_".join([f.split('.')[0] for f in files])

    def funct(folder, files, format):
        f = lambda x : x._test_iter(folder, files, format)
        f.__doc__ = "Checking files %s of format %s"%(name, format)
        return f

    def dbtest(folder, files, format):
        f = lambda x : x._test_index_db(folder, files, format)
        f.__doc__ = "Checking lazy db dict for %s, format %s"%(name, format)
        return f

    setattr(TestMultipleFormats, "test_nuc_%s"%name,
        funct(folder, files, format))
    setattr(TestMultipleFormats, "test_indexdb_%s"%name,
        dbtest(folder, files, format))
    del funct
    del dbtest

#data sets for lazy loading only
for (folder, files, format) in test_only_lazy:
    name = "_".join([f.split('.')[0] for f in files])

    def lazyonly(folder, files, format):
        f = lambda x : x._test_lazy_db_only(folder, files, format)
        f.__doc__ = "Checking files %s of format %s"%(name, format)
        return f

    setattr(TestMultipleFormats, "test_lazy_%s"%name,
        lazyonly(folder, files, format))
    del lazyonly



if __name__ == "__main__":
    #Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
