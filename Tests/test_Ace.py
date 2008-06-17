import os
import sys
import unittest

from Bio.Sequencing import Ace

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [AceTestOne, AceTestTwo, AceTestThree]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class AceTestOne(unittest.TestCase):
    def setUp(self):
        self.handle = open("Ace/contig1.ace")

    def tearDown(self):
        self.handle.close()

    def contig_summary(self,c,x):
        """Print a summary of a contig."""
        print '\n--------------------'
        print 'Contig no.',x,'has',len(c.reads),'reads.'
        print 'CO:',c.name,c.nbases,c.nreads,c.nsegments,c.uorc
        print c.sequence[:10],'...',c.sequence[len(c.sequence)/2-5:len(c.sequence)/2+5],'...',c.sequence[-10:]
        print c.quality[:10],'...',c.quality[len(c.quality)/2-5:len(c.quality)/2+5],'...',c.quality[-10:]
        print '# AF:',len(c.af),'; # BS:',len(c.bs)
        for i in (len(c.af)/2,-1):
            print 'AF(%d) %s %s %d' % (i,c.af[i].name,c.af[i].coru,c.af[i].padded_start),';',
        print
        if c.bs :
            for i in (len(c.bs)/2,-1):
                print 'BS(%d) %s %d %d' % (i,c.bs[i].name,c.bs[i].padded_start,c.bs[i].padded_end),';',
            print
        else :
            print 'BS: none'
        print 'CT:',
        if not c.ct:
            print 'none'
        else:
            print 
            for i in c.ct:
                print i.name, i.tag_type, i.program, i.padded_start, i.padded_end, i.date, i.info
        print 'WA:',
        if not c.wa:
            print 'none'
        else:
            print 
            for i in c.wa:
                print i.tag_type, i.program, i.date, i.info
        for r in c.reads:
            print 'RD:',r.rd.name,r.rd.padded_bases,r.rd.info_items,r.rd.read_tags,
            print r.rd.sequence[:10],'...',r.rd.sequence[len(r.rd.sequence)/2-5:len(r.rd.sequence)/2+5],'...',r.rd.sequence[-10:]
            print 'QA:',r.qa.qual_clipping_start,r.qa.qual_clipping_end,r.qa.align_clipping_start,r.qa.align_clipping_end
            if r.ds:
                print 'DS:','CHROMAT_FILE:',r.ds.chromat_file,'PHD_FILE:',r.ds.phd_file,'TIME:',r.ds.time,
                print 'CHEM:',r.ds.chem,'DYE:',r.ds.dye,'TEMPLATE:',r.ds.template,'DIRECTION:',r.ds.direction
            else:
                print 'DS: none'
            print 'RT:',
            if not r.rt:
                print 'none'
            else:
                print
                for i in r.rt:
                    print i.name, i.tag_type, i.program, i.padded_start, i.padded_end, i.date
            print 'WR:',
            if not r.wr:
                print 'none'
            else:
                print
                for i in r.wr:
                    print i.name, i.aligned, i.program, i.date

    def t_check_ACEParser(self):
        """Test to check that ACEParser can parse the whole file into one record."""
        aceparser=Ace.ACEParser()
        record=aceparser.parse(self.handle)
        x=0
        print '\nInput file: %s' % self.handle.name
        print '\nContigs:',record.ncontigs,'; Reads:',record.nreads
        print 'all WA:',
        if not record.wa:
            print 'none'
        else:
            print 
            for i in record.wa:
                print i.tag_type, i.program, i.date, i.info
        for c in record.contigs:
            self.contig_summary(c,x)
            x+=1
    
    def t_check_record_parser(self):
        """Test to check that record parser parses each contig into a record."""
        recparser=Ace.RecordParser()
        it=Ace.Iterator(self.handle,recparser)
        x=0
        print '\nInput file: %s' % self.handle.name
        while 1:
            r=it.next()
            if not r:
                break
            self.contig_summary(r,x)
            x+=1


class AceTestTwo(AceTestOne) :
    """Test parsing example output from CAP3.

    The sample input file seq.cap.ace was downloaded from:
    http://genome.cs.mtu.edu/cap/data/seq.cap.ace
    """
    def setUp(self):
        self.handle = open("Ace/seq.cap.ace")

class AceTestThree(AceTestOne) :
    """Test parsing example ACE input file for CONSED.

    The sample input file was downloaded from:
    http://bozeman.mbt.washington.edu/consed/distributions/README.16.0.txt
    """
    def setUp(self):
        self.handle = open("Ace/consed_sample.ace")


if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))

