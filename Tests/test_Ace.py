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
    tests = [AceTestOne]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class AceTestOne(unittest.TestCase):
    def setUp(self):
        self.handle = open("Ace/contig1.ace")

    def tearDown(self):
        self.handle.close()

    def t_check_ACEParser(self):
        """Test to check that ACEParser can parse the whole file into one record. 
        """
        aceparser=Ace.ACEParser()
        record=aceparser.parse(self.handle)
        record.sort()
        x=0
        print 'Contigs:',record.contigs,'; Reads:',record.reads
        for r in record.records:
            x+=1
            print '--------------------'
            print 'Contig no.',x,'has',len(r.rd),'reads.'
            print 'CO:',r.contig_name,r.bases,r.reads,r.segments,r.uorc
            print r.sequence[:10],'...',r.sequence[len(r.sequence)/2-5:len(r.sequence)/2+5],'...',r.sequence[-10:]
            print r.quality[:10],'...',r.quality[len(r.quality)/2-5:len(r.quality)/2+5],'...',r.quality[-10:]
            print '# AF:',len(r.af),'; # BS:',len(r.bs)
            for i in (0,len(r.af)/2,-1):
                print 'AF(%d) %s %s %d' % (i,r.af[i].name,r.af[i].coru,r.af[i].padded_start),';',
            print
            for i in (0,len(r.bs)/2,-1):
                print 'BS(%d) %s %d %d' % (i,r.bs[i].name,r.bs[i].padded_start,r.bs[i].padded_end),';',
            print
            if len(r.rd)==len(r.qa) and len(r.ds)==len(r.rd):
                print 'No. of reads in RD,QA,DS matches.'
            else:
                print 'No. of reads in RD,QA,DS does not match.'
            for i,(rd,qa,ds) in enumerate(zip(r.rd,r.qa,r.ds)):
                print 'RD:',rd.name,rd.padded_bases,rd.info_items,rd.read_tags,
                print rd.sequence[:10],'...',rd.sequence[len(rd.sequence)/2-5:len(rd.sequence)/2+5],'...',rd.sequence[-10:]
                print 'QA:',qa.qual_clipping_start,qa.qual_clipping_end,qa.align_clipping_start,qa.align_clipping_end
                print 'DS:','CHROMAT_FILE:',ds.chromat_file,'PHD_FILE:',ds.phd_file,'TIME:',ds.time,
                print 'CHEM:',ds.chem,'DYE:',ds.dye,'TEMPLATE:',ds.template,'DIRECTION:',ds.direction
            print 'RT:',
            if r.rt==[]:
                print 'none'
            else:
                print
                for i in r.rt:
                    print i.name, i.tag_type, i.program, i.padded_start, i.padded_end, i.date
            print 'WR:',
            if r.wr==[]:
                print 'none'
            else:
                print
                for i in r.wr:
                    print i.name, i.aligned, i.program, i.date
            print 'CT:',
            if r.ct==[]:
                print 'none'
            else:
                print 
                for i in r.ct:
                    print i.name, i.tag_type, i.program, i.padded_start, i.padded_end, i.date, i.info
            print 'WA:',
            if r.wa==[]:
                print 'none'
            else:
                print 
                for i in r.wa:
                    print i.tag_type, i.program, i.date, i.info
        print 'Unrelated tags'
        print 'RT:',
        if record.rt==[]:
            print 'none'
        else:
            print
            for i in record.rt:
                print i.name, i.tag_type, i.program, i.padded_start, i.padded_end, i.date
        print 'WR:',
        if record.wr==[]:
            print 'none'
        else:
            print 
            for i in record.wr:
                print i.name, i.aligned, i.program, i.date
        print 'CT:',
        if record.ct==[]:
            print 'none'
        else:
            print
            for i in record.ct:
                print i.name, i.tag_type, i.program, i.padded_start, i.padded_end, i.date, i.info
        print 'WA:',
        if record.wa==[]:
            print 'none'
        else:
            print 
            for i in record.wa:
                print i.tag_type, i.program, i.date, i.info

    def t_check_record_parser(self):
        """Test to check that record parser parses each contig into a record.
        """
        recparser = Ace.RecordParser()
        it = Ace.Iterator(self.handle,recparser)
        x=0
        while 1:
            r=it.next()
            if not r:
                break
            x+=1
            print '\n--------------------'
            print 'Contig no.',x,'has',len(r.rd),'reads.'
            print 'CO:',r.contig_name,r.bases,r.reads,r.segments,r.uorc
            print r.sequence[:10],'...',r.sequence[len(r.sequence)/2-5:len(r.sequence)/2+5],'...',r.sequence[-10:]
            print r.quality[:10],'...',r.quality[len(r.quality)/2-5:len(r.quality)/2+5],'...',r.quality[-10:]
            print '# AF:',len(r.af),'; # BS:',len(r.bs)
            for i in (0,len(r.af)/2,-1):
                print 'AF(%d) %s %s %d' % (i,r.af[i].name,r.af[i].coru,r.af[i].padded_start),';',
            print
            for i in (0,len(r.bs)/2,-1):
                print 'BS(%d) %s %d %d' % (i,r.bs[i].name,r.bs[i].padded_start,r.bs[i].padded_end),';',
            print
            if len(r.rd)==len(r.qa) and len(r.ds)==len(r.rd):
                print 'No. of reads in RD,QA,DS matches.'
            else:
                print 'No. of reads in RD,QA,DS does not match.'
            for i,(rd,qa,ds) in enumerate(zip(r.rd,r.qa,r.ds)):
                print 'RD:',rd.name,rd.padded_bases,rd.info_items,rd.read_tags,
                print rd.sequence[:10],'...',rd.sequence[len(rd.sequence)/2-5:len(rd.sequence)/2+5],'...',rd.sequence[-10:]
                print 'QA:',qa.qual_clipping_start,qa.qual_clipping_end,qa.align_clipping_start,qa.align_clipping_end
                print 'DS:','CHROMAT_FILE:',ds.chromat_file,'PHD_FILE:',ds.phd_file,'TIME:',ds.time,
                print 'CHEM:',ds.chem,'DYE:',ds.dye,'TEMPLATE:',ds.template,'DIRECTION:',ds.direction
            print 'RT:',
            if r.rt==[]:
                print 'none'
            else:
                print 
                for i in r.rt:
                    print i.name, i.tag_type, i.program, i.padded_start, i.padded_end, i.date
            print 'WR:',
            if r.wr==[]:
                print 'none'
            else:
                print
                for i in r.wr:
                    print i.name, i.aligned, i.program, i.date
            print 'CT:',
            if r.ct==[]:
                print 'none'
            else:
                print 
                for i in r.ct:
                    print i.name, i.tag_type, i.program, i.padded_start, i.padded_end, i.date, i.info
            print 'WA:',
            if r.wa==[]:
                print 'none'
            else:
                print 
                for i in r.wa:
                    print i.tag_type, i.program, i.date, i.info

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))

