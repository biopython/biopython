import unittest

from Bio.Sequencing import Phd


class PhdTestOne(unittest.TestCase):
    def setUp(self):
        self.handle = open("Phd/phd1")

    def tearDown(self):
        self.handle.close()

    def test_check_record_parser(self):
        """Test to check that record parser parses all records of a contig.
        """
        rs = Phd.parse(self.handle)
        x=0
        print
        for r in rs:
            x+=1
            print '---- Record',x,'----'
            print 'Filename:',r.file_name
            comments=r.comments.items()
            comments.sort()
            print '\n'.join([c[0]+': '+str(c[1]) for c in comments])
            allsites=['-'.join(s) for s in r.sites]
            print allsites[:10],'...',allsites[len(allsites)/2-5:len(allsites)/2+5],'...',allsites[-10:]
            print r.seq.tostring()[:10],r.seq.tostring()[-10:]
            print r.seq_trimmed.tostring()[:10],r.seq_trimmed.tostring()[-10:]
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
