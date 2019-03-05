import unittest

from Bio import SeqIO

class TestNibReader(unittest.TestCase):

    def test_read_bigendian(self):
        handle = open('Nib/test_bigendian.nib', 'rb')
        records = SeqIO.parse(handle, 'nib')
        record = next(records)
        handle.close()
        self.assertEqual(str(record.seq),
                         'ACGTAAACCGTACCCGTANANCANNNNACNANNANCN')

    def test_read_littleendian(self):
        handle = open('Nib/test_littleendian.nib', 'rb')
        records = SeqIO.parse(handle, 'nib')
        record = next(records)
        handle.close()
        self.assertEqual(str(record.seq),
                         'ACGTAAACCGTACCCGTANANCANNNNACNANNANCN')


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

