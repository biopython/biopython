import unittest, sys, urllib2

import EUtils
from EUtils import ThinClient, DBIdsClient, Datatypes
import support_test


picklestore = None
PICKLESTORE_FILENAME = "test_DBIdsClient.pickle"

class TestSimple(unittest.TestCase):
    def setup(self):
        self.client = DBIdsClient.DBIdsClient(
            picklestore.client(1))

    def test1(self):
        self.setup()
        results = self.client.search("dalke", retmax = 50)
        if len(results) < 30:
            raise AssertionError("Should be at least 30 hits")

# Add more tests ....        


def main():
    global picklestore

    picklestore = support_test.UsePickleStore(PICKLESTORE_FILENAME)
    
    if len(sys.argv) == 2:
        if sys.argv[1] == "--use-live":
            picklestore = support_test.NoPickleStore()
            del sys.argv[1]
        elif sys.argv[1] == "--create-pickle":
            picklestore = support_test.CreatePickleStore(PICKLESTORE_FILENAME)
            del sys.argv[1]
        elif sys.argv[1] == "--use-pickle":
            picklestore = support_test.UsePickleStore(PICKLESTORE_FILENAME)
            del sys.argv[1]
        elif sys.argv[1] == "--help":
            print """\
Usage: %s [--use-live|--create-pickle|--use-pickle|--help]

If you specify nothing, this runs the unittest against the golden pickle file.
""" % (sys.argv[0],)
            sys.exit(0)
    try:
        unittest.main()
    finally:
        picklestore.done()
        
if __name__ == "__main__":
    #ThinClient.DUMP_URL = 1
    main()
