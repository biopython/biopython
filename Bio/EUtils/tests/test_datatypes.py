# Testing properties of datatypes

import unittest
from EUtils import Datatypes

class TestDBIds(unittest.TestCase):
    def testSubtraction(self):
        dbids1 = Datatypes.DBIds("protein", list("ABEFS"))
        dbids2 = Datatypes.DBIds("protein", list("AEF"))
        diff = dbids1-dbids2
        self.assertEquals(diff, Datatypes.DBIds('protein', ['B', 'S']))

if __name__ == "__main__":
    unittest.main()
