
import unittest

def suite():
    modules_to_test = ('testCla', 'testDes', 'testHie', 'testScop',
                       'testResidues', 'testRaf') 
    alltests = unittest.TestSuite()
    for module in map(__import__, modules_to_test):
        alltests.addTest(unittest.findTestCases(module))
    return alltests

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
