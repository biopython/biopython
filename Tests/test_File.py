# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import with_statement

import os.path
import unittest
import shutil
from StringIO import StringIO
import tempfile


from Bio import File



data = """This
is
a multi-line
file"""



### StringHandle

h = File.StringHandle(data)
print repr(h.readline())  # 'This'
print len(h.readlines())  # 3
print repr(h.readline())  # ''
h.close()



### UndoHandle

h = File.UndoHandle(File.StringHandle(data))

print h.readline()   # 'This'
print h.peekline()   # 'is'
print h.readline()   # 'is'
h.saveline("saved")
print h.peekline()   # 'saved'
h.saveline("another")
print h.readline()   # 'another'
print h.readline()   # 'saved'

# Test readlines after saveline
h.saveline("saved again")
lines = h.readlines()
print repr(lines[0])   # 'saved again'
print repr(lines[1])   # 'a multi-line'
print repr(lines[2])   # 'file'

# should be empty now
print repr(h.readline())       # ''

h.saveline("save after empty")
print h.readline()             # 'save after empty'
print repr(h.readline())       # ''

# test read method
h = File.UndoHandle(File.StringHandle("some text"))
h.saveline("more text")
print h.read()                 # 'more textsome text'

class AsHandleTestCase(unittest.TestCase):

    def setUp(self):
        # Create a directory to work in
        self.temp_dir = tempfile.mkdtemp(prefix='biopython-test')

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def _path(self, *args):
        return os.path.join(self.temp_dir, *args)

    def test_handle(self):
        "Test as_handle with a file-like object argument"
        p = self._path('test_file.fasta')
        with open(p, 'wb') as fp:
            with File.as_handle(fp) as handle:
                self.assertEqual(fp, handle, "as_handle should "
                        "return argument when given a file-like object")
                self.assertFalse(handle.closed)

            self.assertFalse(handle.closed,
                    "Exiting as_handle given a file-like object should not "
                    "close the file")

    def test_path(self):
        "Test as_handle with a path argument"
        p = self._path('test_file.fasta')
        mode = 'wb'
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(p, handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    def test_stringio(self):
        s = StringIO()
        with File.as_handle(s) as handle:
            self.assertEqual(s, handle)
