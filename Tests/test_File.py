# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import string
from TestSupport import verbose, TestFailed
from Bio import File



### StringHandle

if verbose:
    print "Running tests on StringHandle"

data = """This
is
a multi-line
file"""

h = File.StringHandle(data)
try:
    assert string.rstrip(h.readline()) == 'This', "readline"
    assert len(h.readlines()) == 3, "readlines"
    assert h.readline() == ''
    h.close()
except Exception, x:
    raise TestFailed, "StringHandle (%s)" % x



### UndoHandle

if verbose:
    print "Running tests on UndoHandle"

data = """This
is
a multi-line
file"""

# os.pipe is not available on MS-DOS.  If DOS compatibility turns
# out to be important, we may have to save the data to a temporary
# file and make a real file handle.
r, w = os.pipe()
os.fdopen(w, 'w').write(data)
h = File.UndoHandle(os.fdopen(r, 'r'))

try:
    assert string.rstrip(h.readline()) == 'This'
    assert string.rstrip(h.peekline()) == 'is'
    assert string.rstrip(h.readline()) == 'is'
    h.saveline("saved")
    assert h.peekline() == 'saved'
    h.saveline("another")
    assert h.readline() == 'another'
    assert h.readline() == 'saved'

    # Test readlines after saveline
    h.saveline("saved again")
    lines = h.readlines()
    assert string.strip(lines[0]) == 'saved again'
    assert string.strip(lines[1]) == 'a multi-line'
    assert string.strip(lines[2]) == 'file'
    
    # should be empty now
    assert h.readline() == ''
    
    h.saveline("save after empty")
    assert h.readline() == 'save after empty'
    assert h.readline() == ''
    h.close()  # should pass this to the original handle
except Exception, x:
    raise TestFailed, "UndoHandle (%s)" % x
