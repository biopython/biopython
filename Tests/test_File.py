# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import string
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
