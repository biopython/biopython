# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import string
from TestSupport import verbose, TestFailed
from Bio import ParserSupport



### TaggingConsumer

if verbose:
    print "Running tests on TaggingConsumer"

class TestHandle:
    def __init__(self):
        self.written = []
    def write(self, s):
        self.written.append(s)

h = TestHandle()
tc = ParserSupport.TaggingConsumer(handle=h, colwidth=5)
tc.start_section()
tc.test1('myline')
tc.end_section()
try:
    assert h.written[0] == '***** start_section\n'
    assert h.written[1] == 'test1: myline\n'
    assert h.written[2] == '***** end_section\n'
except:
    raise TestFailed, "TaggingConsumer"



### OopsHandle

if verbose:
    print "Running tests on OopsHandle"

data = """This
is
a multi-line
file"""

# os.pipe is not available on MS-DOS.  If DOS compatibility turns
# out to be important, we may have to save the data to a temporary
# file and make a real file handle.
r, w = os.pipe()
os.fdopen(w, 'w').write(data)
h = ParserSupport.OopsHandle(os.fdopen(r, 'r'))

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
except:
    raise TestFailed, "OopsHandle"
    


### read_and_call

if verbose:
    print "Running tests on read_and_call"

data = """>gi|132871|sp|P19947|RL30_BACSU 50S RIBOSOMAL PROTEIN L30 (BL27)
MAKLEITLKRSVIGRPEDQRVTVRTLGLKKTNQTVVHEDNAAIRGMINKVSHLVSVKEQ
>gi|132679|sp|P19946|RL15_BACSU 50S RIBOSOMAL PROTEIN L15
MKLHELKPSEGSRKTRNRVGRGIGSGNGKTAGKGHKGQNARSGGGVRPGFEGGQMPLFQRLPKRGFTNIN
RKEYAVVNLDKLNGFAEGTEVTPELLLETGVISKLNAGVKILGNGKLEKKLTVKANKFSASAKEAVEAAG
GTAEVI"""

r, w = os.pipe()
os.fdopen(w, 'w').write(data)
h = ParserSupport.OopsHandle(os.fdopen(r, 'r'))

rac = ParserSupport.read_and_call
lines = []
def m(line):
    lines.append(line)
try:
    rac(h, m)
    assert lines[-1][:10] == '>gi|132871'
    rac(h, m, start='MAKLE', end='KEQ', contains='SVIG')

    # These should be errors.  If they're not, then complain.
    try: rac(h, m, start='foobar')
    except SyntaxError: pass
    else: assert 0
    try: rac(h, m, end='foobar')
    except SyntaxError: pass
    else: assert 0
    try: rac(h, m, contains='foobar')
    except SyntaxError: pass
    else: assert 0
    try: rac(h, m, blank=1)
    except SyntaxError: pass
    else: assert 0
except:
    raise TestFailed, "read_and_call"
        



### attempt_read_and_call

if verbose:
    print "Running tests on attempt_read_and_call"

data = """>gi|132871|sp|P19947|RL30_BACSU 50S RIBOSOMAL PROTEIN L30 (BL27)
MAKLEITLKRSVIGRPEDQRVTVRTLGLKKTNQTVVHEDNAAIRGMINKVSHLVSVKEQ
>gi|132679|sp|P19946|RL15_BACSU 50S RIBOSOMAL PROTEIN L15
MKLHELKPSEGSRKTRNRVGRGIGSGNGKTAGKGHKGQNARSGGGVRPGFEGGQMPLFQRLPKRGFTNIN
RKEYAVVNLDKLNGFAEGTEVTPELLLETGVISKLNAGVKILGNGKLEKKLTVKANKFSASAKEAVEAAG
GTAEVI"""

r, w = os.pipe()
os.fdopen(w, 'w').write(data)
h = ParserSupport.OopsHandle(os.fdopen(r, 'r'))

arac = ParserSupport.attempt_read_and_call
lines = []
def m(line):
    lines.append(line)

try:
    assert arac(h, m, contains="RIBOSOMAL PROTEIN")
    assert not arac(h, m, start="foobar") # make sure it undo's right
    assert not arac(h, m, blank=1)
    assert arac(h, m, end="LVSVKEQ")
except:
    raise TestFailed, "attempt_read_and_call"
