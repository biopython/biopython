# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import string
import sys
from Bio import File
from Bio import ParserSupport

# pyUnit

def pb(b):
    if b:
        return 1
    return 0

### TaggingConsumer

print "Running tests on TaggingConsumer"

class TestHandle:
    def write(self, s):
        print s
        
h = TestHandle()
tc = ParserSupport.TaggingConsumer(handle=h, colwidth=5)
tc.start_section()  # '***** start_section\n'
tc.test1('myline')  # 'test1: myline\n'
tc.end_section()    # '***** end_section\n'




### is_blank_line

print "Running tests on is_blank_line"

is_blank_line = lambda *args, **keywds: \
                pb(ParserSupport.is_blank_line(*args, **keywds))

print is_blank_line('\n')                              # 1
print is_blank_line('\r\n')                            # 1
print is_blank_line('\r')                              # 1
print is_blank_line('')                                # 1
print is_blank_line('', allow_spaces=1)                # 1
print is_blank_line('', allow_spaces=0)                # 1
print is_blank_line(string.whitespace, allow_spaces=1) # 1
print is_blank_line('hello')                           # 0
print is_blank_line('hello', allow_spaces=1)           # 0
print is_blank_line('hello', allow_spaces=0)           # 0
print is_blank_line(string.whitespace, allow_spaces=0) # 0


### safe_readline

print "Running tests on safe_readline"

data = """This
file"""

h = File.UndoHandle(File.StringHandle(data))

safe_readline = ParserSupport.safe_readline
print safe_readline(h)    # "This"
print safe_readline(h)    # "file"
try: safe_readline(h)
except ValueError: print "correctly failed"
else: print "ERROR, should have failed"


### safe_peekline

print "Running tests on safe_peekline"
safe_peekline = ParserSupport.safe_peekline

data = """This
file"""

h = File.UndoHandle(File.StringHandle(data))

print safe_peekline(h) # "This"
h.readline()
print safe_peekline(h) # "file"
h.readline()
try: safe_peekline(h)
except ValueError: print "correctly failed"
else: print "ERROR, should have failed"
h.saveline('hello')
print safe_peekline(h) # 'hello'


### read_and_call

print "Running tests on read_and_call"

data = """>gi|132871|sp|P19947|RL30_BACSU 50S RIBOSOMAL PROTEIN L30 (BL27)
MAKLEITLKRSVIGRPEDQRVTVRTLGLKKTNQTVVHEDNAAIRGMINKVSHLVSVKEQ
>gi|132679|sp|P19946|RL15_BACSU 50S RIBOSOMAL PROTEIN L15
MKLHELKPSEGSRKTRNRVGRGIGSGNGKTAGKGHKGQNARSGGGVRPGFEGGQMPLFQRLPKRGFTNIN
RKEYAVVNLDKLNGFAEGTEVTPELLLETGVISKLNAGVKILGNGKLEKKLTVKANKFSASAKEAVEAAG
GTAEVI


"""

h = File.UndoHandle(File.StringHandle(data))

rac = ParserSupport.read_and_call
lines = []
def m(line):
    lines.append(line)
    
rac(h, m)
print lines[-1][:10]   # '>gi|132871'
rac(h, m, start='MAKLE', end='KEQ', contains='SVIG')
rac(h, m, blank=0)

# These should be errors.  If they're not, then complain.
try: rac(h, m, blank=1)
except ValueError: print "correctly failed"
else: print "ERROR, should have failed"
try: rac(h, m, start='foobar')
except ValueError: print "correctly failed"
else: print "ERROR, should have failed"
try: rac(h, m, end='foobar')
except ValueError: print "correctly failed"
else: print "ERROR, should have failed"
try: rac(h, m, contains='foobar')
except ValueError: print "correctly failed"
else: print "ERROR, should have failed"
try: rac(h, m, blank=0)
except ValueError: print "correctly failed"
else: print "ERROR, should have failed"
        



### attempt_read_and_call

print "Running tests on attempt_read_and_call"

data = """>gi|132871|sp|P19947|RL30_BACSU 50S RIBOSOMAL PROTEIN L30 (BL27)
MAKLEITLKRSVIGRPEDQRVTVRTLGLKKTNQTVVHEDNAAIRGMINKVSHLVSVKEQ
>gi|132679|sp|P19946|RL15_BACSU 50S RIBOSOMAL PROTEIN L15
MKLHELKPSEGSRKTRNRVGRGIGSGNGKTAGKGHKGQNARSGGGVRPGFEGGQMPLFQRLPKRGFTNIN
RKEYAVVNLDKLNGFAEGTEVTPELLLETGVISKLNAGVKILGNGKLEKKLTVKANKFSASAKEAVEAAG
GTAEVI"""

h = File.UndoHandle(File.StringHandle(data))

arac = lambda *args, **keywds: \
       pb(ParserSupport.attempt_read_and_call(*args, **keywds))
lines = []
def m(line):
    lines.append(line)

print arac(h, m, contains="RIBOSOMAL PROTEIN")   # 1
print arac(h, m, start="foobar")                 # 0
print arac(h, m, blank=1)                        # 0
print arac(h, m, end="LVSVKEQ")                  # 1

