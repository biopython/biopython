# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import string
import sys
from Bio import File
from Bio import ParserSupport

# pyUnit
import unittest

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

is_blank_line = ParserSupport.is_blank_line
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
except SyntaxError: print "correctly failed"
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
except SyntaxError: print "correctly failed"
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
except SyntaxError: print "correctly failed"
else: print "ERROR, should have failed"
try: rac(h, m, start='foobar')
except SyntaxError: print "correctly failed"
else: print "ERROR, should have failed"
try: rac(h, m, end='foobar')
except SyntaxError: print "correctly failed"
else: print "ERROR, should have failed"
try: rac(h, m, contains='foobar')
except SyntaxError: print "correctly failed"
else: print "ERROR, should have failed"
try: rac(h, m, blank=0)
except SyntaxError: print "correctly failed"
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

arac = ParserSupport.attempt_read_and_call
lines = []
def m(line):
    lines.append(line)

print arac(h, m, contains="RIBOSOMAL PROTEIN")   # 1
print arac(h, m, start="foobar")                 # 0
print arac(h, m, blank=1)                        # 0
print arac(h, m, end="LVSVKEQ")                  # 1

# --- EventGenerator
print "Running tests on EventGenerator"

class TestConsumer:
    """Collect information from callbacks from the EventGenerator.
    """
    def __init__(self):
        self.info = {}

    def main_tag(self, content):
        raise AssertionError("We should never call this.")

    def single_tag(self, content):
        self.info['single_tag'] = content

    def multiple_tags(self, content):
        self.info['multiple_tags'] = content

def strip_and_combine(line_list):
    """Function to be sure the optional finalizer works.
    """
    stripped_line_list = map(string.strip, line_list)
    return string.join(stripped_line_list, ',')

class EventGeneratorTest(unittest.TestCase):
    """Test the EventGenerator class to be sure callbacks are correct.
    """
    def setUp(self):
        self.interest_tags = ['single_tag', 'multiple_tags']

        import Martel

        test_format = \
           Martel.Group("main_tag",
                        Martel.Group("single_tag", Martel.Str("Spam")) +
                        Martel.Rep(Martel.Group("multiple_tags",
                                                Martel.Str("  Lots of Spam"))))
        self.test_parser = test_format.make_parser()

    def t_basic_callback(self):
        """Test the basic callback functionality.
        """
        consumer = TestConsumer()
        event_gen = ParserSupport.EventGenerator(consumer,
                                                 self.interest_tags)

        self.test_parser.setContentHandler(event_gen)
        self.test_parser.parseString("Spam" + "  Lots of Spam" +
                                     "  Lots of Spam")
        
        assert consumer.info["single_tag"] == ["Spam"], \
               "Single tag parsing failed: %s" % consumer.info["single_tag"]
        assert consumer.info["multiple_tags"] == \
               ["  Lots of Spam", "  Lots of Spam"], \
               "Multiple tag parsing failed: %s" % \
               consumer.info["multiple_tags"]

    def t_finalizer_callback(self):
        """Test the ability to pass a finalizer to the consumer.
        """
        consumer = TestConsumer()
        event_gen = ParserSupport.EventGenerator(consumer, self.interest_tags,
                                                 strip_and_combine)
        
        self.test_parser.setContentHandler(event_gen)
        self.test_parser.parseString("Spam" + "  Lots of Spam" +
                                     "  Lots of Spam")
        
        assert consumer.info["single_tag"] == "Spam", \
               "Single tag parsing failed: %s" % consumer.info["single_tag"]
        assert consumer.info["multiple_tags"] == "Lots of Spam,Lots of Spam", \
               "Multiple tag parsing failed: %s" % \
               consumer.info["multiple_tags"]

    def t_exempt_finalizer_callback(self):
        """Test the ability to exempt tags from being processed by a finalizer.
        """
        consumer = TestConsumer()
        event_gen = ParserSupport.EventGenerator(consumer, self.interest_tags,
                                                 strip_and_combine,
                                                 ["multiple_tags"])
        self.test_parser.setContentHandler(event_gen)
        self.test_parser.parseString("Spam" + "  Lots of Spam" +
                                     "  Lots of Spam")

        assert consumer.info["single_tag"] == "Spam", \
               "Single tag parsing failed: %s" % consumer.info["single_tag"]
        assert consumer.info["multiple_tags"] == \
               ["  Lots of Spam", "  Lots of Spam"], \
               "Multiple tag parsing failed: %s" % \
               consumer.info["multiple_tags"]
        

all_tests = [EventGeneratorTest]

runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)

test_loader = unittest.defaultTestLoader
test_loader.testMethodPrefix = 't_'

for cur_test in all_tests:
    test_suite = test_loader.loadTestsFromTestCase(cur_test)
    runner.run(test_suite)

        
        
        


