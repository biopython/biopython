#!/usr/bin/env python

# To do:
# - Let user specify the parser class on the command line.
# - Let user specify a sequence file to BLAST on the net.
# - Script should help debug connection to NCBI website.

import os, sys
import re
import getopt
import traceback

from Bio import ParserSupport
from Bio.Blast import NCBIStandalone, NCBIWWW

CONTEXT = 5   # show 5 lines of context around the error in the format file

USAGE = """%s [-h] [-v] [-p] [-n] [-o] <testfile>

This script helps diagnose problems with the BLAST parser.

OPTIONS:

-h    Show this help file.

-v    Verbose output.

-p    <testfile> is a protein file.

-n    <testfile> is a nucleotide file.

-o    <testfile> is a BLAST output file.

""" % sys.argv[0]

class DebuggingConsumer:
    def __init__(self, decorated=None):
        self.linenum = 0
        if decorated is None:
            decorated = ParserSupport.AbstractConsumer()
        self.decorated = decorated
        self._prev_attr = None
    def _decorated_section(self):
        getattr(self.decorated, self._prev_attr)()
    def _decorated(self, data):
        getattr(self.decorated, self._prev_attr)(data)
        self.linenum += 1
    def __getattr__(self, attr):
        self._prev_attr = attr
        if attr.startswith('start_') or attr.startswith('end_'):
            return self._decorated_section
        else:
            return self._decorated

def chomp(line):
    return re.sub(r"[\r\n]*$", "", line)

def choose_parser(outfile):
    data = open(outfile).read()
    ldata = data.lower()
    if ldata.find("<html>") >= 0 or ldata.find("<pre>") >= 0:
        return NCBIWWW.BlastParser
    if ldata.find("results from round") >= 0 or ldata.find("converged!") >= 0:
        return NCBIStandalone.PSIBlastParser
    return NCBIStandalone.BlastParser

def test_blast_output(outfile):
    # Try to auto-detect the format
    if 1:
        print "No parser specified.  I'll try to choose one for you based"
        print "on the format of the output file."
        print
        
        parser_class = choose_parser(outfile)
        print "It looks like you have given output that should be parsed"
        print "with %s.%s.  If I'm wrong, you can select the correct parser" %\
              (parser_class.__module__, parser_class.__name__)
        print "on the command line of this script (NOT IMPLEMENTED YET)."
    else:
        raise NotImplementedError
        parser_class = NCBIWWW.BlastParser
        print "Using %s to parse the file." % parser_class.__name__
    print

    scanner_class = parser_class()._scanner.__class__
    consumer_class = parser_class()._consumer.__class__

    #parser_class()._scanner.feed(
    #    open(outfile), ParserSupport.TaggingConsumer())
    print "I'm going to run the data through the parser to see what happens..."
    parser = parser_class()
    try:
        rec = parser.parse_file(outfile)
    except KeyboardInterrupt, SystemExit:
        raise
    except Exception, x:
        exception_info = str(x)
        print "Dang, the parsing failed."
    else:
        print "Parsing succeeded, no problems detected."
        print "However, you should check to make sure the following scanner"
        print "trace looks reasonable."
        print
        parser_class()._scanner.feed(
            open(outfile), ParserSupport.TaggingConsumer())
        return 0
    print

    print "Alright.  Let me try and figure out where in the parser the"
    print "problem occurred..."
    etype, value, tb = sys.exc_info()
    ftb = traceback.extract_tb(tb)
    ftb.reverse()
    class_found = None
    for err_file, err_line, err_function, err_text in ftb:
        if hasattr(consumer_class, err_function):
            class_found = consumer_class
            break
        elif hasattr(scanner_class, err_function):
            class_found = scanner_class
            break
    if class_found is None:
        print "Sorry, I could not pinpoint the error to the parser."
        print "There's nothing more I can tell you."
        print "Here's the traceback:"
        traceback.print_exception(etype, value, tb)
        return 1
    else:
        print "I found the problem in %s.%s.%s, line %d:" % \
              (class_found.__module__, class_found.__name__,
               err_function, err_line)
        print "    %s" % err_text
        print "This output caused an %s to be raised with the" % etype
        print "information %r." % exception_info
    print

    print "Let me find the line in the file that triggers the problem..."
    parser = parser_class()
    scanner, consumer = parser._scanner, parser._consumer
    consumer = DebuggingConsumer(consumer)
    try:
        scanner.feed(open(outfile), consumer)
    except etype, x:
        pass
    else:
        print "Odd, the exception disappeared!  What happened?"
        return 3
    print "It's caused by line %d:" % consumer.linenum
    lines = open(outfile).readlines()
    start, end = consumer.linenum-CONTEXT, consumer.linenum+CONTEXT+1
    if start < 0:
        start = 0
    if end > len(lines):
        end = len(lines)
    ndigits = len(str(end))
    for linenum in range(start, end):
        line = chomp(lines[linenum])
        if linenum == consumer.linenum:
            prefix = '*'
        else:
            prefix = ' '
        
        s = "%s%*d %s" % (prefix, ndigits, linenum, line)
        s = s[:80]
        print s
    print

    if class_found == scanner_class:
        print "Problems in %s are most likely caused by changed formats." % \
              class_found.__name__
        print "You can start to fix this by going to line %d in module %s." % \
              (err_line, class_found.__module__)
        print "Perhaps the scanner needs to be made more lenient by accepting"
        print "the changed format?"
        print

        if VERBOSITY <= 0:
            print "For more help, you can run this script in verbose mode"
            print "to see detailed information about how the scanner"
            print "identifies each line."
        else:
            print "OK, let's see what the scanner's doing!"
            print
            print "*"*20 + " BEGIN SCANNER TRACE " + "*"*20
            try:
                parser_class()._scanner.feed(
                    open(outfile), ParserSupport.TaggingConsumer())
            except etype, x:
                pass
            print "*"*20 + " END SCANNER TRACE " + "*"*20
        print
            
    elif class_found == consumer_class:
        print "Problems in %s can be caused by two things:" % \
              class_found.__name__
        print "    - The format of the line parsed by '%s' changed." % \
              err_function
        print "    - The scanner misidentified the line."
        print "Check to make sure '%s' should parse the line:" % \
              err_function
        s = "    %s" % chomp(lines[consumer.linenum])
        s = s[:80]
        print s
        print "If so, debug %s.%s.  Otherwise, debug %s." % \
              (class_found.__name__, err_function, scanner_class.__name__)
    

VERBOSITY = 0
if __name__ == '__main__':
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "hpnov")
    except getopt.error, x:
        print >>sys.stderr, x
        sys.exit(-1)
    if len(args) != 1:
        print >>sys.stderr, USAGE
        sys.exit(-1)
    TESTFILE, = args
    if not os.path.exists(TESTFILE):
        print >>sys.stderr, "I could not find file: %s" % TESTFILE
        sys.exit(-1)

    PROTEIN = NUCLEOTIDE = OUTPUT = None
    for opt, arg in optlist:
        if opt == '-h':
            print USAGE
            sys.exit(0)
        elif opt == '-p':
            PROTEIN = 1
        elif opt == '-n':
            NUCLEOTIDE = 1
        elif opt == '-o':
            OUTPUT = 1
        elif opt == '-v':
            VERBOSITY += 1

    if len([x for x in (PROTEIN, NUCLEOTIDE, OUTPUT) if x is not None]) != 1:
        OUTPUT = 1
        #print >>sys.stderr, "Exactly one of -p, -n, or -o should be specified."
        #sys.exit(-1)
    if PROTEIN or NUCLEOTIDE:
        print >>sys.stderr, "-p and -n not implemented yet"
        sys.exit(-1)
    test_blast_output(TESTFILE)
        
