import string, sys, cgi, time

from xml.sax import handler

class CheckGood(handler.ContentHandler, handler.ErrorHandler):
    def __init__(self):
        handler.ContentHandler.__init__(self)
        self.good_parse = 0
    def startDocument(self):
        self.good_parse = 0
    def endDocument(self):
        self.good_parse = 1


class Storage:
    def __init__(self):
        self.data = []
    def __getitem__(self, want_name):
        for name, exp, s in self.data:
            if name == want_name:
                x = Storage()
                x.add_test(name, exp, s)
                return x
        raise KeyError, want_name
    def add_test(self, name, exp, s):
        self.data.append( (name, exp, s) )

    def add_test_lines(self, name, exp, s):
        if s[-1:] == "\n":
            s = s[:-1]
        lines = string.split(s, "\n")
        for i in range(len(lines)):
            self.add_test(name + " %d" % (i+1), exp, lines[i] + "\n")

    def test(self, show_timings = 0, debug_level = 0):
        good = CheckGood()
        for name, exp, s in self.data:
            t1 = time.time()
            parser = exp.make_parser(debug_level = debug_level)
            t2 = time.time()
            parser.setContentHandler(good)
            try:
                t3 = t4 = time.time()
                parser.parseString(s)
                t4 = time.time()
            except KeyboardInterrupt:
                raise
            except:
                print "Cannot parse!!"
                print s
                print repr(str(exp))
                print "-" * 70
                raise
            else:
                if show_timings:
                    # Why are the times sometimes inverted???
                    t = "gen %.4f run %.4f gen takes %3.1f%%" % \
                        (t2-t1, t4-t3, 100*(t2-t1) / (t4-t3+t2-t1))
                    print "Good parse", name, t
                else:
                    print "Good parse", name

    def dump(self, debug_level = 0):
        h = Dump(sys.stdout)
        for name, exp, s in self.data:
            print "****************************", name, "********************"
            parser = exp.make_parser(debug_level = debug_level)
            parser.setContentHandler(h)
            parser.setErrorHandler(h)
            parser.parseString(s)


class Dump(handler.ContentHandler, handler.ErrorHandler):
    def __init__(self, outfile = None):
        handler.ContentHandler.__init__(self)
        if outfile is None:
            outfile = sys.stdout
        self.write = outfile.write
    def startElement(self, name, attrs):
        self.write('<%s>' % name)
    def characters(self, content):
        self.write(cgi.escape(content))
    def endElement(self, name):
        self.write('</%s>' % name)

    def startDocument(self):
        self.write("-------> Start\n")
    def endDocument(self):
        self.write("\n-------> End\n")

    def error(self, exc):
        self.write("\n-------> error " + str(exc) + "\n")
    def fatalError(self, exc):
        self.write("\n-------> fatal error - " + str(exc) + "\n")

def test_file(format, infile, outfile = sys.stdout):
    h = Dump(outfile)
    
    parser = format.make_parser()
    parser.setContentHandler(h)
    parser.setErrorHandler(h)
    parser.parseFile(infile)

def test_string(format, str):
    h = Dump(sys.stdout)
    
    parser = format.make_parser()
    parser.setContentHandler(h)
    parser.setErrorHandler(h)
    parser.parseString(str)

def time_file(format, infile, debug_level = 0):
    good = CheckGood()
    parser = format.make_parser(debug_level)
    parser.setContentHandler(good)
    t1 = time.time()
    parser.parseFile(infile)
    t2 = time.time()
    print "Total time", t2-t1
