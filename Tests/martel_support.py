import string, sys, cgi

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

    def test(self):
        good = CheckGood()
        for name, exp, s in self.data:
            parser = exp.make_parser()
            parser.setContentHandler(good)
            try:
                parser.parseString(s)
            except KeyboardInterrupt:
                raise
            except:
                pass
            if not good.good_parse:
                print "Cannot parse", name
                print s
                print repr(str(exp))
                print "-" * 70
                raise AssertionError, "cannot parse"
            else:
                print "Good parse", name

    def dump(self):
        h = Dump(sys.stdout)
        for name, exp, s in self.data:
            print "****************************", name, "********************"
            parser = exp.make_parser()
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

def test_file(format, infile):
    h = Dump(sys.stdout)
    
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
