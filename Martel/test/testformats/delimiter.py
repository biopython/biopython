"""Example of using Martel on a simple delimited file

"""
import Martel
from Martel import RecordReader

def delimiter(delim):
    assert len(delim) == 1, \
           "delimiter can only be a single character long, not %s" % repr(delim)
    assert delim not in "\n\r", "Cannot use %s as a delimiter" % repr(delim)
    
    field = Martel.Group("field", Martel.Rep(Martel.AnyBut(delim + "\r\n")))

    line = field + Martel.Rep(Martel.Str(delim) + field) + Martel.AnyEol()
    record = Martel.Group("record", line)

    format = Martel.ParseRecords("delimited", {}, record,
                                 RecordReader.CountLines, (1,))
    return format

tabformat = delimiter("\t")
spaceformat = delimiter(" ")
colonformat = delimiter(":")
commaformat = delimiter(",")

if __name__ == "__main__":
    from xml.sax import saxutils
    parser = colonformat.make_parser()
    parser.setContentHandler(saxutils.XMLGenerator())
    parser.parseFile(open("/etc/passwd"))
