# Ensure that attributes really work

from xml.sax import handler
import Martel
from Martel import RecordReader

class GrabElements(handler.ContentHandler, handler.ErrorHandler):
    def startDocument(self):
        self.elements = []
    def startElement(self, name, attrs):
        self.elements.append( (name, attrs) )

##def fixup_dict(d):
##    return d
##    dict = {}
##    for k, v in d.items():
##        dict[k] = v
##    return dict

def get_element(format, text = "X"):
    parser = format.make_parser()
    grab = GrabElements()
    parser.setContentHandler(grab)
    parser.parseString(text)
    assert len(grab.elements) == 1, len(grab.elements)
    return grab.elements[0]

def cmp_dicts(d1, d2):
    items1 = d1.items()
    items1.sort()

    items2 = d2.items()
    items2.sort()

    return cmp(items1, items2)

def check_dicts(d1, d2):
    if cmp_dicts(d1, d2) != 0:
        raise AssertionError("No match: %s, %s" % (d1, d2))

def check_element( (name1, dict1), (name2, dict2) ):
    assert name1 == name2, (name1, name2)
    check_dicts(dict1, dict2)


def test_none():
    ele = get_element(Martel.Group("spam", Martel.Str("X")))
    check_element(ele, ("spam", {}))

    ele = get_element(Martel.Group("spam", Martel.Str("X"), {}))
    check_element(ele, ("spam", {}))

def test_single():
    ele = get_element(
        Martel.Group("spam", Martel.Str("X"), {"format": "swissprot"}))
    check_element(ele, ("spam", {"format": "swissprot"}))

def test_multi():
    ele = get_element(Martel.Re("(?P<qwe?def=%7Edalke&a=1&b=2&cc=33>X)"))
    check_element(ele, ("qwe", {"def": "~dalke",
                                "a": "1",
                                "b": "2",
                                "cc": "33"}))
def test_escape():
    name, attrs = get_element(Martel.Re("(?P<qwe?a=%7E>X)"))
    check_dicts({"a": "~"}, attrs)

    name, attrs = get_element(Martel.Re("(?P<qwe?a=%7e>X)"))
    check_dicts({"a": "~"}, attrs)

    name, attrs = get_element(Martel.Re("(?P<qwe?a=>X)"))
    check_dicts({"a": ""}, attrs)

    name, attrs = get_element(Martel.Re("(?P<qwe?>X)"))
    check_dicts({}, attrs)

    name, attrs = get_element(Martel.Re("(?P<qwe?a=%48%65%6c%6c%6f>X)"))
    check_dicts({"a": "Hello"}, attrs)

    name, attrs = get_element(Martel.Re("(?P<qwe?a=%7e%7E&b=%7e%7E>X)"))
    check_dicts({"a": "~~", "b": "~~"}, attrs)

    name, attrs = get_element(Martel.Re("(?P<qwe?a=%7e%7E&b=>X)"))
    check_dicts({"a": "~~", "b": ""}, attrs)

def test_same_tag():
    format = Martel.Re("(?P<char?type=a>a+)(?P<char?type=b>b+)")
    parser = format.make_parser()
    grab = GrabElements()
    parser.setContentHandler(grab)
    parser.parseString("aaabb")
    assert len(grab.elements) == 2, len(grab.elements)

    check_element(grab.elements[0], ("char", {"type": "a"}))
    check_element(grab.elements[1], ("char", {"type": "b"}))
    

def test_record_parser():
    format = Martel.Re("(?P<term?field=first>...)"
                       "(?P<term?field=second>...)"
                       "(?P<last>.)\R")
    format = Martel.ParseRecords("all", {"author": "guido"}, format,
                                 RecordReader.CountLines, (1,) )
    parser = format.make_parser()
    grab = GrabElements()
    parser.setContentHandler(grab)
    parser.parseString("aaabbbZ\ncccdddZ\n")
    elements = grab.elements
    assert len(elements) == 7
    check_element(elements[0], ("all", {"author": "guido"}))
    check_element(elements[1], ("term", {"field": "first"}))
    check_element(elements[2], ("term", {"field": "second"}))
    check_element(elements[3], ("last", {}))
    check_element(elements[4], ("term", {"field": "first"}))
    check_element(elements[5], ("term", {"field": "second"}))
    check_element(elements[6], ("last", {}))
    

def test_header_footer_parser():
    # Check that I can pass same tag names in the header, record and
    # footer but not have them collide.

    header_format = Martel.Re("(?P<term?pos=header>a+)\R")
    record_format = Martel.Re("(?P<term?pos=body>b+)\R")
    footer_format = Martel.Re("(?P<term?pos=footer>c+)\R")

    format = Martel.HeaderFooter(
        "all", {"state": "New Mexico"},
        header_format, RecordReader.CountLines, (1,),
        record_format, RecordReader.CountLines, (1,),
        footer_format, RecordReader.CountLines, (1,),
        )

    parser = format.make_parser()
    grab = GrabElements()
    parser.setContentHandler(grab)
    parser.parseString("a\nbb\nbb\nccc\n")
    elements = grab.elements
    assert len(elements) == 5, len(elements)
    check_element(elements[0], ("all", {"state": "New Mexico"}))
    check_element(elements[1], ("term", {"pos": "header"}))
    check_element(elements[2], ("term", {"pos": "body"}))
    check_element(elements[3], ("term", {"pos": "body"}))
    check_element(elements[4], ("term", {"pos": "footer"}))
    

def test():
    test_none()
    test_single()
    test_multi()
    test_escape()
    test_same_tag()
    test_record_parser()
    test_header_footer_parser()
    
if __name__ == "__main__":
    test()
