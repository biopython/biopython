# tests of the iterator interface

import Martel
from Martel.formats import swissprot38
from xml.sax import handler

import test_swissprot38

text = test_swissprot38.record1 + test_swissprot38.record2


def test1():
    # Does it read all of the records?
    iterator = swissprot38.format.make_iterator("swissprot38_record")
    stream = iterator.iterateString(text, handler.ContentHandler)
    i = 0
    while 1:
        x = stream.next()
        if x is None:
            break
        i = i + 1
    assert i == 2, i
            
def test2():
    # Is is reading a record at a time?
    iterator = swissprot38.format.make_iterator("swissprot38_record")
    stream = iterator.iterateString(test_swissprot38.record1 + "X" * 100,
                                    handler.ContentHandler)
    x = stream.next()
    assert x is not None, x
    try:
        x = stream.next()
    except KeyboardInterrupt:
        raise
    except:
        pass
    else:
        raise AssertionError, "should not allow X's"

def test3():
    # test the non-record reader parser
    iterator = swissprot38.format.expression.make_iterator("swissprot38_record")
    stream = iterator.iterateString(text, handler.ContentHandler)
    i = 0
    while 1:
        x = stream.next()
        if x is None:
            break
        i = i + 1
    assert i == 2, i
    
    
def test():
    print "Iteration through all records"
    test1()
    print "Iteration a record at a time"
    test2()
    print "Test non-record-based parsers"
    test3()

if __name__ == "__main__":
    test()
