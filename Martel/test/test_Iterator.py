# tests of the iterator interface

import os
from testformats import swissprot38
import Martel
from Martel import RecordReader, Parser
from xml.sax import handler

import test_swissprot38

text = test_swissprot38.record1 + test_swissprot38.record2


def test1():
    # Does it read all of the records?
    iterator = swissprot38.format.make_iterator("swissprot38_record")
    stream = iterator.iterateString(text, handler.ContentHandler())
    i = 0
    while 1:
        try:
            x = stream.next()
        except StopIteration:  # for IterParser objects
            x = None
        if x is None:
            break
        i = i + 1
    assert i == 2, i
            
def test2():
    # Is is reading a record at a time?
    iterator = swissprot38.format.make_iterator("swissprot38_record")
    stream = iterator.iterateString(test_swissprot38.record1 + "X" * 100,
                                    handler.ContentHandler())
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
    stream = iterator.iterateString(text, handler.ContentHandler())
    i = 0
    while 1:
        x = stream.next()
        if x is None:
            break
        i = i + 1
    assert i == 2, i

def test4():
    # Make sure the default returns LAX items
    exp = Martel.Re("(?P<term>(?P<a>a+)(?P<b>b+))+")
    x = exp.make_iterator("term").iterateString("aabbabaaaabb")
    term = x.next()
    assert len(term["a"]) == 1 and term["a"][0] == "aa", term["a"]
    assert len(term["b"]) == 1 and term["b"][0] == "bb", term["b"]
    term = x.next()
    assert len(term["a"]) == 1 and term["a"][0] == "a", term["a"]
    assert len(term["b"]) == 1 and term["b"][0] == "b", term["b"]
    term = x.next()
    assert len(term["a"]) == 1 and term["a"][0] == "aaaa", term["a"]
    assert len(term["b"]) == 1 and term["b"][0] == "bb", term["b"]
    term = x.next()
    assert term is None, "Did not stop correctly"

def test5():
    # Does 'iter' work?
    try:
        iter
    except NameError:
        print "Test skipped - missing 'iter' builtin from Python 2.2."
        return
    exp = Martel.Re("(?P<term>(?P<a>a+)(?P<b>b+))+")
    x = exp.make_iterator("term")
    it = iter(x.iterateString("aabbabaaaabb"))
    term = it.next()
    assert len(term["a"]) == 1 and term["a"][0] == "aa", term["a"]
    assert len(term["b"]) == 1 and term["b"][0] == "bb", term["b"]
    term = it.next()
    assert len(term["a"]) == 1 and term["a"][0] == "a", term["a"]
    assert len(term["b"]) == 1 and term["b"][0] == "b", term["b"]
    term = it.next()
    assert len(term["a"]) == 1 and term["a"][0] == "aaaa", term["a"]
    assert len(term["b"]) == 1 and term["b"][0] == "bb", term["b"]
    try:
        it.next()
        raise AssertionError("Did not stop correctly")
    except StopIteration:
        pass

def test_header_footer1():
    exp = Martel.HeaderFooter("dataset", {},
                              Martel.Re("header\R"),
                              RecordReader.CountLines, (1,),
                              
                              Martel.Re("a(?P<b>b*)a\R"),
                              RecordReader.CountLines, (1,),
                              
                              Martel.Re("end\R"),
                              RecordReader.CountLines, (1,),)
    lines = [
        "header",
        "aa",
        "aba",
        "abba",
        "end"
        ]
    text = "\n".join(lines) + "\n"

    i = 0
    for info in exp.make_iterator("b").iterateString(text):
        assert len(info["b"]) == 1
        assert len(info["b"][0]) == i, (info["b"][0], i)
        i = i +1
    assert i == 3, i
                        

def test_header_footer2():
    exp = Martel.HeaderFooter("dataset", {},
                              None,
                              None, None,
                              
                              Martel.Re("a(?P<b>b*)a\R"),
                              RecordReader.CountLines, (1,),
                              
                              Martel.Re("end\R"),
                              RecordReader.CountLines, (1,),)
    lines = [
        "aa",
        "aba",
        "abba",
        "end"
        ]
    text = "\n".join(lines) + "\n"

    i = 0
    for info in exp.make_iterator("b").iterateString(text):
        assert len(info["b"]) == 1
        assert len(info["b"][0]) == i, (info["b"][0], i)
        i = i +1
    assert i == 3, i
                        

def test_header_footer3():
    exp = Martel.HeaderFooter("dataset", {},
                              None,
                              None, None,
                              
                              Martel.Re("a(?P<b>b*)a\R"),
                              RecordReader.CountLines, (1,),
                              
                              None,
                              None, None)
    lines = [
        "aa",
        "aba",
        "abba",
        ]
    text = "\n".join(lines) + "\n"

    i = 0
    for info in exp.make_iterator("b").iterateString(text):
        assert len(info["b"]) == 1
        assert len(info["b"][0]) == i, (info["b"][0], i)
        i = i +1
    assert i == 3, i
                        

def test_header_footer4():
    # Test that the errors are correct
    exp = Martel.HeaderFooter("dataset", {},
                              Martel.Re("header\R"),
                              RecordReader.CountLines, (1,),
                              
                              Martel.Re("a(?P<b>b*)a\R"),
                              RecordReader.CountLines, (1,),
                              
                              Martel.Re("end\R"),
                              RecordReader.CountLines, (1,),)
    lines = [
        "HEADER",
        "aa",
        "aba",
        "abba",
        "end"
        ]
    text = "\n".join(lines) + "\n"

    try:
        for info in exp.make_iterator("b").iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 0
                        

def test_header_footer5():
    exp = Martel.HeaderFooter("dataset", {},
                              None,
                              None, None,
                              
                              Martel.Re("a(?P<b>b*)a\R"),
                              RecordReader.CountLines, (1,),
                              
                              Martel.Re("end\R"),
                              RecordReader.CountLines, (1,),)
    lines = [
        "aa",
        "aba",
        "abba",
        "END"
        ]
    text = "\n".join(lines) + "\n"

    try:
        for info in exp.make_iterator("b").iterateString(text):
            assert len(info["b"]) == 1
    except Parser.ParserPositionException, exc:
        assert exc.pos == 12, exc.pos
                        

def test_header_footer6():
    exp = Martel.HeaderFooter("dataset", {},
                              None,
                              None, None,
                              
                              Martel.Re("a(?P<b>b*)a\R"),
                              RecordReader.CountLines, (1,),
                              
                              None,
                              None, None)
    lines = [
        "aA",
        "aBbbba",
        "abba",
        ]
    text = "\n".join(lines) + "\n"

    try:
        for info in exp.make_iterator("b").iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 1, exc.pos

def test_header_footer7():
    exp = Martel.HeaderFooter("dataset", {},
                              None,
                              None, None,
                              
                              Martel.Re("a(?P<b>b*)a\R"),
                              RecordReader.CountLines, (1,),
                              
                              None,
                              None, None)
    lines = [
        "aa",
        "aBbbba",
        "abba",
        ]
    text = "\n".join(lines) + "\n"

    try:
        for info in exp.make_iterator("b").iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 4, exc.pos
                        

    
def test():
    # print "Iteration through all records"
    test1()
    # print "Iteration a record at a time"
    test2()
    # print "Test non-record-based parsers"
    test3()
    # print "Test return of LAX objects"
    test4()
    # print "Test Python 2.2 iterators"
    test5()

    test_header_footer1()
    test_header_footer2()
    test_header_footer3()
    test_header_footer4()
    test_header_footer5()
    test_header_footer6()
    test_header_footer7()

if __name__ == "__main__":
    test()
