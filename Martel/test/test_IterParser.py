import Martel
from Martel import IterParser, RecordReader, LAX, Parser

# Used a lot
def gen_iterator():
    return IterParser.IterHeaderFooter(
        Martel.Re(r"a*\R").make_parser(),
        RecordReader.CountLines,
        (1,),

        Martel.Group("spam", Martel.Re(r"b*\Rc*\R")).make_parser(debug_level = 1),
        RecordReader.CountLines,
        (2,),

        Martel.Re(r"d*\R").make_parser(),
        RecordReader.CountLines,
        (1,),

        "spam")

def test_hf1():
    ip = gen_iterator()
    
    lines = ["aaaaaaaaa",
             "b",
             "c",
             "bb",
             "cc", 
             "bbb",
             "ccc",
             "d"]
    text = "\n".join(lines) + "\n"
             

    i = 1
    for x in ip.iterateString(text):
        assert x["spam"][0] == "b" * i + "\n" + "c" * i + "\n"
        i = i + 1

def test_hf2():
    ip = IterParser.IterHeaderFooter(
        None,
        None,
        None,

        Martel.Group("spam", Martel.Re(r"b*\Rc*\R")).make_parser(),
        RecordReader.CountLines,
        (2,),

        Martel.Re(r"d*\R").make_parser(),
        RecordReader.CountLines,
        (1,),

        "spam")

    lines = ["b",
             "c",
             "bb",
             "cc", 
             "bbb",
             "ccc",
             "d"]
    text = "\n".join(lines) + "\n"

    i = 1
    for x in ip.iterateString(text):
        assert x["spam"][0] == "b" * i + "\n" + "c" * i + "\n"
        i = i + 1

def test_hf3():
    ip = IterParser.IterHeaderFooter(
        None,
        None,
        None,

        Martel.Group("spam", Martel.Re(r"b*\Rc*\R")).make_parser(),
        RecordReader.CountLines,
        (2,),

        None,
        None,
        None,

        "spam")

    lines = ["b",
             "c",
             "bb",
             "cc", 
             "bbb",
             "ccc",
             ]
    text = "\n".join(lines) + "\n"

    i = 1
    for x in ip.iterateString(text):
        assert x["spam"][0] == "b" * i + "\n" + "c" * i + "\n"
        i = i + 1

def test_hf4():
    ip = IterParser.IterHeaderFooter(
        Martel.Re(r"a*\R").make_parser(),
        RecordReader.CountLines,
        (1,),

        Martel.Group("spam", Martel.Re(r"b*\Rc*\R")).make_parser(),
        RecordReader.CountLines,
        (2,),

        None,
        None,
        None,

        "spam")

    lines = ["aaaaaaaaa",
             "b",
             "c",
             "bb",
             "cc", 
             "bbb",
             "ccc",
             ]
    text = "\n".join(lines) + "\n"

    i = 1
    for x in ip.iterateString(text):
        assert x["spam"][0] == "b" * i + "\n" + "c" * i + "\n"
        i = i + 1


def test_hf5():
    # This has a syntax error in the header
    ip = gen_iterator()

    lines = ["aaaaAaaaa",
             "b",
             "c",
             "bb",
             "cc", 
             "bbb",
             "ccc",
             "d"]
    text = "\n".join(lines) + "\n"
             

    try:
        for x in ip.iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 4

def test_hf6():
    # This has a syntax error in the first line of the first record
    ip = gen_iterator()

    lines = ["aaaaaaaaa",
             "-",
             "c",
             "bb",
             "cc", 
             "bbb",
             "ccc",
             "d"]
    text = "\n".join(lines) + "\n"
             

    try:
        for x in ip.iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 10

def test_hf7():
    # This has a syntax error in the first line of the second record
    ip = gen_iterator()

    lines = ["aaaaaaaaa",
             "b",
             "c",
             "b-",
             "cc", 
             "bbb",
             "ccc",
             "d"]
    text = "\n".join(lines) + "\n"
             

    try:
        for x in ip.iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 14, exc.pos

def test_hf8():
    # This has a syntax error in the footer
    ip = gen_iterator()

    lines = ["aaaaaaaaa",
             "b",
             "c",
             "bb",
             "cc", 
             "bbb",
             "ccc",
             "d-"]
    text = "\n".join(lines) + "\n"
             

    try:
        for x in ip.iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 29, exc.pos


### the simple record iterator
def test_ri1():
    ip = IterParser.IterRecords(
        Martel.Group("spam", Martel.Re(r"b*\Rc*\R")).make_parser(),
        RecordReader.CountLines,
        (2,),
        
        "spam")

    lines = ["b",
             "c",
             "bb",
             "cc",
             "bbb",
             "ccc",
             ]
    text = "\n".join(lines) + "\n"
    i = 1
    for x in ip.iterateString(text):
        assert x["spam"][0] == "b" * i + "\n" + "c" * i + "\n"
        i = i + 1


def test_ri2():
    # error in the first record
    ip = IterParser.IterRecords(
        Martel.Group("spam", Martel.Re(r"b*\Rc*\R")).make_parser(debug_level = 1),
        RecordReader.CountLines,
        (2,),
        
        "spam")

    lines = ["b",
             "-",
             "bb",
             "cc",
             "bbb",
             "ccc",
             ]
    text = "\n".join(lines) + "\n"
    try:
        for x in ip.iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 2, exc.pos

def test_ri3():
    # error in the second record
    ip = IterParser.IterRecords(
        Martel.Group("spam", Martel.Re(r"b*\Rc*\R")).make_parser(debug_level = 1),
        RecordReader.CountLines,
        (2,),
        
        "spam")

    lines = ["b",
             "c",
             "b-",
             "cc",
             "bbb",
             "ccc",
             ]
    text = "\n".join(lines) + "\n"
    try:
        for x in ip.iterateString(text):
            pass
    except Parser.ParserPositionException, exc:
        assert exc.pos == 5, exc.pos

def test_make_iterparsers1():
    exp = Martel.ParseRecords("dataset", {},
                              Martel.Group("spam", Martel.Re(r"a*\R")),
                              RecordReader.CountLines, (1,))
    iterator = exp.make_iterator("spam")
    assert isinstance(iterator, IterParser.IterRecords)
    lines = []
    for i in range(0, 10):
        lines.append("a" * i + "\n")
    text = "".join(lines)

    i = 0
    for rec in iterator.iterateString(text):
        assert len(rec["spam"][0][:-1]) == i, (i, rec["spam"][0])
        i = i + 1
    assert i == 10


def test_make_iterparsers2():
    exp = Martel.HeaderFooter("dataset", {},
                              Martel.Group("header", Martel.Re(r"(a*\R)*")),
                              RecordReader.Until, ("b",),
                              
                              Martel.Group("record", Martel.Re(r"(b*\R)*")),
                              RecordReader.Until, ("c",),
                              
                              Martel.Group("footer", Martel.Re(r"(c*\R)*")),
                              RecordReader.Everything, (),)
    
    iterator = exp.make_iterator("record")
    assert isinstance(iterator, IterParser.IterHeaderFooter), iterator
    lines = ["a"
             "aa",
             "aaaaaaa",
             "b",
             "bb",
             "bbbb",
             "bbbbbbbb",
             "bbbbbbbbbbbbbbbb",
             "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
             "cccc",
             "cc",
             "c",
             ]
    
    text = "\n".join(lines) + "\n"

    i = 0
    for rec in iterator.iterateString(text):
        i = i + 1
    assert i == 1, i

def test_missing_end1():
    exp = Martel.ParseRecords("dataset", {},
                              Martel.Group("record", Martel.Re(r"(b*\R)*a\R")),
                              RecordReader.EndsWith, ("a",))
    lines = [
        "bbb",
        "bb",
        "a",
        "b",
        "a",
        "a",
        ]
    text = "\n".join(lines) + "\n"

    iterator = exp.make_iterator("record")
    # This should work
    for x in iterator.iterateString(text):
        pass

    # This should not
    lines.append("c")
    text = "\n".join(lines) + "\n"
    try:
        for x in iterator.iterateString(text):
            pass
        raise AssertionError
    except Parser.ParserPositionException, exc:
        assert exc.pos == 15, exc.pos

def test_missing_end2():
    # Same as the test_missing_end1 but using HeaderFooter
    exp = Martel.HeaderFooter("dataset", {},
                              None, None, None,
                              Martel.Group("record", Martel.Re(r"(b*\R)*a\R")),
                              RecordReader.EndsWith, ("a",),
                              None, None, None
                              )
    lines = [
        "bbb",
        "bb",
        "a",
        "b",
        "a",
        "a",
        ]
    text = "\n".join(lines) + "\n"

    iterator = exp.make_iterator("record")
    # This should work
    for x in iterator.iterateString(text):
        pass

    # This should not
    lines.append("c")
    text = "\n".join(lines) + "\n"
    try:
        for x in iterator.iterateString(text):
            pass
        raise AssertionError
    except Parser.ParserPositionException, exc:
        assert exc.pos == 15, exc.pos

def test_missing_end3():
    # This one is missing the footer
    exp = Martel.HeaderFooter("dataset", {},
                              None, None, None,
                              
                              Martel.Group("record", Martel.Re(r"(b*\R)*a\R")),
                              RecordReader.EndsWith, ("a",),

                              Martel.Group("footer", Martel.Re(r"c\R")),
                              RecordReader.CountLines, (1,)
                              )
    lines = [
        "bbb",
        "bb",
        "a",
        "b",
        "a",
        "a",
        "c",  # This will be removed for the test
        ]
    text = "\n".join(lines) + "\n"

    iterator = exp.make_iterator("record")
    # This should work
    for x in iterator.iterateString(text):
        pass

    # This should not
    lines.pop()
    text = "\n".join(lines) + "\n"
    try:
        for x in iterator.iterateString(text):
            pass
        raise AssertionError
    except Parser.ParserPositionException, exc:
        assert exc.pos == 15, exc.pos


def test():
    test_hf1()
    test_hf2()
    test_hf3()
    test_hf4()
    test_hf5()
    test_hf6()
    test_hf7()
    test_hf8()

    test_ri1()
    test_ri2()
    test_ri3()

    test_make_iterparsers1()
    test_make_iterparsers2()
    
    test_missing_end1()
    test_missing_end2()
    test_missing_end3()
    
if __name__ == "__main__":
    test()
