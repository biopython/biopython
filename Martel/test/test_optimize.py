from xml.sax import handler

from Martel import *
from Martel import Parser, optimize, Expression

def test1():
    format = Group(None, Group("test", Str("A")))
    assert format.name is None
    format2 = optimize.optimize(format)
    assert format.name is None    # must not affect original expressions
    assert format2.name == "test" # got rid of top level
    assert format2.expression.string == "A"

def test2():
    format = Group("test", Group(None, Str("A")))
    format = optimize.optimize(format)
    assert format.name == "test"
    assert format.expression.string == "A"

def test3():
    format = optimize.optimize(Rep1(Group(None, Str("A"))))
    assert format.expression.string == "A"

def test4():
    format = optimize.optimize(Group(None, Group(None, Group(None, Str("A")))))
    assert format.string == "A"

def test5():
    format = optimize.optimize(Group(None, Str("A")) + \
                               Group("B", Str("B")) + \
                               Group(None, Str("C")))
    assert format.expressions[0].string == "A"
    assert format.expressions[1].expression.string == "B"
    assert format.expressions[2].string == "C"

def test6():
    format = optimize.optimize(Str("A") + Str("B"))
    assert format.string == "AB"

def test7():
    format = optimize.optimize(Str("A") + Str("B") + Str("C") + Str("D"))
    assert format.string == "ABCD"
    
def test8():
    format = optimize.optimize((Str("A") + Str("B")) + Str("C"))
    assert format.string == "ABC"

def test9():
    # Shows that subexpressions are combined first
    format = optimize.optimize(Str("A") + (Str("B") + Str("C")))
    assert format.string == "ABC"

def test10():
    format = optimize.optimize(Any("A") + Str("B"))
    assert format.string == "AB"

def test11():
    format = optimize.optimize(AnyBut("A") + Str("B"))
    assert format.expressions[1].string == "B"

def test12():
    format = optimize.optimize(Str("A") + Expression.Any("B"))
    assert format.string == "AB"

def test13():
    format = optimize.optimize(Expression.Any("A") + \
                               Expression.Any("B") + \
                               Expression.Any("C") + \
                               ((Expression.Any("D") + \
                                 Expression.Any("E")) + \
                                Expression.Any("F")))
    assert format.string == "ABCDEF"

def test14():
    format = optimize.optimize(Str("A") + Expression.Any("B", invert = 1))
    assert isinstance(format, Expression.Seq)

def test15():
    format = optimize.optimize(Str("A") + Any("BC"))
    assert isinstance(format, Expression.Seq)

def test16():
    format = optimize.optimize(Str("A") + Expression.Literal("B"))
    assert format.string == "AB"

def test17():
    format = optimize.optimize(Str("A") + Expression.Literal("B", invert = 1))
    assert isinstance(format, Expression.Seq)

def test18():
    format = optimize.optimize(
        Group(None, Group(None, Str("A") + Str("B")) + Str("C")) + \
        Group(None, Str("D") + Expression.Any("E") + \
              Group(None, Str("FG") + Str("H")) + Group(None, Str("I"))))
    assert format.string == "ABCDEFGHI"

def test19():
    exp = Re("This is a test")
    assert isinstance(exp, Expression.Seq)
    assert exp.expressions[0].char == "T"
    assert exp.expressions[-1].char == "t"
    format = optimize.optimize(exp)
    assert format.string == "This is a test"

class GetErrorPos(handler.ContentHandler, handler.ErrorHandler):
    def __init__(self):
        handler.ContentHandler.__init__(self)
        self.pos = None
    def startDocument(self):
        self.pos = None
    def error(self, exc):
        if isinstance(exc, Parser.ParserPositionException):
            self.pos = exc.pos
    fatalError = error


def check_error_pos(text, exp):
    parser = exp.make_parser()
    err = GetErrorPos()
    parser.setErrorHandler(err)
    parser.parseString(text)
    return err.pos
                              

def test20():
    test_format = Group("full_format",
                        Rep(Group("a_record",
                                  Str("A") + Opt(Str("B")) + Str("\n"))) + \
                        Str("Z"))

    text = "A\nAC\nAB\nZ"
    pos = check_error_pos(text, test_format)
    assert pos == 0
    
    format = select_names(test_format, ())
    format = optimize.optimize(format)
    pos = check_error_pos(text, format)
    assert pos == 2

def test():
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    test9()
    test10()
    test11()
    test12()
    test13()
    test14()
    test15()
    test16()
    test17()
    test18()
    test19()
    test20()

if __name__ == "__main__":
    test()
