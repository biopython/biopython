# regression tests for the Expression module

from Martel import Re
from Martel.Expression import *

from Martel.Expression import _minimize_any_range, _minimize_escape_range, \
     _minimize_escape_char, _verify_name


def compare(s1, s2):
    assert s1 == s2, "%s != %s" % (repr(s1), repr(s2))

def test__add__():
    compare(str(Literal("a") + Literal("b") + Literal("c")), "abc")

def test__or__():
    compare(str(Literal("a") | Literal("b") | Literal("c")), "a|b|c")

    # compose a couple of trees while I'm here
    compare(str(Literal("a") | (Literal("b") + Literal("c"))), "a|bc")
    compare(str(Literal("a") + (Literal("b") | Literal("c"))), "a(b|c)")

def test_any():
    compare(str(Any("abcdef")), "[a-f]")
    compare(str(Any("135")), "[135]")
    compare(str(Any("0123456789")), r"[\d]")

    compare(str(Any("abcdef", 1)), "[^a-f]")
    compare(str(Any("135", 1)), "[^135]")
    compare(str(Any("0123456789", 1)), r"[^\d]")

def test_assert():
    compare(str(Assert(Literal("a"))), "(?=a)")
    compare(str(Assert(Literal("a"), 1)), "(?!a)")

def test_at_beginning():
    compare(str(AtBeginning()), "^")

def test_at_end():
    compare(str(AtEnd()), "$")

def test_dot():
    compare(str(Dot()), ".")

def test_group():
    compare(str(Group(None, Literal("a"))), "(a)")
    compare(str(Group("name1", Literal("a"))), "(?P<name1>a)")

    compare(str(Group(None, Str("abcdef"))), "(abcdef)")
    compare(str(Group("name1", Str("abcdef"))), "(?P<name1>abcdef)")
    try:
        Group("illegal name", Literal("a"))
    except AssertionError:
        pass
    else:
        raise AssertionError, "wasn't catching illegal name"

def test_groupref():
    compare(str(GroupRef("name1")), "(?P=name1)")
    compare(str(GroupRef("N:a:m-e.1")), "(?P=N:a:m-e.1)")

def test_literal():
    compare(str(Literal("x")), "x")
    compare(str(Literal("y")), "y")

def test_max_repeat():
    compare(str(MaxRepeat(Literal("a"))), "a*")
    compare(str(MaxRepeat(Literal("a"), 1)), "a+")
    compare(str(MaxRepeat(Literal("a"), 0, 1)), "a?")
    compare(str(MaxRepeat(Literal("a"), 1, 1)), "a")
    compare(str(MaxRepeat(Literal("a"), 5, 5)), "a{5}")
    compare(str(MaxRepeat(Literal("a"), 0, 2)), "a{,2}")
    compare(str(MaxRepeat(Literal("a"), 8)), "a{8,}")
    compare(str(MaxRepeat(Literal("a"), 3, 9)), "a{3,9}")

    compare(str(MaxRepeat(Literal("a"), "z", "z")), "a{z}")

    # Not that these two are supported, but I expect they will be someday
    compare(str(MaxRepeat(Literal("a"), "x", 9)), "a{x,9}")
    compare(str(MaxRepeat(Literal("a"), "x", "y")), "a{x,y}")

    compare(str(MaxRepeat(Str("test"), 5, 6)), "(test){5,6}")
    compare(str(MaxRepeat( Literal("a") | Str("test"))), "(a|test)*")
    compare(str(MaxRepeat( Literal("a") + Str(" test"))), "(a test)*")
            
def test_null_op():
    compare(str(NullOp()), "")
    compare(str(NullOp() + Str("TeSt")), str(Str("TeSt")))

def test_str():
    compare(str(Str("Andrew Dalke")), "Andrew Dalke")

def test_alt():
    compare(str(Alt( (Any("abc", 1), Dot(), Literal("Q"), Str("hello")) )),
            "[^a-c]|.|Q|hello")

    # Check that an Alt of Alts is handled correctly
    L = Literal
    alt1 = Alt( (L("a"), L("b"), L("c")) )
    alt2 = Alt( (L("x"), L("y"), L("z")) )

    alt = Alt( (alt1, alt2, alt1, alt2) )

    compare(str(alt), "a|b|c|x|y|z|a|b|c|x|y|z")
    

def test_seq():
    compare(str(Seq( (Any("abc", 1), Dot(), Literal("Q"), Str("hello")) )),
            "[^a-c].Qhello")

    # check that the Alt is correctly grouped
    L = Literal
    alt = Alt( (L("a"), L("b"), L("c")) )
    
    compare(str(Seq( (alt, alt) )), "(a|b|c)(a|b|c)")
    compare(str(Seq( (alt, L("x"), alt) )), "(a|b|c)x(a|b|c)")

    seq = Seq( (L("a"), L("b"), L("c")) )
    compare(str(Seq( (seq, seq) )), "abcabc")

    compare(str(Seq( (seq, alt, seq, alt) )), "abc(a|b|c)abc(a|b|c)")
        

def test_nocase():
    compare(str(NoCase(Literal("A"))), "[Aa]")
    compare(str(NoCase(Literal("9.8"))), "9\\.8")
    compare(str(NoCase(Str("A"))), "[Aa]")
    compare(str(NoCase(Str("AB"))), "[Aa][Bb]")
    compare(str(NoCase(Re("[ab]"))), "[ABab]")
    compare(str(NoCase(Re("[abC]"))), "[A-Ca-c]")
    compare(str(NoCase(Any("1AbC9"))), "[19A-Ca-c]")
    compare(str(NoCase(Str("age = 10"))), "[Aa][Gg][Ee] = 10")
    compare(str(NoCase(Alt( (Str("A")|Str("B")|Str("C"),)))), "[Aa]|[Bb]|[Cc]")

def test_nodes():
    test__add__()
    test__or__()
    test_any()
    test_assert()
    test_at_beginning()
    test_at_end()
    test_dot()
    test_group()
    test_groupref()
    test_literal()
    test_max_repeat()
    test_null_op()
    test_str()
    test_alt()
    test_seq()
    test_nocase()

    

def test_minimize():
    data = (
        ("a", "a"),
        ("ab", "ab"),
        ("abc", "a-c"),
        ("abcdef", "a-f"),
        ("abcj", "a-cj"),
        ("abcjk", "a-cjk"),
        ("abcjkl", "a-cj-l"),
        ("0123456789", r"\d"),
        (string.lowercase, r"a-z"),
        (string.lowercase + string.uppercase, r"A-Za-z"),
        (string.lowercase + string.uppercase + string.digits, r"\dA-Za-z"),
        ("\007\010", r"\a\b"),
        ("1357", "1357"),
        ("", ""),
        ("$", "$"),
        ("\\", "\\\\"),
        ("-0123456789", r"\d-"),
        ("-", "-"),
        ("^", r"\^"),
        ("^-", r"\^-"),
        ("[]", r"\[\]"),
        ("\001", "\\\001"),
        ("\001\002\003\004\005\006\007", "\\\001-\\a"),
        )
    for input, output in data:
        x = _minimize_any_range(input)
        assert x == output, "%s : %s != %s" % (repr(input), repr(output), repr(x))

def test_group_names():
    a = Literal('a')
    b = Literal('b')
    g1 = Group("foo", a)
    g2 = Group("bar", b)
    g = Group("spam", Alt([g1, g2]))
    assert str(g) == "(?P<spam>(?P<foo>a)|(?P<bar>b))", str(g)
    assert g.group_names() == ("spam", "foo", "bar"), g.group_names()

    h = g.copy()
    assert str(h) == "(?P<spam>(?P<foo>a)|(?P<bar>b))", str(h)
    
    g2.name = "baz"
    assert g.group_names() == ("spam", "foo", "baz"), g.group_names()

    assert str(h) == "(?P<spam>(?P<foo>a)|(?P<bar>b))", str(h)

    g._select_names(["foo", "baz"])
    assert g.group_names() == ("foo", "baz"), g.group_names()
        

def test_valid_names():
    good_names = (
        "x",
        "name",
        "spec:header",
        "xmlns:xlink",
        ":",
        ":something",
        ":012",
        "_",
        "_9",
        "_-_",
        "A-",
        "this-field",
        ":::",
        "MixedCase",
        "UPPER",
        "this._is_.valid",
        "x0")
    for name in good_names:
        _verify_name(name)
        
    bad_names = (
        "0",
        "\b",
        "-",
        "A B",
        "AB ",
        " AB",
        ".invalid",
        )
    for name in bad_names:
        try:
            _verify_name(name)
        except AssertionError:
            pass
        else:
            raise AssertionError, "Should not have allowed %s" % repr(name)

def test():
    test_valid_names()
    test_nodes()
    test_minimize()
    test_group_names()

if __name__ == "__main__":
    test()
