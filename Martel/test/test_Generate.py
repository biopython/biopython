# Test that the parser generation works.

from Martel.Generate import *
from Martel.Generate import _generate

from Martel import convert_re, Parser
import re

from xml.sax import handler


def test():
    class _Test(handler.ContentHandler, handler.ErrorHandler):
        def __init__(self):
            handler.ContentHandler.__init__(self)
            self.good_parse = 0
        def startDocument(self):
            self.good_parse = 1
        def fatalError(self, exc):
            if isinstance(exc, Parser.ParserPositionException):
                # Called when there aren't enough characters
                self.good_parse = 0

        def error(self, exc):
            # shouldn't be called with these parsers
            raise exc

    cb = _Test()
    
    patterns = (
        ("a", ("a",), ("A", "", "Z")),
        ("[a-z]", ("a", "b", "q"), ("A", "-")),
        ("[^abc]", ("A", "d", "f"), ("a", "b", "c", "ab", "")),
        ("a+", ("a", "aaa"), ("A", "baa")),
        ("a*", ("", "a", "aaa"), ()),
        ("\\]", ("]",), ("a",)),
        ("a*$", ("a", "aaa"), ("A", "baa", "aaaaab")),
        ("(ab|ac)", ("ab", "ac"), ("aa", "A", "a", "cb")),
        ("(ab|ac)*$", ("", "ab", "ac", "abacabac", "ababab"),
                      ("aa", "A", "a", "cb", "ababababaca")),
        ("ab{3}$", ("abbb",), ("abb", "bbb", "abbbb")),
        ("ab{3,}$", ("abbb", "abbbb", "abbbbbbbbb"), ("abb", "bbb", "abbbc")),
        ("ab{3,}", ("abbb", "abbbb", "abbbbbbbbb"), ("abb", "bbb")),
        ("abc$|abcd|bc|d", ("abc", "abcd", "bc", "d"),
                           ("xabc", "ab", "a", "", "abce")),
        ("^a.*", ("a", "aa"), ("b", "ba", "c", "")),
        ("^[^b]+", ("a", "aa", "c"), ("b", "ba", "")),
        ("a(?!b).b?", ("aa", "ac", "aab"), ("a", "ab", "abc")),
        ("a(?=[bc])..", ("abx", "acx", "aba"), ("ac", "ab", "adb")),

        ("ab?[bc]?", ("a", "ab", "abb", "ac"), ("", "cab", "x")),
        ("ab{2,4}c?", ("abb", "abbb", "abbbb", "abbbbc"),
         ("ab", "abc", "xabbb")),
        ("ab{2,4}$", ("abb", "abbb", "abbbb"),
         ("ab", "abc", "xabbb", "abbbbc", "abbbbbb")),
        ("ab{2,4}cd?", ("abbc", "abbbc", "abbbbc", "abbbbcd"),
         ("abc", "abbbbbc", "abcbbc")),
        ("ab?c", ("ac", "abc"), ("abb", "abbc", "abbbc")),

        (r"\R", ("\n", "\r", "\r\n"), (" ", "\r\r", "\n\n", "\r\n ")),
        (r"a\Rb\R", ("a\nb\n", "a\rb\r", "a\r\nb\r\n", "a\rb\r\n"),
         ("ab", "a", "a\n\nb\n", "a\nb", "a\r\nb")),
        (r"ID [^\R]+\R", ("ID A123\n", "ID A123\r", "ID A123\r\n"),
         ("ID A123\n\n", "ID A123\r\r", "ID A123", "ID \n")),
        )
    for re_pat, good_list, bad_list in patterns:
        print "Testing", repr(re_pat)
        tree = convert_re.make_expression(re_pat)
        exp = Parser.Parser(tuple(_generate(tree, {}, 0)))
        exp.setContentHandler(cb)
        exp.setErrorHandler(cb)
        if string.find(re_pat, r"\R") == -1:  # \R is a Martel-specific flag
            pat = re.compile(re_pat)
        else:
            pat = None
        pat2 = re.compile(str(tree))
        
        for word in good_list:
            exp.parseString(word)

            if pat is not None:
                m = pat.match(word)
                assert m, "Re problem recognizing " + repr(word)
                assert m.end() == len(word), "Did not parse all of %s: %d" % \
                       (repr(word), m.end())

            m = pat2.match(word)
            assert m, "created Re problem recognizing " + repr(word)
            assert m.end() == len(word), "Did not parse all of created %s: %d"\
                   % (repr(word), m.end())
            
            assert cb.good_parse, "Problem not recognizing " + repr(word)
            
        for word in bad_list:
            exp.parseString(word)

            if pat is not None:
                m = pat.match(word)
                assert not m or m.end() != len(word), \
                       "Re should not recognize " + repr(word)

            m = pat2.match(word)
            assert not m or m.end() != len(word), \
                   "created Re should not recognize " + repr(word)
                   
            assert not cb.good_parse, \
                   "Should not recognize %s\ntagtable is %s" % \
                   (repr(word), repr(exp.tagtable))

if __name__ == "__main__":
    test()
    
