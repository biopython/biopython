from Martel.convert_re import *

def test():
    data = (
        ("a", "a"),
        ("ab", "ab"),
        ("a|b", "[ab]"),
        ("[ab]", "[ab]"),
        ("a|b|c|d", "[a-d]"),
        ("ab|bc", "ab|bc"),
        ("a*", "a*"),
        ("a+", "a+"),
        ("a?", "a?"),
        ("a{3}", "a{3}"),
        ("a{3,8}", "a{3,8}"),
        ("a{3,3}", "a{3}"),
        ("(?:test1)", "(test1)"),
        ("(?P<foo>test1)", "(?P<foo>test1)"),
        (".", "."),
        ("((?P<x>..?)|(?P<y>..*))", "((?P<x>..?)|(?P<y>..*))"),
        ("\w", "[\\dA-Z_a-z]"),
        ("\s", "[\\t-\\r ]"),
        ("\\d", "[\\d]"),
        ("[0123456789]", "[\\d]"),
        ("(0|1|2|3|4|5|6|7|8|9)", "([\\d])"),
        (r"This is (?!not\.)nothing\.", r"This is (?!not\.)nothing\."),
        (r"This is (?=not\.)nothing\.", r"This is (?=not\.)nothing\."),
        ("^start", "^start"),
        ("not^start", "not^start"),
        ("end$", "end$"),
        ("end$not", "end$not"),
        ("[a-z]", "[a-z]"),
        ("[^a-z]", "[^a-z]"),
        ("[^b]", "[^b]"),
        ("[a-zA-Z]", "[A-Za-z]"),
        ("[ababab]", "[ab]"),
        ("a{foo,bar}", "a{foo,bar}"),
        ("a{foo,foo}", "a{foo}"),
        ("[\\]]", "\\]"),   # it special cases a single character in []s
        ("[A\\]]", "[A\\]]"),
        )
    for input, output in data:
        result = str(make_expression(input))
        assert result == output, "input %s : expected %s but got %s" % \
               (repr(input), repr(output), repr(result))

if __name__ == "__main__":
    test()
