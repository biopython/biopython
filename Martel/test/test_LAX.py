import string

import Martel
from Martel import LAX

def test1():
    fields = ( ["Andrew", "Dalke", "12"],
               ["Liz", "Nelson",  "22"],
               ["Mandrake", "Moose",  "23"],
               ["Lisa", "Marie",  "91"], )
    text = ""
    for line in fields:
        text = text + string.join(line, " ") + "\n"
    
    format = Martel.Rep1(
                 Martel.Group("line",
                     Martel.Word("name", {"type": "first"}) + \
                     Martel.Spaces() + \
                     Martel.Word("name", {"type": "last"}) + \
                     Martel.Spaces() + \
                     Martel.Integer("age") + \
                     Martel.AnyEol()
                              ))
    iterator = format.make_iterator("line")
    i = 0
    for record in iterator.iterateString(text, LAX.LAX()):
        assert record["name"] == fields[i][:2], (record["name"], fields[i][:2])
        assert record["age"] == fields[i][2:3], (record["age"], fields[i][2:3])
        i = i + 1

    i = 0
    for record in iterator.iterateString(text, LAX.LAXAttrs()):
        assert [x[0] for x in record["name"]] == fields[i][:2], \
                ([x[0] for x in record["name"]], fields[i][:2])
        assert [x[0] for x in record["age"]] == fields[i][2:3], \
               ([x[0] for x in record["age"]], fields[i][2:3])
        assert record["name"][0][1]["type"] == "first"
        assert record["name"][1][1]["type"] == "last"
        assert record["age"][0][1].keys() == []
        i = i + 1

def test_filter():
    #  8 stretches of "a"s
    # 10 stretches of "b"s
    #  4 stretches of "c"s
    data = "ababcbaaaababbbabccbaabcabcba"
    format = Martel.Re("((?P<a>a+)|(?P<b>b+)|(?P<c>c+))+")
    parser = format.make_parser()
    lax = LAX.LAX(["b", "c"])
    parser.setContentHandler(lax)
    parser.parseString(data)
    assert lax.has_key("a") == 0
    assert len(lax["b"]) == 10
    assert len(lax["c"]) == 4
    

def test():
    test1()
    test_filter()

if __name__ == "__main__":
    test()
