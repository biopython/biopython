import os, string
import Martel
from Martel import RecordReader
from cStringIO import StringIO

def test_count(reader, check_remainder = 1):
    count = 0
    while 1:
        record = reader.next()
        if record is None:
            if check_remainder:
                infile, text = reader.remainder()
                if text:
                    raise AssertionError, text
                line = infile.readline()
                if line:
                    raise AssertionError, repr(line)
            return count
        #print repr(record)
        #print string.split(record, "\n")[0]
        count = count + 1

def test_fasta():
    s = """\
>seq1
SEQUENCE

>seq2
MULTI
LINE
SEQUENCE

>seq3

>seq4

"""
    reader = RecordReader.StartsWith(StringIO(s), ">")
    assert test_count(reader) == 4

    reader = RecordReader.StartsWith(StringIO(""), ">")
    assert test_count(reader) == 0

    reader = RecordReader.StartsWith(StringIO(">seq\n"), ">")
    assert test_count(reader) == 1

sp_sample = os.path.join("samples", "sample.swissprot")

def test_start():
    reader = RecordReader.StartsWith(open(sp_sample), "ID")
    assert test_count(reader, 0) == 8

def test_start_lines():
    lookahead = open(sp_sample).read()
    reader = RecordReader.StartsWith(StringIO(""), "ID", lookahead = lookahead)
    assert test_count(reader) == 8
    
def test_end():
    reader = RecordReader.EndsWith(open(sp_sample), "//\n")
    assert test_count(reader) == 8

def test_end_lines():
    lookahead = open(sp_sample).read()
    reader = RecordReader.EndsWith(StringIO(""), "//\n", lookahead = lookahead)
    assert test_count(reader) == 8


def test_until():
    s = "a\nID\nb\nID\n"
    reader = RecordReader.Until(StringIO(s), "ID")
    assert test_count(reader, check_remainder = 0) == 1
    assert reader.remainder()[1]

def test_until_lines():
    lookahead = "a\nID\nb\nID\n"
    reader = RecordReader.Until(StringIO(""), "ID", lookahead = lookahead)
    assert test_count(reader, check_remainder = 0) == 1
    assert reader.remainder()[1]

def test_count_lines():
    s = "1\n2\n3\n4\n5\n6\n7\n8\n"
    reader = RecordReader.CountLines(StringIO(s), 2)
    assert test_count(reader) == 4

def test_count_lines_lines():
    lookahead = "1\n2\n3\n4\n5\n6\n7\n8\n"
    reader = RecordReader.CountLines(StringIO(""), 2, lookahead = lookahead)
    assert test_count(reader) == 4

def test_nothing():
    s = "1\n2\n3\n4\n5\n6\n7\n8\n"
    infile = StringIO(s)
    reader = RecordReader.Nothing(infile)
    assert test_count(reader, check_remainder = 0) == 0
    assert infile.readline() == "1\n"

def test_nothing_lines():
    lookahead = "1\n2\n3\n4\n5\n6\n7\n8\n"
    reader = RecordReader.Nothing(StringIO(""), lookahead = lookahead)
    assert test_count(reader, check_remainder = 0) == 0
    file, result = reader.remainder()
    assert result == lookahead, (result, lookahead)


def test_everything():
    s = "1\n2\n3\n4\n5\n6\n7\n8\n"
    infile = StringIO(s)
    reader = RecordReader.Everything(infile)
    assert test_count(reader) == 1
    assert not infile.readline()

def test_everything_lines():
    lookahead = "1\n2\n3\n4\n5\n6\n7\n8\n"
    reader = RecordReader.Everything(StringIO(""), lookahead = lookahead)
    assert test_count(reader) == 1

def test():
    test_start()
    test_start_lines()
    test_end()
    test_end_lines()
    test_until()
    test_until_lines()
    test_count_lines()
    test_count_lines_lines()
    test_nothing()
    test_nothing_lines()
    test_everything()
    test_everything_lines()

    test_fasta()

    
if __name__ == "__main__":
    test()
    print "All tests passed."
