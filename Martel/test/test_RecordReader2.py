# Still more tests of the RecordReader.  This one stresses the ability
# to find newlines and pass in lookahead text

from cStringIO import StringIO
import string
from Martel import RecordReader
from mx import TextTools as TT

def count_records(reader):
    i = 0
    #print "Testing new reader", reader
    while 1:
        x = reader.next()
        #print "Read", repr(x)
        if x is None:
            break
        i = i + 1
    return i

# Notice how this omits the final newline?
def normalize(s):
    return string.join(TT.splitlines(s), "\n")
    


data1 = """\
ID   Q1234
DE   Some protein
SQ   blah
     ABCDE FGHIJ
//
ID   Q2345
DE   ID
SQ   lahb
 ID  just checking
     BCDEF GHIJA
//
ID   Q3456
DE   me proteinSo
SQ   ahbl
     CDEFG HIJID
//
"""
data1 = normalize(data1)

### StartsWith

def test_startswith_generic():
    to_eol_data = (("A\n", 1),
                   ("AA\n", 1),
                   ("A\nB\n", 1),
                   ("A\nB\nB\n", 1),
                   ("A\nB\nB\nA\n", 2),
                   ("A\nA\nB\nA\n", 3),
                   ("A\n A\nA\n", 2),
                   ("A A A A A A A A A A A A A A\nB\nA\n", 2),
                   )
    for s, expected in to_eol_data:
        reader = RecordReader.StartsWith(StringIO(s), "A")
        count = count_records(reader)
        assert count == expected, (s, expected, count) # skips to EOL
                   

def test_startswith_SP():
    # Check using a SWISS-PROT-like format
    for ending in ("\n", "\r", "\r\n"):
        s = string.replace(data1, "\n", ending)
        for final in ("", ending):
            d = s + final
            
            for i in range(5, 20):
                infile = StringIO(d)
                reader = RecordReader.StartsWith(infile, "ID", i)
                count = count_records(reader)
                assert count == 3, (ending, final, i, count, d)

            for i in range(6, 20):
                infile = StringIO(d)
                reader = RecordReader.StartsWith(infile, "ID ", i)
                count = count_records(reader)
                assert count == 3, (ending, final, i, count, d)

def test_startswith_exhaustive(ending):
    # Exhaustive test of the various combinations.  Should catch most
    # edge conditions.
    for base in ("A" + ending, "A" + ending + "BA" + ending):
        for repeat in range(0, 15):
            s = base * repeat
            infile = StringIO(s)
            #for marker in ("A", "A\n"):  # Don't use; bug when using "\n"
            for marker in ("A",):
                for look in range(5):
                    lookahead = base * look
                    for readhint in range(4, 10):
                            infile.seek(0)
                            reader = RecordReader.StartsWith(\
                                infile, marker, readhint, lookahead)
                            count = count_records(reader)
                            assert count == repeat + look, \
                                   (count, ending, base, repeat, marker,
                                    look, readhint)
                    infile.seek(0)
                    reader = RecordReader.StartsWith(infile, marker,
                                                     lookahead = lookahead)
                    count = count_records(reader)
                    assert count == repeat + look, \
                           (count, ending, repeat, marker, look)

def test_startswith_remainder():
    # Make sure the remainder method works
    for repeat in range(20):
        vals = map(lambda x: "A\n%d\n" % x, range(repeat))
        data = string.join(vals, "")
        infile = StringIO(data)
        for look in range(10) + range(10, len(data), 5):
            infile.seek(look)
            lookahead = data[:look]
            reader = RecordReader.StartsWith(infile, "A",
                                             lookahead = lookahead)
            all = ""
            while 1:
                file, lh = reader.remainder()
                pos = file.tell()
                rest = file.read()
                assert all + lh + rest == data, (all, lh, rest, data)
                file.seek(pos)
                assert data.startswith(all), (data, all)
                record = reader.next()
                if record is None:
                    break
                all = all + record
            assert all == data, (all, data)
            
                
def test_startswith_errors():
    # Check the failure cases.  Actually, there's only one.

    # Doesn't start with A
    for s in ("B", "B\n", " A\n", " A", "B\nA\n"):
        try:
            infile = StringIO(s)
            # The current implementation will fail here, but the
            # interface spec allows the error to be unreported
            # until reading the record.
            reader = RecordReader.StartsWith(infile, "A")
            rec = reader.next()
            raise AssertionError, "should not allow %r" % s
        except RecordReader.ReaderError:
            pass
        else:
            raise AssertionError, "should not get here"

def test_startswith():
    print "Testing StartsWith"

    print " ... generic"
    test_startswith_generic()

    print " ... newline variations"
    test_startswith_SP()

    for ending in ("\n", "\r", "\r\n"):
        print " ... exhaustive testing against %r" % ending
        test_startswith_exhaustive(ending)

    print " ... remainder"
    test_startswith_remainder()

    print " ... format errors"
    test_startswith_errors()

### EndsWith

def test_endswith_generic():
    to_eol_data = (("A\n", 1),
                   ("AA\n", 1),
                   ("B\nA\n", 1),
                   ("B\nB\nA\n", 1),
                   ("A\nA\n", 2),
                   ("A A\nA\n", 2),  # this changes with an "A\n" reader
                   ("A", 1),
                   ("A\nA A\nA\n", 3), # this changes with an "A\n" reader
                   ("A\nA A\nA", 3), # this changes with an "A\n" reader
                   )
    for s, expected in to_eol_data:
        reader = RecordReader.EndsWith(StringIO(s), "A")
        count = count_records(reader)
        assert count == expected, (s, expected, count)  # skips to EOL


    newline_data = (("A\n", 1),
                    #("AA\n", 1),  # not legal
                    ("B\nA\n", 1),
                    ("B\nB\nA\n", 1),
                    ("A\nA\n", 2),
                    ("A A\nA\n", 1),  # this changed with an "A\n" reader
                    ("A", 1),
                    ("A\nA A\nA\n", 2), # this changed with an "A\n" reader
                    ("A\nA A\nA", 2), # this changed with an "A\n" reader
                    )
    for s, expected in newline_data:
        reader = RecordReader.EndsWith(StringIO(s), "A\n")
        count = count_records(reader)
        assert count == expected, (s, expected, count)  # expects newline

        
def test_endswith_SP():
    # Check using a SWISS-PROT-like format
    for ending in ("\n", "\r", "\r\n"):
        s = string.replace(data1, "\n", ending)
        for final in ("", ending):
            d = s + final
            
            loop = 0
            for i in range(5, 20):
                infile = StringIO(d)
                reader = RecordReader.EndsWith(infile, "//", i)
                count = count_records(reader)
                assert count == 3, (ending, final, i, count, d)

            for i in range(5, 20):
                infile = StringIO(d)
                reader = RecordReader.EndsWith(infile, "//\n", i)
                count = count_records(reader)
                assert count == 3, (ending, final, i, count, d)


def test_endswith_exhaustive(ending):
    # Exhaustive test of the various combinations.  Should catch most
    # edge conditions.
    for base in ("A" + ending, "BA" + ending + "A" + ending):
        for repeat in range(0, 15):
            s = base * repeat
            infile = StringIO(s)
            for marker in ("A", "A\n"):
                for look in range(5):
                    lookahead = base * look
                    for readhint in range(4, 10):
                            infile.seek(0)
                            reader = RecordReader.EndsWith(\
                                infile, marker, readhint, lookahead)
                            count = count_records(reader)
                            assert count == repeat + look, \
                                   (count, ending, base, repeat, marker,
                                    look, readhint)
                    infile.seek(0)
                    reader = RecordReader.EndsWith(infile, marker,
                                                   lookahead = lookahead)
                    count = count_records(reader)
                    assert count == repeat + look, \
                           (count, ending, repeat, marker, look)

def test_endswith_remainder():
    # Make sure the remainder method works
    for repeat in range(20):
        vals = map(lambda x: "%d\nA\n" % x, range(repeat))
        data = string.join(vals, "")
        infile = StringIO(data)
        for look in range(10) + range(10, len(data), 5):
            infile.seek(look)
            lookahead = data[:look]
            reader = RecordReader.EndsWith(infile, "A",
                                           lookahead = lookahead)
            all = ""
            while 1:
                file, lh = reader.remainder()
                pos = file.tell()
                rest = file.read()
                assert all + lh + rest == data, (all, lh, rest, data)
                file.seek(pos)
                assert data.startswith(all), (data, all)
                record = reader.next()
                if record is None:
                    break
                all = all + record
            assert all == data, (all, data)

def test_endswith_errors():
    # Check the failure cases.

    # Could no record at all
    # Could be some records followed by an incomplete record
    # Could be a line which partially matches the data
    for s in ("B", "B\n", "A\nB\n", "A\nB\nA\nB\n", "A\nB\nA\n ", "AA",
              "AA\n", "A\nB\nA X\n"):
        has_error = 0
        infile = StringIO(s)
        try:
            reader = RecordReader.EndsWith(infile, "A\n")
        except RecordReader.ReaderError:
            has_error = 1

        if not has_error:
            while not has_error:
                try:
                    rec = reader.next()
                except RecordReader.ReaderError:
                    has_error = 1
                if not has_error and rec is None:
                    break
        if not has_error:
            raise AssertionError, "should not get here with %r" % s

    # Could no record at all
    # Could be some records followed by an incomplete record
    # *Allowed* to read rest of line
    for s in ("B", "B\n", "A\nB\n", "A\nB\nA\nB\n", "A\nB\nA\n "):
        has_error = 0
        infile = StringIO(s)
        try:
            reader = RecordReader.EndsWith(infile, "A")
        except RecordReader.ReaderError:
            has_error = 1

        if not has_error:
            while not has_error:
                try:
                    rec = reader.next()
                except RecordReader.ReaderError:
                    has_error = 1
                if not has_error and rec is None:
                    break
        if not has_error:
            raise AssertionError, "should not get here with %r" % s
    


def test_endswith():
    print "Testing EndsWith"

    print " ... generic"
    test_endswith_generic()

    print " ... newline variations"
    test_endswith_SP()

    for ending in ("\n", "\r", "\r\n"):
        print " ... exhaustive testing against %r" % ending
        test_endswith_exhaustive(ending)
        
    print " ... remainder"
    test_endswith_remainder()

    print " ... format errors"
    test_endswith_errors()

### Until

# Don't need to do that much testing since the code is built on
# top of the StartsWith reader, which has already been tested.
def test_until():
    # Can only read at most one record
    print "Testing Until"
    test_data = ("A",
                 "A\n",
                 "A\nB",
                 "A\nB\n",
                 "A\nBCDE\nQWE\nTRYU\nA\n",
                 "A\nA\nA\nA\nA\nA\nA\n",
                 "AB",
                 "AB\n",
                 "AB\nAC",
                 "AB\nAC\n",
                 )

    for ending in ("\n", "\r", "\r\n"):
        for pre in ("", "B\n", "BA\n", "CA\nBA\n"):
            pre_nl = string.replace(pre, "\n", ending)
            for look in (0, 1, 3):
                for text in test_data:
                    text_nl = string.replace(text, "\n", ending)
                    s = pre_nl + text_nl
                    reader = RecordReader.Until(StringIO(s[look:]),
                                                "A",
                                                lookahead = s[:look],
                                                sizehint = 4)
                    found_record = None
                    while 1:
                        rec = reader.next()
                        if rec is None:
                            break
                        assert not found_record, \
                               "Already found %r but also found %r in %r" % \
                               (found_record, rec, s)
                        found_record = rec

                    assert found_record == pre_nl, \
                           "Expecting record %r, found %r in %r" % \
                           (pre_nl, found_record, s)
                    infile, remainder = reader.remainder()
                    remainder = remainder + infile.read()
                    assert remainder == text_nl, \
                           "Expecting remainder %r, found %r in %r" % \
                           (text_nl, remainder, s)

### CountLines

def test_count_lines():
    # Create a set of 'i' lines and read 'count' lines at a time.
    # Either count divides i or it doesn't.  If it does, the reader
    # should go to completion.  If it does not, the reader should
    # have a remainder whose size can be verified.
    
    print "Testing CountLines"
    for ending in ("\n", "\r", "\r\n"):
        print " ... exhaustive testing against %r" % ending
        s = ""
        for i in range(25):
            for count in range(1,i+1):
                for look in (0, 2, 5):
                    rep, final = divmod(i, count)
                    reader = RecordReader.CountLines(StringIO(s[look:]),
                                                     count,
                                                     lookahead = s[:look],
                                                     sizehint = 1)
                    all = ""
                    while rep > 0:
                        rec = reader.next()
                        lines = string.split(rec, ending)
                        assert len(lines)-1 == count, \
                               "Expecting %d lines, got %d in %r using %r" % \
                               (count, len(lines)-1, rec, ending)
                        all = all + rec
                        rep = rep - 1

                    if final == 0:
                        # Reader should be at the end of input
                        rec = reader.next()
                        assert rec is None, \
                               "Should be at end of reader, got %r" % rec
                        rec = reader.next()
                        assert rec is None, \
                               "data after end of reader, got %r" % rec
                    else:
                        # There should be a remainder of size i % count lines
                        infile, remainder = reader.remainder()
                        text = remainder + infile.read()
                        all = all + text
                        lines = string.split(text, ending)
                        assert len(lines)-1 == final, \
                               "Expecting %d final lines, got %d" % \
                               (final, len(lines)-1)
                        try:
                            rec = reader.next()
                            raise AssertionError, \
                                  "Got unexpected final record, %r" % rec
                        except RecordReader.ReaderError:
                            pass
                    assert all == s, \
                           "record data %r doesn't rebuild input %r" % \
                           (all, s)
            s = s + str(i) + ending
    
### Nothing

def test_nothing():
    print "Testing Nothing"
    
    s = "This is a test.\nThis is only a test.\nHad this been an actual...\n"
    for ending in ("\n", "\r", "\r\n"):
        data = string.replace(s, "\n", ending)
        for look in (0, 1, 2, 5):
            reader = RecordReader.Nothing(StringIO(data[look:]),
                                          sizehint = 1,
                                          lookahead = data[:look])
            rec = reader.next()
            assert rec is None, "should be empty, not %r" % rec
            rec = reader.next()
            assert rec is None, "2nd time should also be empty, not %r" % rec
            
            infile, remainder = reader.remainder()
            remainder = remainder + infile.read()
            assert remainder == data, "Why %r when input was %r?" % \
                   (remainder, data)

### Everything

def test_everything():
    print "Testing Everything"

    s = "This is a test.\nThis is only a test.\nHad this been an actual...\n"
    for ending in ("\n", "\r", "\r\n"):
        data = string.replace(s, "\n", ending)
        for look in (0, 1, 2, 5):
            reader = RecordReader.Everything(StringIO(data[look:]),
                                             sizehint = 1,
                                             lookahead = data[:look])
            rec = reader.next()
            assert rec == data, "Record %r is not same as input %r" % \
                   (rec, data)
            infile, remainder = reader.remainder()
            remainder = remainder + infile.read()
            assert not remainder, "Why is there a remainder of %r?" % \
                   remainder

            rec = reader.next()
            assert rec is None, "Expecting None after final read, got %r" % \
                   rec
            rec = reader.next()
            assert rec is None, "Expecting None (again), got %r" % rec
            


### test driver
    
def test():
    test_startswith()
    test_endswith()
    test_until()
    test_count_lines()
    test_nothing()
    test_everything()
    
if __name__ == "__main__":
    test()
