# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

# The existing parsers are in-memory.  For large data files, like
# swissprot, that requires too much memory.

# On the other hand, the records aren't all that large (there's just a
# lot of them.)  This module has readers which scan the input to get a
# record as a string.

import string
from mx import TextTools as TT

SIZEHINT = 100000

class ReaderError(TypeError):
    pass

class RecordReader:
    def __init__(self, infile):
        self.infile = infile
    def next(self):
        raise NotImplementedError
    def remainder(self):
        raise NotImplementedError
    
def _startswith_tagtable_rest_of_line(text):
    return (
        # Ensure the text starts with the given word
        ("begin", TT.Word, text, TT.MatchFail, +1),

        # Read to the end of line
        (None, TT.AllInSet, TT.invset('\r\n'), +1, +1),

        # Read the end of line
        (None, TT.Is, '\n', +1, +4),  # matches '\n' or
        (None, TT.Is, '\r', +2, +1),  # '\r' followed by
        (None, TT.Is, '\n', +2, +2),  # optional '\n'

        # Check if EOF (allow EOF if no EOL found)
        (None, TT.EOF, TT.Here, +1, TT.MatchOk),

        # Not EOF, so look for the next line starting with text
        ("begin", TT.Word, text, +1, -5),

        # Not what I am looking for, so read to the end of line
        (None, TT.AllInSet, TT.invset('\r\n'), +1, +1),

        # Read the end of line then test the next line
        (None, TT.Is, '\n', +1, -2),  # '\n'
        (None, TT.Is, '\r', +2, +1),  # '\r' followed by
        (None, TT.Is, '\n', -4, -4),  # optional '\n'
        # Allow termination at EOF
        (None, TT.EOF, TT.Here, TT.MatchFail, TT.MatchOk),    
        )

def _startswith_tagtable_newline(text):
    return (
        # Ensure the text starts with the given word ...
        ("begin", TT.Word, text, TT.MatchFail, +1),

        # ... followed by the end of line
        (None, TT.Is, '\n', +1, +4),  # matches '\n' or
        (None, TT.Is, '\r', +2, +1),  # '\r' followed by
        (None, TT.Is, '\n', +2, +2),  # optional '\n'

        # Check if EOF instead of a newline (allow EOF if found)
        # Otherwise, this means the line starts with the text but
        # doesn't have a successive newline.
        # XXX BUG! When looking for "A\n" should not fail on "AA\n"!
        (None, TT.EOF, TT.Here, TT.MatchFail, TT.MatchOk),

        # Look for the next line starting with text
        ("begin", TT.Word, text, +1, -4),

        # Not what I am looking for, so read to the end of line
        (None, TT.AllInSet, TT.invset('\r\n'), +1, +1),

        # Read the end of line then test the next line
        (None, TT.Is, '\n', +1, -2),  # '\n'
        (None, TT.Is, '\r', +2, +1),  # '\r' followed by
        (None, TT.Is, '\n', -4, -4),  # optional '\n'
        # Allow termination at EOF
        (None, TT.EOF, TT.Here, TT.MatchFail, TT.MatchOk),    
        )
        

def _find_begin_positions(text, tagtable):
    success, tags, pos = TT.tag(text, tagtable)
    # print "XXX", success, tags, pos, len(text)
    if not success:
        raise ReaderError("invalid format starting with %s" % repr(text[:50]))
    if pos != len(text):
        raise ReaderError, \
            "could not parse to end of text (ended at %d of %d)" % \
            (pos, len(text))
    return [tag[1] for tag in tags]
                        

class StartsWith(RecordReader):
    def __init__(self, infile, text, sizehint = SIZEHINT, lookahead = ""):
        RecordReader.__init__(self, infile)
        self.text = text
        self.sizehint = sizehint

        pos = string.find(text, "\n")
        if pos != -1:
            if pos != len(text)-1:
                raise AssertionError, "'\\n' can only exist at the end of the string"
            text = text[:-1]
            has_newline = 1
        else:
            has_newline = 0
        assert len(text), "StartsWith text size is too short"
        assert len(text) < sizehint - 2, \
               "StartsWith text size larger than sizehint allows"


        if has_newline:
            raise NotImplementedError, "there's a bug in the '\\n' option"
            self.tagtable = _startswith_tagtable_newline(text)
        else:
            self.tagtable = _startswith_tagtable_rest_of_line(text)

        self.lookahead = lookahead

        # Start parsing here.  This guarantees the first line is in
        # the right format.
        if len(self.lookahead) < len(text) + 2:
            self.lookahead += infile.read(sizehint)
        if self.lookahead:
                self.positions = _find_begin_positions(self.lookahead,
                                                       self.tagtable)
        else:
            self.positions = [0]
        self.index = 0
        
    def next(self):
        # Are any precomputed positions remaining?
        if self.index + 1 < len(self.positions):
            # Yes, so return the text in the range
            s = self.lookahead[self.positions[self.index]:
                               self.positions[self.index+1]]
            self.index += 1
            return s

        # The last record may be incomplete, so reset the
        # lookahead to be just that text
        self.lookahead = self.lookahead[self.positions[-1]:]

        # Read past at least the start of the second record or to the
        # end of file.
        positions = [self.positions[-1]]
        while 1:
            data = self.infile.read(self.sizehint)
            if not data:
                break
            self.lookahead = self.lookahead + data
            positions = _find_begin_positions(self.lookahead, self.tagtable)
            if len(positions) > 1:
                break
        if len(positions) > 1:
            self.positions = positions
            self.index = 1
            return self.lookahead[positions[0]:positions[1]]
        elif not self.lookahead:
            # No data (either empty file or at EOF)
            self.positions = [0]
            self.index = 0
            return None
        else:
            # Read to the end of file and it's the last record
            assert len(positions) == 1
            self.positions = [0]
            self.index = 0
            s = self.lookahead
            self.lookahead = ""
            return s

    def remainder(self):
        return self.infile, self.lookahead[self.positions[self.index]:]

def _endswith_tagtable_newline(text):
    return (
        # Is the current line the end of record marker?
        (None, TT.Word, text, +6, +1),
 
        # Make sure it ends the line
        ("end", TT.Is, '\n', +1, -1),  # matches '\n'
        (None, TT.Is, '\r', +4, +1),
        ("end", TT.Is, '\n', +1, -3),
        (None, TT.Skip, -1, +1, +1),
        ("end", TT.Skip, +1, -5, -5),
 
        # Not the end of record marker, so read to the end of line
        (None, TT.AllInSet, TT.invset('\r\n'), +1, +1),
 
        # Check if EOF
        (None, TT.EOF, TT.Here, +1, TT.MatchOk),
 
        # Not EOF, so scarf any newlines
        (None, TT.AllInSet, TT.set('\r\n'), TT.MatchFail, -8),
        )

def _endswith_tagtable_rest_of_line(text):
    return (
        # Is the current line the end of record marker?
        (None, TT.Word, text, +8, +1),

        # Read whatever else is on that line (could be nothing)
        (None, TT.AllInSet, TT.invset('\r\n'), +1, +1),
 
        # Get the end of line
        ("end", TT.Is, '\n', +1, -2),  # matches '\n'
        (None, TT.Is, '\r', +4, +1),
        ("end", TT.Is, '\n', +1, -4),
        (None, TT.Skip, -1, +1, +1),
        ("end", TT.Skip, +1, -6, -6),
 
        # Check if EOF (only tests when the end of record line has no \n)
        # Only time this should fail is with a bug in TT.
        ("end", TT.EOF, TT.Here, TT.MatchFail, TT.MatchOk),
        
        # Not the end of record marker, so read to the end of line
        (None, TT.AllInSet, TT.invset('\r\n'), +1, +1),
 
        # Check if EOF
        (None, TT.EOF, TT.Here, +1, TT.MatchOk),
 
        # Not EOF, so scarf any newlines and try again
        (None, TT.AllInSet, TT.set('\r\n'), TT.MatchFail, -10),
        )
        

def _find_end_positions(text, tagtable):
    success, tags, pos = TT.tag(text, tagtable)
    #print "XXX", success, tags, pos, len(text), repr(text)
    if not success:
        raise ReaderError("invalid format starting with %s" % repr(text[:50]))
    if pos != len(text):
        raise ReaderError, \
            "could not parse to end of text (ended at %d of %d)" % \
            (pos, len(text))
    return [tag[2] for tag in tags]

class EndsWith(RecordReader):
    def __init__(self, infile, text, sizehint = SIZEHINT, lookahead = ""):
        RecordReader.__init__(self, infile)
        self.text = text
        self.sizehint = sizehint
        
        pos = string.find(text, "\n")
        if pos != -1:
            if pos != len(text)-1:
                raise AssertionError, "'\\n' can only exist at the end of the string"
            text = text[:-1]
            has_newline = 1
        else:
            has_newline = 0
        assert len(text) < sizehint - 2, \
               "EndsWith text size larger than sizehint allows"

        if has_newline:
            self.tagtable = _endswith_tagtable_newline(text)
        else:
            self.tagtable = _endswith_tagtable_rest_of_line(text)


        self.lookahead = lookahead
        self.positions = []
        self.index = 0
        self.pos = 0

    def next(self):
        # Are any precomputed positions remaining?
        if self.index < len(self.positions):
            # Yes, so return that text
            newpos = self.positions[self.index]
            s = self.lookahead[self.pos:newpos]
            self.pos = newpos
            self.index = self.index + 1
            return s

        # No available information, so use what remains to seed the
        # next level.
        lookahead = self.lookahead[self.pos:]

        data = ""
        positions = []
        # Add text until I've found a record or there is no more data
        while 1:
            data = self.infile.read(self.sizehint)
            if not data:
                if not positions:
                    positions = _find_end_positions(lookahead, self.tagtable)
                break
            lookahead = lookahead + data
            positions = _find_end_positions(lookahead, self.tagtable)
            if len(positions) > 1:
                del positions[-1]
                break

        self.lookahead = lookahead
        self.positions = positions
        
        if positions:
            self.index = 1
            self.pos = positions[0]
            return lookahead[:positions[0]]
        elif not lookahead:
            # No data (either empty file or at EOF)
            self.pos = 0
            self.index = 0
            return None

        # This is likely an unterminated record.  However, it could be
        # that there is no terminal end-of-line character so check for
        # that case.
        if lookahead[-1:] not in "\r\n":
            special_case = lookahead + "\n"
            positions = _find_end_positions(special_case, self.tagtable)
            if positions:
                assert len(positions) == 1, "this case should not occur"
                pos = positions[0]
                assert pos == len(special_case), "wrong sizes: %d and %d" % \
                       (pos, len(special_case))
                self.lookahead = ""
                self.positions = []
                self.pos = 0
                self.index = 0
                return lookahead

        # Really could not find a terminator
        self.index = 0
        self.pos = 0
        raise ReaderError("Last record not terminated: at %s ..." %
                          repr(self.lookahead[:50]))

    def remainder(self):
        return self.infile, self.lookahead[self.pos:]



class Until(RecordReader):
    def __init__(self, infile, text, sizehint = SIZEHINT, lookahead = ""):
        RecordReader.__init__(self, infile)
        self.text = text
        self.lookahead = lookahead
        self.sizehint = sizehint
        self.found = 0

        if text[-1] == "\n":
            raise NotImplementedError, "Until reader does not support '\\n'"
        if "\n" in text:
            raise AssertionError, "'\\n' can only exist at the end of the string"

    def next(self):
        if self.found:
            return None

        # Use the StartsWith reader to get to the end of this record.
        # Need to fake the first line..
        fake = self.text + "\n"
        reader = StartsWith(self.infile, self.text, self.sizehint,
                            fake + self.lookahead)
        rec = reader.next()
        rec = rec[len(fake):]  # remove the fake data
        self.infile, self.lookahead = reader.remainder()
        self.found = 1
        return rec

    def remainder(self):
        return self.infile, self.lookahead

# Tag the last byte of every newline
_tag_lines_tagtable = (
    # Skip non-newline characters
    (None, TT.AllInSet, TT.invset('\r\n'), +1, +1),

    # Check if newline
    ("newline", TT.Is, '\n', +1, -1),  # can be '\n'
    (None, TT.Is, '\r', +3, +1),       # or start a '\r' followed by ..
    ("newline", TT.Is, '\n', +1, -3),  #  .. an optional '\n'
    ("newline", TT.Skip, 0, -4, -4),   # get here with just an '\r'
    (None, TT.EOF, TT.Here, -5, TT.MatchOk),  # stop at end of text
    )


class CountLines(RecordReader):
    """Read a specified (fixed) number of lines"""
    def __init__(self, infile, count, sizehint = SIZEHINT, lookahead = ""):
        assert count > 0, "CountLines reader must read at least one line"
        assert lookahead > 0, "Must read at least a character at a time"
        assert sizehint > 0, "sizehint must be positive"
        RecordReader.__init__(self, infile)
        self.sizehint = sizehint
        self.lookahead = lookahead
        self.count = count
        self.pos = 0
        self.positions = []
        self.index = 0

    def next(self):
        if self.index + self.count < len(self.positions):
            self.index = self.index + self.count
            endpos = self.positions[self.index-1]
            s = self.lookahead[self.pos:endpos]
            self.pos = endpos
            return s
            
        lookahead = self.lookahead[self.pos:]
        while 1:
            positions = _find_end_positions(lookahead, _tag_lines_tagtable)
            if len(positions) > self.count:
                # Last line may be incomplete, as with "\r" of "\r\n"
                # Hmm, is this *really* needed?  Doesn't hurt. XXX
                del positions[-1]
                break
            data = self.infile.read(self.sizehint)
            if not data:
                break
            lookahead = lookahead + data

        self.lookahead = lookahead
        self.pos = 0
        self.positions = positions

        if not lookahead:
            return None

        if len(positions) >= self.count:
            self.index = self.count
            endpos = self.positions[self.count-1]
            s = lookahead[0:endpos]
            self.pos = endpos
            return s

        # Commented out to require final newline
        # Don't allow that case since it is more frequent that the line
        # count is wrong.  (Could change in future releases.)
##        elif len(positions) == self.count - 1 and not data:
##            # This was the last record, and it has no trailing newline
##            s = self.lookahead
##            self.lookahead = ""
##            self.positions = []
##            self.index = 0
##            return s
        raise ReaderError, \
              "Only found %d lines, expecting %d (starting with %s ...)" % \
              (len(positions), self.count, repr(lookahead[:20]))

    def remainder(self):
        return self.infile, self.lookahead[self.pos:]

class Nothing(RecordReader):
    """Reads nothing"""
    def __init__(self, infile, sizehint = SIZEHINT, lookahead = ""):
        RecordReader.__init__(self, infile)
        self.lookahead = lookahead

    def next(self):
        return None
    
    def remainder(self):
        return self.infile, self.lookahead

class Everything(RecordReader):
    """Reads everything"""
    def __init__(self, infile, sizehint = SIZEHINT, lookahead = ""):
        RecordReader.__init__(self, infile)
        self.lookahead = lookahead
        self.found = 0

    def next(self):
        if self.found:
            return None
        s = self.lookahead + self.infile.read()
        self.lookahead = ""
        self.found = 1
        return s

    def remainder(self):
        return self.infile, self.lookahead
