"""Wrapper around a file handle to allow reseeks to the beginning.

This will only seek back to the initial position, and only handles the
'read' method, which is all that Martel uses.
"""


from cStringIO import StringIO

class ReseekFile:
    def __init__(self, file):
        self.file = file
        self.buffer_file = StringIO()
        self.at_beginning = 1
        try:
            self.beginning = file.tell()
        except (IOError, AttributeError):
            self.beginning = 0
        self._use_buffer = 1
        
    def seek(self, offset, whence = 0):
        if whence != 0:
            raise TypeError("Unexpected whence value of %s; expecting 0" % \
                            (whence,))
        if offset != self.beginning:
            raise TypeError("Unexpected offset value of %r; expecting '%s'" % \
                             (offset, self.beginning))
        self.buffer_file.seek(0)
        self.at_beginning = 1
        
    def tell(self):
        if not self.at_beginning:
            raise TypeError("ReseekFile cannot tell except at the beginning of file")
        return self.beginning

    def _read(self, size):
        if size < 0:
            y = self.file.read()
            z = self.buffer_file.read() + y
            self.buffer_file.write(y)
            return z
        if size == 0:
            return ""
        x = self.buffer_file.read(size)
        if len(x) < size:
            y = self.file.read(size - len(x))
            self.buffer_file.write(y)
            return x + y
        return x
        
    def read(self, size = -1):
        x = self._read(size)
        if self.at_beginning and x:
            self.at_beginning = 0
        self._check_no_buffer()
        return x

    def _check_no_buffer(self):
        if self._use_buffer == 0 and self.buffer_file.tell() == \
                                        len(self.buffer_file.getvalue()):
            self.seek = self.file.seek
            self.tell = self.file.tell
            self.read = self.file.read
            del self.buffer_file

    def nobuffer(self):
        self._use_buffer = 0

def prepare_input_source(source):
    from xml.sax import saxutils
    source = saxutils.prepare_input_source(source)
    # Is this correct?  Don't know - don't have Unicode exprerience
    f = source.getCharacterStream() or source.getByteStream()
    try:
        f.tell()
    except (AttributeError, IOError):
        f = ReseekFile.ReseekFile(f)
        source.setByteStream(f)
        source.setCharacterStream(None)
    return source

def test_reads(test_s, file, seek0):
    assert file.read(2) == "Th"
    assert file.read(3) == "is "
    assert file.read(4) == "is a"
    assert file.read(0) == ""
    assert file.read(0) == ""
    assert file.read(6) == " test."
    file.seek(seek0)
    assert file.read(2) == "Th"
    assert file.read(3) == "is "
    assert file.read(4) == "is a"
    assert file.read(0) == ""
    assert file.read(0) == ""
    assert file.read(6) == " test."
    assert file.read(1) == "\n"
    assert file.read(5) == "12345"
    assert file.read() == "67890\n"
    file.seek(seek0)
    assert file.read() == test_s
    file.seek(seek0)

    
def test():
    s = "This is a test.\n1234567890\n"
    file = StringIO(s)
    # Test with a normal file
    x = file.tell()
    test_reads(s, file, x)
    test_reads(s, file, x)

    # Test with a ReseekFile wrapper
    rf = ReseekFile(file)
    y = rf.tell()
    rf.seek(y)
    test_reads(s, rf, y)
    assert rf.read() == s
    assert rf.read() == ""

    # Make sure the tell offset is correct (may not be 0)
    file = StringIO("X" + s)
    file.read(1)
    rf = ReseekFile(file)
    y = rf.tell()
    test_reads(s, rf, y)
    rf.seek(y)
    test_reads(s, rf, y)
    assert rf.read() == s
    assert rf.read() == ""

    # Test the ability to turn off buffering and have changes
    # propogate correctly
    file = StringIO("X" + s)
    file.read(1)
    rf = ReseekFile(file)
    y = rf.tell()
    assert y == 1
    rf.read(1000)
    rf.seek(y)
    rf.nobuffer()
    assert rf.tell() == y
    test_reads(s, rf, y)
    rf.seek(y)
    test_reads(s, rf, y)
    assert rf.read() == s
    assert rf.read() == ""

    # turn off buffering after partial reads
    file = StringIO("X" + s)
    file.read(1)
    rf = ReseekFile(file)
    y = rf.tell()
    rf.read(5)
    rf.seek(y)
    rf.nobuffer()
    assert rf.read() == s

    file = StringIO("X" + s)
    file.read(1)
    rf = ReseekFile(file)
    y = rf.tell()
    t = rf.read(5)
    rf.seek(y)
    rf.nobuffer()
    assert rf.read(5) == t

    file = StringIO("X" + s)
    file.read(1)
    rf = ReseekFile(file)
    y = rf.tell()
    t = rf.read(5)
    assert t == s[:5]
    rf.seek(y)
    rf.nobuffer()
    assert rf.read(8) == s[:8]

    file = StringIO("X" + s)
    file.read(1)
    rf = ReseekFile(file)
    y = rf.tell()
    t = rf.read(5)
    assert t == s[:5]
    rf.nobuffer()
    assert rf.read(8) == s[5:5+8]

    # Should only do this test on Unix systems
    import os
    infile = os.popen("echo HELLO_THERE")
    infile.read(1)
    rf = ReseekFile(infile)
    y = rf.tell()
    assert rf.read(1) == "E"
    assert rf.read(2) == "LL"
    rf.seek(y)
    assert rf.read(4) == "ELLO"
    rf.seek(y)
    assert rf.read(1) == "E"
    rf.nobuffer()
    assert rf.read(1) == "L"
    assert rf.read(4) == "LO_T"
    assert rf.read(4) == "HERE"
    try:
        rf.seek(y)
        raise AssertionError("Cannot seek here!")
    except IOError:
        pass
    try:
        rf.tell()
        raise AssertionError("Cannot tell here!")
    except IOError:
        pass


if __name__ == "__main__":
    test()
