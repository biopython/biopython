import compression

class Location:
    """Handle for a record (use 'text' to get the record's text)"""
    def __init__(self, namespace, name, filename, startpos, length):
        self.namespace = namespace
        self.name = name
        self.filename = filename
        self.startpos = startpos
        self.length = length
    def __repr__(self):
        return "Location(namespace = %r, name = %r, filename = %r, startpos = %r, length = %r)" % (self.namespace, self.name, self.filename, self.startpos, self.length)
    def __str__(self):
        return "Location(%s:%s at %s: %s, %s)" % \
               (self.namespace, self.name,
                self.filename,self.startpos, self.length)
    def __getattr__(self, key):
        if key == "text":
            infile = compression.open_file(self.filename)
            if hasattr(infile, "seek"):
                infile.seek(self.startpos)
                return infile.read(self.length)
            # read 1MB chunks at a time
            CHUNKSIZE = 1000000
            count = 0
            while count + CHUNKSIZE < self.startpos:
                infile.read(CHUNKSIZE)
                count += CHUNKSIZE
            infile.read(self.startpos - count)
            return infile.read(self.length)
        elif key == "__members__":
            return ["text"]
        raise AttributeError(key)

