import compression

class Location:
    def __init__(self, namespace, name, filename, startpos, length):
        self.namespace = namespace
        self.name = name
        self.filename = filename
        self.startpos = startpos
        self.length = length
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
        raise AttributeError(key)

