import commands, os

_uncompress_table = {
    ".bz": "bzip2",
    ".BZ": "bzip2",
    ".gz": "gzip",
    ".GZ": "gzip",
    ".Z": "compress",
    }

def open_file(filename, mode = "rb"):
    ext = os.path.splitext(filename)[1]
    type = _uncompress_table.get(ext)
    if type is None:
        return open(filename, mode)
    if type == "gzip":
        import gzip
        gzip.open(filename, mode)
    if type == "bzip2":
        cmd = "bzcat --decompress"
        cmd += commands.mkarg(filename)
        return os.popen(cmd, mode)
    if type == "compress":
        cmd = "zcat -d"
        cmd += commands.mkarg(filename)
        return os.popen(cmd, mode)
    raise AssertionError("What's a %r?" % type)
            

