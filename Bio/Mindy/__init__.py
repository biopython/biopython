import os, sys

_open = open  # rename for internal use -- gets redefined below

def open(dbname, mode = "r"):
    text = _open(os.path.join(dbname, "config.dat"), "rb").read()
    line = text.split("\n")[0]
    if line == "index\tBerkeleyDB/1":
        import BerkeleyDB
        return BerkeleyDB.open(dbname, mode)
    elif line == "index\tflat/1":
        import FlatDB
        return FlatDB.open(dbname, mode)

    raise TypeError("Unknown index type: %r" % (line,))
    

def main():
    from Bio import Std
    import XPath
    import FlatDB
    XPath.xpath_index(
        #dbname = "sprot_flat",
        dbname = "sprot_small",
        filenames = ["/home/dalke/ftps/swissprot/smaller_sprot38.dat",
        #filenames = ["/home/dalke/ftps/swissprot/sprot38.dat",
                     ],
        primary_namespace = "entry",
        extract_info = [
        ("entry", "//entry_name"),
        ("accession", "//%s[@type='accession']" % (Std.dbid.tag,)),
        ],
        #creator_factory = FlatDB.CreateFlatDB,
        )


if __name__ == "__main__":
    main()
