import os, sys

_open = open  # rename for internal use -- gets redefined below

def create_berkeleydb(dbname, unique, data_fields, format = "sequence"):
    import BerkeleyDB
    return BerkeleyDB.CreateBerkeleyDB(dbname, unique, data_fields, format)

def create_flatdb():
    import FlatDB
    return FlatDB.CreateFlatDB(dbname, unique, data_fields, format)

def open(dbname):
    text = _open(os.path.join(dbname, "config.dat"), "rb").read()
    line = text.split("\n")[0]
    if line == "index\tBerkeleyDB/1":
        import BerkeleyDB
        return BerkeleyDB.open(dbname)
    elif line == "index\tflat/1":
        import FlatDB
        return FlatDB.open(dbname)

    raise TypeError("Unknown index type: %r" % (line,))
    

def main():
    from Bio import Std
    import XPath
    import FlatDB
    XPath.xpath_index(
        #dbname = "sprot_flat",
        dbname = "sprot_fullbdb",
        #filenames = ["/home/dalke/ftps/swissprot/smaller_sprot38.dat",
        filenames = ["/home/dalke/ftps/swissprot/sprot38.dat",
                     ],
        unique = "entry",
        extract_info = [
        ("entry", "//entry_name"),
        ("accession", "//%s[@type='accession']" % (Std.dbid.tag,)),
        ],
        #creator_class = FlatDB.CreateFlatDB,
        )


if __name__ == "__main__":
    main()
