import os, sys

_open = open  # rename for internal use -- gets redefined below

def create_berkeleydb(dbname, unique, data_fields, format = "sequence"):
    import BerkeleyDB
    return BerkeleyDB.CreateBerkeleyDB(dbname, unique, data_fields, format)

def create_flatdb():
    import FlatDB
    return FlatDB.CreateFlatDB(dbname, unique, data_fields, format)

def open(dbname):
    line = _open(os.path.join(dbname, "BIOINDEX.dat")).readline()
    if line == "index\tBerkeleyDB/1\n":
        import BerkeleyDB
        return BerkeleyDB.open(dbname)
    elif line == "index\tflat/1\n":
        import FlatDB
        return FlatDB.open(dbname)

    raise TypeError("Unknown index type: %r" % (dbname,))
    

def main():
    from Bio import Std
    import XPath
    import FlatDB
    XPath.xpath_index(
        dbname = "sprot",
        filenames = ["/home/dalke/ftps/swissprot/smaller_sprot38.dat",
                     ],
        unique = "entry",
        extract_info = [
        ("entry", "//entry_name"),
        ("accession", "//%s[@type='accession']" % (Std.dbid.tag,)),
        ],
        creator_class = FlatDB.CreateFlatDB,
        )


if __name__ == "__main__":
    main()
