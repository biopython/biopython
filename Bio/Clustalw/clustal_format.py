"""clustal_format.py

A parser to read information from a clustal formatted file (*.aln).

This uses Andrew Dalke's Martel to do the parsing dirty work for me.
So all we need to do here is set up a big ol' regular expression to
let Martel know what the file looks like."""
# standard library
import sys

# Martel stuff
try:
    import Martel
except ImportError:
    print "Ooops, Clustalw parsing requires Martel, available from:"
    print "http://www.biopython.org/~dalke/Martel/"
    sys.exit(0)
    

# define everything we will parse at a ton of regular expressions with
# specific callbacks
version = Martel.Group("version",
                       Martel.Re("\d.\d\d?"))

header = Martel.Group("header",
                     Martel.Str("CLUSTAL ") +
                     Martel.Re(".+") +
                     Martel.MaxRepeat(Martel.AnyEol(), 0, 3))

seq_id = Martel.Group("seq_id",
                      Martel.Re("[-a-zA-Z:;^_'\",\+\#\|\[\]\(\)\/\.\d]+"))

# space between the sequence and id
seq_space = Martel.Group("seq_space",
                         Martel.Re("[ ]+"))

seq_info = Martel.Group("seq_info",
                        Martel.Re("[-a-zA-Z.]+"))

# you can output an optional number to tell you where you are in the sequence
# we need to swallow this up if it is here
seq_num = Martel.Group("seq_num",
                       Martel.Re("[ ]+") +
                       Martel.Re("[\d]+"))

seq_line = Martel.Group("seq_line", seq_id + seq_space + seq_info +
                        Martel.Opt(seq_num) +
                        Martel.Str("\n"))

match_stars = Martel.Group("match_stars",
                           Martel.Re("[ :\.\*]+") +
                           Martel.Opt(Martel.AnyEol()))

# separator between blocks
new_block = Martel.Group("new_block",
                         Martel.AnyEol())

block_info = Martel.Group("block_info",
                          Martel.Rep1(seq_line) +
                          Martel.Opt(match_stars) +
                          Martel.Rep(new_block))


# define the format we can import to parse clustal files, one header
# plus multiple lines of alignments
format = Martel.Group("clustalx",
                      header +
                      Martel.Rep1(block_info))
