from Martel import *

# BLASTN 2.0a19MP-WashU [05-Feb-1998] [Build decunix3.2 01:53:29 05-Feb-1998]

blastn_prog = Group("prog_name", Str("BLASTN")) + \
              Str(" ") + \
              ToSep("prog_version", " ") + \
              ToEol("prog_descr")

skip_to_query_info = SkipLinesUntil(Str("Query="))

# This is tricky because the description can include a statement like
#  '(492 letters)' at just the right place to confuse things.

# Query= YFL037W YFL037W, Chr VI from 56336 to 57709
#          (1374 letters)

# Query=  YAL001C YAL001C, Chr I from 147596 to 147665, and 147756 to 151168,
#     reverse complement
#         (3483 letters)

# Query= gi|2959129|gb|AA864816|AA864816 oh38a08.s1 NCI_CGAP_Kid6 Homo
# sapiens cDNA clone IMAGE:1460054 3' similar to SW:MY16_MOUSE P17564
# MYELOID DIFFERENTIATION PRIMARY RESPONSE PROTEIN MYD116. ;.
#          (492 letters)
query_info = Str("Query=") + Opt(Spaces()) + ToEol("query_descr") + \
             Rep(AssertNot(Spaces() + Str("(") + Digits() + \
                           Str(" letters)") + AnyEol() + AnyEol()) + \
                 Opt(Spaces()) + ToEol("query_descr")) + \
             Spaces() + Str("(") + Digits("query_size") + Str(" letters)") + \
             AnyEol()


# This requires there to be at least one comma
_comma_digits = Re(r"\d(\d\d?)?(,\d\d\d)+")
def Number(name = None):
    if name is None:
        return _comma_digits | Digits()
    return Group(name, _comma_digits, {"type": "comma"}) | \
           Group(name, Digits(), {"type": "digits"})

# Database:  YeastORF-P
#            6183 sequences; 2,900,438 total letters.

# Database:  Non-redundant GenBank CDS translations+PDB+SwissProt+SPupdate+PIR
#            307,320 sequences; 92,696,426 total letters.

# Database: mgpep
#            468 sequences; 170,400 total letters

# Database:  plantPept
#            6418 sequences; 2,370,771 total letters.

# Database: Non-redundant GenBank CDS
# translations+PDB+SwissProt+SPupdate+PIR
#            301,269 sequences; 90,873,415 total letters


query_database = Str("Database:  ") + ToEol("query_database") + \
                 Spaces() + Number("num_sequences") + \
                 Str(" sequences;") + Spaces() + \
                 Number("total_letters") + \
                 Spaces() + Str("total letters.") + AnyEol()

# Skip lines up until the report
blah1 = SkipLinesUntil(Str(">"))

# This is tricky since it's (theoretical) possible the description
# could have a "Length =" at exactly the right spot to confuse things.
# So stop before getting to the line with a 'Length =' where the
# line afterwards is blank.
hit_descr = Str(">") + ToEol("hit_descr") + \
            Rep(AssertNot(Spaces() + Str("Length = ") + Number() + \
                          AnyEol() + AnyEol()) + \
            Spaces() + ToEol("hit_descr"))
hit_length = Spaces() + Str("Length = ") + Number("hit_length") + \
             AnyEol() + AnyEol()

##hit_info_1 = Str(" Score = ") + Integer("hit_score") + Str(" ( ") + \
##             Float("hit_num_bits") + Str(" bits, Expect = ") + \
##             Float("hit_expect") + Str(", Sum P(2) = ") + Float("hit_P2") + \
##             AnyEol()

hit_info_1 = Str(" Score = ") + Digits("hit_score") + Str(" (") + \
             Opt(Spaces()) + Float("hit_num_bits") + \
             Str(" bits), Expect = ") + \
             Float("hit_expect") + Str(", Sum P(4) = ") + Float("hit_P4") + \
             AnyEol()

hit_info_2 = Str(" Identities = ") + Integer("hit_identities") + \
             Str("/") + Integer("hit_total_size") + \
             ToSep(sep = ",") + Str(" Positives = ") + \
             Integer("hit_positives") + Str("/") + \
             ToSep(sep = "=") + Spaces() + \
             Word("hit_query_strand") + Str(" / ") + \
             Word("hit_target_strand") + \
             ToEol()
            
query_line = Str("Query:") + Spaces() + Digits("hsp_query_start") + \
             Spaces() + ToSep("hsp_query_seq", " ") + \
             Digits("hsp_query_end") + AnyEol()

query_line = Str("Query:") + Spaces() + Digits("hsp_query_start") + \
             Spaces() + ToSep("hsp_query_seq", " ") + \
             Digits("hsp_query_end") + AnyEol()
alignment_line = Str("      ") + ToEol("hsp_alignment_seq")
sbjct_line = Str("Sbjct:") + Spaces() + Digits("hsp_subject_start") + \
             Spaces() + ToSep("hsp_subject_seq", " ") + \
             Digits("hsp_subject_end") + AnyEol()
