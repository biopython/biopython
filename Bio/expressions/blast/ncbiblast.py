from Martel import *
from Bio import Std

spaces_line = Opt(Spaces()) + AnyEol()

blast_filter = Str("BLASTN", "BLASTP", "BLASTX", "TBLASTX", "TBLASTN")
ncbi_blastn_filter = Re(r"BLASTN 2.\d+\.\d+ ")
ncbi_blastp_filter = Re(r"BLASTP 2.\d+\.\d+ ")
ncbi_blastx_filter = Re(r"BLASTX 2.\d+\.\d+ ")
ncbi_tblastn_filter = Re(r"TBLASTN 2.\d+\.\d+ ")
ncbi_tblastx_filter = Re(r"TBLASTX 2.\d+\.\d+ ")


# BLASTN 2.0.11 [Jan-20-2000]
blastn_header = (Std.application_name(Str("BLASTN"), {"app": "blastn"}) + 
                 Str(" ") +
                 Std.application_version(UntilSep(sep=" ")) +
                 ToEol())

# BLASTP 2.0.11 [Jan-20-2000]
blastp_header = (Std.application_name(Str("BLASTP"), {"app": "blastp"}) + 
                Str(" ") +
                Std.application_version(UntilSep(sep=" ")) +
                ToEol())

# BLASTX 2.0.11 [Jan-20-2000]
blastx_header = (Std.application_name(Str("BLASTX"), {"app": "blastx"}) + 
                 Str(" ") +
                 Std.application_version(UntilSep(sep=" ")) +
                ToEol())

# TBLASTN 2.0.11 [Jan-20-2000]
tblastn_header = (Std.application_name(Str("TBLASTN"), {"app": "tblastn"}) + 
                  Str(" ") +
                  Std.application_version(UntilSep(sep=" ")) +
                  ToEol())

# TBLASTX 2.0.11 [Jan-20-2000]
tblastx_header = (Std.application_name(Str("TBLASTX"), {"app": "tblastx"}) + 
                  Str(" ") +
                  Std.application_version(UntilSep(sep=" ")) +
                  ToEol())

# This is tricky because the description can include a statement like
#  '(492 letters)' at just the right place to confuse things.
query_letters = (Spaces() +
                 Str("(") +
                 Std.query_size(Digits()) +
                 Str(" letters)") +
                 AnyEol() +
                 Assert(AnyEol())
                 )

query_descr = (Str("Query=") +
               Opt(Spaces()) +
               Std.query_description(UntilEol()) +
               AnyEol() +
               Rep(AssertNot(query_letters) +
                   Std.query_description(UntilEol()) +
                   AnyEol())
               )

query_descr = Std.query_description_block(query_descr)

# This requires there to be at least one comma
_comma_digits = Re(r"\d(\d\d?)?(,\d\d\d)+")
def Number(name = None, attrs = {}):
    if name is None:
        assert not attrs
        return _comma_digits | Digits()
    return Group(name, _comma_digits, attrs) | \
           Digits(name, attrs)


query_database = (Str("Database:") + Opt(Spaces()) +
                  Std.database_name(UntilEol()) + AnyEol() +
                  Spaces() +
                  Std.database_num_sequences(Number(),
                                             {"bioformat:encoding": "comma"}) +
                  Str(" sequences;") + Spaces() + 
                  Std.database_num_letters(Number(),
                                           {"bioformat:encoding": "comma"}) +
                  Spaces() + Str("total letters") + ToEol())


#                                                                    Score     E
# Sequences producing significant alignments:                        (bits)  Value
# 
# U51677 Human non-histone chromatin protein HMG1 (HMG1) gene, co...  4129  0.0
# L38477 Mus musculus (clone Clebp-1) high mobility group 1 prote...   353  7e-95

to_table = (SkipLinesUntil(Str("Sequences producing significant")) +
            ToEol() +
            AnyEol())

table_entry = Std.search_table_entry(
    (AssertNot(AnyEol()) +
     Std.search_table_description(Re(r".{66}")) + Spaces() +
     Std.search_table_value(Float(), {"name": "score",
                                      "bioformat:type": "float"}) +
     Spaces() +
     Std.search_table_value(Float(), {"name": "evalue",
                                      "bioformat:type": "float"}) +
     AnyEol()))

table = Std.search_table(Rep(table_entry))


header = Std.search_header(
    SkipLinesUntil(Str("Query=")) +
    query_descr +
    query_letters +
    AnyEol() +
    query_database +
    Opt(to_table + table) + 
    SkipLinesUntil(Str(">"))
    )


# >U51677 Human non-histone chromatin protein HMG1 (HMG1) gene, complete
#             cds.
#             Length = 2575
#  
#  Score = 4129 bits (2083), Expect = 0.0
#  Identities = 2167/2209 (98%)
#  Strand = Plus / Plus

# This is tricky since it's (theoretical) possible the description
# could have a "Length =" at exactly the right spot to confuse things.
# So stop before getting to the line with a 'Length =' where the
# line afterwards is blank.
hit_length = (Spaces() + Str("Length = ") +
              Std.hit_length(Digits()) + AnyEol() +
              Opt(Spaces()) + AnyEol())

hit_descr = (Str(">") + Std.hit_description(UntilEol()) + AnyEol() +
             Rep(AssertNot(hit_length) +
                 Std.hit_description(UntilEol()) + AnyEol())
             )


num_bits = Std.hsp_value(Float(), {"name": "bits",
                                   "bioformat:type": "float",
                                   })

expect = Std.hsp_value(Float(), {"name": "expect",
                                 "bioformat:type": "float"})

num_identical = Std.hsp_value(Digits(), {"name": "identical",
                                          "bioformat:type": "int"})

hsp_length = Std.hsp_value(Digits(), {"name": "length",
                                      "bioformat:type": "int"})

num_positives = Std.hsp_value(Digits(), {"name": "positives",
                                         "bioformat:type": "int"})

# What's after the bits?
score = (Re(r" *Score = *") + num_bits +
         Re(r" bits \((?P<XXX1>\d+)\), Expect = *") +
         expect + AnyEol())

# XXX What's the expect_num mean?
tblastn_score = (Re(r" *Score = *") + num_bits +
                 Re(r" bits \((?P<XXX1>\d+)\), "
                    r"Expect(\((?P<expect_num>\d+)\))? = *") +
                 expect + AnyEol())
    
blastn_identical = (Re(r" *Identities = *") +
                    num_identical +
                    Str("/") +
                    hsp_length +
                    ToEol())

#  Identities = 154/168 (91%), Positives = 154/168 (91%)
blastp_identical = (Re(r" *Identities = *") +
                     num_identical +
                     Str("/") +
                     hsp_length +
                     ToSep(sep = ",") +
                     Re(" Positives = *") +
                     num_positives +
                     ToEol())

#  Frame = +1
# can be +1, +2, +3, -1, -2, -3
frame = (Str(" Frame = ") + Std.hsp_frame(UntilEol(), {"which": "query"}) +
         AnyEol())
tblastx_frame = (Str(" Frame = ") +
                 Std.hsp_frame(UntilSep(sep = " "), {"which": "query"}) + 
                 Str(" / ") +
                 Std.hsp_frame(UntilEol(), {"which": "subject"}) +
                 AnyEol())

# Strand = Plus / Plus
strand = (Re(r" *Strand = ") +
          (Std.hsp_strand(Str("Plus"), {"strand": "+1", "which": "query"}) | 
           Std.hsp_strand(Str("Minus"), {"strand": "-1", "which": "query"})) +
          Str(" / ") +
          (Std.hsp_strand(Str("Plus"), {"strand": "+1", "which": "subject"}) | 
           Std.hsp_strand(Str("Minus"), {"strand": "-1", "which": "subject"}))+
          AnyEol())

query_line = Std.hsp_seqalign_line(
    Std.hsp_seqalign_leader(
      Std.hsp_seqalign_name(Str("Query")) +
      Re(": *") +
      Std.hsp_seqalign_start(Digits()) +
      Re(" *")
    ) +
    Std.hsp_seqalign_seq(UntilSep(sep = " ")) +
    Re(" *") +
    Std.hsp_seqalign_end(Digits()) +
    AnyEol(),
    {"which": "query"})

homology_line = Std.hsp_seqalign_line(
    Std.hsp_seqalign_seq(Str("      ") +
                         UntilEol()) +
    AnyEol(),
    {"which": "homology"})

subject_line = Std.hsp_seqalign_line(
    Std.hsp_seqalign_leader(
      Std.hsp_seqalign_name(Str("Sbjct")) +
      Re(": *") +
      Std.hsp_seqalign_start(Digits()) +
      Re(" *")
    ) +
    Std.hsp_seqalign_seq(UntilSep(sep = " ")) +
    Re(" *") +
    Std.hsp_seqalign_end(Digits()) +
    AnyEol(),
    {"which": "subject"})

alignment = Std.hsp_seqalign(query_line + homology_line + subject_line +
                             Rep1(spaces_line))

blastn_hsp = Std.hsp(score + blastn_identical + strand +
                     spaces_line + spaces_line +
                     Rep(alignment))
blastp_hsp = Std.hsp(score + blastp_identical + spaces_line +
                     Rep(alignment))
blastx_hsp = Std.hsp(score + blastp_identical + frame + spaces_line +
                     Rep(alignment))
tblastn_hsp = Std.hsp(tblastn_score + blastp_identical + frame + spaces_line +
                      Rep(alignment))
tblastx_hsp = Std.hsp(tblastn_score + blastp_identical + tblastx_frame +
                      spaces_line + spaces_line +
                      Rep(alignment))

blastn_hit = Std.hit(hit_descr + hit_length + Rep(blastn_hsp))
blastp_hit = Std.hit(hit_descr + hit_length + Rep(blastp_hsp))
blastx_hit = Std.hit(hit_descr + hit_length + Rep(blastx_hsp))
tblastn_hit = Std.hit(hit_descr + hit_length + Rep(tblastn_hsp))
tblastx_hit = Std.hit(hit_descr + hit_length + Rep(tblastx_hsp))


# Already know the database

# parse the name: value lines
def _nv(s, expr):
    return Str(s) + expr + AnyEol()
def float_stat(name):
    return Std.search_statistic(Float(), {"name": name,
                                          "bioformat:type": "float"})
def _num_stat(name, signed, comma):
    d = {"name": name,
         "bioformat:type": "int"}
    if signed:
        exp = Integer()
    elif comma:
        exp = Number()
        d["bioformat:encoding"] = "comma"
    else:
        exp = Digits()
    return exp, d

def int_stat(name, signed = 0, comma = 0):
    exp, d = _num_stat(name, signed, comma)
    return Std.search_statistic(exp, d)

def int_param(name, signed = 0, comma = 0):
    exp, d = _num_stat(name, signed, comma)
    return Std.search_parameter(exp, d)

def float_param(name):
    return Std.search_parameter(Float(), {"name": name,
                                          "bioformat:type": "float"})


db_stats = (Str("  Database:") + ToEol() +
            Str("    Posted date:") + ToEol() +
            Str("  Number of letters") + ToEol() +
            Str("  Number of sequences") + ToEol() +
            spaces_line)

ungapped_lambda_stats = (
    Str("Lambda     K      H") + AnyEol() +
    ( Spaces() + float_stat("lambda") + 
      Spaces() + float_stat("k") + 
      Spaces() + float_stat("h") +
      ToEol() )
)
gapped_lambda_stats = (
    Str("Gapped") + AnyEol() + 
    Str("Lambda     K      H") + AnyEol() + 
    ( Spaces() + float_stat("gapped_lambda") + 
      Spaces() + float_stat("gapped_k") + 
      Spaces() + float_stat("gapped_h") +
      ToEol() )
    )
lambda_stats = (ungapped_lambda_stats +
                AnyEol() +
                gapped_lambda_stats)

matrix_stats = (
    SkipLinesUntil(Str("Matrix:")) +
    # Matrix: blastn matrix:1 -3
    _nv("Matrix: ", Std.search_statistic(UntilEol(), {"name": "matrix"}))
    )

gap_penalties_stats = (
    # Gap Penalties: Existence: 5, Extension: 2
    _nv("Gap Penalties: Existence: ",
        int_param("existence_gap_penalty") + 
        Str(", Extension: ") +
        int_param("extension_gap_penalty"))
    )
    
generic_info1 = (
    # Number of Hits to DB: 1123218
    # !! I've seen negative numbers here !!
    # Number of Hits to DB: -331804600
    _nv("Number of Hits to DB: ", int_stat("num_hits_to_db", signed = 1)) +

    # Number of Sequences: 442729
    _nv("Number of Sequences: ", int_stat("num_sequences")) +

    # Number of extensions: 1123218
    _nv("Number of extensions: ", int_stat("num_extensions")) +

    # Number of successful extensions: 86816
    _nv("Number of successful extensions: ",
        int_stat("num_successful_extensions")) +

    # Number of sequences better than 10.0: 106
    _nv("Number of sequences better than 10.0: ",
        int_stat("num_sequences_better_than_10"))
    )

generic_info2 = (
    # length of query: 2575
    _nv("length of query: ", int_stat("query_length")) +

    # length of database: 675,252,082
    _nv("length of database: ", int_stat("database_length", comma = 1)) +
          
    # effective HSP length: 21
    _nv("effective HSP length: ", int_stat("effective_hsp_length"))+

    # effective length of query: 2554
    _nv("effective length of query: ", int_stat("effective_query_length")) +
                                                   
    # effective length of database: 665,954,773
    _nv("effective length of database: ",
        int_stat("effective_database_length", comma = 1)) +
    
    # effective search space: 1700848490242
    _nv("effective search space: ", int_stat("effective_search_space")) +
          
    # effective search space used: 1700848490242
    _nv("effective search space used: ",
        int_stat("effective_search_space_used"))
    )

hsp_info = (
    # Number of HSP's better than 10.0 without gapping: 136
    _nv("Number of HSP's better than 10.0 without gapping: ",
        int_stat("num_hsps_better_than_10_without_gapping")) +

    # Number of HSP's successfully gapped in prelim test: 14
    _nv("Number of HSP's successfully gapped in prelim test: ",
        int_stat("num_hsps_successfully_gapped_in_prelim_test")) +

    # Number of HSP's that attempted gapping in prelim test: 880
    _nv("Number of HSP's that attempted gapping in prelim test: ",
        int_stat("num_hsps_attempted_gapping_in_prelim_test")) +

    # Number of HSP's gapped (non-prelim): 291
    _nv("Number of HSP's gapped (non-prelim): ",
        int_stat("num_hsps_gapped"))
    )

# frameshift window, decay const: 50,  0.1
frameshift_info = _nv("frameshift window, decay const: ",
                      int_stat("frameshift_window") +
                      Re(r", *") +
                      float_stat("frameshift_decat_const"))

def _bit_info(s, name):
    return _nv(s,
               int_stat(name) + 
               Spaces() + Str("(") + Opt(Spaces()) +
               float_stat(name + "_bits") +
               Str(" bits)"))

            # T: 0
t_info = _nv("T: ", int_stat("t"))

            # A: 0
a_info = _nv("A: ", int_stat("1"))

          # X1: 6 (11.9 bits)
x1_info = _bit_info("X1: ", "x1")
          
          # X2: 10 (19.8 bits)
x2_info = _bit_info("X2: ", "x2")

          # X3: 64 (24.9 bits)
x3_info = _bit_info("X3: ", "x3")

          # S1: 12 (24.3 bits)
s1_info = _bit_info("S1: ", "s1")

          # S2: 19 (38.2 bits)
s2_info = _bit_info("S2: ", "s2")
    

blastn_ending = (db_stats +
                 lambda_stats + 
                 matrix_stats +
                 gap_penalties_stats +
                 generic_info1 +
                 generic_info2 +
                 t_info +
                 a_info +
                 x1_info +
                 x2_info +
                 s1_info +
                 s2_info
                 )

blastp_ending = (db_stats +
                 lambda_stats + 
                 matrix_stats +
                 gap_penalties_stats +
                 generic_info1 +
                 hsp_info +
                 generic_info2 + 
                 t_info +
                 a_info +
                 x1_info +
                 x2_info +
                 x3_info +
                 s1_info +
                 s2_info
                 )

blastx_ending = (db_stats +
                 lambda_stats + 
                 matrix_stats +
                 gap_penalties_stats +
                 generic_info1 +
                 hsp_info +
                 generic_info2 +
                 frameshift_info +
                 t_info +
                 a_info +
                 x1_info +
                 x2_info +
                 x3_info +
                 s1_info +
                 s2_info
                 )
tblastn_ending = blastx_ending
tblastx_ending = (db_stats +
                 ungapped_lambda_stats + 
                 matrix_stats +
                 generic_info1 +
                 generic_info2 +
                 frameshift_info +
                 t_info +
                 a_info +
                 x1_info +
                 x2_info +
                 s1_info +
                 s2_info
                 )

blastn = Std.record(blastn_header + header + Rep(blastn_hit) + blastn_ending)
blastp = Std.record(blastp_header + header + Rep(blastp_hit) + blastp_ending)
blastx = Std.record(blastx_header + header + Rep(blastx_hit) + blastx_ending)
tblastn = Std.record(tblastn_header + header + Rep(tblastn_hit) + tblastn_ending)
tblastx = Std.record(tblastx_header + header + Rep(tblastx_hit) + tblastx_ending)
