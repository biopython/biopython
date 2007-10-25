import warnings
warnings.warn("Bio.expressions was deprecated, as it does not work with recent versions of mxTextTools. If you want to continue to use this module, please get in contact with the Biopython developers at biopython-dev@biopython.org to avoid permanent removal of this module from Biopython", DeprecationWarning)


from Martel import *
from Martel import Time
from Bio import Std
import ncbiblast

wublast_version = Re("2\.[^-]+-WashU")

# BLASTN 2.0a19MP-WashU [05-Feb-1998] [Build decunix3.2 01:53:29 05-Feb-1998]
# BLASTP 2.0a19MP-WashU [05-Feb-1998] [Build sol2.5-ultra 01:47:24 05-Feb-1998]
blastn_appheader = replace_groups(ncbiblast.blastn_appheader,
                                  [("appversion", wublast_version)])
blastp_appheader = replace_groups(ncbiblast.blastp_appheader,
                                  [("appversion", wublast_version)])
blastx_appheader = replace_groups(ncbiblast.blastx_appheader,
                                  [("appversion", wublast_version)])
tblastn_appheader = replace_groups(ncbiblast.tblastn_appheader,
                                  [("appversion", wublast_version)])
tblastx_appheader = replace_groups(ncbiblast.tblastx_appheader,
                                   [("appversion", wublast_version)])

spaces_line = ncbiblast.spaces_line
Number = ncbiblast.Number

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

query_database = (Str("Database:") + Opt(Spaces()) +
                  Std.database_name(UntilEol()) + AnyEol() +
                  Spaces() +
                  Std.database_num_sequences(Number(),
                                           {"bioformat:decode": "int.comma"}) +
                  Str(" sequences;") + Spaces() + 
                  Std.database_num_letters(Number(),
                                           {"bioformat:decode": "int.comma"}) +
                  Spaces() + Str("total letters.") + ToEol())
#                        notice the "."  -----^^^


#                                                                      Smallest
#                                                                        Sum
#                                                               High  Probability
# Sequences producing High-scoring Segment Pairs:              Score  P(N)      N
#  
# sp|P09429|HMG1_HUMAN HIGH MOBILITY GROUP PROTEIN HMG1 (HM...  1144  8.6e-117  1
# sp|P10103|HMG1_BOVIN HIGH MOBILITY GROUP PROTEIN HMG1 (HM...  1140  2.3e-116  1

to_table = (SkipLinesUntil(Str("Sequences producing High-scoring Segment")) +
            ToEol() +
            AnyEol())

table_entry = (AssertNot(AnyEol()) +
               Std.search_table_description(Re(r".{60}")) + Spaces() +
               Std.search_table_value(Float(), {"name": "score",
                                                "bioformat:decode": "float"}) +
               Spaces() +
               Std.search_table_value(Float(), {"name": "P",
                                                "bioformat:decode": "float"}) +
               Spaces() +
               Std.search_table_value(Digits(), {"name": "N",
                                              "bioformat:decode": "float"}) +
               AnyEol())
     
header = replace_groups(ncbiblast.header,
                        (("query_database", query_database),
                         ("to_table", to_table),
                         ("table_entry", table_entry)))

# >sp|P09429|HMG1_HUMAN HIGH MOBILITY GROUP PROTEIN HMG1 (HMG-1).
#             Length = 214
#  
#  Score = 1144 (402.7 bits), Expect = 8.6e-117, P = 8.6e-117
#  Identities = 214/214 (100%), Positives = 214/214 (100%)
#  
# Query:     1 GKGDPKKPRGKMSSYAFFVQTCREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFE 60
#              GKGDPKKPRGKMSSYAFFVQTCREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFE
# Sbjct:     1 GKGDPKKPRGKMSSYAFFVQTCREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFE 60

hit_descr = ncbiblast.hit_descr
hit_length = ncbiblast.hit_length


#  Score = 1144 (402.7 bits), Expect = 8.6e-117, P = 8.6e-117
#  Score = 79 (27.8 bits), Expect = 0.0088, Sum P(2) = 0.0088
#  Score = 65 (22.9 bits), Expect = 3.1, Sum P(3) = 0.96

blastp_score = (Re(r" *Score = *") +
                Std.hsp_value(Digits(), {"name": "bits",
                                         "bioformat:decode": "int",
                                         }) +
                # XXX Skipping the "402.7 bits" field, what's it called?
                ToSep(sep = "=") +
                Spaces() +
                Std.hsp_value(Float(), {"name": "expect",
                                         "bioformat:decode": "float",
                                         }) +
                # What should I do with the (2)?, or the (3)?
                (Str(", P =") | Re(r", Sum P\(\d+\) =")) + 
                Spaces() +
                Std.hsp_value(Float(), {"name": "P",
                                         "bioformat:decode": "float",
                                         }) +
                ToEol()
                )

#  Identities = 214/214 (100%), Positives = 214/214 (100%)
identities = (Re(r" *Identities = ") +
              Std.hsp_value(Digits(), {"name": "identical",
                                       "bioformat:decode": "int",
                                       }) +
              Str("/") +
              
              Std.hsp_value(Digits(), {"name": "length",
                                       "bioformat:decode": "int",
                                       }) +
              ToSep(sep = ",") +
              Str(" Positives =") +
              Spaces() +
              Std.hsp_value(Digits(), {"name": "positives",
                                       "bioformat:decode": "int",
                                       }) +
              ToEol())

blastp_hsp_header = blastp_score + identities + spaces_line


######################### blastn_hit


# >U51677 Human non-histone chromatin protein HMG1 (HMG1) gene, complete cds.
#         Length = 2575
#  
#   Plus Strand HSPs:
#  
#  Score = 12875 (1931.8 bits), Expect = 0.0, P = 0.0
#  Identities = 2575/2575 (100%), Positives = 2575/2575 (100%), Strand = Plus / Plus
#  
# Query:     1 ATGGGCAAAGGAGATCCTAAGAAGCCGAGAGGCAAAATGTCATCATATGCATTTTTTGTG 60
#              ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Sbjct:     1 ATGGGCAAAGGAGATCCTAAGAAGCCGAGAGGCAAAATGTCATCATATGCATTTTTTGTG 60
# ....

# Don't need to parse this since it comes from the identities line
blastn_strand =Re(r" *(Plus|Minus) Strand HSPs:\R")

blastn_score = blastp_score

blastn_identities = (
    Re(r" *Identities = ") +
    Std.hsp_value(Digits(), {"name": "identical",
                             "bioformat:decode": "int",
                             }) +
    Str("/") +
    
    Std.hsp_value(Digits(), {"name": "length",
                             "bioformat:decode": "int",
                             }) +
    ToSep(sep = ",") +
    Str(" Positives =") +
    Spaces() +
    Std.hsp_value(Digits(), {"name": "positives",
                             "bioformat:decode": "int",
                             }) +
    ToSep(sep = ",") + Str(" ") +
    ncbiblast.strand
    # that 'strand' already includes a newline
    )


blastn_hsp_header = (Opt(blastn_strand + spaces_line) +
                     blastn_score + blastn_identities + spaces_line)

######################### blastx

#                                                                      Smallest
#                                                                        Sum
#                                                      Reading  High  Probability
# Sequences producing High-scoring Segment Pairs:        Frame Score  P(N)      N
#  
# sp|P09429|HMG1_HUMAN HIGH MOBILITY GROUP PROTEIN HMG1 ... +2   311  5.0e-108  4
# sp|P10103|HMG1_BOVIN HIGH MOBILITY GROUP PROTEIN HMG1 ... +2   311  1.3e-107  4
# sp|P07155|HMG1_MOUSE HIGH MOBILITY GROUP PROTEIN HMG1 ... +2   311  2.7e-107  4                                        

blastx_table_entry = (
    AssertNot(AnyEol()) +
    Std.search_table_description(Re(r".{57}")) + Spaces() +
    Std.search_table_value(Integer(), {"name": "frame",
                                       "bioformat:decode": "int"}) +
    Spaces() +
    Std.search_table_value(Float(), {"name": "score",
                                     "bioformat:decode": "float"}) +
    Spaces() +
    Std.search_table_value(Float(), {"name": "P",
                                     "bioformat:decode": "float"}) +
    Spaces() +
    Std.search_table_value(Digits(), {"name": "N",
                                      "bioformat:decode": "int"}) +
    AnyEol())


blastx_to_database = (
 Str("  Translating both strands of query sequence in all 6 reading frames") +
 AnyEol() +
 AnyEol())

blastx_identities = (
    Re(r" *Identities = ") +
    Std.hsp_value(Digits(), {"name": "identical",
                             "bioformat:decode": "int",
                             }) +
    Str("/") +
    
    Std.hsp_value(Digits(), {"name": "length",
                             "bioformat:decode": "int",
                             }) +
    ToSep(sep = ",") +
    Str(" Positives =") +
    Spaces() +
    Std.hsp_value(Digits(), {"name": "positives",
                             "bioformat:decode": "int",
                             }) +
    ToSep(sep = "=") + Str(" ") +
    Std.hsp_frame(Integer(), {"which": "query"}) +
    AnyEol()
    )

blastx_score = blastp_score
blastx_strand = blastn_strand

blastx_hsp_header = (Opt(blastx_strand + spaces_line) +
                     blastx_score + blastx_identities + spaces_line)


#############################
######################### tblastn
tblastn_strand =Re(r" *(Plus|Minus) Strand HSPs:\R")

tblastn_score = blastp_score

tblastn_identities = (
    Re(r" *Identities = ") +
    Std.hsp_value(Digits(), {"name": "identical",
                             "bioformat:decode": "int",
                             }) +
    Str("/") +
    
    Std.hsp_value(Digits(), {"name": "length",
                             "bioformat:decode": "int",
                             }) +
    ToSep(sep = ",") +
    Str(" Positives =") +
    Spaces() +
    Std.hsp_value(Digits(), {"name": "positives",
                             "bioformat:decode": "int",
                             }) +
    ToSep(sep = ",") + 
    ncbiblast.frame
    # that 'strand' already includes a newline
    )


tblastn_hsp_header = (Opt(tblastn_strand + spaces_line) +
                     tblastn_score + tblastn_identities + spaces_line)



#############################
# May have to skip the WARNING section
ending_start = (Opt(Str("WARNING") + ToEol() +
                    Rep(AssertNot(Str(">")) +
                        AssertNot(Str("Parameters")) +
                        ToEol())) +
                Str("Parameters:") + ToEol())

parameters = SkipLinesUntil(Str("Statistics:"))

_nv = ncbiblast._nv
int_stat = ncbiblast.int_stat

timestamp = Time.make_expression(
    "%(Mon) %(Jan) %(day) %(hour):%(minute):%(second) %(YYYY)")

statistics = (
    Str("Statistics:") + AnyEol() +
    spaces_line +
        #  Database:  /dbEST1/wublast/database/swissprot/sprot.fa
    _nv("  Database: ",
        Std.search_statistic(UntilEol(), {"name": "database"})) +
    
        #    Title:  swissprot
    _nv("    Title: ",
        Std.search_statistic(UntilEol(), {"name": "database_title"})) +
    
        #    Release date:  unknown
    _nv("    Release date: ",
        Std.search_statistic(UntilEol(), {"name": "release_date"})) +
    
        #    Posted date:  11:37 AM GMT Feb 21, 2000
    _nv("    Posted date: ",
        Std.search_statistic(UntilEol(), {"name": "posted_date"})) +
    
        #    Format:  BLAST
    _nv("    Format: ",
        Std.search_statistic(UntilEol(), {"name": "format"})) +
    
        #  # of letters in database:  26,335,942
    _nv("  # of letters in database:  ",
        int_stat("num_letters_in_database", comma = 1)) +
    
        #  # of sequences in database:  72,615
    _nv("  # of sequences in database:  ",
        int_stat("num_sequences_in_database", comma = 1)) +

        #  # of database sequences satisfying E:  1119
    _nv("  # of database sequences satisfying E:  ",
        int_stat("num_sequences_satisying_E")) +

        #  No. of states in DFA:  549 (54 KB)
    _nv("  No. of states in DFA: ",
        Std.search_statistic(UntilEol(), {"name": "num_dfa_states"})) +
        
        #  Total size of DFA:  117 KB (128 KB)
    _nv("  Total size of DFA: ",
        Std.search_statistic(UntilEol(), {"name": "total_dfa_size"})) +

        #  Time to generate neighborhood:  0.01u 0.00s 0.01t  Elapsed: 00:00:00
    _nv("  Time to generate neighborhood: ",
        Std.search_statistic(UntilEol(),
                             {"name": "neighborhood_generation_time"})) +

        #  No. of threads or processors used:  2
    Opt(
    	_nv("  No. of threads or processors used:  ",
        	int_stat("num_threads"))) +

        #  Search cpu time:  38.04u 0.15s 38.19t  Elapsed: 00:00:33
    _nv("  Search cpu time: ",
        Std.search_statistic(UntilEol(), {"name": "search_cpu_time"})) +
        
        #  Total cpu time:  38.96u 0.29s 39.25t  Elapsed: 00:00:34
    _nv("  Total cpu time: ",
        Std.search_statistic(UntilEol(), {"name": "total_time"})) +

        #  Start:  Mon Feb 21 14:33:13 2000   End:  Mon Feb 21 14:33:47 2000
    _nv("  Start:  ",
        Std.search_statistic(timestamp, {"name": "start_time"}) +
        Spaces() + Str("End:  ") +
        Std.search_statistic(timestamp, {"name": "end_time"})) +
    Opt(spaces_line +
            #WARNINGS ISSUED:  2
        _nv("WARNINGS ISSUED:  ",
            int_stat("num_warnings")))
    )

ending = ending_start + parameters + statistics
blastp_ending = ending_start + parameters + statistics
blastn_ending = ending_start + parameters + statistics
blastx_ending = ending_start + parameters + statistics

format = replace_groups(ncbiblast.format,
                        (("to_table", to_table),
                         ("table_entry", table_entry),
                         ("ending", ending)))

blastp = replace_groups(format,
                        (("appheader", blastp_appheader),
                         ("hsp_header", blastp_hsp_header)))

blastn = replace_groups(format,
                        (("appheader", blastn_appheader),
                         ("hsp_header", blastn_hsp_header)))

blastx = replace_groups(format,
                        (("appheader", blastx_appheader),
                         ("table_entry", blastx_table_entry),
                         ("to_database", blastx_to_database),
                         ("hsp_header", blastx_hsp_header)))
tblastn = replace_groups(format,
                        (("appheader", tblastn_appheader),
                         ("table_entry", blastx_table_entry),
                         ("hsp_header", tblastn_hsp_header)))
			 

def main(outf):
    from xml.sax import saxutils
    parser = blastn.make_parser(debug_level = 0)
    parser.setContentHandler(saxutils.XMLGenerator())
    parser.parse(outf)
    #parser.parse("wu-sh_blastp.out")
    
if __name__ == "__main__":
    main()
