"""Parser for the SWISS-PROT 38 format.

You probably want to use the variables 'record' (for a single record)
and 'format' (for a set of records).

"""
import Martel
from Martel import RecordReader

def Simple(tag, tag_data):
    return Martel.Group(tag, Martel.Str(tag + "   ") +
                        Martel.Group(tag_data, Martel.Re("[^\R]*")) +
                        Martel.AnyEol()
                        )
#--- ID

ID = Martel.Group("ID", Martel.Re(
    r"ID   (?P<entry_name>\w+) +(?P<data_class_table>\w+); +" \
    r"(?P<molecule_type>\w+); +(?P<sequence_length>\d+) AA\.\R"
    ))

#--- AC

AC = Martel.Group("AC", Martel.Re(
    r"AC   (?P<ac_number>\w+);( (?P<ac_number>\w+);)*\R"
    ))
AC_block = Martel.Group("AC_block", Martel.Rep1(AC))


#--- DT

DT_created = Martel.Group("DT_created", Martel.Re(
    r"DT   (?P<day>\d\d)-(?P<month>...)-(?P<year>\d{4}) \(Rel. "\
    r"(?P<release>\d\d), Created\)\R"
    ))
DT_seq_update = Martel.Group("DT_seq_update", Martel.Re(
    r"DT   (?P<day>\d\d)-(?P<month>...)-(?P<year>\d{4}) \(Rel. "\
    r"(?P<release>\d\d), Last sequence update\)\R"
    ))
DT_ann_update = Martel.Group("DT_ann_update", Martel.Re(
    r"DT   (?P<day>\d\d)-(?P<month>...)-(?P<year>\d{4}) \(Rel. "\
    r"(?P<release>\d\d), Last annotation update\)\R"
    ))

#--- DE

# Only the last DE is supposed to have a ".", but I don't see why *I*
# need to enforce it.
DE = Simple("DE", "description")
DE_block = Martel.Group("DE_block", Martel.Rep1(DE))


#--- GN

GN = Simple("GN", "gene_names")
GN_block = Martel.Group("GN_block", Martel.Rep1(GN))

#--- OS

OS = Simple("OS", "organism_species")
OS_block = Martel.Group("OS_block", Martel.Rep1(OS))



#--- OG

OG = Simple("OG", "organelle")
OG_block = Martel.Group("OG_block", Martel.Rep1(OG))


#--- OC

OC = Simple("OC", "organism_classification")
OC_block = Martel.Group("OC_block", Martel.Rep1(OC))

############ Reference section

#--- RN

# occurs once
RN = Martel.Group("RN", Martel.Re("RN   \[(?P<reference_number>\d+)]\R"))

#--- RP

# occurs once
RP = Simple("RP",  "reference_position")


#--- RC

# 0 or more
RC = Simple("RC",  "reference_comment")
RC_block = Martel.Group("RC_block", Martel.Rep1(RC))

#--- RX

# 0 or 1
RX = Martel.Group("RX",
                  Martel.Re("RX   (?P<bibliographic_database_name>\w+); " \
                            "(?P<bibliographic_identifier>\d+)\.\R"))

#--- RA

# 1 or more
RA = Simple("RA",  "reference_author")
RA_block = Martel.Group("RA_block", Martel.Rep1(RA))


#--- RT

# 0 or more
RT = Simple("RT",  "reference_title")
RT_block = Martel.Group("RT_block", Martel.Rep1(RT))


#--- RL

# 1 or more

RL = Simple("RL",  "reference_location")
RL_block = Martel.Group("RL_block", Martel.Rep1(RL))

reference = Martel.Group("reference",
                         RN + \
                         RP + \
                         Martel.Opt(RC_block) + \
                         Martel.Opt(RX) + \
                         RA_block + \
                         Martel.Opt(RT_block) + \
                         RL_block
                         )


############

#--- CC

CC_begin = Martel.Group("CC",
                        Martel.Re("CC   -!- ") + \
                        Martel.ToEol("comment_text"))
CC =       Martel.Group("CC",
                        Martel.Re("CC       ") + \
                        Martel.ToEol("comment_text"))

single_comment = Martel.Group("comment", 
                              CC_begin +
                              Martel.Rep(CC)
                              )


CC_copyright_begin = Martel.Group("CC_copyright_begin",
                                  Martel.Re("CC   -+\R"))
CC_copyright = Martel.Group("CC_copyright",
                            Martel.Re("CC   (?!-+\R)") + \
                            Martel.ToEol("copyright"))
CC_copyright_end =   Martel.Group("CC_copyright_end",
                                  Martel.Re("CC   -+\R"))

# From N33_HUMAN
bogus_DR_group = Martel.Group("bogus_DR_block",
                       Martel.Re(r"(?P<DR>DR   (?P<database_identifier>MIM); " \
                                 r"(?P<primary_identifier>601385); " \
                                 r"(?P<secondary_identifier>-).\R)")
                       )


comment = Martel.Group("comment_block",
                Martel.Rep(single_comment) + \
                Martel.Opt(bogus_DR_group) + \
                Martel.Opt(CC_copyright_begin + \
                           Martel.Rep(CC_copyright) + \
                           CC_copyright_end \
                           )
                       )

#--- DR

# The ([^.\R]|(?!.\R)\.)+) is because of things like
#   DR   MGD; MGI:95401; EPB4.1.
# where I need to scan up to the last "."

DR_general = Martel.Re("(?P<database_identifier>[^;]+);" \
                       "(?P<primary_identifier>[^;]+); " \
                       "(?P<secondary_identifier>([^.\R]|(?!.\R)\.)+)")

DR_prosite = Martel.Re("(?P<database_identifier>(PROSITE|PFAM)); " \
                       "(?P<primary_identifier>[^;]+); " \
                       "(?P<secondary_identifier>[^;]+); " \
                       "(?P<status_identifier>[^.]+)")

DR_embl = Martel.Re("(?P<database_identifier>EMBL); " \
                    "(?P<primary_identifier>[^;]+); " \
                    "(?P<secondary_identifier>[^;]+); " \
                    "(?P<status_identifier>[^.]+)")

DR = Martel.Group("DR", Martel.Str("DR   ") + \
           Martel.Group("database_reference",
                        Martel.Alt(DR_embl, DR_prosite, DR_general)) + \
           Martel.Str(".") + Martel.AnyEol())

DR_block = Martel.Group("DR_block", Martel.Rep1(DR))



#--- KW

KW = Simple("KW", "keyword")
KW_block = Martel.Group("KW_block", Martel.Rep1(KW))


#--- FT

# FT   DOMAIN       77     88       ASP/GLU-RICH (ACIDIC).
# 123456789012345678901234567890123456789012345678901234567890123456789012345
#          1         2         3         4         5         6         7
# FT   ........ ...... ......       .........................................
# FT   12345678 123456 123456       12345678901234567890123456789012345678901
# FT   .{8}     .{6}   .{6}         [^\R]*
#              1      1      1234567

# "FT   " + ".{8}" + " " + ".{6}" + " " + ".{6}" + "       " + "[^\R]*" + "\R"
# "FT   .{8} .{6} .{6}       [^\R]*\R"

FT_range = Martel.Group("FT",
                        Martel.Re("FT   (?P<ft_name>.{8}) " \
                                  "(?P<ft_from>.{6}) (?P<ft_to>.{6})" \
                                  "(       (?P<ft_description>[^\R]*))?\R")
                        )
FT_continuation = Martel.Group("FT_continuation",
                        Martel.Re("FT                                " \
                                  "(?P<ft_description>[^\R]*)\R")
                        )
FT = Martel.Group("feature", FT_range + Martel.Rep(FT_continuation))

feature_block = Martel.Group("feature_block", Martel.Rep1(FT))


#--- SQ

# SQ   SEQUENCE  XXXX AA; XXXXX MW;  XXXXX CRC32;
# (Those X's don't really indicate the size)

SQ = Martel.Group("SQ",
   Martel.Re("SQ   SEQUENCE +(?P<sequence_length>\d+) AA;" \
             " +(?P<molecular_weight>\d+) MW;" \
             " +(?P<crc32>\w+) CRC32;\R")
                  )
SQ_data = Martel.Group("SQ_data",
                       Martel.Re("     (?P<sequence>[^\R]*)\R"))

sequence = Martel.Group("sequence_block", Martel.Group("SQ_data_block",
                                                 SQ + Martel.Rep(SQ_data)))

#--- //

end = Martel.Group("END", Martel.Str("//") + Martel.AnyEol())

####################### put it all together

record = Martel.Group("swissprot38_record", \
    ID + \
    AC + \
    DT_created + \
    DT_seq_update + \
    DT_ann_update + \
    Martel.Opt(DE_block) + \
    Martel.Opt(GN_block) + \
    Martel.Opt(OS_block) + \
    Martel.Opt(OG_block) + \
    Martel.Opt(OC_block) + \
    Martel.Group("reference_block", Martel.Rep(reference)) + \
    comment + \
    Martel.Opt(DR_block) + \
    Martel.Opt(KW_block) + \
    Martel.Opt(feature_block) + \
    sequence + \
    end
    )

format_expression = Martel.Group("swissprot38", Martel.Rep1(record))

format = Martel.ParseRecords("swissprot38", {}, record,
                             RecordReader.EndsWith, ("//\n",) )
