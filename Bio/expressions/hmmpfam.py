"""Martel expression for the hmmpfam database search program in hmmer.

This has been tested with version 2.2g. I should also make it work with
2.1.1 output.

XXX This isn't completely finished and doesn't do everything quite right.
The main two problems being that it isn't well tested for a wide variety of
outputs and that the family line is not parsed into it's respective parts 
(see multitude of comments on this below).
"""
from Martel import *
from Martel import RecordReader
from Bio import Std

# -- header
# hmmpfam - search one or more sequences against HMM database
program_description = (Std.application_name(Str("hmmpfam")) + 
                       ToEol())

# HMMER 2.2g (August 2001)
program_version = (Str("HMMER ") +
                   Std.application_version(Re(r"\d\.\d\w") |
                                           Re(r"\d\.\d\.\d")) +
                   ToEol())

# Copyright (C) 1992-2001 HHMI/Washington University School of Medicine
# Freely distributed under the GNU General Public License (GPL)

copyright = (ToEol() +
             ToEol())

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HMM file:                 /data/hmm/Pfam_fs
# Sequence file:            hmmpfam_test.fasta
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         
files = (ToEol() +
         Str("HMM file:") + Spaces() +
         Std.database_name(UntilEol()) + AnyEol() +
         Str("Sequence file:") + Spaces() +
         Group("inputfile_name", UntilEol()) + AnyEol() +
         ToEol())

header = Std.search_header(program_description + program_version +
                           copyright + files)

# -- record

#
# Query sequence: AT1G01040
# Accession:      [none]
# Description:    Test Sequence for hmmpfam

sequence_info = (ToEol() + # blank line
                 Str("Query sequence:") + Spaces() +
                 Group("query_name", UntilEol()) + AnyEol() +
                 Str("Accession:") + Spaces() +
                 Group("query_accession", UntilEol()) + AnyEol() +
                 Str("Description:") + Spaces() +
                 Std.query_description(UntilEol()) + AnyEol())

# generic name for models
# most model names are something like PFwhatever, but we also have some
# different kinds like AAA and Sigma54_activat
model_name = Re(r"[\w-]+")

#
# Scores for sequence family classification (score includes all domains):
# Model          Description                              Score    E-value  N 
# --------       -----------                              -----    ------- ---

family_header = (ToEol() + # blank line
                 Str("Scores") + ToEol() +
                 Str("Model") + ToEol() +
                 Str("-----") + ToEol())

# family lines can either be a hit:
# PF00636        RNase3 domain                            354.9   1.1e-118   2
# or no hit:
#         [no hits above thresholds]

no_hit_line = (Spaces() + Str("[no hits above thresholds]") + AnyEol())

family_hit_line = (Group("family_model", model_name) + Spaces() +
                   # XXX This is a pain in the ass, so I gave up
                   # We need to match descriptions but stop when we get to the
                   # score value.
                   # This is very difficult since the description can contain
                   # letters, numbers, (, ), /, - and periods. The numbers and
                   # periods leave me unable to think of a way to distinguish
                   # a "description word" from a score word.
                   # You also can't use spaces as descriptions can have 2
                   # spaces in them (not just one separating every word) and
                   # the description to score value can also be separated by
                   # only 2 spaces.
                   # I can't for the life of me figure this so I just
                   # eat the rest of the line for right now and downstream
                   # apps will have to get the information out from there.
                    ToEol("family_information"))
                   #Group("family_description",
                   #    Rep1(Re(r"[a-zA-Z0-9_()/-]+") +
                   #         Alt(Str("  "), Str(" ")))) +
                   # Alt(Spaces()) +
                   # Float("family_score") + Spaces() +
                   # Float("family_evalue") + Spaces() +
                   # Integer("family_n") + AnyEol())

#
# Parsed for domains:
# Model          Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
# --------       ------- ----- -----    ----- -----      -----  -------
domain_header = (ToEol() + # blank line
                 Str("Parsed for domains:") + AnyEol() +
                 Str("Model") + ToEol() +
                 Str("-----") + ToEol())

# domain lines can also be either a hit:
# PF00270          1/1     239   382 ..     1   141 [.    30.3  2.2e-09
# or not hit:
#         [no hits above thresholds]

# I have no idea what these things exactly mean, but parse 'em anyways
symbol_forward = Str(".") | Str("[")
symbol_reverse = Str(".") | Str("]")
match_symbols = symbol_forward + symbol_reverse

domain_hit_line = (Group("domain_model", model_name) + Spaces() +
                   Group("domain_domain", Integer() + Str("/") + Integer()) +
                   Spaces() +
                   Integer("domain_seq-f") + Spaces() +
                   Integer("domain_seq-r") + Spaces() +
                   Group("domain_seq_symbols", match_symbols) + Spaces() +
                   Integer("domain_hmm-f") + Spaces() +
                   Integer("domain_hmm-t") + Spaces() +
                   Group("domain_hmm_symbols", match_symbols) + Spaces() +
                   Float("domain_score") + Spaces() +
                   Float("domain_evalue") + AnyEol())

# 
# Alignments of top-scoring domains:
alignment_header = (ToEol() + # empty line
                    Str("Alignments of top-scoring domains:") +
                    AnyEol())

# PF00271: domain 1 of 1, from 686 to 767: score 58.6, E = 8.9e-17
#                    *->ell.epgikvarlhG.....glsqeeReeilekFrngkskvLvaTdv
#                       el ++  i++a++ G+++++++++++ ++++ kFr+g + +LvaT+v
#    AT1G01040   686    ELPsLSFIRCASMIGhnnsqEMKSSQMQDTISKFRDGHVTLLVATSV 732  
# 
#                    aarGlDipdvnlVInydlpwnpesyiQRiGRaGRaG<-*
#                    a++GlDi ++n+V  +dl +++  yiQ +GR +R     
#    AT1G01040   733 AEEGLDIRQCNVVMRFDLAKTVLAYIQSRGR-ARKP    767
#
domain_align_header = (Group("dalign_name", model_name) +
                       Str(": domain ") +
                       Integer("dalign_of_domain") +
                       Str(" of ") +
                       Integer("dalign_total_domain") +
                       Str(", from ") +
                       Integer("dalign_domain_start") +
                       Str(" to ") +
                       Integer("dalign_domain_end") +
                       Str(": score ") +
                       Float("dalign_score") +
                       Str(", E = ") +
                       Float("dalign_evalue") +
                       AnyEol())

# occasionally we have a line like:
#                 RF xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# preceeding the alignment. I'm not sure what this means.
rf_line = (Spaces() + Str("RF ") + ToEol())

domain_align_top = (Spaces() +
                    UntilEol("dalign_match_top") + AnyEol())

domain_align_middle = (Spaces() +
                       UntilEol("dalign_match_middle") + AnyEol())

domain_align_bottom = (Spaces() +
                       ToSep("dalign_query_name", " ") + Spaces() +
                       # some match ends can have:
                       # AT1G01230     -   -
                       # instead of an actual line, so there are fixes
                       # in here for that messed up case
                       Alt(Integer("dalign_query_start") + Spaces() +
                           Group("dalign_match_bottom",
                               Re("[\w\-]+")) + Spaces() +
                           Integer("dalign_query_end") +
                           Spaces() + AnyEol(),
                           Str("- ") + ToEol()) +
                       ToEol()) # blank line

domain_alignment = (domain_align_header +
                    Rep1(Opt(rf_line) +
                         domain_align_top +
                         domain_align_middle +
                         domain_align_bottom))

# //
record_end = Str("//") + AnyEol()


record = Std.record(Rep(sequence_info +
                        family_header +
                        (no_hit_line | Rep1(family_hit_line)) +
                        domain_header +
                        (no_hit_line | Rep1(domain_hit_line)) +
                        alignment_header +
                        (no_hit_line | Rep1(domain_alignment)) +
                        record_end
                    ))

format = HeaderFooter("hmmpfam", {},
                      header, RecordReader.CountLines, (8,),
                      record, RecordReader.EndsWith, ("//\n",),
                      None, None, None)
