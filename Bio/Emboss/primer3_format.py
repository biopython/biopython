"""Martel definitions for the output files produced by primer3.
"""
import Martel

any_space = Martel.Re("[ ]+")
blank_line = Martel.AnyEol()

comment_line = Martel.Str("#") + Martel.ToEol()

# comments and blank lines in the file
comments = Martel.Group("comments",
                        blank_line +
                        comment_line +
                        blank_line +
                        comment_line)

#   1 PRODUCT SIZE: 289
product_size = Martel.Group("product_size",
                 Martel.Re("[\d]+"))
start_primer = Martel.Group("start_primer",
                            any_space + Martel.Re("[\d]+") +
                            Martel.Str(" PRODUCT SIZE: "))
primer_start_line = Martel.Group("primer_start_line",
                      start_primer +
                      product_size + Martel.AnyEol())

# a blank line that signifies a new primer is coming up
single_primer_line = Martel.Group("single_primer_line",
                                  blank_line)
                        
#      FORWARD PRIMER    1725   20  59.96  55.00  AGGGAAGGGATGCTAGGTGT
primer_space = Martel.Str(" " * 5)

any_integer = Martel.Re("[\d]+")
any_float = Martel.Re("[\d\.]+")
sequence = Martel.Re("[GATCN]+")

forward_primer_start = Martel.Group("forward_start",
                                   any_integer)
forward_primer_length = Martel.Group("forward_length",
                                     any_integer)
forward_primer_tm = Martel.Group("forward_tm",
                                 any_float)
forward_primer_gc = Martel.Group("forward_gc",
                                 any_float)
forward_primer_seq = Martel.Group("forward_seq",
                                  sequence)

forward_line = Martel.Group("forward_line", 
                 primer_space + Martel.Str("FORWARD PRIMER") +
                 any_space + forward_primer_start + any_space +
                 forward_primer_length + any_space + forward_primer_tm +
                 any_space + forward_primer_gc + any_space + 
                 forward_primer_seq + Martel.AnyEol())

#      REVERSE PRIMER    1994   20  59.99  55.00  AGAAGCACACCTCTCCCTGA
reverse_primer_start = Martel.Group("reverse_start",
                                   any_integer)
reverse_primer_length = Martel.Group("reverse_length",
                                     any_integer)
reverse_primer_tm = Martel.Group("reverse_tm",
                                 any_float)
reverse_primer_gc = Martel.Group("reverse_gc",
                                 any_float)
reverse_primer_seq = Martel.Group("reverse_seq",
                                  sequence)
reverse_line = Martel.Group("reverse_line", 
                 primer_space + Martel.Str("REVERSE PRIMER") +
                 any_space + reverse_primer_start + any_space +
                 reverse_primer_length + any_space + reverse_primer_tm +
                 any_space + reverse_primer_gc + any_space + 
                 reverse_primer_seq + Martel.AnyEol())

#     INTERNAL OLIGO     197   20  59.96  40.00  AATTGCACATAACGGATGCA

internal_oligo_start = Martel.Group("internal_start",
                                   any_integer)
internal_oligo_length = Martel.Group("internal_length",
                                     any_integer)
internal_oligo_tm = Martel.Group("internal_tm",
                                 any_float)
internal_oligo_gc = Martel.Group("internal_gc",
                                 any_float)
internal_oligo_seq = Martel.Group("internal_seq",
                                  sequence)
internal_line = Martel.Group("internal_line", 
                 primer_space + Martel.Str("INTERNAL OLIGO") +
                 any_space + internal_oligo_start + any_space +
                 internal_oligo_length + any_space + internal_oligo_tm +
                 any_space + internal_oligo_gc + any_space + 
                 internal_oligo_seq + Martel.AnyEol())


# XXX This record definition is ugly. But it works :-)
record = Martel.Group("primer3_record",
                      comments + \
                        Martel.Alt(
                         # case 1. primer file with nothing
                         Martel.Str("\n" * 3) +
                         Martel.Opt(Martel.Str("\n" * 4)),
                         # case 2. some primers have been picked
                         Martel.Rep(
                           # case 2a. we are designing a primer pair 
                           Martel.Alt(blank_line + primer_start_line,
                           # case 2b. we are designing a single primer
                                      single_primer_line) +
                           # case 2a. both primer pairs
                           Martel.Alt(forward_line + blank_line + 
                                      reverse_line + blank_line,
                           # case 2b1. Reverse primer
                                      reverse_line + blank_line,
                           # case 2b2, Forward primer
                                      forward_line + blank_line,
                           # case 2b3, Internal oligo 
                                      internal_line + blank_line)) +
                           blank_line + blank_line + Martel.Rep(blank_line)))
                        
                          
