"""Martel format for primersearch output files,
"""
import Martel

blank_line = Martel.AnyEol()

# Primer name D1S2660
primer_name = Martel.Group("primer_name",
                           Martel.ToEol())
primer_name_line = Martel.Str("Primer name ") + primer_name

# Amplimer 1
amplifier = Martel.Group("amplifier",
                         Martel.Re("[\d]+"))
amplimer_line = Martel.Str("Amplimer ") + amplifier + Martel.AnyEol()

# Sequence: AC074298 AC074298
# Telomere associated sequence for Arabidopsis thaliana TEL1N
# CCGGTTTCTCTGGTTGAAAA hits forward strand at 114 with 0 mismatches
# TCACATTCCCAAATGTAGATCG hits reverse strand at [114] with 0 mismatches
seq_indent = Martel.Str("\t")

sequence_id = Martel.Group("sequence_id",
                           Martel.ToEol())
sequence_descr = Martel.Group("sequence_descr",
                              Martel.ToEol())
sequence_info = sequence_id + sequence_descr
forward_strand_info = Martel.Group("forward_strand_info",
                                   Martel.ToEol())
reverse_strand_info = Martel.Group("reverse_strand_info",
                                   Martel.ToEol())
amplifier_sequence = Martel.Group("amplifier_sequence",
                                  sequence_info + 
                                  forward_strand_info +
                                  reverse_strand_info)
amplifier_sequence_lines = seq_indent + Martel.Str("Sequence: ") + \
                           amplifier_sequence

amplifier_length = Martel.Group("amplifier_length", 
                                Martel.Re("[\d]+"))
amplifier_length_line = seq_indent + Martel.Str("Amplimer length: ") + \
                        amplifier_length + Martel.Str(" bp") + \
                        Martel.AnyEol()

record = Martel.Group("primersearch_record",
                      Martel.Rep1(blank_line + primer_name_line +
                        Martel.Rep(amplimer_line +
                                   amplifier_sequence_lines +
                                   amplifier_length_line)))

