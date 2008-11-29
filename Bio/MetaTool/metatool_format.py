# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read MetaTool output files.

This is a huge regular regular expression for MetaTool 3.5 output, built using
the 'regular expressiona on steroids' capabilities of Martel.

http://www2.bioinf.mdc-berlin.de/metabolic/metatool/


This helps us have endlines be consistent across platforms.

"""

# Martel
from Martel import Opt, Alt, Digits, Integer, Group, Str, MaxRepeat
from Martel import Any, AnyBut, RepN, Rep, Rep1, ToEol, AnyEol
from Martel import Expression
from Martel import RecordReader

blank = ' '
tab = '\t'
blank_space = MaxRepeat( Any( blank + tab), 1, 80 )
optional_blank_space = Rep( Any( blank + tab ) )
white_space = " \t" + chr( 10 ) + chr( 13 )
blank_line = optional_blank_space + AnyEol()
lower_case_letter = Group( "lower_case_letter", Any( "abcdefghijklmnopqrstuvwxyz" ) )
digits = "0123456789"

enzyme = Group( "enzyme", optional_blank_space + Digits() +
    optional_blank_space + Str( ':' ) + ToEol() )
reaction = Group( "reaction", optional_blank_space + Digits() +
    optional_blank_space + Str( ":" ) + ToEol() )
not_found_line = Group( "not_found_line", optional_blank_space + Str( "- not found -" ) +
    ToEol() )

enzymes_header = Group( "enzymes_header", optional_blank_space + Str( "enzymes" ) +
     ToEol() )
enzymes_list = Group( "enzymes_list", Alt( Rep1( enzyme ), \
    not_found_line ) )
enzymes_block = Group( "enzymes_block", enzymes_header + Rep( blank_line ) +
    enzymes_list )

reactions_header = Group( "reactions_header", optional_blank_space +
    Str( "overall reaction" ) + ToEol() )
reactions_list = Group( "reactions_list", Alt( Rep1( reaction ), \
    not_found_line ) )
reactions_block = Group( "reactions_block", reactions_header + Rep( blank_line ) +
    reactions_list )

rev = Group( "rev", Opt( lower_case_letter ) )
version = Group( "version", Digits( "version_major") + Any( "." ) +
    Digits( "version_minor") + rev )
metatool_tag = Str( "METATOOL OUTPUT" )
metatool_line = Group( "metatool_line", metatool_tag + blank_space +
    Str( "Version" ) + blank_space + version + ToEol() )

input_file_tag = Str( "INPUT FILE:" )
input_file_line = Group( "input_file_line", input_file_tag + blank_space +
    ToEol( "input_file_name" ) )

metabolite_count_tag = Str( "INTERNAL METABOLITES:" )
metabolite_count_line = Group( "metabolite_count_line",  metabolite_count_tag +
    blank_space + Digits( "num_int_metabolites" ) + ToEol() )

reaction_count_tag = Str( "REACTIONS:" )
reaction_count_line = Group( "reaction_count_line", reaction_count_tag + blank_space +
    Digits( "num_reactions" ) + ToEol() )

type_metabolite = Group( "type_metabolite", Alt( Str( "int" ), \
    Str( "external" ) ) )
metabolite_info = Group( "metabolite_info", optional_blank_space +
    Digits() + blank_space + type_metabolite + blank_space +
#    Integer() + blank_space + Rep1( lower_case_letter ) +
    Rep1( AnyBut( white_space ) ) )
metabolite_line = Group( "metabolite_line", metabolite_info + ToEol() )
metabolites_summary = Group( "metabolites_summary", optional_blank_space + Digits() +
    blank_space + Str( "metabolites" ) + ToEol() )
metabolites_block = Group( "metabolites_block", Rep1( metabolite_line ) +
    metabolites_summary + Rep( blank_line ) )

graph_structure_heading = Group( "graph_structure_heading", optional_blank_space +
    Str( "edges" ) + blank_space + Str( "frequency of nodes" ) + ToEol() )
graph_structure_line = Group( "graph_structure_line", optional_blank_space +
    Digits( "edge_count" ) + blank_space + Digits( "num_nodes" ) + ToEol() )
graph_structure_block =  Group( "graph_structure_block", \
    graph_structure_heading + Rep( blank_line ) +
    Rep1( graph_structure_line ) + Rep( blank_line ) )

sum_is_constant_line = Group( "sum_is_constant_line", optional_blank_space +
    Digits() + optional_blank_space + Any( ":" ) + optional_blank_space +
    Rep1( AnyBut( white_space ) ) +
    Rep( blank_space + Any( "+" ) + blank_space + Rep1( AnyBut( white_space ) ) ) +
    optional_blank_space + Str( "=" ) + ToEol() )
sum_is_constant_block = Group( "sum_is_constant_block", Rep( sum_is_constant_line ) )


stoichiometric_tag = Group( "stoichiometric_tag", Str( "STOICHIOMETRIC MATRIX" ) )
stoichiometric_line = Group( "stoichiometric_line", stoichiometric_tag +
    ToEol() )

not_balanced_tag = Group( "not_balanced_tag", Str( "NOT BALANCED INTERNAL METABOLITES" ) )
not_balanced_line = Group( "not_balanced_line", not_balanced_tag +
    ToEol() )

subsets_tag = Group( "subsets_tag", Str( "SUBSETS OF REACTIONS" ) )
subsets_line = Group( "subsets_line", \
     subsets_tag + ToEol() )

reduced_system_tag = Group( "reduced_system_tag", Str( "REDUCED SYSTEM" ) )
reduced_system_line = Group( "reduced_system_line", reduced_system_tag +
    Rep1(  AnyBut( digits ) ) + Digits( "branch_points" ) +
    Rep1( AnyBut( digits ) ) + Digits() + ToEol() )

kernel_tag = Group( "kernel_tag", Str( "KERNEL" ) )
kernel_line = Group( "kernel_line", kernel_tag + ToEol() )

convex_basis_tag = Group( "convex_basis_tag", Str( "CONVEX BASIS" ) )
convex_basis_line = Group( "convex_basis_line", convex_basis_tag +
    ToEol() )

conservation_relations_tag = Group( "conservation_relations_tag", \
    Str( "CONSERVATION RELATIONS" ) )
conservation_relations_line = Group( "conservation_relations_line", \
    conservation_relations_tag + ToEol() )

elementary_modes_tag = Group( "elementary_modes_tag", \
    Str( "ELEMENTARY MODES" ) )
elementary_modes_line = Group( "elementary_modes_line", \
    elementary_modes_tag + ToEol() )

num_rows = Group( "num_rows", Digits() )
num_cols = Group( "num_cols", Digits() )
matrix_header = Group( "matrix_header", optional_blank_space +
    Str( "matrix dimension" ) + blank_space  + Any( "r" ) +
    num_rows + blank_space +  Any( "x" ) + blank_space +
    Any( "c" ) + num_cols + optional_blank_space + AnyEol() )
matrix_element = Group( "matrix_element", Integer() )
matrix_row = Group( "matrix_row", MaxRepeat( optional_blank_space + matrix_element, \
    "num_cols", "num_cols" ) + ToEol() )
matrix = Group( "matrix", MaxRepeat( matrix_row, "num_rows", "num_rows" ) )

matrix_block = Group( "matrix_block", matrix_header + matrix )
irreversible_vector = Group( "irreversible_vector", \
    MaxRepeat( blank_space + matrix_element, "num_cols", "num_cols" ) + 
    ToEol() )

little_gap = Str( " " )
big_gap = Alt( Str( "\t" ), MaxRepeat( Str( " " ), 2, 80 ) )
unbalanced_metabolite = Group( "unbalanced_metabolite", \
    Rep1( AnyBut( white_space ) ) + Opt( little_gap +
    Rep1( AnyBut( white_space ) ) ) )
not_balanced_data = Group( "not_balanced_data", optional_blank_space +
    unbalanced_metabolite + Rep( big_gap + unbalanced_metabolite ) + ToEol() )

metabolite_roles_heading = Group( "metabolite_roles_heading", \
    Str( "->" ) + ToEol() )
metabolite_role_cols = Group( "metabolite_role_cols", \
    optional_blank_space + Str( "met" ) + blank_space + Str( "cons" ) +
    blank_space + Str( "built" ) +
    blank_space + Str( "reactions" ) + ToEol() )
branch_metabolite = Group( "branch_metabolite", optional_blank_space +
    Rep1( AnyBut( white_space ) ) + blank_space +
    RepN( Digits() + blank_space, 3 ) + Rep1( Any( "ir" ) ) + ToEol() )
non_branch_metabolite = Group( "non_branch_metabolite", optional_blank_space +
    Rep1( AnyBut( white_space ) ) + blank_space +
    RepN( Digits() + blank_space, 3 ) + Rep1( Any( "ir" ) ) + ToEol() )
branch_metabolite_block = Group( "branch_metabolite_block", \
    metabolite_roles_heading +
    metabolite_role_cols + Rep( branch_metabolite ) )
non_branch_metabolite_block = Group( "non_branch_metabolite_block", \
    metabolite_roles_heading +
    metabolite_role_cols + Rep( non_branch_metabolite ) )

end_stoichiometric = Group( "end_stochiometric", \
    Rep( Expression.Assert( not_balanced_tag, 1 ) +
    Expression.Assert( kernel_tag, 1 ) + ToEol() ) )
end_not_balanced = Group( "end_not_balanced", \
    Rep( Expression.Assert( kernel_tag, 1 ) + ToEol() ) )
end_kernel = Group( "end_kernel", \
    Rep( Expression.Assert( subsets_tag, 1 ) + ToEol() ) )
end_subsets = Group( "end_subsets", \
    Rep( Expression.Assert( reduced_system_tag, 1 ) + ToEol() ) )
end_reduced_system = Group( "end_reduced_system", \
    Rep( Expression.Assert( convex_basis_tag, 1 ) + ToEol() ) )
end_convex_basis = Group( "end_convex_basis", \
    Rep( Expression.Assert( conservation_relations_tag, 1 ) + ToEol() ) )
end_conservation_relations = Group( "end_conservation_relations", \
    Rep( Expression.Assert( elementary_modes_tag, 1 ) + ToEol() ) )
end_elementary_modes = Group( "end_elementary_modes", Rep( ToEol() ) )
#    Rep1( AnyBut( '.') ) + Str( "." ) )

input_file_block = Group( "input_file_block", input_file_line +
    Rep( blank_line ) )
metatool_block = Group( "metatool_block", metatool_line + Rep1( blank_line ) )

metabolite_count_block = Group( "metabolite_count_block", \
    metabolite_count_line + Rep( blank_line ) )
reaction_count_block = Group( "reaction_count_block", reaction_count_line +
    Rep( blank_line ) + metabolites_block + Rep( blank_line ) +
    graph_structure_block + Rep( blank_line ) )
stoichiometric_block = Group( "stoichiometric_block", stoichiometric_line +
    Rep( blank_line ) + matrix_block + ToEol() + irreversible_vector +
    end_stoichiometric )
not_balanced_block = Group( "not_balanced_block", not_balanced_line +
    Rep( blank_line ) + not_balanced_data + Rep( blank_line ) )
kernel_block = Group( "kernel_block", kernel_line + Rep( blank_line ) +
    matrix_block + ToEol() + Rep( blank_line ) + enzymes_block +
    Rep( blank_line ) + reactions_block + end_kernel )
subsets_block = Group( "subsets_block", subsets_line + Rep( blank_line ) +
    matrix_block + ToEol() + Rep( blank_line ) + enzymes_block +
    Rep( blank_line ) + reactions_block + end_subsets )
reduced_system_block = Group( "reduced_system_block", reduced_system_line +
    Rep( blank_line ) + matrix_block + ToEol() + irreversible_vector +
    Rep( blank_line ) + branch_metabolite_block + Rep( blank_line ) +
    non_branch_metabolite_block + end_reduced_system )
convex_basis_block = Group( "convex_basis_block", convex_basis_line +
    Rep( blank_line ) + matrix_block + Opt( ToEol() ) + Rep( blank_line ) +
    enzymes_block + Rep( blank_line ) + reactions_block + end_convex_basis )
conservation_relations_block = Group( "conservation_relations_block", \
    conservation_relations_line + Rep( blank_line ) + matrix_block +
    Rep( blank_line ) + sum_is_constant_block +
    end_conservation_relations )
elementary_modes_block = Group( "elementary_modes_block", elementary_modes_line +
    Rep( blank_line ) + matrix_block + Opt( ToEol() ) + Rep( blank_line ) +
    enzymes_block + Rep( blank_line ) + reactions_block + end_elementary_modes )


metatool_record = Group( "metatool_record", metatool_block + input_file_block +
   metabolite_count_block + reaction_count_block + stoichiometric_block +
    Opt( not_balanced_block ) + kernel_block + subsets_block +
    reduced_system_block + convex_basis_block + conservation_relations_block +
    elementary_modes_block )
