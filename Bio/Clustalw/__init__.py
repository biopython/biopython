# Clustalw modules

"""Clustalw.py

A set of classes to interact with the multiple alignment command
line program Clustalw (which is the command line version of the current
graphical Clustalx aligment program).

This requires Clustalw1.8.1 available from:

ftp://ftp-igbmc.u-strasbg.fr/pub/ClustalW/.

Older versions may work, but haven't been tested yet.

functions:
o parse_file
o do_alignment

classes:
o ClustalAlignment
o _AlignCreator
o MultipleAlignCL"""

__all__ = [
    'clustal_format',
    ]

# standard library
import os
import string

# biopython
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Alphabet import IUPAC
import clustal_format
from Bio.Align.Generic import Alignment

# PyXML package
from xml.sax import saxutils
from xml.sax import handler

def parse_file(file_name, alphabet = IUPAC.unambiguous_dna):
    """Parse the given file into a clustal aligment object.
    
    Arguments:
    o file_name - The name of the file to parse.
    o alphabet - The type of alphabet to use for the alignment sequences.
    This should correspond to the type of information contained in the file.
    Defaults to be unambiguous_dna sequence.
    """ 
    align_handler = _AlignCreator(Alphabet.Gapped(alphabet))

    parser = clustal_format.format.make_parser()
    parser.setContentHandler(align_handler)
    parser.setErrorHandler(handler.ErrorHandler())

    to_parse = open(file_name, 'r')
    parser.parseFile(to_parse)
    to_parse.close()

    return align_handler.align

def do_alignment(command_line, alphabet=None):
    """Perform an alignment with the given command line.

    Arguments:
    o command_line - A command line object that can give out
    the command line we will input into clustalw.
    o alphabet - the alphabet to use in the created alignment. If not
    specified IUPAC.unambiguous_dna and IUPAC.protein will be used for
    dna and protein alignment respectively.
    
    Returns:
    o A clustal alignment object corresponding to the created alignment.
    If the alignment type was not a clustal object, None is returned.
    """
    run_clust = os.popen(str(command_line))
    value = run_clust.close()

    # check the return value for errors, as on 1.81 the return value
    # from Clustalw is actually helpful for figuring out errors
    # 1 => bad command line option
    if value == 1:
        raise ValueError("Bad command line option in the command: %s"
                         % str(command_line))
    # 2 => can't open sequence file
    elif value == 2:
        raise IOError("Cannot open sequence file %s"
                      % command_line.sequence_file)
    # 3 => wrong format in sequence file
    elif value == 3:
        raise IOError("Sequence file %s has an invalid format."
                      % command_line.sequence_file)
    # 4 => sequence file only has one sequence
    elif value == 4:
        raise IOError("Sequence file %s has only one sequence present."
                      % command_line.sequence_file)

    # if an output file was specified, we need to grab it
    if command_line.output_file:
        out_file = command_line.output_file
    else:
        out_file = os.path.splitext(command_line.sequence_file)[0] + '.aln'

    # if we can't deal with the format, just return None
    if command_line.output_type and command_line.output_type != 'CLUSTAL':
        return None
    # otherwise parse it into a ClustalAlignment object
    else:
        if not alphabet:
            alphabet = (IUPAC.unambiguous_dna, IUPAC.protein)[
                command_line.type == 'PROTEIN']
            
        # check if the outfile exists before parsing
        if not(os.path.exists(out_file)):
            raise IOError("Output .aln file %s not produced, commandline: %s"
                          % (out_file, command_line))
            
        return parse_file(out_file, alphabet)


class ClustalAlignment(Alignment):
    """Work with the clustal aligment format.

    This format is the default output from clustal -- these files normally
    have an extension of .aln.
    """
    # the default version to use if one isn't set
    DEFAULT_VERSION = '1.81'
    
    def __init__(self, alphabet = Alphabet.Gapped(IUPAC.ambiguous_dna)):
        Alignment.__init__(self, alphabet)

        # represent all of those stars in the aln output format
        self._star_info = ''
        
        self._version = ''

    def __str__(self):
        """Print out the alignment so it looks pretty.

        The output produced from this should also be formatted in valid
        clustal format.
        """
        # if the version isn't set, we need to use the default
        if self._version == '':
            self._version = self.DEFAULT_VERSION
        
        output = "CLUSTAL X (%s) multiple sequence alignment\n\n\n" % \
                 self._version

        cur_char = 0
        max_length = len(self._records[0].seq)

        # keep displaying sequences until we reach the end
        while cur_char != max_length:
            # calculate the number of sequences to show, which will
            # be less if we are at the end of the sequence
            if (cur_char + 50) > max_length:
                show_num = max_length - cur_char
            else:
                show_num = 50

            # go through all of the records and print out the sequences
            # when we output, we do a nice 80 column output, although this
            # may result in truncation of the ids.
            for record in self._records:
                line = record.description[0:30].ljust(36)
                line = line + record.seq.data[cur_char:(cur_char + show_num)]
                
                output = output + line + "\n"

            # now we need to print out the star info, if we've got it
            if self._star_info != '':
                output = output + (" " * 36) + \
                     self._star_info[cur_char:(cur_char + show_num)] + "\n"

            output = output + "\n"
            cur_char = cur_char + show_num

        # have a extra newline, so strip two off and add one before returning
        return string.rstrip(output) + "\n"
            

    def _add_star_info(self, stars):
        """Add all of the stars, which indicate consensus sequence.
        """
        self._star_info = stars

    def _add_version(self, version):
        """Add the version information about the clustal file being read.
        """
        self._version = version

class _AlignCreator(handler.ContentHandler):
    """Handler to create a ClustalAlignment object from clustal file info.

    This handler is used to accept events coming from a Martel parsing
    stream, and acts like a normal SAX handler.

    After parsing, the alignment object created is available as the
    align attribute of the class.
    """
    def __init__(self, alphabet):
        """Create a new handler ready to deal with output from Martel parsing.

        Arguments:
        o alphabet - The alphabet to create all of the new sequences with.
        """
        self.align = ClustalAlignment(alphabet)

        # store sequence info in a dictionary
        self.all_info = {}
        self.all_keys = []

        # the current id we are working with
        self.cur_id = None

        # info so we know how big the ids and sequences are
        self.id_size = 0
        self.space_size = 0
        self.seq_size = 0
        
        # flags so we can keep track of where we are during the parse
        self.in_version = 0
        self.in_stars = 0
        self.in_seq_id = 0
        self.in_space = 0
        self.in_seq = 0
        self.all_star_info = ''

    def startElement(self, name, attrs):
        """Check the various tags for the info we are interested in."""
        if name == "version":
            self.in_version = 1
            self.version_info = ''
        elif name == "seq_id":
            self.in_seq_id = 1
            self.seq_id_info = ''
        elif name == "seq_space":
            self.in_space = 1
            self.space_info = ''
        elif name == "seq_info":
            self.in_seq = 1
            self.seq_info = ''
        elif name == "match_stars":
            self.in_stars = 1
            self.star_info = ''

    def characters(self, content):
        if self.in_version:
            self.version_info = self.version_info + content
        elif self.in_seq_id:
            self.seq_id_info = self.seq_id_info + content
        elif self.in_space:
            self.space_info = self.space_info + content
        elif self.in_seq:
            self.seq_info = self.seq_info + content
        elif self.in_stars:
            self.star_info = self.star_info + content

    def endElement(self, name):
        if name == "version":
            self.in_version = 0
            self.align._add_version(string.strip(self.version_info))
        elif name == "seq_id":
            self.in_seq_id = 0
            self.id_size = len(self.seq_id_info)
            self.cur_id = self.seq_id_info
        elif name == "seq_space":
            self.in_space = 0
            self.space_size = len(self.space_info)
        elif name == "seq_info":
            self.in_seq = 0
            self.seq_size = len(self.seq_info)

            # if the id is already there, add the sequence info
            if self.cur_id in self.all_info.keys():
                self.all_info[self.cur_id] = self.all_info[self.cur_id] + \
                                             self.seq_info
            else:
                self.all_info[self.cur_id] = self.seq_info
                self.all_keys.append(self.cur_id)
                
        elif name == "match_stars":
            id_length = self.id_size + self.space_size
            line_length = id_length + self.seq_size
            
            self.all_star_info = self.all_star_info + \
                            self.star_info[id_length:line_length]
                                                                           
    def endDocument(self):
        # when we are done parsing add all of the info we need
        self.align._add_star_info(self.all_star_info)

        for id in self.all_keys:
            self.align.add_sequence(id, self.all_info[id])
        
            
      
class MultipleAlignCL:
    """Represent a clustalw multiple alignment command line.

    This is meant to make it easy to code the command line options you
    want to submit to clustalw.

    Clustalw has a ton of options and things to do but this is set up to
    represent a clustalw mutliple alignment.

    Warning: I don't use all of these options personally, so if you find
    one to be broken for any reason, please let us know!
    """
    # set the valid options for different parameters
    OUTPUT_TYPES = ['GCG', 'GDE', 'PHYLIP', 'PIR', 'NEXUS']
    OUTPUT_ORDER = ['INPUT', 'ALIGNED']
    OUTPUT_CASE = ['LOWER', 'UPPER']
    OUTPUT_SEQNOS = ['OFF', 'ON']
    RESIDUE_TYPES = ['PROTEIN', 'DNA']
    PROTEIN_MATRIX = ['BLOSUM', 'PAM', 'GONNET', 'ID']
    DNA_MATRIX = ['IUB', 'CLUSTALW']

    def __init__(self, sequence_file, command = 'clustalw'):
        """Initialize some general parameters that can be set as attributes.

        Arguments:
        o sequence_file - The file to read the sequences for alignment from.
        o command - The command used to run clustalw. This defaults to
        just 'clustalw' (ie. assumes you have it on your path somewhere).
        
        General attributes that can be set:
        o is_quick - if set as 1, will use a fast algorithm to create
        the alignment guide tree.
        o allow_negative - allow negative values in the alignment matrix.

        Multiple alignment attributes that can be set as attributes:
        o gap_open_pen - Gap opening penalty
        o gap_ext_pen - Gap extension penalty
        o is_no_end_pen - A flag as to whether or not there should be a gap
        separation penalty for the ends.
        o gap_sep_range - The gap separation penalty range.
        o is_no_pgap - A flag to turn off residue specific gaps
        o is_no_hgap - A flag to turn off hydrophilic gaps
        o h_gap_residues - A list of residues to count a hydrophilic
        o max_div - A percent identity to use for delay (? - I don't undertand
        this!)
        o trans_weight - The weight to use for transitions
        """
        self.sequence_file = sequence_file
        self.command = command
        
        self.is_quick = None
        self.allow_negative = None

        self.gap_open_pen = None
        self.gap_ext_pen = None
        self.is_no_end_pen = None
        self.gap_sep_range = None
        self.is_no_pgap = None
        self.is_no_hgap = None
        self.h_gap_residues = []
        self.max_div = None
        self.trans_weight = None

        # other attributes that should be set via various functions
        # 1. output parameters
        self.output_file = None
        self.output_type = None
        self.output_order = None
        self.change_case = None
        self.add_seqnos = None

        # 2. a guide tree to use
        self.guide_tree = None
        self.new_tree = None

        # 3. matrices
        self.protein_matrix = None
        self.dna_matrix = None

        # 4. type of residues
        self.type = None

    def __str__(self):
        """Write out the command line as a string."""
        cline = self.command + " " + self.sequence_file

        # general options
        if self.type:
            cline = cline + " -TYPE=" + self.type
        if self.is_quick == 1:
            cline = cline + " -INTERACTIVE"
        if self.allow_negative == 1:
            cline = cline + " -NEGATIVE"

        # output options
        if self.output_file:
            cline = cline + " -OUTFILE=" + self.output_file
        if self.output_type:
            cline = cline + " -OUTPUT=" + self.output_type
        if self.output_order:
            cline = cline + " -OUTORDER=" + self.output_order
        if self.change_case:
            cline = cline + " -CASE=" + self.change_case
        if self.add_seqnos:
            cline = cline + " -SEQNOS=" + self.add_seqnos
        if self.new_tree:
            # clustal does not work if -align is written -ALIGN
            cline = cline + " -NEWTREE=" + self.new_tree + " -align"

        # multiple alignment options
        if self.guide_tree:
            cline = cline + " -USETREE=" + self.guide_tree
        if self.protein_matrix:
            cline = cline + " -MATRIX=" + self.protein_matrix
        if self.dna_matrix:
            cline = cline + " -DNAMATRIX=" + self.dna_matrix
        if self.gap_open_pen:
            cline = cline + " -GAPOPEN=" + self.gap_open_pen
        if self.gap_ext_pen:
            cline = cline + " -GAPEXT=" + self.gap_ext_pen
        if self.is_no_end_pen == 1:
            cline = cline + " -ENDGAPS"
        if self.gap_sep_range:
            cline = cline + " -GAPDIST=" + self.gap_sep_range
        if self.is_no_pgap == 1:
            cline = cline + " -NOPGAP"
        if self.is_no_hgap == 1:
            cline = cline + " -NOHGAP"
        if len(self.h_gap_residues) != 0:
            # stick the list of residues together as one big list o' residues
            residue_list = ''
            for residue in self.h_gap_residues:
                residue_list = residue_list + residue
            cline = cline + " -HGAPRESIDUES=" + residue_list
        if self.max_div:
            cline = cline + " -MAXDIV=" + self.max_div
        if self.trans_weight:
            cline = cline + " -TRANSWEIGHT=" + self.trans_weight

        return cline

    def set_output(self, output_file, output_type = None, output_order = None,
                   change_case = None, add_seqnos = None):
        """Set the output parameters for the command line.
        """
        self.output_file = output_file

        if output_type:
            output_type = string.upper(output_type)
            if output_type not in self.OUTPUT_TYPES:
                raise ValueError("Invalid output type %s. Valid choices are %s"
                                 % (output_type, self.OUTPUT_TYPES))
            else:
                self.output_type = output_type

        if output_order:
            output_order = string.upper(output_order)
            if output_order not in self.OUTPUT_ORDER:
                raise ValueError("Invalid output order %s. Valid choices are %s"
                                 % (output_order, self.OUTPUT_ORDER))
            else:
                self.output_order = output_order

        if change_case:
            change_case = string.upper(change_case)
            if output_type != "GDE":
                raise ValueError("Change case only valid for GDE output.")
            elif change_case not in self.CHANGE_CASE:
                raise ValueError("Invalid change case %s. Valid choices are %s"
                                 % (change_case, self.CHANGE_CASE))
            else:
                self.change_case = change_case

        if add_seqnos:
            add_seqnos = string.upper(add_seqnos)
            if output_type:
                raise ValueError("Add SeqNos only valid for CLUSTAL output.")
            elif add_seqnos not in self.OUTPUT_SEQNOS:
                raise ValueError("Invalid seqnos option %s. Valid choices: %s"
                                 % (add_seqnos, self.OUTPUT_SEQNOS))
            else:
                self.add_seqnos = add_seqnos

    def set_guide_tree(self, tree_file):
        """Provide a file to use as the guide tree for alignment.

        Raises:
        o IOError - If the tree_file doesn't exist."""
        if not(os.path.exists(tree_file)):
            raise IOError("Could not find the guide tree file %s." %
                          tree_file)
        else:
            self.guide_tree = tree_file

    def set_new_guide_tree(self, tree_file):
        """Set the name of the guide tree file generated in the alignment.
        """
        self.new_tree = tree_file
        
    def set_protein_matrix(self, protein_matrix):
        """Set the type of protein matrix to use.

        Protein matrix can be either one of the defined types (blosum, pam,
        gonnet or id) or a file with your own defined matrix.
        """
        if string.upper(protein_matrix) in self.PROTEIN_MATRIX:
            self.protein_matrix = string.upper(protein_matrix)
        elif os.path.exists(protein_matrix):
            self.protein_matrix = protein_matrix
        else:
            raise ValueError("Invalid matrix %s. Options are %s or a file." %
                             (string.upper(protein_matrix),
                              self.PROTEIN_MATRIX))

    def set_dna_matrix(self, dna_matrix):
        """Set the type of DNA matrix to use.

        The dna_matrix can either be one of the defined types (iub or clustalw)
        or a file with the matrix to use."""
        if string.upper(dna_matrix) in self.DNA_MATRIX:
            self.dna_matrix = string.upper(dna_matrix)
        elif os.path.exists(dna_matrix):
            self.dna_matrix = dna_matrix
        else:
            raise ValueError("Invalid matrix %s. Options are %s or a file." %
                             (dna_matrix, self.DNA_MATRIX))

    def set_type(self, residue_type):
        """Set the type of residues within the file.

        Clustal tries to guess whether the info is protein or DNA based on
        the number of GATCs, but this can be wrong if you have a messed up
        protein or DNA you are working with, so this allows you to set it
        explicitly.
        """
        residue_type = string.upper(residue_type)
        if residue_type in self.RESIDUE_TYPES:
            self.type = residue_type
        else:
            raise ValueError("Invalid residue type %s. Valid choices are %s"
                             % (residue_type, self.RESIDUE_TYPES))

