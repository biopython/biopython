# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for calling ClustalW and parsing its output (DEPRECATED).

This module has been superseded by the Bio.AlignIO framework for
alignment parsing, and the ClustalW command line wrapper in
Bio.Align.Applications for calling the tool. These are both described
in the current version of the Biopython Tutorial and Cookbook.
This means Bio.Clustalw is now deprecated and likely to be
removed in future releases of Biopython.

A set of classes to interact with the multiple alignment command
line program clustalw. 

Clustalw is the command line version of the graphical Clustalx 
aligment program.

This requires clustalw available from:

ftp://ftp-igbmc.u-strasbg.fr/pub/ClustalW/.

functions:
o read
o parse_file
o do_alignment

classes:
o ClustalAlignment
o MultipleAlignCL"""

import Bio
import warnings
warnings.warn("Bio.Clustalw is deprecated. Please use the Bio.AlignIO framework for alignment parsing, and the ClustalW command line wrapper in Bio.Align.Applications for calling the tool. These are both described in the current version of the Biopython Tutorial and Cookbook.", Bio.BiopythonDeprecationWarning)

# standard library
import os
import sys
import subprocess

# biopython
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align.Generic import Alignment
from Bio.Application import _escape_filename

def parse_file(file_name, alphabet = IUPAC.unambiguous_dna, debug_level = 0):
    """Parse the given file into a clustal aligment object (OBSOLETE).
    
    Arguments:
    o file_name - The name of the file to parse.
    o alphabet - The type of alphabet to use for the alignment sequences.
    This should correspond to the type of information contained in the file.
    Defaults to be unambiguous_dna sequence.

    There is a deprecated optional argument debug_level which has no effect.

    This function is obsolete, and any new code should call Bio.AlignIO
    instead. For example using Bio.Clustalw, you might have:

    >>> from Bio import Clustalw
    >>> from Bio import Alphabet
    >>> filename = "Clustalw/protein.aln"
    >>> alpha = Alphabet.Gapped(Alphabet.generic_protein)
    >>> align = Clustalw.parse_file(filename, alphabet=alpha)
    >>> print align.get_alignment_length()
    411
    >>> clustalw_string = str(align)

    This becomes:

    >>> from Bio import AlignIO
    >>> from Bio import Alphabet
    >>> filename = "Clustalw/protein.aln"
    >>> alpha = Alphabet.Gapped(Alphabet.generic_protein)
    >>> align = AlignIO.read(open(filename), "clustal", alphabet=alpha)
    >>> print align.get_alignment_length()
    411
    >>> assert clustalw_string == align.format("clustal")
    """ 

    import warnings
    warnings.warn("This function is obsolete, and any new code should call Bio.AlignIO instead.", PendingDeprecationWarning)
    # Avoid code duplication by calling Bio.AlignIO to do this for us.
    handle = open(file_name, 'r')
    from Bio import AlignIO
    generic_alignment = AlignIO.read(handle, "clustal")
    handle.close()

    #Force this generic alignment into a ClustalAlignment... nasty hack
    if isinstance(alphabet, Alphabet.Gapped):
        alpha = alphabet
    else:
        alpha = Alphabet.Gapped(alphabet)
    clustal_alignment = ClustalAlignment(alpha)
    clustal_alignment._records = generic_alignment._records
    for record in clustal_alignment._records:
        record.seq.alphabet = alpha

    try:
        clustal_alignment._version = generic_alignment._version
    except AttributeError:
        #Missing the version, could be a 3rd party tool's output
        pass

    try :       
        clustal_alignment._star_info = generic_alignment._star_info
    except AttributeError:
        #Missing the consensus, again, this is not always present
        pass

    return clustal_alignment

def do_alignment(command_line, alphabet=None):
    """Perform an alignment with the given command line (OBSOLETE).
    
    Arguments:
    o command_line - A command line object that can give out
    the command line we will input into clustalw.
    o alphabet - the alphabet to use in the created alignment. If not
    specified IUPAC.unambiguous_dna and IUPAC.protein will be used for
    dna and protein alignment respectively.
    
    Returns:
    o A clustal alignment object corresponding to the created alignment.
    If the alignment type was not a clustal object, None is returned.

    This function (and the associated command line object) are now obsolete.
    Please use the Bio.Align.Applications.ClustalwCommandline wrapper with
    the Python subprocess module (and Bio.AlignIO for parsing) as described
    in the tutorial.
    """
    import warnings
    warnings.warn("This function (and the associated command line object) are now obsolete. Please use the Bio.Align.Applications.ClustalwCommandline wrapper with the Python subprocess module (and Bio.AlignIO for parsing) as described in the tutorial.", PendingDeprecationWarning)
    #We don't need to supply any piped input, but we setup the
    #standard input pipe anyway as a work around for a python
    #bug if this is called from a Windows GUI program.  For
    #details, see http://bugs.python.org/issue1124861
    child_process = subprocess.Popen(str(command_line),
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True,
                                     shell=(sys.platform!="win32")
                                     )
    #Use .communicate as can get deadlocks with .wait(), see Bug 2804
    child_process.communicate() #ignore the stdout and strerr data
    value = child_process.returncode
    child_process.stdin.close()
    child_process.stdout.close()
    child_process.stderr.close()
    
    # check the return value for errors, as on 1.81 the return value
    # from Clustalw is actually helpful for figuring out errors
    # TODO - Update this for new error codes using in clustalw 2

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
    """Work with the clustal aligment format (OBSOLETE).

    This format is the default output from clustal -- these files normally
    have an extension of .aln.

    This obsolete alignment object is a subclass of the more general alignment
    object used in Bio.AlignIO. The old practical difference is here str(align)
    would give the alignment as a string in clustal format, whereas in general
    you must do align.format("clustal"), which supports other formats too.
    """
    # the default version to use if one isn't set
    import warnings
    warnings.warn("This class is obsolete.", PendingDeprecationWarning)
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
        return self.format("clustal")

class MultipleAlignCL:
    """Represent a clustalw multiple alignment command line (OBSOLETE).

    This command line wrapper is considerd obsolete. Please use the replacement
    Bio.Align.Applications.ClustalwCommandline wrapper instead, which uses the
    standardised Bio.Application style interface. This is described in the
    tutorial, with examples using ClustalW.
    """
    import warnings
    warnings.warn("This command line wrapper is considerd obsolete. Please use the replacement Bio.Align.Applications.ClustalwCommandline wrapper instead, which uses the standardised Bio.Application style interface. This is described in the tutorial, with examples using ClustalW.", PendingDeprecationWarning)
    # set the valid options for different parameters
    OUTPUT_TYPES = ['GCG', 'GDE', 'PHYLIP', 'PIR', 'NEXUS', 'FASTA']
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

        #On Linux with clustalw 1.83, you can do:
        #clustalw input.faa
        #clustalw /full/path/input.faa
        #clustalw -INFILE=input.faa
        #clustalw -INFILE=/full/path/input.faa
        #
        #Note these fail (using DOS style slashes):
        #
        #clustalw /INFILE=input.faa
        #clustalw /INFILE=/full/path/input.faa
        #
        #On Windows XP with clustalw.exe 1.83, these work at
        #the command prompt:
        #
        #clustalw.exe input.faa
        #clustalw.exe /INFILE=input.faa
        #clustalw.exe /INFILE="input.faa"
        #clustalw.exe /INFILE="with space.faa"
        #clustalw.exe /INFILE=C:\full\path\input.faa
        #clustalw.exe /INFILE="C:\full path\with spaces.faa"
        #
        #Sadly these fail:
        #clustalw.exe "input.faa"
        #clustalw.exe "with space.faa"
        #clustalw.exe C:\full\path\input.faa
        #clustalw.exe "C:\full path\with spaces.faa"
        #
        #Testing today (using a different binary of clustalw.exe 1.83),
        #using -INFILE as follows seems to work.  However I had once noted:
        #These also fail but a minus/dash does seem to
        #work with other options (!):
        #clustalw.exe -INFILE=input.faa
        #clustalw.exe -INFILE=C:\full\path\input.faa
        #
        #Also these fail:
        #clustalw.exe "/INFILE=input.faa"
        #clustalw.exe "/INFILE=C:\full\path\input.faa"
        #
        #Thanks to Emanuel Hey for flagging this on the mailing list.
        #
        #In addition, both self.command and self.sequence_file
        #may contain spaces, so should be quoted. But clustalw
        #is fussy.
        cline = _escape_filename(self.command)
        cline += ' -INFILE=%s' % _escape_filename(self.sequence_file)

        # general options
        if self.type:
            cline += " -TYPE=%s" % self.type
        if self.is_quick == 1:
            #Some versions of clustalw are case sensitive,
            #and require -quicktree rather than -QUICKTREE
            cline += " -quicktree"
        if self.allow_negative == 1:
            cline += " -NEGATIVE"

        # output options
        if self.output_file:
            cline += " -OUTFILE=%s" % _escape_filename(self.output_file)
        if self.output_type:
            cline += " -OUTPUT=%s" % self.output_type
        if self.output_order:
            cline += " -OUTORDER=%s" % self.output_order
        if self.change_case:
            cline += " -CASE=%s" % self.change_case
        if self.add_seqnos:
            cline += " -SEQNOS=%s" % self.add_seqnos
        if self.new_tree:
            # clustal does not work if -align is written -ALIGN
            cline += " -NEWTREE=%s -align" % _escape_filename(self.new_tree)

        # multiple alignment options
        if self.guide_tree:
            cline += " -USETREE=%s" % _escape_filename(self.guide_tree)
        if self.protein_matrix:
            cline += " -MATRIX=%s" % self.protein_matrix
        if self.dna_matrix:
            cline += " -DNAMATRIX=%s" % self.dna_matrix
        if self.gap_open_pen:
            cline += " -GAPOPEN=%s" % self.gap_open_pen
        if self.gap_ext_pen:
            cline += " -GAPEXT=%s" % self.gap_ext_pen
        if self.is_no_end_pen == 1:
            cline += " -ENDGAPS"
        if self.gap_sep_range:
            cline += " -GAPDIST=%s" % self.gap_sep_range
        if self.is_no_pgap == 1:
            cline += " -NOPGAP"
        if self.is_no_hgap == 1:
            cline += " -NOHGAP"
        if len(self.h_gap_residues) != 0:
            # stick the list of residues together as one big list o' residues
            residue_list = ''
            for residue in self.h_gap_residues:
                residue_list = residue_list + residue
            cline += " -HGAPRESIDUES=%s" % residue_list
        if self.max_div:
            cline += " -MAXDIV=%s" % self.max_div
        if self.trans_weight:
            cline += " -TRANSWEIGHT=%s" % self.trans_weight

        return cline

    def set_output(self, output_file, output_type = None, output_order = None,
                   change_case = None, add_seqnos = None):
        """Set the output parameters for the command line.
        """
        self.output_file = output_file

        if output_type:
            output_type = output_type.upper()
            if output_type not in self.OUTPUT_TYPES:
                raise ValueError("Invalid output type %s. Valid choices are %s"
                                 % (output_type, self.OUTPUT_TYPES))
            else:
                self.output_type = output_type

        if output_order:
            output_order = output_order.upper()
            if output_order not in self.OUTPUT_ORDER:
                raise ValueError("Invalid output order %s. Valid choices are %s"
                                 % (output_order, self.OUTPUT_ORDER))
            else:
                self.output_order = output_order

        if change_case:
            change_case = change_case.upper()
            if output_type != "GDE":
                raise ValueError("Change case only valid for GDE output.")
            elif change_case not in self.CHANGE_CASE:
                raise ValueError("Invalid change case %s. Valid choices are %s"
                                 % (change_case, self.CHANGE_CASE))
            else:
                self.change_case = change_case

        if add_seqnos:
            add_seqnos = add_seqnos.upper()
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
        if protein_matrix.upper() in self.PROTEIN_MATRIX:
            self.protein_matrix = protein_matrix.upper()
        elif os.path.exists(protein_matrix):
            self.protein_matrix = protein_matrix
        else:
            raise ValueError("Invalid matrix %s. Options are %s or a file." %
                             (protein_matrix.upper(), self.PROTEIN_MATRIX))

    def set_dna_matrix(self, dna_matrix):
        """Set the type of DNA matrix to use.

        The dna_matrix can either be one of the defined types (iub or clustalw)
        or a file with the matrix to use."""
        if dna_matrix.upper() in self.DNA_MATRIX:
            self.dna_matrix = dna_matrix.upper()
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
        residue_type = residue_type.upper()
        if residue_type in self.RESIDUE_TYPES:
            self.type = residue_type
        else:
            raise ValueError("Invalid residue type %s. Valid choices are %s"
                             % (residue_type, self.RESIDUE_TYPES))

def _test():
    """Run the Bio.Clustalw module's doctests (PRIVATE).

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    _test()
