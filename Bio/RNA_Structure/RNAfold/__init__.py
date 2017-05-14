# -*- coding: utf-8 -*-
# Copyright 2017 by Joanna Zbijewska, Agata Gruszczyńska, Michał Karlicki.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that:
# 1. RMSD scripts and license is available separately.
#    We added it in file: calculate_rmsd and README2.
#
# 2. RNAstructure license is available separately.
#    Please consult rna.urmc.rochester.edu .
#
#example of usage:
#
#cline = RNAfoldCommandLine(infile="~/Desktop/jakis_plik.fasta",energyModel=0, #dangles=2,maxBPspan=1)
#print("About to run: %s" % cline)
#std_output, err_output = cline()
#print(std_output)


from __future__ import print_function
import re
from Bio.Application import _Option, _Switch, AbstractCommandline


"""
example of usage:

cline = RNAfoldCommandLine(infile="~/Desktop/jakis_plik.fasta",energyModel=0, dangles=2,maxBPspan=1)
print("About to run: %s" % cline)
std_output, err_output = cline()
print(std_output)

"""


class _RNAfoldMinimalCommandLine(AbstractCommandline):
    """Base Commandline object for RNAfold wrapper.
    It deals with shared options:
    - help              Print help and exit
    - detailed-help     Print help, including all details and hidden options, and exit
    - full-help         Print help, including hidden options, and exit
    - version           Print version and exit
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [

            _Switch(["-h", "help"],
                    """help
                    Print help and exit"""),

            _Switch(["--detailed-help", "detailed_help"],
                    """detailed_help
                    Print help, including all details and hidden options, and exit"""),

            _Switch(["--full-help", "full_help"],
                    """full_help
                    Print help, including hidden options, and exit"""),

            _Switch(["-V", "--version", "version"],
                    """version
                    Print version and exit"""),
            ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

class RNAfoldCommandLine(_RNAfoldMinimalCommandLine):
    """

    """
    def __init__(self, cmd="RNAfold", **kwargs):
        #assert cmd is not None
        self.parameters = [
#General options
    #Command line options which alter the general behavior of this program
            _Option(["-v", "--verbose", "verbose"],
                    """verbose
                    Be verbose. (default=off)"""),
            _Option(["-i", "--infile", "infile"], #infile=<filename>
                    """infile
                    Read a file instead of reading from stdin.
                    The default behavior of RNAfold is to read input from
                    stdin. Using this parameter the user can specify an input
                    file name where data is read from.""",filename=True,equate=False),

            _Option(["-o", "--outfile", "outfile"], #outfile=<filename>
                    """outfile
                    Print output to file instead of stdout
                    Specifying a file name/prefix will print the output into a file
                    instead of stdout. If a FASTA header precedes the input sequence,
                    the appropriate fasta ID is used as infix for the file name.
                    Each generated file will be suffixed with the file extension ’.fold’.
                    If a file with the same filename already exists, any output of the program
                    will be appended to it.""", filename=True, equate=False),
            _Option(["-t", "--layout-type", "layout_type"], #layout-type=INT
                    """layout_type
                    Choose the layout algorithm. Simple radial layout if 0, or naview if 1
                    (default='1')""",equate=False),
            _Option(["--noPS", "noPS"],
                    """noPS
                    Do not produce postscript drawing of the mfe structure.
                    (default=off)""",equate=False),
            _Option(["--noconv", "noconv"],
                    """noconv
                    Do not automatically substitute nucleotide "T" with "U"
                    (default=off)"""),
            _Option(["--auto-id", "auto_id"],
                    """auto_id
                    Automatically generate an ID for each sequence.
                    (default=off)

                    The default mode of RNAfold is to automatically determine an ID from
                    the input sequence data if the input file format allows to do that.
                    Sequence IDs are usually given in the FASTA header of input sequences.
                    If this flag is active, RNAfold ignores any IDs retrieved from the input
                    and automatically generates an ID for each sequence. This ID consists
                    of a prefix and an increasing number. This flag can also be used to add
                    a FASTA header to the output even if the input has none.""",equate=False),
            _Option(["--id-prefix", "id_prefix"], #id-prefix=prefix
                    """id_prefix
                    Prefix for automatically generated IDs (as used in output file names)
                    (default=‘sequence’)

                    If this parameter is set, each sequence will be prefixed with the provided
                    string. Hence, the output files will obey the following naming scheme:
                    "prefix_xxxx_ss.ps" (secondary structure plot), "prefix_xxxx_dp.ps" (dot−plot),
                    "prefix_xxxx_dp2.ps" (stack probabilities), etc. where xxxx is the sequence
                    number.
                    Note: Setting this parameter implies −−auto−id.""",equate=False),
            _Option(["--id-digits","id_digits"], #digits=INT
                    """id_digits
                    Specify the number of digits of the counter in automatically generated
                    alignment IDs.

                    (default=‘4’)

                    When alignments IDs are automatically generated, they receive
                    an increasing number, starting with 1. This number will always
                    be left−padded by leading zeros, such that the number takes up
                    a certain width. Using this parameter, the width can be specified
                    to the users need. We allow numbers in the range [1:18].
                    This option implies auto−id.""",equate=False),
            _Option(["--id-start", "id_start"], #id-start=LONG
                    """id_start
                    Specify the first number in automatically generated alignment IDs.

                    (default=‘1’)

                    When sequence IDs are automatically generated, they receive an increasing
                    number, usually starting with 1. Using this parameter, the first number
                    can be specified to the users requirements. Note: negative numbers are
                    not allowed. Note: Setting this parameter implies to ignore any IDs
                    retrieved from the input data, i.e. it activates the auto−id flag.""",equate=False),
#Structure Constrains:
    #Command line options to interact with the structure constraints feature of this program
            _Option(["--maxBPspan", "maxBPspan"], #maxBPspan=INT
                    """maxBPspan
                    Set the maximum base pair span
                    (default='-1')""",equate=False),
            _Option(["-C", "--constraint", "constraint"], #constraint=<filename>
                    """constraint
                    Calculate structures subject to constraints.
                    (default='')

                    The program reads first the sequence, then a string containing
                    constraints on the structure encoded with the symbols:
                    . (no constraint for this base)
                    | (the corresponding base has to be paired)
                    x (the base is unpaired)
                    < (base i is paired with a base j>i)
                    > (base i is paired with a base j<i)
                    and matching brackets ( ) (base i pairs base j)

                    With the exception of "|", constraints will disallow all pairs
                    conflicting with the constraint. This is usually sufficient
                    to enforce the constraint, but occasionally a base may stay unpaired
                    in spite of constraints. PF folding ignores constraints of type "|".
                    """,equate=False, filename=True),
            _Switch(["--batch", "batch"],
                    """batch
                    Use constraints for multiple sequences.
                    (default=off)
                    Usually, constraints provided from input file only apply to
                    a single input sequence. Therefore, RNAfold will stop its computation
                    and quit after the first input sequence was processed. Using this switch,
                    RNAfold processes multiple input sequences and applies the same provided
                    constraints to each of them."""),
            _Switch(["--canonicalBPonly", "canonicalBPonly"],
                    """canonicalBPonly
                    Remove non−canonical base pairs from the structure constraint
                    (default=off"""),
            _Option(["--enforceConstraint", "enforceConstraint"],
                    """enforceConstraint
                    Enforce base pairs given by round brackets ( ) in structure constraint
                    (default=off)""",equate=False),
            _Option(["--shape", "shape"], #shape=<filename>
                    """shape
                    Use SHAPE reactivity data to guide structure predictions""",equate=False, filename=True),
            _Option(["--shapeMethod", "shapeMethod"], #shapeMethod=[D/Z/W]+[optional parameters]
                    """shapeMethod
                    Select method to incorporate SHAPE reactivity data.

                    (default='D')

                    The following methods can be used to convert SHAPE reactivities
                    into pseudo energy contributions.

                    ’D’: Convert by using a linear equation according to Deigan et al 2009.
                    The calculated pseudo energies will be applied for every nucleotide
                    involved in a stacked pair. This method is recognized by a capital ’D’
                    in the provided parameter, i.e.: shapeMethod="D" is the default setting.
                    The slope ’m’ and the intercept ’b’ can be set to a non−default value
                    if necessary, otherwise m=1.8 and b=−0.6. To alter these parameters,
                    e.g. m=1.9 and b=−0.7, use a parameter string like this:shapeMethod="Dm1.9b−0.7".
                    You may also provide only one of the two parameters like: shapeMethod="Dm1.9"
                    or shapeMethod="Db−0.7".

                    ’Z’: Convert SHAPE reactivities to pseudo energies according to Zarringhalam et al 2012.
                    SHAPE reactivities will be converted to pairing probabilities by using
                    linear mapping. Aberration from the observed pairing probabilities will
                    be penalized during the folding recursion. The magnitude of the penalties
                    can affected by adjusting the factor beta (e.g. shapeMethod="Zb0.8").

                    ’W’: Apply a given vector of perturbation energies to unpaired nucleotides
                    according to Washietl et al 2012. Perturbation vectors can be calculated by
                    using RNApvmin.""",equate=False),
            _Option(["--shapeConversion", "shapeConversion"], #shapeConversion=M/C/S/L/O + [optional parameters]
                    """shapeConversion
                    Select method to convert SHAPE reactivities to pairing probabilities.

                    (default=‘O’)

                    This parameter is useful when dealing with the SHAPE incorporation according to
                    Zarringhalam et al. The following methods can be used to convert SHAPE reactivities
                    into the probability for a certain nucleotide to be unpaired.

                    ’M’: Use linear mapping according to Zarringhalam et al.

                    ’C’: Use a cutoff−approach to divide into paired and unpaired nucleotides (e.g. "C0.25")

                    ’S’: Skip the normalizing step since the input data already represents probabilities
                    for being unpaired rather than raw reactivity values

                    ’L’: Use a linear model to convert the reactivity into a probability for being
                    unpaired (e.g. "Ls0.68i0.2" to use a slope of 0.68 and an intercept of 0.2)

                    ’O’: Use a linear model to convert the log of the reactivity into a probability
                    for being unpaired (e.g. "Os1.6i−2.29" to use a slope of 1.6 and an intercept of −2.29)
                    """,equate=False),
            _Option(["--motif", "motif"], #motif=SEQUENCE,STRUCTURE,ENERGY
                    """motif
                    Specify stabilizing effect of ligand binding to a particular
                    sequence/structure motif.

                    Some ligands binding to RNAs require and/or induce particular sequence
                    and structure motifs. For instance they bind to an interior loop,
                    or small hairpin loop. If for such cases a binding free energy is known,
                    the binding and therefore stabilizing effect of the ligand can be included
                    in the folding recursions. Interior loop motifs are specified as concatenations
                    of 5’ and 3’ motif, separated by an ’&’ character.

                    Energy contributions must be specified in kcal/mol.

                    See the manpage for an example usage of this option."""),
            _Option(["--commands", "commands"], #commands=<filename>
                    """commands
                    Read additional commands from file
                    Commands include hard and soft constraints, but also structure motifs in hairpin
                    and interior loops that need to be treeted differently. Furthermore,
                    commands can be set for unstructured and structured domains.""",equate=False, filename=True),
#Algorithms:
    #Select additional algorithms which should be included in the calculations.
    #The Minimum free energy (MFE) and a structure representative are calculated in any case.
            _Option(["-p","--partfunc", "partfunc"], #partfunc=INT
                    """partfunc
                    Calculate the partition function and base pairing probability matrix.

                    (default=‘1’)

                    In addition to the MFE structure we print a coarse representation of the pair
                    probabilities in form of a pseudo bracket notation followed by the ensemble
                    free energy. This notation makes use of the letters " . , | { } ( ) " denoting
                    bases that are essentially unpaired, weakly paired, strongly paired
                    without preference, weakly upstream (downstream) paired, or strongly
                    up− (down−)stream paired bases, respectively. On the next line the centroid
                    structure as derived from the pair probabilities together with its free energy
                    and distance to the ensemble is shown. Finally it prints the frequency
                    of the mfe structure, and the structural diversity (mean distance between
                    the structures in the ensemble). See the description of pf_fold()
                    and mean_bp_dist() and centroid() in the RNAlib documentation for details.
                    Note that unless you also specify −d2 or −d0, the partition function
                    and mfe calculations will use a slightly different energy model.
                    See the discussion of dangling end options below.

                    An additionally passed value to this option changes the behavior of partition
                    function calculation: −p0 Calculate the partition function but not the pair
                    probabilities, saving about 50%  in runtime. This prints the ensemble free
                    energy −kT ln(Z). −p2 Compute stack probabilities, i.e. the probability
                    that a pair (i,j) and the immediately interior pair (i+1,j−1) are formed
                    simultaneously in addition to pair probabilities. A second postscript
                    dot plot called "name_dp2.ps", or "dot2.ps" (if the sequence does not have
                    a name), is produced that contains pair probabilities in the upper right half
                    and stack probabilities in the lower left.""",equate=False),
            _Option(["--MEA", "MEA"], #MEA=gamma #dziala
                    """MEA
                    Calculate an MEA (maximum expected accuracy) structure, where the expected
                    accuracy is computed from the pair probabilities: each base pair (i,j) gets
                    a score 2*gamma*p_ij and the score of an unpaired base is given by
                    the probability of not forming a pair.

                    (default=‘1.’)

                    The parameter gamma tunes the importance of correctly predicted pairs versus
                    unpaired bases. Thus, for small values of gamma the MEA structure will contain
                    only pairs with very high probability. Using MEA implies p for computing
                    the pair probabilities.""",),
            _Option(["-S","--pfScale", "pfScale"], #pfScale=scaling factor
                    """pfScale
                    In the calculation of the pf use scale*mfe as an estimate for the ensemble
                    free energy (used to avoid overflows).

                    The default is 1.07, useful values are 1.0 to 1.2. Occasionally needed for long
                    sequences. You can also recompile the program to use double precision."""),
            _Switch(["-c","--circ", "circ"],
                    """circ
                    Assume a circular (instead of linear) RNA molecule.
                    (default=off)"""),
            _Switch(["--ImFeelingLucky", "ImFeelingLucky"],
                    """ImFeelingLucky
                    Return exactly one stochastically backtracked structure

                    (default=off)

                    This function computes the partition function and returns exactly one secondary
                    structure stochastically sampled from the Boltzmann equilibrium according to its
                    probability in the ensemble"""),
            _Option(["--bppmTreshold", "bppmTreshold"], #bppmTreshold=<value>
                    """bppmTreshold
                    Set the threshold for base pair probabilities included in the postscript
                    output

                    (default=‘1e−5’)

                    By setting the threshold the base pair probabilities that are included
                    in the output can be varied. By default only those exceeding 1e−5
                    in probability will be shown as squares in the dot plot.
                    Changing the threshold to any other value allows for increase or decrease
                    of data."""),
            _Option(["-g","--gquad", "gquad"],
                    """gquad
                    Incoorporate G−Quadruplex formation into the structure prediction algorithm.
                    (default=off)""",equate=False),
#Model Details
            _Option(["-T","--temp", "temp"], #temp=DOUBLE
                    """temp
                    Rescale energy parameters to a temperature of temp C.
                    Default is 37C.""",equate=False),
            #_Option(["-4","--noTetra", "noTetra"],
                    #"""noTetra
                    #Do not include special tabulated stabilizing energies for tri−,
                    #tetra− and hexaloop hairpins.
                    #(default=off)
                    #Mostly for testing."""),
            _Option(["-d","--dangles", "dangles"], #Działa
                    """dangles
                    How to treat "dangling end" energies for bases adjacent to helices
                    in free ends and multi−loops

                    (default=‘2’)

                    With −d1 only unpaired bases can participate in at most one dangling end,
                    this is the default for mfe folding but unsupported for the partition function
                    folding.

                    With −d2 this check is ignored, dangling energies will be added for the bases
                    adjacent to a helix on both sides in any case; this is the default for partition
                    function folding (partfunc). The option dangles0 ignores dangling ends altogether
                    (mostly for debugging). With dangles3 mfe folding will allow coaxial
                    stacking of adjacent helices in multi−loops. At the moment the implementation
                    will not allow coaxial stacking of the two interior pairs in a loop of degree 3
                    and works only for mfe folding.

                    Note that by default (as well as with dangles1 and dangles3) pf and mfe folding treat
                    dangling ends differently. Use dangles2 in addition to partfunc to ensure that both
                    algorithms use the same energy model.""",non_space=True),
            _Switch(["--noLP", "noLP"],
                    """noLP
                    Produce structures without lonely pairs (helices of length 1).

                    (default=off)

                    For partition function folding this only disallows pairs that can only occur isolated.
                    Other pairs may still occasionally occur as helices of length 1."""),
            _Switch(["--noGU", "noGU"],
                    """noGU
                    Do not allow GU pairs
                    (default=off)"""),
            _Switch(["--noClosingGU", "noClosingGU"],
                    """noClosingGU
                    Do not allow GU pairs at the end of helices
                    (default=off)"""),
            _Option(["-P", "--paramFile", "paramFile"], #paramFile=paramfile
                    """paramFile

                    Read energy parameters from paramfile, instead of using the default parameter set.

                    A sample parameter file should accompany your distribution. See the RNAlib
                    documentation for details on the file format.""",equate=False, filename=True),
            _Option(["--nsp", "nsp"], #nsp=STRING
                    """nsp
                    Allow other pairs in addition to the usual AU,GC,and GU pairs.

                    Its argument is a comma separated list of additionally allowed pairs.
                    If the first character is a "−" then AB will imply that AB and BA are allowed
                    pairs. e.g. RNAfold nsp −GA will allow GA and AG pairs. Nonstandard pairs
                    are given 0 stacking energy.""",equate=False),
            _Option(["-e","--energyModel", "energyModel"], #energyModel=INT
                    """energyModel
                    Rarely used option to fold sequences from the artificial ABCD... alphabet,
                    where A pairs B, C−D etc. Use the energy parameters for GC (e 1) or AU (e 2)
                    pairs.""",equate=False),
            _Option(["--betaScale", "betaScale"], #betaScale=DOUBLE
                    """betaScale

                    Set the scaling of the Boltzmann factors

                    (default=‘1.’)

                    The argument provided with this option enables to scale the thermodynamic
                    temperature used in the Boltzmann factors independently from the temperature
                    used to scale the individual energy contributions of the loop types.
                    The Boltzmann factors then become exp(−dG/(kT*betaScale)) where k
                    is the Boltzmann constant, dG the free energy contribution of the state
                    and T the absolute temperature."""),

            ]
        _RNAfoldMinimalCommandLine.__init__(self, cmd, **kwargs)

def sequence(x):
    x = x.replace(" ","\n")
    return(x.split("\n"))[1]
def dot_parenthesis(x):
    x = x.replace(" ","\n")
    return(x.split("\n"))[2]
def energy(x):
    x = x.replace(" ","\n")
    return(x.split("\n"))[3]
