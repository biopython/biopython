"""Generic.py

Contains classes to deal with generic sequence alignment stuff not
specific to a particular program or format.

classes:
o Alignment"""

# biopython
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Alphabet import IUPAC

class Alignment:
    """Represent a set of alignments.

    This is a base class to represent alignments, which should be subclassed
    to deal with an alignment in a specific format.
    """
    def __init__(self, alphabet):
        """Initialize a new Alignment object.

        Arguments:
        o alphabet - The alphabet to use for the sequence objects that are
        created. This alphabet must be a gapped type.
        """
        self._alphabet = alphabet
        # hold everything at a list of seq record objects
        self._records = []

    def get_all_seqs(self):
        """Return all of the sequences involved in the alignment.

        The return value is a list of SeqRecord objects.
        """
        return self._records 

    def get_alignment_length(self):
        """Return the maximum length of the alignment.

        All objects in the alignment should (hopefully) have the same
        length. This function will go through and find this length
        by finding the maximum length of sequences in the alignment.
        """
        max_length = 0

        for record in self._records:
            if len(record.seq) > max_length:
                max_length = len(record.seq)

        return max_length

    def add_sequence(self, descriptor, sequence):
        """Add a sequence to the alignment.

        This doesn't do any kind of alignment, it just adds in the sequence
        object, which is assumed to be prealigned with the existing
        sequences.

        Arguments:
        o descriptor - The descriptive id of the sequence being added.
        o sequence - A string with sequence info.
        """
        new_seq = Seq(sequence, self._alphabet)
        new_record = SeqRecord(new_seq, description = descriptor)
        self._records.append(new_record)

    def dumb_consensus(self, threshold = .3, ambiguous = "N",
                       consensus_alpha = None):
        """Output a fast consensus sequence of the alignment.

        This doesn't do anything fancy at all. It will just go through the
        sequence residue by residue and count up the number of each type
        of residue (ie. A or G or T or C for DNA) in all sequences in the
        alignment. If the percentage of the most common residue type is
        greater then the passed threshold, then we will add that residue type,
        otherwise an ambiguous character will be added.

        This could be made a lot fancier (ie. to take a substitution matrix
        into account), but it just meant for a quick and dirty consensus.

        Arguments:
        o threshold - The threshold value that is required to add a particular
        atom.
        o ambiguous - The ambiguous character to be added when the threshold is
        not reached.
        o consensus_alpha - The alphabet to return for the consensus sequence.
        If this is None, then we will try to guess the alphabet.
        """
        consensus = ''

        # find the length of the consensus we are creating
        con_len = self.get_alignment_length()

        # go through each seq item
        for n in range(con_len):
            # keep track of the counts of the different atoms we get
            atom_dict = {}
            num_atoms = 0

            for record in self._records:
                # make sure we haven't run past the end of any sequences
                # if they are of different lengths
                if n < len(record.seq):
                    if record.seq[n] != '-' and record.seq[n] != '.':
                        if record.seq[n] not in atom_dict.keys():
                            atom_dict[record.seq[n]] = 1
                        else:
                            atom_dict[record.seq[n]] = \
                              atom_dict[record.seq[n]] + 1

                        num_atoms = num_atoms + 1

            max_atoms = []
            max_size = 0

            for atom in atom_dict.keys():
                if atom_dict[atom] > max_size:
                    max_atoms = [atom]
                    max_size = atom_dict[atom]
                elif atom_dict[atom] == max_size:
                    max_atoms.append(atom)

            if (len(max_atoms) == 1) and ((max_size/num_atoms) >= threshold):
                consensus = consensus + max_atoms[0]
            else:
                consensus = consensus + ambiguous

        # guess the alphabet to use if one isn't specified
        if consensus_alpha is None:
            if not(isinstance(self._records[0].seq.alphabet, Alphabet.Gapped)):
                raise ValueError \
                      ("Non-gapped alphabet found in alignment object.")
        
            if isinstance(self._records[0].seq.alphabet.alphabet,
                          Alphabet.ProteinAlphabet):
                alpha = IUPAC.protein
            elif isinstance(self._records[0].seq.alphabet.alphabet,
                            Alphabet.DNAAlphabet):
                alpha = IUPAC.ambiguous_dna
            elif isinstance(self._records[0].seq.alphabet.alphabet,
                            Alphabet.RNAAlphabet):
                alpha = IUPAC.ambiguous_rna
            else:
                raise ValueError("Could not determine the type of alphabet.")
        else:
            alpha = consensus_alpha
        

        return Seq(consensus, alpha)



