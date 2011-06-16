"""Find and deal with motifs in biological sequence data.

Representing DNA (or RNA or proteins) in a neural network can be difficult
since input sequences can have different lengths. One way to get around
this problem is to deal with sequences by finding common motifs, and counting
the number of times those motifs occur in a sequence. This information can
then be used for creating the neural networks, with occurances of motifs
going into the network instead of raw sequence data.
"""
# biopython
from Bio.Alphabet import _verify_alphabet
from Bio.Seq import Seq

# local modules
from Pattern import PatternRepository

class MotifFinder(object):
    """Find motifs in a set of Sequence Records.
    """
    def __init__(self, alphabet_strict = 1):
        """Initialize a finder to get motifs.

        Arguments:

        o alphabet_strict - Whether or not motifs should be
        restricted to having all of there elements within the alphabet
        of the sequences. This requires that the Sequences have a real
        alphabet, and that all sequences have the same alphabet.
        """
        self.alphabet_strict = alphabet_strict

    def find(self, seq_records, motif_size):
        """Find all motifs of the given size in the passed SeqRecords.

        Arguments:

        o seq_records - A list of SeqRecord objects which the motifs
        will be found from.

        o motif_size - The size of the motifs we want to look for.

        Returns:
        A PatternRepository object that contains all of the motifs (and their
        counts) found in the training sequences).
        """
        motif_info = self._get_motif_dict(seq_records, motif_size)

        return PatternRepository(motif_info)

    def _get_motif_dict(self, seq_records, motif_size):
        """Return a dictionary with information on motifs.

        This internal function essentially does all of the hard work for
        finding motifs, and returns a dictionary containing the found motifs
        and their counts. This is internal so it can be reused by
        find_motif_differences.
        """
        if self.alphabet_strict:
            alphabet = seq_records[0].seq.alphabet
        else:
            alphabet = None

        # loop through all records to find the motifs in the sequences
        all_motifs = {}
        for seq_record in seq_records:
            # if we are working with alphabets, make sure we are consistent
            if alphabet is not None:
                assert seq_record.seq.alphabet == alphabet, \
                       "Working with alphabet %s and got %s" % \
                       (alphabet, seq_record.seq.alphabet)

            # now start finding motifs in the sequence
            for start in range(len(seq_record.seq) - (motif_size - 1)):
                motif = seq_record.seq[start:start + motif_size].tostring()

                # if we are being alphabet strict, make sure the motif
                # falls within the specified alphabet
                if alphabet is not None:
                    motif_seq = Seq(motif, alphabet)
                    if _verify_alphabet(motif_seq):
                        all_motifs = self._add_motif(all_motifs, motif)

                # if we are not being strict, just add the motif
                else:
                    all_motifs = self._add_motif(all_motifs, motif)

        return all_motifs

    def find_differences(self, first_records, second_records, motif_size):
        """Find motifs in two sets of records and return the differences.

        This is used for finding motifs, but instead of just counting up all
        of the motifs in a set of records, this returns the differences
        between two listings of seq_records.

        o first_records, second_records - Two listings of SeqRecord objects
        to have their motifs compared.

        o motif_size - The size of the motifs we are looking for.

        Returns:
        A PatternRepository object that has motifs, but instead of their
        raw counts, this has the counts in the first set of records
        subtracted from the counts in the second set.
        """
        first_motifs = self._get_motif_dict(first_records, motif_size)
        second_motifs = self._get_motif_dict(second_records, motif_size)

        motif_diffs = {}

        # first deal with all of the keys from the first motif
        for cur_key in first_motifs:
            if cur_key in second_motifs:
                motif_diffs[cur_key] = first_motifs[cur_key] - \
                                       second_motifs[cur_key]
            else:
                motif_diffs[cur_key] = first_motifs[cur_key]

        # now see if there are any keys from the second motif
        # that we haven't got yet.
        missing_motifs = list(second_motifs)

        # remove all of the motifs we've already added
        for added_motif in motif_diffs:
            if added_motif in missing_motifs:
                missing_motifs.remove(added_motif)

        # now put in all of the motifs we didn't get
        for cur_key in missing_motifs:
            motif_diffs[cur_key] = 0 - second_motifs[cur_key]

        return PatternRepository(motif_diffs)
                
    def _add_motif(self, motif_dict, motif_to_add):
        """Add a motif to the given dictionary.
        """
        # incrememt the count of the motif if it is already present
        if motif_to_add in motif_dict:
            motif_dict[motif_to_add] += 1
        # otherwise add it to the dictionary
        else:
            motif_dict[motif_to_add] = 1

        return motif_dict
    
class MotifCoder(object):
    """Convert motifs and a sequence into neural network representations.

    This is designed to convert a sequence into a representation that
    can be fed as an input into a neural network. It does this by
    representing a sequence based the motifs present.
    """
    def __init__(self, motifs):
        """Initialize an input producer with motifs to look for.

        Arguments:

        o motifs - A complete list of motifs, in order, that are to be
        searched for in a sequence.
        """
        self._motifs = motifs

        # check to be sure the motifs make sense (all the same size)
        self._motif_size = len(self._motifs[0])
        for motif in self._motifs:
            if len(motif) != self._motif_size:
                raise ValueError("Motif %s given, expected motif size %s"
                                 % (motif, self._motif_size))

    def representation(self, sequence):
        """Represent a sequence as a set of motifs.

        Arguments:

        o sequence - A Bio.Seq object to represent as a motif.

        This converts a sequence into a representation based on the motifs.
        The representation is returned as a list of the relative amount of
        each motif (number of times a motif occured divided by the total
        number of motifs in the sequence). The values in the list correspond
        to the input order of the motifs specified in the initializer.
        """
        # initialize a dictionary to hold the motifs in this sequence
        seq_motifs = {}
        for motif in self._motifs:
            seq_motifs[motif] = 0
        
        # count all of the motifs we are looking for in the sequence
        for start in range(len(sequence) - (self._motif_size - 1)):
            motif = sequence[start:start + self._motif_size].tostring()

            if motif in seq_motifs:
                seq_motifs[motif] += 1

        # normalize the motifs to go between zero and one
        min_count = min(seq_motifs.values())
        max_count = max(seq_motifs.values())

        # as long as we have some motifs present, normalize them
        # otherwise we'll just return 0 for everything 
        if max_count > 0:
            for motif in seq_motifs.keys():
                seq_motifs[motif] = (float(seq_motifs[motif] - min_count)
                                     / float(max_count))

        # return the relative motif counts in the specified order
        motif_amounts = []
        for motif in self._motifs:
            motif_amounts.append(seq_motifs[motif])

        return motif_amounts
        
