"""Find and deal with signatures in biological sequence data.

In addition to representing sequences according to motifs (see Motif.py
for more information), we can also use Signatures, which are conserved
regions that are not necessarily consecutive. This may be useful in
the case of very diverged sequences, where signatures may pick out
important conservation that can't be found by motifs (hopefully!).
"""
# biopython
from Bio.Alphabet import _verify_alphabet
from Bio.Seq import Seq

# local stuff
from Pattern import PatternRepository

class SignatureFinder(object):
    """Find Signatures in a group of sequence records.

    In this simple implementation, signatures are just defined as a
    two motifs separated by a gap. We need something a lot smarter than
    this to find more complicated signatures.
    """
    def __init__(self, alphabet_strict = 1):
        """Initialize a finder to get signatures.

        Arguments:

        o alphabet_strict - Specify whether signatures should be required
        to have all letters in the signature be consistent with the
        alphabet of the original sequence. This requires that all Seqs
        used have a consistent alphabet. This helps protect against getting
        useless signatures full of ambiguity signals.
        """
        self._alphabet_strict = alphabet_strict

    def find(self, seq_records, signature_size, max_gap):
        """Find all signatures in a group of sequences.

        Arguments:

        o seq_records - A list of SeqRecord objects we'll use the sequences
        from to find signatures.

        o signature_size - The size of each half of a signature (ie. if this
        is set at 3, then the signature could be AGC-----GAC)

        o max_gap - The maximum gap size between two parts of a signature.
        """
        sig_info = self._get_signature_dict(seq_records, signature_size,
                                            max_gap)

        return PatternRepository(sig_info)

    def _get_signature_dict(self, seq_records, sig_size, max_gap):
        """Return a dictionary with all signatures and their counts.

        This internal function does all of the hard work for the
        find_signatures function.
        """
        if self._alphabet_strict:
            alphabet = seq_records[0].seq.alphabet
        else:
            alphabet = None

        # loop through all records to find signatures
        all_sigs = {}
        for seq_record in seq_records:
            # if we are working with alphabets, make sure we are consistent
            if alphabet is not None:
                assert seq_record.seq.alphabet == alphabet, \
                       "Working with alphabet %s and got %s" % \
                       (alphabet, seq_record.seq.alphabet)

            # now start finding signatures in the sequence
            largest_sig_size = sig_size * 2 + max_gap
            for start in range(len(seq_record.seq) - (largest_sig_size - 1)):
                # find the first part of the signature
                first_sig = seq_record.seq[start:start + sig_size].tostring()

                # now find all of the second parts of the signature
                for second in range(start + 1, (start + 1) + max_gap):
                    second_sig = seq_record.seq[second: second + sig_size].tostring()

                    # if we are being alphabet strict, make sure both parts
                    # of the sig fall within the specified alphabet
                    if alphabet is not None:
                        first_seq = Seq(first_sig, alphabet)
                        second_seq = Seq(second_sig, alphabet)
                        if _verify_alphabet(first_seq) \
                        and _verify_alphabet(second_seq):
                            all_sigs = self._add_sig(all_sigs,
                                                     (first_sig, second_sig))

                    # if we are not being strict, just add the motif
                    else:
                        all_sigs = self._add_sig(all_sigs,
                                                 (first_sig, second_sig))

        return all_sigs

    def _add_sig(self, sig_dict, sig_to_add):
        """Add a signature to the given dictionary.
        """
        # incrememt the count of the signature if it is already present
        if sig_to_add in sig_dict:
            sig_dict[sig_to_add] += 1
        # otherwise add it to the dictionary
        else:
            sig_dict[sig_to_add] = 1

        return sig_dict

class SignatureCoder(object):
    """Convert a Sequence into its signature representatives.

    This takes a sequence and a set of signatures, and converts the
    sequence into a list of numbers representing the relative amounts
    each signature is seen in the sequence. This allows a sequence to
    serve as input into a neural network.
    """
    def __init__(self, signatures, max_gap):
        """Initialize with the signatures to look for.

        Arguments:

        o signatures - A complete list of signatures, in order, that
        are to be searched for in the sequences. The signatures should
        be represented as a tuple of (first part of the signature,
        second_part of the signature) -- ('GATC', 'GATC').

        o max_gap - The maximum gap we can have between the two
        elements of the signature.
        """
        self._signatures = signatures
        self._max_gap = max_gap

        # check to be sure the signatures are all the same size
        # only do this if we actually have signatures
        if len(self._signatures) > 0:
            first_sig_size = len(self._signatures[0][0])
            second_sig_size = len(self._signatures[0][1])

            assert first_sig_size == second_sig_size, \
                   "Ends of the signature do not match: %s" \
                   % self._signatures[0]

            for sig in self._signatures:
                assert len(sig[0]) == first_sig_size, \
                       "Got first part of signature %s, expected size %s" % \
                       (sig[0], first_sig_size)
                assert len(sig[1]) == second_sig_size, \
                       "Got second part of signature %s, expected size %s" % \
                       (sig[1], second_sig_size)

    def representation(self, sequence):
        """Convert a sequence into a representation of its signatures.

        Arguments:

        o sequence - A Seq object we are going to convert into a set of
        signatures.

        Returns:
        A list of relative signature representations. Each item in the
        list corresponds to the signature passed in to the initializer and
        is the number of times that the signature was found, divided by the
        total number of signatures found in the sequence.
        """
        # check to be sure we have signatures to deal with,
        # otherwise just return an empty list
        if len(self._signatures) == 0:
            return []
        
        # initialize a dictionary to hold the signature counts
        sequence_sigs = {}
        for sig in self._signatures:
            sequence_sigs[sig] = 0

        # get a list of all of the first parts of the signatures
        all_first_sigs = []
        for sig_start, sig_end in self._signatures:
            all_first_sigs.append(sig_start)
        
        # count all of the signatures we are looking for in the sequence
        sig_size = len(self._signatures[0][0])
        smallest_sig_size = sig_size * 2

        for start in range(len(sequence) - (smallest_sig_size - 1)):
            # if the first part matches any of the signatures we are looking
            # for, then expand out to look for the second part
            first_sig = sequence[start:start + sig_size].tostring()
            if first_sig in all_first_sigs:
                for second in range(start + sig_size,
                                    (start + sig_size + 1) + self._max_gap):
                    second_sig = sequence[second:second + sig_size].tostring()

                    # if we find the motif, increase the counts for it
                    if (first_sig, second_sig) in sequence_sigs:
                        sequence_sigs[(first_sig, second_sig)] += 1

        # -- normalize the signature info to go between zero and one
        min_count = min(sequence_sigs.values())
        max_count = max(sequence_sigs.values())

        # as long as we have some signatures present, normalize them
        # otherwise we'll just return 0 for everything 
        if max_count > 0:
            for sig in sequence_sigs:
                sequence_sigs[sig] = (float(sequence_sigs[sig] - min_count)
                                      / float(max_count))

        # return the relative signature info in the specified order
        sig_amounts = []
        for sig in self._signatures:
            sig_amounts.append(sequence_sigs[sig])

        return sig_amounts
        
