# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Generic functionality useful for all gene representations.

This module contains classes which can be used for all the different
types of patterns available for representing gene information (ie. motifs,
signatures and schemas). These are the general classes which should be
handle any of the different specific patterns.
"""
# standard library
import random

# biopython
from Bio.Alphabet import _verify_alphabet
from Bio.Seq import Seq, MutableSeq

__docformat__ = "restructuredtext en"

class PatternIO(object):
    """Allow reading and writing of patterns to files.

    This just defines a simple persistance class for patterns, making
    it easy to write them to a file and read 'em back.
    """
    def __init__(self, alphabet=None):
        """Intialize the reader and writer class.

        Arguments:

        o alphabet - An optional argument specifying the alphabet
        which patterns should follow. If an alphabet is set it'll be used
        to verify that all patterns follow it.

        Attributes:
        o separator - A character to use in separating items in a signature
        when it is written to a file and read back. This character should
        not be in the possible alphabet of the sequences, or there will
        be trouble.
        """
        self._alphabet = alphabet

        self.separator = ";"

    def write(self, pattern_list, output_handle):
        """Write a list of patterns to the given handle.
        """
        for pattern in pattern_list:
            # deal with signatures, concatentate them with the separator
            if isinstance(pattern, list) or isinstance(pattern, tuple):
                string_pattern = self.separator.join(pattern)
            # deal with the normal cases
            else:
                string_pattern = pattern

            output_handle.write("%s\n" % string_pattern)

    def write_seq(self, seq_pattern_list, output_handle):
        """Convenience function to write Seq objects to a file.

        This can take Seqs and MutableSeqs, and write them to a file
        as strings.
        """
        # convert the seq patterns into just string patterns
        all_patterns = []

        for seq_pattern in seq_pattern_list:
            if isinstance(seq_pattern, MutableSeq):
                seq = seq_pattern.toseq()
                all_patterns.append(str(seq))
            elif isinstance(seq_pattern, Seq):
                all_patterns.append(str(seq_pattern))
            else:
                raise ValueError("Unexpected pattern type %r" % seq_pattern)

        self.write(all_patterns, output_handle)

    def read(self, input_handle):
        """Read patterns from the specified handle.
        """
        all_patterns = []

        while True:
            cur_line = input_handle.readline()

            if not(cur_line):
                break

            cur_pattern = cur_line.rstrip()
            # split up signatures
            if self.separator in cur_pattern:
                cur_pattern = tuple(cur_pattern.split(self.separator))

            if self._alphabet is not None:
                # make single patterns (not signatures) into lists, so we
                # can check signatures and single patterns the same
                if not isinstance(cur_pattern, tuple):
                    test_pattern = [cur_pattern]
                else:
                    test_pattern = cur_pattern
                for pattern_item in test_pattern:
                    pattern_seq = Seq(pattern_item, self._alphabet)
                    if not(_verify_alphabet(pattern_seq)):
                        raise ValueError("Pattern %s not matching alphabet %s"
                                         % (cur_pattern, self._alphabet))

            all_patterns.append(cur_pattern)

        return all_patterns


class PatternRepository(object):
    """This holds a list of specific patterns found in sequences.

    This is designed to be a general holder for a set of patterns and
    should be subclassed for specific implementations (ie. holding Motifs
    or Signatures.
    """
    def __init__(self, pattern_info):
        """Initialize a repository with patterns,

        Arguments:

            - pattern_info - A representation of all of the patterns found in
              a finder search. This should be a dictionary, where the keys
              are patterns, and the values are the number of times a pattern is
              found.

        The patterns are represented interally as a list of two
        tuples, where the first element is the number of times a pattern
        occurs, and the second is the pattern itself. This makes it easy
        to sort the list and return the top N patterns.
        """
        self._pattern_dict = pattern_info

        # create the list representation
        self._pattern_list = []
        for pattern_name in self._pattern_dict:
            self._pattern_list.append((self._pattern_dict[pattern_name],
                                       pattern_name))

        self._pattern_list.sort()
        self._pattern_list.reverse()

    def get_all(self):
        """Retrieve all of the patterns in the repository.
        """
        patterns = []
        for pattern_info in self._pattern_list:
            patterns.append(pattern_info[1])

        return patterns

    def get_random(self, num_patterns):
        """Retrieve the specified number of patterns randomly.

        Randomly selects patterns from the list and returns them.

        Arguments:

        o num_patterns - The total number of patterns to return.
        """
        all_patterns = []

        while len(all_patterns) < num_patterns:
            # pick a pattern, and only add it if it is not already present
            new_pattern_info = random.choice(self._pattern_list)

            if new_pattern_info[1] not in all_patterns:
                all_patterns.append(new_pattern_info[1])

        return all_patterns

    def get_top_percentage(self, percent):
        """Return a percentage of the patterns.

        This returns the top 'percent' percentage of the patterns in the
        repository.
        """
        all_patterns = self.get_all()

        num_to_return = int(len(all_patterns) * percent)

        return all_patterns[:num_to_return]

    def get_top(self, num_patterns):
        """Return the specified number of most frequently occurring patterns

        Arguments:

        o num_patterns - The number of patterns to return.
        """
        all_patterns = []
        for pattern_info in self._pattern_list[:num_patterns]:
            all_patterns.append(pattern_info[1])

        return all_patterns

    def get_differing(self, top_num, bottom_num):
        """Retrieve patterns that are at the extreme ranges.

        This returns both patterns at the top of the list (ie. the same as
        returned by get_top) and at the bottom of the list. This
        is especially useful for patterns that are the differences between
        two sets of patterns.

        Arguments:

        o top_num - The number of patterns to take from the top of the list.

        o bottom_num - The number of patterns to take from the bottom of
        the list.
        """
        all_patterns = []
        # first get from the top of the list
        for pattern_info in self._pattern_list[:top_num]:
            all_patterns.append(pattern_info[1])

        # then from the bottom
        for pattern_info in self._pattern_list[-bottom_num:]:
            all_patterns.append(pattern_info[1])

        return all_patterns

    def remove_polyA(self, at_percentage=.9):
        """Remove patterns which are likely due to polyA tails from the lists.

        This is just a helper function to remove pattenrs which are likely
        just due to polyA tails, and thus are not really great motifs.
        This will also get rid of stuff like ATATAT, which might be a
        useful motif, so use at your own discretion.

        XXX Could we write a more general function, based on info content
        or something like that?

        Arguments:

        o at_percentage - The percentage of A and T residues in a pattern
        that qualifies it for being removed.
        """
        remove_list = []
        # find all of the really AT rich patterns
        for pattern_info in self._pattern_list:
            pattern_at = float(pattern_info[1].count('A') + pattern_info[1].count('T')) / len(pattern_info[1])
            if pattern_at > at_percentage:
                remove_list.append(pattern_info)

        # now remove them from the master list
        for to_remove in remove_list:
            self._pattern_list.remove(to_remove)

    def count(self, pattern):
        """Return the number of times the specified pattern is found.
        """
        try:
            return self._pattern_dict[pattern]
        except KeyError:
            return 0
