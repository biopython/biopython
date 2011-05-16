"""Extract information from alignment objects.

In order to try and avoid huge alignment objects with tons of functions,
functions which return summary type information about alignments should
be put into classes in this module.

classes:
o SummaryInfo
o PSSM
"""

# standard library
import math
import sys

# biopython modules
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SubsMat import FreqTable

# Expected random distributions for 20-letter protein, and
# for 4-letter nucleotide alphabets
Protein20Random = 0.05
Nucleotide4Random = 0.25
class SummaryInfo(object):
    """Calculate summary info about the alignment.

    This class should be used to caclculate information summarizing the
    results of an alignment. This may either be straight consensus info
    or more complicated things.
    """
    def __init__(self, alignment):
        """Initialize with the alignment to calculate information on.
           ic_vector attribute. A dictionary. Keys: column numbers. Values:
        """
        self.alignment = alignment
        self.ic_vector = {}

    def dumb_consensus(self, threshold = .7, ambiguous = "X",
                       consensus_alpha = None, require_multiple = 0):
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
        o require_multiple - If set as 1, this will require that more than
        1 sequence be part of an alignment to put it in the consensus (ie.
        not just 1 sequence and gaps).
        """
        # Iddo Friedberg, 1-JUL-2004: changed ambiguous default to "X"
        consensus = ''

        # find the length of the consensus we are creating
        con_len = self.alignment.get_alignment_length()

        # go through each seq item
        for n in range(con_len):
            # keep track of the counts of the different atoms we get
            atom_dict = {}
            num_atoms = 0

            for record in self.alignment._records:
                # make sure we haven't run past the end of any sequences
                # if they are of different lengths
                if n < len(record.seq):
                    if record.seq[n] != '-' and record.seq[n] != '.':
                        if record.seq[n] not in atom_dict:
                            atom_dict[record.seq[n]] = 1
                        else:
                            atom_dict[record.seq[n]] += 1

                        num_atoms = num_atoms + 1

            max_atoms = []
            max_size = 0

            for atom in atom_dict:
                if atom_dict[atom] > max_size:
                    max_atoms = [atom]
                    max_size = atom_dict[atom]
                elif atom_dict[atom] == max_size:
                    max_atoms.append(atom)

            if require_multiple and num_atoms == 1:
                consensus += ambiguous
            elif (len(max_atoms) == 1) and ((float(max_size)/float(num_atoms))
                                         >= threshold):
                consensus += max_atoms[0]
            else:
                consensus += ambiguous

        # we need to guess a consensus alphabet if one isn't specified
        if consensus_alpha is None:
            consensus_alpha = self._guess_consensus_alphabet(ambiguous)

        return Seq(consensus, consensus_alpha)

    def gap_consensus(self, threshold = .7, ambiguous = "X",
                       consensus_alpha = None, require_multiple = 0):
        """Same as dumb_consensus(), but allows gap on the output.

        Things to do: Let the user define that with only one gap, the result
        character in consensus is gap. Let the user select gap character, now
        it takes the same is input.
        """
        # Iddo Friedberg, 1-JUL-2004: changed ambiguous default to "X"
        consensus = ''

        # find the length of the consensus we are creating
        con_len = self.alignment.get_alignment_length()

        # go through each seq item
        for n in range(con_len):
            # keep track of the counts of the different atoms we get
            atom_dict = {}
            num_atoms = 0

            for record in self.alignment._records:
                # make sure we haven't run past the end of any sequences
                # if they are of different lengths
                if n < len(record.seq):
                    if record.seq[n] not in atom_dict:
                        atom_dict[record.seq[n]] = 1
                    else:
                        atom_dict[record.seq[n]] += 1

                    num_atoms += 1

            max_atoms = []
            max_size = 0

            for atom in atom_dict:
                if atom_dict[atom] > max_size:
                    max_atoms = [atom]
                    max_size = atom_dict[atom]
                elif atom_dict[atom] == max_size:
                    max_atoms.append(atom)

            if require_multiple and num_atoms == 1:
                consensus += ambiguous
            elif (len(max_atoms) == 1) and ((float(max_size)/float(num_atoms))
                                         >= threshold):
                consensus += max_atoms[0]
            else:
                consensus += ambiguous

        # we need to guess a consensus alphabet if one isn't specified
        if consensus_alpha is None:
            #TODO - Should we make this into a Gapped alphabet?
            consensus_alpha = self._guess_consensus_alphabet(ambiguous)

        return Seq(consensus, consensus_alpha)
          
    def _guess_consensus_alphabet(self, ambiguous):
        """Pick an (ungapped) alphabet for an alignment consesus sequence.

        This just looks at the sequences we have, checks their type, and
        returns as appropriate type which seems to make sense with the
        sequences we've got.
        """
        #Start with the (un-gapped version of) the alignment alphabet
        a = Alphabet._get_base_alphabet(self.alignment._alphabet)

        #Now check its compatible with all the rest of the sequences
        for record in self.alignment:
            #Get the (un-gapped version of) the sequence's alphabet
            alt =  Alphabet._get_base_alphabet(record.seq.alphabet)
            if not isinstance(alt, a.__class__):
                raise ValueError \
                ("Alignment contains a sequence with an incompatible alphabet.")

        #Check the ambiguous character we are going to use in the consensus
        #is in the alphabet's list of valid letters (if defined).
        if hasattr(a, "letters") and a.letters is not None \
        and ambiguous not in a.letters:
            #We'll need to pick a more generic alphabet...
            if isinstance(a, IUPAC.IUPACUnambiguousDNA):
                if ambiguous in IUPAC.IUPACUnambiguousDNA().letters:
                    a = IUPAC.IUPACUnambiguousDNA()
                else:
                    a = Alphabet.generic_dna
            elif isinstance(a, IUPAC.IUPACUnambiguousRNA):
                if ambiguous in IUPAC.IUPACUnambiguousRNA().letters:
                    a = IUPAC.IUPACUnambiguousRNA()
                else:
                    a = Alphabet.generic_rna
            elif isinstance(a, IUPAC.IUPACProtein):
                if ambiguous in IUPAC.ExtendedIUPACProtein().letters:
                    a = IUPAC.ExtendedIUPACProtein()
                else:
                    a = Alphabet.generic_protein
            else:
                a = Alphabet.single_letter_alphabet
        return a

    def replacement_dictionary(self, skip_chars = []):
        """Generate a replacement dictionary to plug into a substitution matrix
        
        This should look at an alignment, and be able to generate the number
        of substitutions of different residues for each other in the
        aligned object.

        Will then return a dictionary with this information:
        {('A', 'C') : 10, ('C', 'A') : 12, ('G', 'C') : 15 ....}

        This also treats weighted sequences. The following example shows how
        we calculate the replacement dictionary. Given the following
        multiple sequence alignments:

        GTATC  0.5
        AT--C  0.8
        CTGTC  1.0

        For the first column we have:
        ('A', 'G') : 0.5 * 0.8 = 0.4
        ('C', 'G') : 0.5 * 1.0 = 0.5
        ('A', 'C') : 0.8 * 1.0 = 0.8

        We then continue this for all of the columns in the alignment, summing
        the information for each substitution in each column, until we end
        up with the replacement dictionary.

        Arguments:
        o skip_chars - A list of characters to skip when creating the dictionary.
        For instance, you might have Xs (screened stuff) or Ns, and not want
        to include the ambiguity characters in the dictionary.
        """
        # get a starting dictionary based on the alphabet of the alignment
        rep_dict, skip_items = self._get_base_replacements(skip_chars)

        # iterate through each record
        for rec_num1 in range(len(self.alignment._records)):
            # iterate through each record from one beyond the current record
            # to the end of the list of records
            for rec_num2 in range(rec_num1 + 1, len(self.alignment._records)):
                # for each pair of records, compare the sequences and add
                # the pertinent info to the dictionary
                rep_dict = self._pair_replacement(
                    self.alignment._records[rec_num1].seq,
                    self.alignment._records[rec_num2].seq,
                    self.alignment._records[rec_num1].annotations.get('weight',1.0),
                    self.alignment._records[rec_num2].annotations.get('weight',1.0),
                    rep_dict, skip_items)

        return rep_dict

    def _pair_replacement(self, seq1, seq2, weight1, weight2,
                          start_dict, ignore_chars):
        """Compare two sequences and generate info on the replacements seen.

        Arguments:
        o seq1, seq2 - The two sequences to compare.
        o weight1, weight2 - The relative weights of seq1 and seq2.
        o start_dict - The dictionary containing the starting replacement
        info that we will modify.
        o ignore_chars - A list of characters to ignore when calculating
        replacements (ie. '-').

        Returns:
        o A replacment dictionary which is modified from initial_dict with
        the information from the sequence comparison.
        """
        # loop through each residue in the sequences
        for residue_num in range(len(seq1)):
            residue1 = seq1[residue_num]
            try:
                residue2 = seq2[residue_num]
            # if seq2 is shorter, then we just stop looking at replacements
            # and return the information
            except IndexError:
                return start_dict

            # if the two residues are characters we want to count
            if (residue1 not in ignore_chars) and (residue2 not in ignore_chars):
                try:
                    # add info about the replacement to the dictionary,
                    # modified by the sequence weights
                    start_dict[(residue1, residue2)] += weight1 * weight2
                                         
                # if we get a key error, then we've got a problem with alphabets
                except KeyError:
                    raise ValueError("Residues %s, %s not found in alphabet %s"
                                     % (residue1, residue2,
                                        self.alignment._alphabet))

        return start_dict


    def _get_all_letters(self):
        """Returns a string containing the expected letters in the alignment."""
        all_letters = self.alignment._alphabet.letters
        if all_letters is None \
        or (isinstance(self.alignment._alphabet, Alphabet.Gapped) \
        and all_letters == self.alignment._alphabet.gap_char):
            #We are dealing with a generic alphabet class where the
            #letters are not defined!  We must build a list of the
            #letters used...
            set_letters = set()
            for record in self.alignment:
                #Note the built in set does not have a union_update
                #which was provided by the sets module's Set
                set_letters = set_letters.union(record.seq)
            list_letters = list(set_letters)
            list_letters.sort()
            all_letters = "".join(list_letters)
        return all_letters

    def _get_base_replacements(self, skip_items = []):
        """Get a zeroed dictonary of all possible letter combinations.

        This looks at the type of alphabet and gets the letters for it.
        It then creates a dictionary with all possible combinations of these
        letters as keys (ie. ('A', 'G')) and sets the values as zero.

        Returns:
        o The base dictionary created
        o A list of alphabet items to skip when filling the dictionary.Right
        now the only thing I can imagine in this list is gap characters, but
        maybe X's or something else might be useful later. This will also
        include any characters that are specified to be skipped.
        """
        base_dictionary = {}
        all_letters = self._get_all_letters()

        # if we have a gapped alphabet we need to find the gap character
        # and drop it out
        if isinstance(self.alignment._alphabet, Alphabet.Gapped):
            skip_items.append(self.alignment._alphabet.gap_char)
            all_letters = all_letters.replace(self.alignment._alphabet.gap_char,'')

        # now create the dictionary
        for first_letter in all_letters:
            for second_letter in all_letters:
                if (first_letter not in skip_items and
                    second_letter not in skip_items):
                    base_dictionary[(first_letter, second_letter)] = 0

        return base_dictionary, skip_items


    def pos_specific_score_matrix(self, axis_seq = None,
                                  chars_to_ignore = []):
        """Create a position specific score matrix object for the alignment.

        This creates a position specific score matrix (pssm) which is an
        alternative method to look at a consensus sequence.

        Arguments:
        o chars_to_ignore - A listing of all characters not to include in
        the pssm.  If the alignment alphabet declares a gap character,
        then it will be excluded automatically.
        o axis_seq - An optional argument specifying the sequence to
        put on the axis of the PSSM. This should be a Seq object. If nothing
        is specified, the consensus sequence, calculated with default
        parameters, will be used.

        Returns:
        o A PSSM (position specific score matrix) object.
        """
        # determine all of the letters we have to deal with
        all_letters = self._get_all_letters()
        assert all_letters

        if not isinstance(chars_to_ignore, list):
            raise TypeError("chars_to_ignore should be a list.")

        # if we have a gap char, add it to stuff to ignore
        if isinstance(self.alignment._alphabet, Alphabet.Gapped):
            chars_to_ignore.append(self.alignment._alphabet.gap_char)
        
        for char in chars_to_ignore:
            all_letters = all_letters.replace(char, '')

        if axis_seq:
            left_seq = axis_seq
            assert len(axis_seq) == self.alignment.get_alignment_length()
        else:
            left_seq = self.dumb_consensus()

        pssm_info = []
        # now start looping through all of the sequences and getting info
        for residue_num in range(len(left_seq)):
            score_dict = self._get_base_letters(all_letters)
            for record in self.alignment._records:
                try:
                    this_residue = record.seq[residue_num]
                # if we hit an index error we've run out of sequence and
                # should not add new residues
                except IndexError:
                    this_residue = None
                    
                if this_residue and this_residue not in chars_to_ignore:
                    weight = record.annotations.get('weight', 1.0)
                    try:
                        score_dict[this_residue] += weight
                    # if we get a KeyError then we have an alphabet problem
                    except KeyError:
                        raise ValueError("Residue %s not found in alphabet %s"
                                     % (this_residue,
                                        self.alignment._alphabet))

            pssm_info.append((left_seq[residue_num],
                              score_dict))


        return PSSM(pssm_info)
                    
    def _get_base_letters(self, letters):
        """Create a zeroed dictionary with all of the specified letters.
        """
        base_info = {}
        for letter in letters:
            base_info[letter] = 0

        return base_info

    def information_content(self, start = 0,
                            end = None,
                            e_freq_table = None, log_base = 2,
                            chars_to_ignore = []):
        """Calculate the information content for each residue along an alignment.

        Arguments:
        o start, end - The starting an ending points to calculate the
        information content. These points should be relative to the first
        sequence in the alignment, starting at zero (ie. even if the 'real'
        first position in the seq is 203 in the initial sequence, for
        the info content, we need to use zero). This defaults to the entire
        length of the first sequence.
        o e_freq_table - A FreqTable object specifying the expected frequencies
        for each letter in the alphabet we are using (e.g. {'G' : 0.4,
        'C' : 0.4, 'T' : 0.1, 'A' : 0.1}). Gap characters should not be
        included, since these should not have expected frequencies.
        o log_base - The base of the logathrim to use in calculating the
        information content. This defaults to 2 so the info is in bits.
        o chars_to_ignore - A listing of characterw which should be ignored
        in calculating the info content.

        Returns:
        o A number representing the info content for the specified region.

        Please see the Biopython manual for more information on how information
        content is calculated.
        """
        # if no end was specified, then we default to the end of the sequence
        if end is None:
            end = len(self.alignment._records[0].seq)
        
        if start < 0 or end > len(self.alignment._records[0].seq):
            raise ValueError \
                  ("Start (%s) and end (%s) are not in the range %s to %s"
                   % (start, end, 0, len(self.alignment._records[0].seq)))
        # determine random expected frequencies, if necessary
        random_expected = None
        if not e_freq_table:
            #TODO - What about ambiguous alphabets?
            base_alpha = Alphabet._get_base_alphabet(self.alignment._alphabet)
            if isinstance(base_alpha, Alphabet.ProteinAlphabet):
                random_expected = Protein20Random
            elif isinstance(base_alpha, Alphabet.NucleotideAlphabet):
                random_expected = Nucleotide4Random
            else:
                errstr = "Error in alphabet: not Nucleotide or Protein, "
                errstr += "supply expected frequencies"
                raise ValueError(errstr)
            del base_alpha
        elif not isinstance(e_freq_table, FreqTable.FreqTable):
            raise ValueError("e_freq_table should be a FreqTable object")
            

        # determine all of the letters we have to deal with
        all_letters = self._get_all_letters()
        for char in chars_to_ignore:
            all_letters = all_letters.replace(char, '')

        info_content = {}
        for residue_num in range(start, end):
            freq_dict = self._get_letter_freqs(residue_num,
                                               self.alignment._records,
                                               all_letters, chars_to_ignore)
            # print freq_dict,
            column_score = self._get_column_info_content(freq_dict,
                                                         e_freq_table,
                                                         log_base,
                                                         random_expected)

            info_content[residue_num] = column_score
        # sum up the score
        total_info = sum(info_content.itervalues())
        # fill in the ic_vector member: holds IC for each column
        for i in info_content:
            self.ic_vector[i] = info_content[i]
        return total_info

    def _get_letter_freqs(self, residue_num, all_records, letters, to_ignore):
        """Determine the frequency of specific letters in the alignment.

        Arguments:
        o residue_num - The number of the column we are getting frequencies
        from.
        o all_records - All of the SeqRecords in the alignment.
        o letters - The letters we are interested in getting the frequency
        for.
        o to_ignore - Letters we are specifically supposed to ignore.

        This will calculate the frequencies of each of the specified letters
        in the alignment at the given frequency, and return this as a
        dictionary where the keys are the letters and the values are the
        frequencies.
        """
        freq_info = self._get_base_letters(letters)

        total_count = 0
        # collect the count info into the dictionary for all the records
        for record in all_records:
            try:
                if record.seq[residue_num] not in to_ignore:
                    weight = record.annotations.get('weight',1.0)
                    freq_info[record.seq[residue_num]] += weight
                    total_count += weight
            # getting a key error means we've got a problem with the alphabet 
            except KeyError:
                raise ValueError("Residue %s not found in alphabet %s"
                                 % (record.seq[residue_num],
                                    self.alignment._alphabet))

        if total_count == 0:
            # This column must be entirely ignored characters
            for letter in freq_info:
                assert freq_info[letter] == 0
                #TODO - Map this to NA or NaN?
        else:
            # now convert the counts into frequencies
            for letter in freq_info:
                freq_info[letter] = freq_info[letter] / total_count

        return freq_info
            
    def _get_column_info_content(self, obs_freq, e_freq_table, log_base,
                                 random_expected):
        """Calculate the information content for a column.

        Arguments:
        o obs_freq - The frequencies observed for each letter in the column.
        o e_freq_table - An optional argument specifying the expected
        frequencies for each letter. This is a SubsMat.FreqTable instance.
        o log_base - The base of the logathrim to use in calculating the
        info content.
        """
        try:
            gap_char = self.alignment._alphabet.gap_char
        except AttributeError:
            #The alphabet doesn't declare a gap - there could be none
            #in the sequence... or just a vague alphabet.
            gap_char = "-" #Safe?
            
        if e_freq_table:
            if not isinstance(e_freq_table, FreqTable.FreqTable):
                raise ValueError("e_freq_table should be a FreqTable object")
            # check the expected freq information to make sure it is good
            for key in obs_freq:
                if (key != gap_char and key not in e_freq_table):
                    raise ValueError("Expected frequency letters %s "
                                     "do not match observed %s" \
                                     % (e_freq_table.keys(),
                                        obs_freq.keys() - [gap_char]))
        
        total_info = 0.0

        for letter in obs_freq:
            inner_log = 0.0
            # if we have expected frequencies, modify the log value by them
            # gap characters do not have expected frequencies, so they
            # should just be the observed frequency.
            if letter != gap_char:
                if e_freq_table:
                    inner_log = obs_freq[letter] / e_freq_table[letter]
                else:
                    inner_log = obs_freq[letter] / random_expected
            # if the observed frequency is zero, we don't add any info to the
            # total information content
            if inner_log > 0:
                letter_info = (obs_freq[letter] * 
                               math.log(inner_log) / math.log(log_base))
                total_info += letter_info
        return total_info 

    def get_column(self,col):
        return self.alignment.get_column(col)

class PSSM(object):
    """Represent a position specific score matrix.

    This class is meant to make it easy to access the info within a PSSM
    and also make it easy to print out the information in a nice table.

    Let's say you had an alignment like this:
    GTATC
    AT--C
    CTGTC

    The position specific score matrix (when printed) looks like:

      G A T C
    G 1 1 0 1
    T 0 0 3 0
    A 1 1 0 0
    T 0 0 2 0
    C 0 0 0 3

    You can access a single element of the PSSM using the following:

    your_pssm[sequence_number][residue_count_name]

    For instance, to get the 'T' residue for the second element in the
    above alignment you would need to do:

    your_pssm[1]['T']
    """
    def __init__(self, pssm):
        """Initialize with pssm data to represent.

        The pssm passed should be a list with the following structure:

        list[0] - The letter of the residue being represented (for instance,
        from the example above, the first few list[0]s would be GTAT...
        list[1] - A dictionary with the letter substitutions and counts.
        """
        self.pssm = pssm

    def __getitem__(self, pos):
        return self.pssm[pos][1]

    def __str__(self):
        out = " "
        all_residues = self.pssm[0][1].keys()
        all_residues.sort()
        
        # first print out the top header
        for res in all_residues:
            out += "   %s" % res
        out += "\n"

        # for each item, write out the substitutions
        for item in self.pssm:
            out += "%s " % item[0]
            for res in all_residues:
                out += " %.1f" % item[1][res]

            out += "\n"
        return out

    def get_residue(self, pos):
        """Return the residue letter at the specified position.
        """
        return self.pssm[pos][0]


def print_info_content(summary_info,fout=None,rep_record=0):
    """ Three column output: position, aa in representative sequence,
        ic_vector value"""
    fout = fout or sys.stdout
    if not summary_info.ic_vector:
        summary_info.information_content()
    rep_sequence = summary_info.alignment._records[rep_record].seq
    positions = summary_info.ic_vector.keys()
    positions.sort()
    for pos in positions:
        fout.write("%d %s %.3f\n" % (pos, rep_sequence[pos],
                   summary_info.ic_vector[pos]))

if __name__ == "__main__":
    print "Quick test"
    from Bio import AlignIO
    from Bio.Align.Generic import Alignment

    filename = "../../Tests/GFF/multi.fna"
    format = "fasta"
    expected = FreqTable.FreqTable({"A":0.25,"G":0.25,"T":0.25,"C":0.25},
                                   FreqTable.FREQ,
                                   IUPAC.unambiguous_dna)

    alignment = AlignIO.read(open(filename), format)
    for record in alignment:
        print record.seq.tostring()
    print "="*alignment.get_alignment_length()
    
    summary = SummaryInfo(alignment)
    consensus = summary.dumb_consensus(ambiguous="N")
    print consensus
    consensus = summary.gap_consensus(ambiguous="N")
    print consensus
    print
    print summary.pos_specific_score_matrix(chars_to_ignore=['-'],
                                            axis_seq=consensus)
    print
    #Have a generic alphabet, without a declared gap char, so must tell
    #provide the frequencies and chars to ignore explicitly.
    print summary.information_content(e_freq_table=expected,
                                      chars_to_ignore=['-'])
    print
    print "Trying a protein sequence with gaps and stops"

    alpha = Alphabet.HasStopCodon(Alphabet.Gapped(Alphabet.generic_protein, "-"), "*")
    a = Alignment(alpha)
    a.add_sequence("ID001", "MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-")
    a.add_sequence("ID002", "MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*")
    a.add_sequence("ID003", "MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*")
    print a
    print "="*a.get_alignment_length()

    s = SummaryInfo(a)
    c = s.dumb_consensus(ambiguous="X")
    print c
    c = s.gap_consensus(ambiguous="X")
    print c
    print
    print s.pos_specific_score_matrix(chars_to_ignore=['-', '*'], axis_seq=c)

    print s.information_content(chars_to_ignore=['-', '*'])

    
    print "Done"
