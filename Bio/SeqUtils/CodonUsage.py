# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Methods for codon usage calculations."""


from __future__ import print_function

import math
from .CodonUsageIndices import SharpEcoliIndex
from .CodonUsageIndices import ErillLabEcoliIndex

from Bio import SeqIO  # To parse a FASTA file


CodonsDict = {
    'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
    'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
    'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
    'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
    'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
    'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
    'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
    'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
    'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
    'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
    'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
    'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}


# this dictionary shows which codons encode the same AA
SynonymousCodons = {
    'CYS': ['TGT', 'TGC'],
    'ASP': ['GAT', 'GAC'],
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'GLN': ['CAA', 'CAG'],
    'MET': ['ATG'],
    'ASN': ['AAC', 'AAT'],
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
    'LYS': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
    'PHE': ['TTT', 'TTC'],
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
    'ILE': ['ATC', 'ATA', 'ATT'],
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'HIS': ['CAT', 'CAC'],
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'TRP': ['TGG'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
    'GLU': ['GAG', 'GAA'],
    'TYR': ['TAT', 'TAC']}

# DNA bases that can occupy each codon position
CodonBases = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}

class CodonAdaptationIndex(object):
    """A codon adaptation index (CAI) implementation.

    Implements the codon adaptation index (CAI) described by Sharp and
    Li (Nucleic Acids Res. 1987 Feb 11;15(3):1281-95).

    NOTE - This implementation does not currently cope with alternative genetic
    codes: only the synonymous codons in the standard table are considered.
    """

    def __init__(self):
        """Initialize the class."""
        self.index = {}
        self.codon_count = {}

    # use this method with predefined CAI index
    def set_cai_index(self, index):
        """Set up an index to be used when calculating CAI for a gene.

        Just pass a dictionary similar to the SharpEcoliIndex in the
        CodonUsageIndices module.
        """
        self.index = index

    def generate_index(self, fasta_file):
        """Generate a codon usage index from a FASTA file of CDS sequences.

        Takes a location of a Fasta file containing CDS sequences
        (which must all have a whole number of codons) and generates a codon
        usage index.

        RCSU values
        """
        # first make sure we're not overwriting an existing index:
        if self.index != {} or self.codon_count != {}:
            raise ValueError("an index has already been set or a codon count "
                             "has been done. Cannot overwrite either.")

        # count codon occurrences in the file.
        self._count_codons(fasta_file)

        # now to calculate the index we first need to sum the number of times
        # synonymous codons were used all together.
        for aa in SynonymousCodons:
            total = 0.0
            # RCSU values are CodonCount/((1/num of synonymous codons) * sum of
            # all synonymous codons)
            rcsu = []
            codons = SynonymousCodons[aa]

            for codon in codons:
                total += self.codon_count[codon]

            # calculate the RSCU value for each of the codons
            for codon in codons:
                denominator = float(total) / len(codons)
                rcsu.append(self.codon_count[codon] / denominator)

            # now generate the index W=RCSUi/RCSUmax:
            rcsu_max = max(rcsu)
            for codon_index, codon in enumerate(codons):
                self.index[codon] = rcsu[codon_index] / rcsu_max

    def cai_for_gene(self, dna_sequence):
        """Calculate the CAI (float) for the provided DNA sequence (string).

        This method uses the Index (either the one you set or the one you
        generated) and returns the CAI for the DNA sequence.
        """
        cai_value, cai_length = 0, 0

        # if no index is set or generated, the default SharpEcoliIndex will
        # be used.
        if self.index == {}:
            self.set_cai_index(SharpEcoliIndex)

        if dna_sequence.islower():
            dna_sequence = dna_sequence.upper()

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            if codon in self.index:
                # these two codons are always one, exclude them:
                if codon not in ['ATG', 'TGG']:
                    cai_value += math.log(self.index[codon])
                    cai_length += 1
            # some indices may not include stop codons:
            elif codon not in ['TGA', 'TAA', 'TAG']:
                raise TypeError("illegal codon in sequence: %s.\n%s"
                                % (codon, self.index))

        return math.exp(cai_value / (cai_length - 1.0))

    def _count_codons(self, fasta_file):
        with open(fasta_file, 'r') as handle:

            # make the codon dictionary local
            self.codon_count = CodonsDict.copy()

            # iterate over sequence and count all the codons in the FastaFile.
            for cur_record in SeqIO.parse(handle, "fasta"):
                # make sure the sequence is lower case
                if str(cur_record.seq).islower():
                    dna_sequence = str(cur_record.seq).upper()
                else:
                    dna_sequence = str(cur_record.seq)
                for i in range(0, len(dna_sequence), 3):
                    codon = dna_sequence[i:i + 3]
                    if codon in self.codon_count:
                        self.codon_count[codon] += 1
                    else:
                        raise TypeError("illegal codon %s in gene: %s"
                                        % (codon, cur_record.id))

    def print_index(self):
        """Print out the index used.

        This just gives the index when the objects is printed.
        """
        for i in sorted(self.index):
            print("%s\t%.3f" % (i, self.index[i]))
            
            
class normRelativeCodonAdaptationIndex(object): 
    """A normalized Relative Codon Adaptation index implementation. 
       Implements the normalized Relative Codon Adaptation index (nRCA) 
       described by O'Neill, Or and Erill (PLoS One. 2013 Oct 7;8(10):e76177)
       NOTE - This implementation does not currently cope with alternative 
       genetic codes: only the synonymous codons in the standard table are 
       considered.
       
       The nRCA index works similarly to the Codon Adapation Index (CAI) from
       Sharp and Li (1987). It uses a reference set of highly-expressed genes
       provided by the user, and computes the alignment of any candidate coding
       sequence with the codon usage patterns seen in that reference set.
       A known problem of CAI is that it computes the weight (index) for each
       codon as the ratio between the frequency of that codon in the reference 
       set and the largest frequency among its synonymous codons. This 
       implicitly assumes that the background nucleotide distribution is 
       uniform. On mutationally-biased genomes (such as those found in many 
       bacterial clades), this assumption will backfire, because CAI will
       attribute to translational selection the patterns of genome-wide 
       mutational bias observed in the reference set (i.e. a weakly used codon
       will be overvalued by CAI if it is aligned with the overall GC% content
       of the reference set, and vice versa). nRCA removes this bias by 
       computing each codon index as the ration between the observed codon
       frequency in the reference set and its expected frequency (inferred as
       the product of the positional frequencies of each the codon bases).
       This has been shown to improve the correlation of nRCA with gene 
       expression values with respect to CAI (Fox and Erill, DNA Research, 
       17:3, 185-196, 2010).
       Based on the BioPython Module Bio.SeqUtils.CodonUsage code for CAI.
    """
    
    def __init__(self): 
        """Initializes the class""" 
        self.index = {} 
        self.codon_count = {} 
        self.first_pos_count = {}
        self.second_pos_count = {}
        self.third_pos_count = {}
    
    # use this method with predefined CAI index 
    def set_nrca_index(self, index): 
        """Set up an index to be used when calculating nRCA for a gene. 
            Just pass a dictionary similar to the ErillLabEcoliIndex in the 
            CodonUsageIndices module. 
        """ 
        self.index = index 
    
    def generate_index(self, fasta_file): 
        """Generate a codon usage index from a FASTA file of CDS sequences. 
           Takes a location of a FASTA file containing CDS sequences 
           (which must all have a whole number of codons) and generates the 
           nRCA codon usage index. 
        """ 
    
        # first make sure we're not overwriting an existing index: 
        if self.index != {} or self.codon_count != {}: 
            raise ValueError("an index has already been set or a codon count " 
                             "has been done. Cannot overwrite either.") 

        # count codon occurrences in the file, as well as codon position base 
        # counts
        self._codon_count(fasta_file) 
        
        # compute total number of codons
        total_codons = sum(self.codon_count.values()) + 0.0

        #add pseudocoount    
        for codon, count in self.codon_count.iteritems():
            self.codon_count[codon] += 1 / total_codons

        for base, count in self.first_pos_count.iteritems():
            self.first_pos_count[base] += 0.25 / total_codons
            self.second_pos_count[base] += 0.25 / total_codons
            self.third_pos_count[base] += 0.25 / total_codons
            
        # initialize dictionaries for frequencies
        codon_freq = CodonsDict.copy()
        first_pos_freq = CodonBases.copy()
        second_pos_freq = CodonBases.copy()
        third_pos_freq = CodonBases.copy()

        #get relative frequencies for codons and codon positions
        for cod, value in self.codon_count.iteritems():
            codon_freq[cod] = self.codon_count[cod] / total_codons   
        
        for base, value in self.first_pos_count.iteritems():
            first_pos_freq[base] = self.first_pos_count[base] / total_codons
        
        for base, value in self.second_pos_count.iteritems():
            second_pos_freq[base] = self.second_pos_count[base] / total_codons
        
        for base, value in self.third_pos_count.iteritems():
            third_pos_freq[base] = self.third_pos_count[base] / total_codons        
        
        #compute the unnormalized index value
        codon_w = CodonsDict.copy()
        for codon, w in codon_w.iteritems():
            #compute the expected frequency for each codon
            expected_freq = math.exp(math.log(first_pos_freq[codon[0]]) + \
                            math.log(second_pos_freq[codon[1]]) + \
                            math.log(third_pos_freq[codon[2]]))
            #compute the index value
            codon_w[codon] = codon_freq[codon] / expected_freq
                        
        # using the computed codon frequencies, we now normalize for amino 
        # acid usage
        # now to calculate normalized index we first need to sum the number of times 
        # synonymous codons were used all together. 
        for aa in SynonymousCodons:
            codons = SynonymousCodons[aa]
            
            # compute maximum index for codons encoding this amino acid            
            max_codon_w = 0
            for codon in codons:
                if codon_w[codon] > max_codon_w:
                    max_codon_w = codon_w[codon]
            
            # compute the normalized index
            for codon in codons:
                self.index[codon] = codon_w[codon] / max_codon_w
                        
    
    def nrca_for_gene(self, dna_sequence): 
        """Calculate the nRCA (float) for the provided DNA sequence (string). 
           This method uses the index (either the one you set or the one you 
           generated) and returns the nRCA value for the DNA sequence. 
        """ 
    
        # initialize nRCA value and length to compute geometric mean on
        nrca_value = 0
        nrca_length = 0
        
        # if no index is set or generated, the default ErillLabEcoliIndex will 
        # be used. 
        if self.index == {}: 
            self.set_nrca_index(ErillLabEcoliIndex) 
        
        # uppercase DNA sequence
        if dna_sequence.islower(): 
            dna_sequence = dna_sequence.upper() 
        
        # go through each codon in the sequence and add its contribution
        # to the geometric mean (adding in log-space to get a product)
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            if codon in self.index:
                # ATG and TGG codons are unique for their amino acids
                # stop codons (TGA, TAA and TAG) are not used by nRCA,
                # so we exclude them:
                if codon not in ['ATG', 'TGG', 'TGA', 'TAA', 'TAG']:
                    nrca_value += math.log(self.index[codon])
                    nrca_length += 1
                elif codon not in ['ATG', 'TGG', 'TGA', 'TAA', 'TAG']: 
                    raise TypeError("Illegal codon in sequence: %s.\n%s"
                    % (codon, self.index))
        
        # compute nRCA for sequence (geometric mean) in log-space
        return math.exp(nrca_value / (nrca_length - 1.0))
    
    def _codon_count(self, fasta_file): 
        with open(fasta_file, 'r') as handle: 
    
            # make the codon dictionary and codon positions local 
            self.codon_count = CodonsDict.copy()
            self.first_pos_count = CodonBases.copy()
            self.second_pos_count = CodonBases.copy()
            self.third_pos_count = CodonBases.copy()
            
            # iterate over sequence and count all the codons in the FASTA file
            for cur_record in SeqIO.parse(handle, "fasta"):
                # make sure the sequence is lower case
                if str(cur_record.seq).islower():
                    dna_sequence = str(cur_record.seq).upper()
                else: 
                    dna_sequence = str(cur_record.seq)
                    for i in range(0, len(dna_sequence), 3):
                        codon = dna_sequence[i:i + 3]
                        if codon in self.codon_count:
                            self.codon_count[codon] += 1
                            self.first_pos_count[codon[0]] += 1
                            self.second_pos_count[codon[1]] += 1
                            self.third_pos_count[codon[2]] += 1
                        else:
                            raise TypeError("Illegal codon %s in gene: %s"
                                            % (codon, cur_record.id))
       

            
    def print_index(self):
        """Print out the index used.
            This just gives the index when the objects is printed.
        """ 
        for i in sorted(self.index):
            print("%s\t%.3f" % (i, self.index[i]))
