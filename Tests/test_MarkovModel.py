# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

try:
    from numpy import asarray
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(\
        "Install NumPy if you want to use Bio.MarkovModel.")

from Bio import MarkovModel

def print_mm(markov_model):
    print "STATES: %s" % ' '.join(markov_model.states)
    print "ALPHABET: %s" % ' '.join(markov_model.alphabet)
    print "INITIAL:"
    for i in range(len(markov_model.p_initial)):
        print "  %s: %.2f" % (
            markov_model.states[i], markov_model.p_initial[i])
    print "TRANSITION:"
    for i in range(len(markov_model.p_transition)):
        x = ["%.2f" % x for x in markov_model.p_transition[i]]
        print "  %s: %s" % (markov_model.states[i], ' '.join(x))
    print "EMISSION:"
    for i in range(len(markov_model.p_emission)):
        x = ["%.2f" % x for x in markov_model.p_emission[i]]
        print "  %s: %s" % (markov_model.states[i], ' '.join(x))



print "TESTING train_visible"
states = ["0", "1", "2", "3"]
alphabet = ["A", "C", "G", "T"]
training_data = [
    ("AACCCGGGTTTTTTT", "001112223333333"),
    ("ACCGTTTTTTT", "01123333333"),
    ("ACGGGTTTTTT", "01222333333"),
    ("ACCGTTTTTTTT", "011233333333"),
    ]
print "Training HMM"
mm = MarkovModel.train_visible(states, alphabet, training_data)
print "Classifying"

#print MarkovModel.find_states(mm, "AACGTT")
#Don't just print this, as the float may have different
#precision on different platforms.  This returns a list
#containing a tuple containing a list (fine), and a float.
states = MarkovModel.find_states(mm, "AACGTT")
for state_list, state_float in states:
    print "State %s, %0.10f" % (repr(state_list), state_float)
print_mm(mm)




print "TESTING baum welch"
states = ["CP", "IP"]
alphabet = ["cola", "ice_t", "lem"]
outputs = [
    (2, 1, 0)
    ]
print "Training HMM"
p_initial = [1.0, 0.0000001]
p_transition = [[0.7, 0.3],
                [0.5, 0.5]]
p_emission = [[0.6, 0.1, 0.3],
              [0.1, 0.7, 0.2]]
N, M = len(states), len(alphabet)
x = MarkovModel._baum_welch(N, M, outputs,
                            p_initial=p_initial,
                            p_transition=p_transition,
                            p_emission=p_emission
                            )
p_initial, p_transition, p_emission = x
mm = MarkovModel.MarkovModel(states, alphabet,
                             p_initial, p_transition, p_emission)
print_mm(mm)


# Test Baum-Welch.  This is hard because it is a non-deterministic
# algorithm.  Each run will result in different states having to
# different emissions.  In order to help this, we need to specify some
# initial probabilities to bias the final results.  This is not
# implemented yet in the MarkovModel module.

## states = [
##     "state0",
##     "state1",
##     "state2",
##     "state3",
##     ]
## alphabet = ["a", "c", "g", "t"]

## training_data = [
##     "aacccgggttttttt",
##     "accgttttttt",
##     "acgggtttttt",
##     "accgtttttttt",
##     "aaccgtttttttt",
##     "aacggttttt",
##     "acccggttttt",
##     "acccggggttt",
##     "aacccggggtttt",
##     "aaccgggtttttt"
##     ]
## print "TRAINING HMM"
## ep = ones((len(states), len(alphabet)))
## hmm = MarkovModel.train_bw(states, alphabet, training_data,
##                            pseudo_emission=ep)
## print "CLASSIFYING"
## states = MarkovModel.find_states(hmm, "aacgtt")
## print states
## print "STATES: %s" % ' '.join(hmm.states)
## print "ALPHABET: %s" % ' '.join(hmm.alphabet)
## print "INITIAL:"
## for i in range(len(hmm.p_initial)):
##     print "  %s: %.2f" % (hmm.states[i], hmm.p_initial[i])
## print "TRANSITION:"
## for i in range(len(hmm.p_transition)):
##     x = ["%.2f" % x for x in hmm.p_transition[i]]
##     print "  %s: %s" % (hmm.states[i], ' '.join(x))
## print "EMISSION:"
## for i in range(len(hmm.p_emission)):
##     x = ["%.2f" % x for x in hmm.p_emission[i]]
##     print "  %s: %s" % (hmm.states[i], ' '.join(x))




# Do some tests from the topcoder competition.

class DNAStrand:
    def mostLikely(self, normal, island, dnastrand):
        states = "NR"
        alphabet = "AGTC"

        normal = [float(x)/100 for x in normal]
        island = [float(x)/100 for x in island]
        
        p_initial = [1.0, 0.0]
        p_initial = asarray(p_initial)

        p_transition = []
        p_transition.append([1.0-normal[-1], normal[-1]])
        p_transition.append([island[-1], 1.0-island[-1]])
        p_transition = asarray(p_transition)
        
        p_emission = []   # 2x4 matrix
        p_emission.append(normal[:4])
        p_emission.append(island[:4])
        p_emission = asarray(p_emission)

        mm = MarkovModel.MarkovModel(
            states, alphabet, p_initial, p_transition, p_emission)

        x = MarkovModel.find_states(mm, dnastrand)
        states, x = x[0]
        return ''.join(states)

ds = DNAStrand()


# NNNN
print ds.mostLikely([30, 20, 30, 20, 10],
                    [10, 40, 10, 40, 20],
                    "TGCC"
                    )

# NNNRRRNNRRNRRN
print ds.mostLikely([4, 14, 62, 20, 44],
                    [39, 15, 4, 42, 25],
                    "CCTGAGTTAGTCGT"
                    )

# NRRRRRRRRRRRNNNNRRRRRRRRR
print ds.mostLikely([45, 36, 6, 13, 25],
                    [24, 18, 12, 46, 25],
                    "CCGTACTTACCCAGGACCGCAGTCC"
                    )

# NRRRRRRRRRR
print ds.mostLikely([75,3,1,21,45],
                    [34,11,39,16,15],
                    "TTAGCAGTGCG"
                    )

# N
print ds.mostLikely([26,37,8,29,16],
                    [31,13,33,23,25],
                    "T"
                    )
