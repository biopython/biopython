from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet

print
print "Testing Seq"
print "==========="

s = Seq.Seq("TCAAAAGGATGCATCATG", IUPAC.unambiguous_dna)

print s.tostring()
print len(s)
print s[0]
print s[-1]
print s[3:5].tostring()

print "Reverse using -1 stride:", repr(s[::-1])

print "Extract every third nucleotide (slicing with stride 3):"
print repr(s[0::3])
print repr(s[1::3])
print repr(s[2::3])

print s.alphabet.letters

t = Seq.Seq("T", IUPAC.unambiguous_dna)
u = s + t
print str(u.alphabet)
print len(u)
assert s.tostring() + "T" == u.tostring()

t = Seq.Seq("T", IUPAC.protein)
try:
    u = s + t
except TypeError:
    print "expected error, and got it"
else:
    print "huh?  ERROR"

t = Seq.Seq("T", IUPAC.ambiguous_dna)
u = s + t
print str(u.alphabet)


from Bio.Seq import MutableSeq
import array

print
print "Testing MutableSeq"
print "=================="

print "Testing creating MutableSeqs in multiple ways"
string_seq = MutableSeq("TCAAAAGGATGCATCATG", IUPAC.ambiguous_dna)
array_seq = MutableSeq(array.array("c", "TCAAAAGGATGCATCATG"),
                       IUPAC.ambiguous_dna)
converted_seq = s.tomutable()

for test_seq in [string_seq]:
    print repr(test_seq)
    print test_seq.tostring()
    print len(test_seq)
    print repr(test_seq.toseq())

    print test_seq[0]
    print repr(test_seq[1:5])
    
    test_seq[1:3] = "GAT"
    print "Set slice with string:", repr(test_seq)
    test_seq[1:3] = test_seq[5:7]
    print "Set slice with MutableSeq:", repr(test_seq)
    test_seq[1:3] = array.array("c", "GAT")
    print "Set slice with array:", repr(test_seq)

    test_seq[3] = "G"
    print "Set item:", repr(test_seq)

    del test_seq[4:5]
    print "Delete slice:", repr(test_seq)
    del test_seq[3]
    print "Delete item:", repr(test_seq)

    test_seq.append("C")
    print "Append:", repr(test_seq)
    test_seq.insert(4, "G")
    print "Insert:", repr(test_seq)

    print "Pop off the last item:", test_seq.pop()

    test_seq.remove("G")
    print "Removed Gs:", repr(test_seq)

    try:
        test_seq.remove("Z")
        raise AssertionError("Did not get expected value error.")
    except ValueError:
        print "Expected value error and got it"

    print "A count:", test_seq.count("A")
    print "A index:", test_seq.index("A")

    test_seq.reverse()
    print "Reversed Seq:", repr(test_seq)

    print "Reverse using -1 stride:", repr(test_seq[::-1])
    

    test_seq.extend("GAT")
    test_seq.extend(MutableSeq("TTT", IUPAC.ambiguous_dna))
    print "Extended Seq:", repr(test_seq)

    del test_seq[4:6:-1]
    print "Delete stride slice:", repr(test_seq)

    print "Extract every third nucleotide (slicing with stride 3):"
    print repr(test_seq[0::3])
    print repr(test_seq[1::3])
    print repr(test_seq[2::3])
    
    print "Setting wobble codon to N (set slice with stride 3):"
    test_seq[2::3] = "N" * len(test_seq[2::3])
    print repr(test_seq)

###########################################################################
from Bio.Data.IUPACData import ambiguous_dna_complement, ambiguous_rna_complement
from Bio.Data.IUPACData import ambiguous_dna_values, ambiguous_rna_values
from sets import Set

print
print "Checking ambiguous complements"
print "=============================="

#See bug 2380, Bio.Nexus was polluting the dictionary.
assert "-" not in ambiguous_dna_values
assert "?" not in ambiguous_dna_values

def complement(sequence) :
    #TODO - Add a complement function to Bio/Seq.py?
    #There is already a complement method on the Seq and MutableSeq objects.
    return Seq.reverse_complement(sequence)[::-1]

print
print "DNA Ambiguity mapping:", ambiguous_dna_values
print "DNA Complement mapping:", ambiguous_dna_complement
for ambig_char, values in ambiguous_dna_values.iteritems() :
    compl_values = complement(values)
    print "%s={%s} --> {%s}=%s" % \
        (ambig_char, values, compl_values, ambiguous_dna_complement[ambig_char])
    assert Set(compl_values) == Set(ambiguous_dna_values[ambiguous_dna_complement[ambig_char]])
    
print
print "RNA Ambiguity mapping:", ambiguous_rna_values
print "RNA Complement mapping:", ambiguous_rna_complement
for ambig_char, values in ambiguous_rna_values.iteritems() :
    compl_values = complement(values).replace("T","U") #need to help as no alphabet
    print "%s={%s} --> {%s}=%s" % \
        (ambig_char, values, compl_values, ambiguous_rna_complement[ambig_char])
    assert Set(compl_values) == Set(ambiguous_rna_values[ambiguous_rna_complement[ambig_char]])

print
print "Reverse complements:"
for sequence in [Seq.Seq("".join(ambiguous_rna_values)),
            Seq.Seq("".join(ambiguous_dna_values)),
            Seq.Seq("".join(ambiguous_rna_values), Alphabet.generic_rna),
            Seq.Seq("".join(ambiguous_dna_values), Alphabet.generic_dna),
            Seq.Seq("".join(ambiguous_rna_values), IUPAC.IUPACAmbiguousRNA()),
            Seq.Seq("".join(ambiguous_dna_values), IUPAC.IUPACAmbiguousDNA()),
            Seq.Seq("AWGAARCKG")]:  # Note no U or T
        print "%s -> %s" \
              % (repr(sequence), repr(Seq.reverse_complement(sequence)))
        assert sequence.tostring() \
           == Seq.reverse_complement(Seq.reverse_complement(sequence)).tostring(), \
           "Dobule reverse complement didn't preserve the sequence!"
print

###########################################################################

test_seqs = [s,t,u,
             Seq.Seq("ATGAAACTG"),
             #TODO - Fix ambiguous translation
             #Seq.Seq("ATGAARCTG"),
             #Seq.Seq("AWGAARCKG"),  # Note no U or T
             #Seq.Seq("".join(ambiguous_rna_values)),
             #Seq.Seq("".join(ambiguous_dna_values)),
             #Seq.Seq("".join(ambiguous_rna_values), Alphabet.generic_rna),
             #Seq.Seq("".join(ambiguous_dna_values), Alphabet.generic_dna),
             #Seq.Seq("".join(ambiguous_rna_values), IUPAC.IUPACAmbiguousDNA()),
             #Seq.Seq("".join(ambiguous_dna_values), IUPAC.IUPACAmbiguousRNA()),
             #Seq.Seq("AWGAARCKG", Alphabet.generic_dna), 
             Seq.Seq("AUGAAACUG", Alphabet.generic_dna), 
             Seq.Seq("ATGAAACTG", Alphabet.generic_rna), 
             Seq.Seq("ATGAAACTG", Alphabet.generic_nucleotide), 
             Seq.Seq("AUGAAACTG", Alphabet.generic_nucleotide), #U and T
             Seq.MutableSeq("ATGAAACTG", Alphabet.generic_rna),
             Seq.Seq("ACTGTCGTCT", Alphabet.generic_protein)]

print
print "Transcribe DNA into RNA"
print "======================="
for nucleotide_seq in test_seqs:
    try :
        print "%s -> %s" \
        % (repr(nucleotide_seq) , repr(Seq.transcribe(nucleotide_seq)))
    except ValueError, e :
        print "%s -> %s" \
        % (repr(nucleotide_seq) , str(e))

print
print "Back-transcribe RNA into DNA"
print "============================"
for nucleotide_seq in test_seqs:
    try :
        print "%s -> %s" \
        % (repr(nucleotide_seq) , repr(Seq.transcribe(nucleotide_seq)))
    except ValueError, e :
        print "%s -> %s" \
        % (repr(nucleotide_seq) , str(e))
        
print
print "Reverse Complement"
print "=================="
for nucleotide_seq in test_seqs:
    try :
	print "%s\n-> %s" \
	% (repr(nucleotide_seq) , repr(Seq.reverse_complement(nucleotide_seq)))
    except ValueError, e :
	print "%s\n-> %s" \
	% (repr(nucleotide_seq) , str(e))
        
print
print "Translating"
print "==========="
for nucleotide_seq in test_seqs:
    try :
	print "%s\n-> %s" \
	% (repr(nucleotide_seq) , repr(Seq.translate(nucleotide_seq)))
    except ValueError, e :
	print "%s\n-> %s" \
	% (repr(nucleotide_seq) , str(e))
    
print
print "Seq's .complement() method"
print "=========================="
for nucleotide_seq in test_seqs:
    if isinstance(nucleotide_seq, Seq.Seq) :
        try :
            print "%s -> %s" \
            % (repr(nucleotide_seq) , repr(nucleotide_seq.complement()))
            assert nucleotide_seq.complement().tostring() \
                == Seq.reverse_complement(nucleotide_seq).tostring()[::-1], \
                "Bio.Seq function and method disagree!"
        except ValueError, e :
            print "%s -> %s" \
            % (repr(nucleotide_seq) , str(e))
        
print
print "Seq's .reverse_complement() method"
print "=================================="
for nucleotide_seq in test_seqs:
    if isinstance(nucleotide_seq, Seq.Seq) :
        try :
            print "%s -> %s" \
            % (repr(nucleotide_seq) , repr(nucleotide_seq.reverse_complement()))
            assert nucleotide_seq.reverse_complement().tostring() \
                == Seq.reverse_complement(nucleotide_seq).tostring(), \
                "Bio.Seq function and method disagree!"
        except ValueError, e :
            print "%s -> %s" \
            % (repr(nucleotide_seq) , str(e))
