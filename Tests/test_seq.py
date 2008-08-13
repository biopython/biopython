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
print
print "Testing Seq addition"
print "===================="
dna = [Seq.Seq("ATCG", IUPAC.ambiguous_dna),
       Seq.Seq("GTCA", Alphabet.generic_dna),
       Seq.MutableSeq("GGTCA", Alphabet.generic_dna),
       Seq.Seq("CTG-CA", Alphabet.Gapped(IUPAC.unambiguous_dna, "-")),
       "TGGTCA"]
rna = [Seq.Seq("AUUUCG", IUPAC.ambiguous_rna),
       Seq.MutableSeq("AUUCG", IUPAC.ambiguous_rna),
       Seq.Seq("UCAG", Alphabet.generic_rna),
       Seq.MutableSeq("UC-AG", Alphabet.Gapped(Alphabet.generic_rna, "-")),
       Seq.Seq("U.CAG", Alphabet.Gapped(Alphabet.generic_rna, ".")),
       "UGCAU"]
nuc = [Seq.Seq("ATCG", Alphabet.generic_nucleotide),"UUUTTTACG"]
protein = [Seq.Seq("ATCGPK", IUPAC.protein),
           Seq.Seq("T.CGPK", Alphabet.Gapped(IUPAC.protein, ".")),
           Seq.Seq("T-CGPK", Alphabet.Gapped(IUPAC.protein, "-")),
           Seq.Seq("MEDG-KRXR*", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "*"), "-")),
           Seq.MutableSeq("ME-K-DRXR*XU", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "*"), "-")),
           Seq.Seq("MEDG-KRXR@", Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.extended_protein, "-"), "*")),
           Seq.Seq("ME-KR@", Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.protein, "-"), "@")),
           Seq.Seq("MEDG.KRXR@", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "@"), ".")),
           "TEDDF"]
for a in dna+rna :
    for b in nuc :
        c=a+b
        assert str(c) == str(a) + str(b)
for a in rna :
    for b in rna :
        try :
            c=a+b
            assert str(c) == str(a) + str(b)
        except ValueError, e :
            print "%s + %s\n-> %s" % (repr(a.alphabet), repr(b.alphabet), str(e))
for a in dna :
    for b in dna :
        try :
            c=a+b
            assert str(c) == str(a) + str(b)
        except ValueError, e :
            print "%s + %s\n-> %s" % (repr(a.alphabet), repr(b.alphabet), str(e))
    for b in rna :
        try :
            c=a+b
            assert (isinstance(a,str) or isinstance(b,str)), \
                   "DNA+RNA addition should fail!"
        except TypeError :
            pass
        try :
            c=b+a
            assert (isinstance(a,str) or isinstance(b,str)), \
                   "RNA+DNA addition should fail!"
        except TypeError :
            pass
for a in protein :
    for b in protein :
        try :
            c=a+b
            assert str(c) == str(a) + str(b)
        except ValueError, e:
            print "%s + %s\n-> %s" % (repr(a.alphabet), repr(b.alphabet), str(e))
    for b in nuc+dna+rna :
        try :
            c=a+b
            assert (isinstance(a,str) or isinstance(b,str)), \
                   "Protein+Nucleotide addition should fail!"
        except TypeError :
            pass
for a in nuc :
    for b in dna+rna+nuc :
        c=a+b
        assert str(c) == str(a) + str(b)
for a in dna+rna+nuc :
    for b in protein :
        try :
            c=a+b
            assert (isinstance(a,str) or isinstance(b,str)), \
                   "Nucleotide+Protein addition should fail!"
        except TypeError :
            pass
del dna, rna, nuc, protein

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
             Seq.Seq("ATGAAACUG", IUPAC.unambiguous_dna), 
             Seq.Seq("ATGAAACUGWN", IUPAC.ambiguous_dna), 
             Seq.Seq("ATGAAACTG", Alphabet.generic_rna), 
             Seq.Seq("ATGAAACTG", IUPAC.unambiguous_rna), 
             Seq.Seq("ATGAAACTGWN", IUPAC.ambiguous_rna), 
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
        % (repr(nucleotide_seq) , repr(Seq.back_transcribe(nucleotide_seq)))
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

misc_stops = "TAATAGTGAAGAAGG"
for nucleotide_seq in [misc_stops, Seq.Seq(misc_stops),
                       Seq.Seq(misc_stops, Alphabet.generic_nucleotide),
                       Seq.Seq(misc_stops, Alphabet.DNAAlphabet()),
                       Seq.Seq(misc_stops, IUPAC.unambiguous_dna)] :
    assert "***RR" == str(Seq.translate(nucleotide_seq))
    assert "***RR" == str(Seq.translate(nucleotide_seq, table=1))
    assert "***RR" == str(Seq.translate(nucleotide_seq, table="SGC0"))
    assert "**W**" == str(Seq.translate(nucleotide_seq, table=2))
    assert "**WRR" == str(Seq.translate(nucleotide_seq, \
                                        table='Yeast Mitochondrial'))
    assert "**WSS" == str(Seq.translate(nucleotide_seq, table=5))
    assert "**WSS" == str(Seq.translate(nucleotide_seq, table=9))
    assert "**CRR" == str(Seq.translate(nucleotide_seq, \
                                        table='Euplotid Nuclear'))
    assert "***RR" == str(Seq.translate(nucleotide_seq, table=11))
    assert "***RR" == str(Seq.translate(nucleotide_seq, table='Bacterial'))
del misc_stops

ambig = Set(IUPAC.IUPACAmbiguousDNA.letters)
for c1 in ambig :
    for c2 in ambig :
        for c3 in ambig :
            values = Set([Seq.translate(a+b+c, table=1) \
                          for a in ambiguous_dna_values[c1] \
                          for b in ambiguous_dna_values[c2] \
                          for c in ambiguous_dna_values[c3]])
            if "*" in values and len(values) > 1 :
                #Translation is expected to fail.
                #TODO - Confirm it fails with a try/except
                continue
            t = Seq.translate(c1+c2+c3)
            if t=="*" :
                assert values == Set(["*"])
            elif t=="X" :
                assert len(values) > 1, \
                    "translate('%s') = '%s' not '%s'" \
                    % (c1+c2+c3, t, ",".join(values))
            elif t=="Z" :
                assert values == Set(("E", "Q"))
            elif t=="B" :
                assert values == Set(["D", "N"])
            elif t=="J" :
                assert values == Set(["L", "I"])
            else :
                assert values == Set(t)
            #TODO - Use the Bio.Data.IUPACData module for the
            #ambiguous protein mappings?
del t,c1,c2,c3,ambig

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
