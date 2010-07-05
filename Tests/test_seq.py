# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import sys
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet
from Bio.Data.IUPACData import ambiguous_dna_complement, ambiguous_rna_complement
from Bio.Data.IUPACData import ambiguous_dna_values, ambiguous_rna_values
from Bio.Data.CodonTable import TranslationError


if sys.version_info[0] == 3:
   array_indicator = "u"
else:
   array_indicator = "c"

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
array_seq = MutableSeq(array.array(array_indicator, "TCAAAAGGATGCATCATG"),
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
    test_seq[1:3] = array.array(array_indicator, "GAT")
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
       Seq.Seq("gtca", Alphabet.generic_dna),
       Seq.MutableSeq("GGTCA", Alphabet.generic_dna),
       Seq.Seq("CTG-CA", Alphabet.Gapped(IUPAC.unambiguous_dna, "-")),
       "TGGTCA"]
rna = [Seq.Seq("AUUUCG", IUPAC.ambiguous_rna),
       Seq.MutableSeq("AUUCG", IUPAC.ambiguous_rna),
       Seq.Seq("uCAg", Alphabet.generic_rna),
       Seq.MutableSeq("UC-AG", Alphabet.Gapped(Alphabet.generic_rna, "-")),
       Seq.Seq("U.CAG", Alphabet.Gapped(Alphabet.generic_rna, ".")),
       "UGCAU"]
nuc = [Seq.Seq("ATCG", Alphabet.generic_nucleotide),"UUUTTTACG"]
protein = [Seq.Seq("ATCGPK", IUPAC.protein),
           Seq.Seq("atcGPK", Alphabet.generic_protein),
           Seq.Seq("T.CGPK", Alphabet.Gapped(IUPAC.protein, ".")),
           Seq.Seq("T-CGPK", Alphabet.Gapped(IUPAC.protein, "-")),
           Seq.Seq("MEDG-KRXR*", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "*"), "-")),
           Seq.MutableSeq("ME-K-DRXR*XU", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "*"), "-")),
           Seq.Seq("MEDG-KRXR@", Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.extended_protein, "-"), "@")),
           Seq.Seq("ME-KR@", Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.protein, "-"), "@")),
           Seq.Seq("MEDG.KRXR@", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "@"), ".")),
           "TEDDF"]
for a in dna+rna:
    for b in nuc:
        c=a+b
        assert str(c) == str(a) + str(b)
for a in rna:
    for b in rna:
        try:
            c=a+b
            assert str(c) == str(a) + str(b)
        except ValueError, e:
            print "%s + %s\n-> %s" % (repr(a.alphabet), repr(b.alphabet), str(e))
for a in dna:
    for b in dna:
        try:
            c=a+b
            assert str(c) == str(a) + str(b)
        except ValueError, e:
            print "%s + %s\n-> %s" % (repr(a.alphabet), repr(b.alphabet), str(e))
    for b in rna:
        try:
            c=a+b
            assert (isinstance(a,str) or isinstance(b,str)), \
                   "DNA+RNA addition should fail!"
        except TypeError:
            pass
        try:
            c=b+a
            assert (isinstance(a,str) or isinstance(b,str)), \
                   "RNA+DNA addition should fail!"
        except TypeError:
            pass
for a in protein:
    for b in protein:
        try:
            c=a+b
            assert str(c) == str(a) + str(b)
        except ValueError, e:
            print "%s + %s\n-> %s" % (repr(a.alphabet), repr(b.alphabet), str(e))
    for b in nuc+dna+rna:
        try:
            c=a+b
            assert (isinstance(a,str) or isinstance(b,str)), \
                   "Protein+Nucleotide addition should fail!"
        except TypeError:
            pass
for a in nuc:
    for b in dna+rna+nuc:
        c=a+b
        assert str(c) == str(a) + str(b)
for a in dna+rna+nuc:
    for b in protein:
        try:
            c=a+b
            assert (isinstance(a,str) or isinstance(b,str)), \
                   "Nucleotide+Protein addition should fail!"
        except TypeError:
            pass

###########################################################################
print
print "Testing Seq string methods"
print "=========================="
for a in dna + rna + nuc + protein:
    if not isinstance(a, Seq.Seq) : continue
    assert a.strip().tostring() == a.tostring().strip()
    assert a.lstrip().tostring() == a.tostring().lstrip()
    assert a.rstrip().tostring() == a.tostring().rstrip()
    assert a.lower().tostring() == a.tostring().lower()
    assert a.upper().tostring() == a.tostring().upper()
    test_chars = ["-", Seq.Seq("-"), Seq.Seq("*"), "-X@"]
    alpha = Alphabet._get_base_alphabet(a.alphabet)
    if isinstance(alpha, Alphabet.DNAAlphabet):
        test_chars.append(Seq.Seq("A", IUPAC.ambiguous_dna))
    if isinstance(alpha, Alphabet.RNAAlphabet):
        test_chars.append(Seq.Seq("A", IUPAC.ambiguous_rna))
    if isinstance(alpha, Alphabet.NucleotideAlphabet):
        test_chars.append(Seq.Seq("A", Alphabet.generic_nucleotide))
    if isinstance(alpha, Alphabet.ProteinAlphabet):
        test_chars.append(Seq.Seq("K", Alphabet.generic_protein))
        test_chars.append(Seq.Seq("K-", Alphabet.Gapped(Alphabet.generic_protein,"-")))
        test_chars.append(Seq.Seq("K@", Alphabet.Gapped(IUPAC.protein,"@")))
        #Setup a clashing alphabet sequence
        b = Seq.Seq("-", Alphabet.generic_nucleotide)
    else:
        b = Seq.Seq("-", Alphabet.generic_protein)
    try:
        print a.strip(b).tostring()
        assert False, "Alphabet should have clashed!"
    except TypeError:
        pass #Good!
            
    for chars in  test_chars:
        str_chars = str(chars)
        assert a.strip(chars).tostring() == a.tostring().strip(str_chars)
        assert a.lstrip(chars).tostring() == a.tostring().lstrip(str_chars)
        assert a.rstrip(chars).tostring() == a.tostring().rstrip(str_chars)
        assert a.find(chars) == a.tostring().find(str_chars)
        assert a.find(chars,2,-2) == a.tostring().find(str_chars,2,-2)
        assert a.rfind(chars) == a.tostring().rfind(str_chars)
        assert a.rfind(chars,2,-2) == a.tostring().rfind(str_chars,2,-2)
        assert a.count(chars) == a.tostring().count(str_chars)
        assert a.count(chars,2,-2) == a.tostring().count(str_chars,2,-2)
        #Now check splits
        assert [x.tostring() for x in a.split(chars)] \
               == a.tostring().split(str(chars))
        assert [x.tostring() for x in a.rsplit(chars)] \
               == a.tostring().rsplit(str(chars))
        for max_sep in [0,1,2,999]:
            assert [x.tostring() for x in a.split(chars, max_sep)] \
                   == a.tostring().split(str(chars), max_sep)
            assert [x.tostring() for x in a.rsplit(chars, max_sep)] \
                   == a.tostring().rsplit(str(chars), max_sep)
del a, alpha, chars, str_chars, test_chars
del dna, rna, nuc, protein
###########################################################################
print
print "Checking ambiguous complements"
print "=============================="

#See bug 2380, Bio.Nexus was polluting the dictionary.
assert "-" not in ambiguous_dna_values
assert "?" not in ambiguous_dna_values

def complement(sequence):
    #TODO - Add a complement function to Bio/Seq.py?
    #There is already a complement method on the Seq and MutableSeq objects.
    return Seq.reverse_complement(sequence)[::-1]

def sorted_dict(d):
    """A sorted repr of a dictionary."""
    return "{%s}" % ", ".join("%s: %s" % (repr(k),repr(v)) \
                              for k,v in sorted(d.iteritems()))

print
print "DNA Ambiguity mapping:", sorted_dict(ambiguous_dna_values)
print "DNA Complement mapping:", sorted_dict(ambiguous_dna_complement)
for ambig_char, values in sorted(ambiguous_dna_values.iteritems()):
    compl_values = complement(values)
    print "%s={%s} --> {%s}=%s" % \
        (ambig_char, values, compl_values, ambiguous_dna_complement[ambig_char])
    assert set(compl_values) == set(ambiguous_dna_values[ambiguous_dna_complement[ambig_char]])
    
print
print "RNA Ambiguity mapping:", sorted_dict(ambiguous_rna_values)
print "RNA Complement mapping:", sorted_dict(ambiguous_rna_complement)
for ambig_char, values in sorted(ambiguous_rna_values.iteritems()):
    compl_values = complement(values).replace("T","U") #need to help as no alphabet
    print "%s={%s} --> {%s}=%s" % \
        (ambig_char, values, compl_values, ambiguous_rna_complement[ambig_char])
    assert set(compl_values) == set(ambiguous_rna_values[ambiguous_rna_complement[ambig_char]])

print
print "Reverse complements:"
for sequence in [Seq.Seq("".join(sorted(ambiguous_rna_values))),
            Seq.Seq("".join(sorted(ambiguous_dna_values))),
            Seq.Seq("".join(sorted(ambiguous_rna_values)), Alphabet.generic_rna),
            Seq.Seq("".join(sorted(ambiguous_dna_values)), Alphabet.generic_dna),
            Seq.Seq("".join(sorted(ambiguous_rna_values)).replace("X",""), IUPAC.IUPACAmbiguousRNA()),
            Seq.Seq("".join(sorted(ambiguous_dna_values)).replace("X",""), IUPAC.IUPACAmbiguousDNA()),
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
             "ATGAAACtg",
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
             Seq.Seq("AUGAAACUG", Alphabet.generic_rna), 
             Seq.Seq("ATGAAACTG", IUPAC.unambiguous_dna), 
             Seq.Seq("ATGAAA-CTG", Alphabet.Gapped(IUPAC.unambiguous_dna)),
             Seq.Seq("ATGAAACTGWN", IUPAC.ambiguous_dna), 
             Seq.Seq("AUGAAACUG", Alphabet.generic_rna), 
             Seq.Seq("AUGAAA==CUG", Alphabet.Gapped(Alphabet.generic_rna,"=")),
             Seq.Seq("AUGAAACUG", IUPAC.unambiguous_rna),
             Seq.Seq("AUGAAACUGWN", IUPAC.ambiguous_rna),
             Seq.Seq("ATGAAACTG", Alphabet.generic_nucleotide),
             Seq.Seq("AUGAAACTG", Alphabet.generic_nucleotide), #U and T
             Seq.MutableSeq("ATGAAACTG", Alphabet.generic_dna),
             Seq.MutableSeq("AUGaaaCUG", IUPAC.unambiguous_rna),
             Seq.Seq("ACTGTCGTCT", Alphabet.generic_protein)]
protein_seqs = [Seq.Seq("ATCGPK", IUPAC.protein),
                Seq.Seq("T.CGPK", Alphabet.Gapped(IUPAC.protein, ".")),
                Seq.Seq("T-CGPK", Alphabet.Gapped(IUPAC.protein, "-")),
                Seq.Seq("MEDG-KRXR*", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "*"), "-")),
                Seq.MutableSeq("ME-K-DRXR*XU", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "*"), "-")),
                Seq.Seq("MEDG-KRXR@", Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.extended_protein, "-"), "@")),
                Seq.Seq("ME-KR@", Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.protein, "-"), "@")),
                Seq.Seq("MEDG.KRXR@", Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "@"), "."))]

#Sanity test on the test sequence alphabets (see also enhancement bug 2597)
for nucleotide_seq in test_seqs:
    if hasattr(nucleotide_seq, "alphabet"):
        if "U" in str(nucleotide_seq).upper():
            assert not isinstance(nucleotide_seq.alphabet, Alphabet.DNAAlphabet)
        if "T" in str(nucleotide_seq).upper():
            assert not isinstance(nucleotide_seq.alphabet, Alphabet.RNAAlphabet)
            

print
print "Transcribe DNA into RNA"
print "======================="
for nucleotide_seq in test_seqs:
    try:
        expected = Seq.transcribe(nucleotide_seq)
        assert str(nucleotide_seq).replace("t","u").replace("T","U") == str(expected)
        print "%s -> %s" \
        % (repr(nucleotide_seq) , repr(expected))
    except ValueError, e:
        expected = None
        print "%s -> %s" \
        % (repr(nucleotide_seq) , str(e))
    #Now test the Seq object's method
    if isinstance(nucleotide_seq, Seq.Seq):
        try:
            assert repr(expected) == repr(nucleotide_seq.transcribe())
        except ValueError:
            assert expected is None

for s in protein_seqs:
    try:
        print Seq.transcribe(s)
        assert False, "Transcription shouldn't work on a protein!"
    except ValueError:
        pass
    if not isinstance(s, Seq.Seq) : continue #Only Seq has this method
    try:
        print s.transcribe()
        assert False, "Transcription shouldn't work on a protein!"
    except ValueError:
        pass

print
print "Back-transcribe RNA into DNA"
print "============================"
for nucleotide_seq in test_seqs:
    try:
        expected = Seq.back_transcribe(nucleotide_seq)
        assert str(nucleotide_seq).replace("u","t").replace("U","T") == str(expected)
        print "%s -> %s" \
        % (repr(nucleotide_seq) , repr(expected))
    except ValueError, e:
        expected = None
        print "%s -> %s" \
        % (repr(nucleotide_seq) , str(e))
    #Now test the Seq object's method
    if isinstance(nucleotide_seq, Seq.Seq):
        try:
            assert repr(expected) == repr(nucleotide_seq.back_transcribe())
        except ValueError:
            assert expected is None
            
for s in protein_seqs:
    try:
        print Seq.back_transcribe(s)
        assert False, "Back transcription shouldn't work on a protein!"
    except ValueError:
        pass
    if not isinstance(s, Seq.Seq) : continue #Only Seq has this method
    try:
        print s.back_transcribe()
        assert False, "Back transcription shouldn't work on a protein!"
    except ValueError:
        pass
        
print
print "Reverse Complement"
print "=================="
for nucleotide_seq in test_seqs:
    try:
        expected = Seq.reverse_complement(nucleotide_seq)
        print "%s\n-> %s" \
        % (repr(nucleotide_seq) , repr(expected))
    except ValueError, e:
        expected = None
        print "%s\n-> %s" \
        % (repr(nucleotide_seq) , str(e))
    #Now test the Seq object's method
    #(The MutualSeq object acts in place)
    if isinstance(nucleotide_seq, Seq.Seq):
        try:
            assert repr(expected) == repr(nucleotide_seq.reverse_complement())
            assert repr(expected[::-1]) == repr(nucleotide_seq.complement())
        except ValueError:
            assert expected is None

for s in protein_seqs:
    try:
        print Seq.reverse_complement(s)
        assert False, "Reverse complement shouldn't work on a protein!"
    except ValueError:
        pass
    #Note that these methods are "in place" for the MutableSeq:
    try:
        print s.complement()
        assert False, "Complement shouldn't work on a protein!"
    except ValueError:
        pass
    try:
        print s.reverse_complement()
        assert False, "Reverse complement shouldn't work on a protein!"
    except ValueError:
        pass
   
print
print "Translating"
print "==========="
for nucleotide_seq in test_seqs:
    try:
        expected = Seq.translate(nucleotide_seq)
        print "%s\n-> %s" \
        % (repr(nucleotide_seq) , repr(expected))
    except (ValueError, TranslationError), e:
        expected = None
        print "%s\n-> %s" \
        % (repr(nucleotide_seq) , str(e))
    #Now test the Seq object's method
    if isinstance(nucleotide_seq, Seq.Seq):
        try:
            assert repr(expected) == repr(nucleotide_seq.translate())
        except (ValueError, TranslationError):
            assert expected is None
    #Now check translate(..., to_stop=True)
    try:
        short = Seq.translate(nucleotide_seq, to_stop=True)
    except (ValueError, TranslationError), e:
        short = None
    if expected is not None:
        assert short is not None
        assert str(short) == str(expected.split("*")[0])
    if isinstance(nucleotide_seq, Seq.Seq):
        try:
            assert repr(short) == repr(nucleotide_seq.translate(to_stop=True))
        except (ValueError, TranslationError):
            assert short is None

for s in protein_seqs:
    try:
        print Seq.translate(s)
        assert False, "Translation shouldn't work on a protein!"
    except ValueError:
        pass
    if not isinstance(s, Seq.Seq) : continue #Only Seq has this method
    try:
        print s.translate()
        assert False, "Translation shouldn't work on a protein!"
    except ValueError:
        pass


misc_stops = "TAATAGTGAAGAAGG"
for nucleotide_seq in [misc_stops, Seq.Seq(misc_stops),
                       Seq.Seq(misc_stops, Alphabet.generic_nucleotide),
                       Seq.Seq(misc_stops, Alphabet.DNAAlphabet()),
                       Seq.Seq(misc_stops, IUPAC.unambiguous_dna)]:
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

for s in protein_seqs:
    try:
        print Seq.translate(s)
        assert False, "Shouldn't work on a protein!"
    except ValueError:
        pass

assert Seq.translate("TAT")=="Y"
assert Seq.translate("TAR")=="*"
assert Seq.translate("TAN")=="X"
assert Seq.translate("NNN")=="X"

assert Seq.translate("TAt")=="Y"
assert Seq.translate("TaR")=="*"
assert Seq.translate("TaN")=="X"
assert Seq.translate("nnN")=="X"

assert Seq.translate("tat")=="Y"
assert Seq.translate("tar")=="*"
assert Seq.translate("tan")=="X"
assert Seq.translate("nnn")=="X"

for codon in ["TA?", "N-N", "AC_", "Ac_"]:
    try:
        print Seq.translate(codon)
        assert "Translating %s should have failed" % repr(codon)
    except TranslationError:
        pass

ambig = set(IUPAC.IUPACAmbiguousDNA.letters)
for c1 in ambig:
    for c2 in ambig:
        for c3 in ambig:
            values = set([Seq.translate(a+b+c, table=1) \
                          for a in ambiguous_dna_values[c1] \
                          for b in ambiguous_dna_values[c2] \
                          for c in ambiguous_dna_values[c3]])
            t = Seq.translate(c1+c2+c3)
            if t=="*":
                assert values == set("*")
            elif t=="X":
                assert len(values) > 1, \
                    "translate('%s') = '%s' not '%s'" \
                    % (c1+c2+c3, t, ",".join(values))
            elif t=="Z":
                assert values == set("EQ")
            elif t=="B":
                assert values == set("DN")
            elif t=="J":
                assert values == set("LI")
            else:
                assert values == set(t)
            #TODO - Use the Bio.Data.IUPACData module for the
            #ambiguous protein mappings?
del t,c1,c2,c3,ambig

print
print "Seq's .complement() method"
print "=========================="
for nucleotide_seq in test_seqs:
    if isinstance(nucleotide_seq, Seq.Seq):
        try:
            print "%s -> %s" \
            % (repr(nucleotide_seq) , repr(nucleotide_seq.complement()))
            assert nucleotide_seq.complement().tostring() \
                == Seq.reverse_complement(nucleotide_seq).tostring()[::-1], \
                "Bio.Seq function and method disagree!"
        except ValueError, e:
            print "%s -> %s" \
            % (repr(nucleotide_seq) , str(e))
        
print
print "Seq's .reverse_complement() method"
print "=================================="
for nucleotide_seq in test_seqs:
    if isinstance(nucleotide_seq, Seq.Seq):
        try:
            print "%s -> %s" \
            % (repr(nucleotide_seq) , repr(nucleotide_seq.reverse_complement()))
            assert nucleotide_seq.reverse_complement().tostring() \
                == Seq.reverse_complement(nucleotide_seq).tostring(), \
                "Bio.Seq function and method disagree!"
        except ValueError, e:
            print "%s -> %s" \
            % (repr(nucleotide_seq) , str(e))
