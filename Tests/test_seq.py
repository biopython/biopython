
from Bio import Seq
from Bio.Alphabet import IUPAC
s = Seq.Seq("TCAAAAGGATGCATCATG", IUPAC.unambiguous_dna)

print s.tostring()
print len(s)
print s[0]
print s[-1]
print s[3:5].tostring()

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

print "\nTesting MutableSeq"

print "Testing creating MutableSeqs in multiple ways"
string_seq = MutableSeq("TCAAAAGGATGCATCATG", IUPAC.ambiguous_dna)
array_seq = MutableSeq(array.array("c", "TCAAAAGGATGCATCATG"),
                       IUPAC.ambiguous_dna)
converted_seq = s.tomutable()

for test_seq in [string_seq]:
    print test_seq
    print test_seq.tostring()
    print len(test_seq)
    print test_seq.toseq()

    print test_seq[0]
    print test_seq[1:5]

    test_seq[1:3] = "GAT"
    print "Set slice with string:", test_seq
    test_seq[1:3] = test_seq[5:7]
    print "Set slice with MutableSeq:", test_seq
    test_seq[1:3] = array.array("c", "GAT")
    print "Set slice with array:", test_seq

    test_seq[3] = "G"
    print "Set item:", test_seq

    del test_seq[4:5]
    print "Delete slice:", test_seq
    del test_seq[3]
    print "Delete item:", test_seq

    test_seq.append("C")
    print "Append:", test_seq
    test_seq.insert(4, "G")
    print "Insert:", test_seq

    print "Pop off the last item:", test_seq.pop()

    test_seq.remove("G")
    print "Removed Gs:", test_seq

    try:
        test_seq.remove("Z")
        raise AssertionError("Did not get expected value error.")
    except ValueError:
        print "Expected value error and got it"

    print "A count:", test_seq.count("A")
    print "A index:", test_seq.index("A")

    test_seq.reverse()
    print "Reversed Seq:", test_seq

    test_seq.extend("GAT")
    test_seq.extend(MutableSeq("TTT", IUPAC.ambiguous_dna))
    print "Extended Seq:", test_seq
