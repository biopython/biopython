
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
