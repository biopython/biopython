from Bio.Prosite import Pattern
from Bio import Seq

unchanged_patterns = (
    ("A.", ("A", "BA"), ("B", "BBP", "")),
    ("x.", ("A", "BA", "B"), ("",)),
    ("{A}.", ("B", "BBP"), ("A", "AAA", "")),
    ("[AB].", ("AB", "PA", "PPAB", "PPABP", "B"), ("C", "CCCCCCFDS")),
    ("{AB}.", ("C", "CDFDC"), ("", "A", "B", "BA", "BAB")),
    ("A-B.", ("AB", "AACAB", "ABAAA"), ("BA", "PPDAAPB", "A", "B")),
    ("A-x.", ("AB", "AACAB", "ABAAA"), ("BA", "A", "B")),
    ("[AB]-[BC].", ("AB", "AC", "BB", "BC"), ("BA", "CA", "A", "C", "PPQ")),
    ("[AB]-[BC]-[CD].", ("ABC", "QACD", "QBBCQ"), ("ABA", "A", "DCB")),
    ("[AB]-[BC]-[CDEFGHIKLMN].", ("ABC", "ACN"), ("AB", "PERPPWPTW")),
    ("{AB}-[BC]-[CD].", ("CCC", "CCD", "QCCCP"), ("ABC",)),
    ("[AB]-{BC}-[CD].", ("AAD", "APD", "QBACP"), ("ABC", "PAAEQ")),
    ("[AB]-[BC]-{CD}.", ("ABA", "ABP", "QBBAP"), ("ABC", "PAAE")),
    ("{AB}-[BC]-{CD}.", ("CCB", "PBP"), ("ABA", "PACR")),
    ("A-B-C-D.", ("QABCDQ",), ("ABCP",)),
    ("<A.", ("ABC", "A"), ("PAB", "BAAAAAA")),
    ("<{A}.", ("PAB", "BAAAAAA"), ("ABC", "A")),
    ("<[AB].", ("A", "B", "AP", "BBQ"), ("QBB", "PA")),
    ("<{AB}.", ("QBB", "PA"), ("A", "B", "AP", "BBQ")),
    ("<A-B.", ("AB", "ABB"), ("BAB", "PABQ")),
    ("<[AB]-[BC].", ("AB", "AC", "BB", "BC"), ("PAB", "QBB", "PPPPPQABBCAC")),
    ("<[AB]-[BC]-[CD].", ("ABC", "ACC"), ("CCD", "Q", "QPPAABC")),
    ("<[AB]-[BC]-[CDEFGHIKLMN].", ("ABC", "BBN"), ("QABC", "PBBN", "PPQ")),
    ("<{AB}-[BC]-[CD].", ("CCC", "QBDQ"), ("ABD", "PDD")),
    ("<[AB]-{BC}-[CD].", ("AEC", "APD"), ("CDD", "AEE", "ACC", "QAEC")),
    ("<[AB]-[BC]-{CD}.", ("ABB",), ("ABC", "QABB")),
    ("<{AB}-[BC]-{CD}.", ("CBB",), ("ABA", "QCDB")),
    ("A>.", ("A", "BA", "AA", "QQQA"), ("AB", "AAAAB")),
    ("{A}>.", ("B", "AB", "AAAQ"), ("A", "BA", "QQQQQQA")),
    ("[AB]>.", ("A", "B", "AB", "QQA"), ("Q", "AQ", "BQ", "ABQ")),
    ("{AB}>.", ("Q", "AQ", "BQ", "ABQ"), ("A", "B", "AB", "QQA")),
    ("A-B>.", ("AB", "QAB"), ("ABA", "ABQ", "BA")),
    ("[AB]-[BC]>.", ("AB", "AC", "BB", "BC"), ("ABQ", "ACQ", "Q", "QQQ")),
    ("[AB]-[BC]-[CD]>.", ("ABC", "BBD", "QABC"), ("ABCQ", "BBDQ", "QQQ")),
    ("[AB]-[BC]-[CDEFGHIKLMN]>.", ("ABN", "QQABN"), ("ABB", "BBB", "ACCD")),
    ("{AB}-[BC]-[CD]>.", ("CCC", "CCD"), ("CCDQ", "Q")),
    ("[AB]-{BC}-[CD]>.", ("AAC", "QBQD"), ("AACQ", "Q")),
    ("[AB]-[BC]-{CD}>.", ("ABB", "QABE"), ("ACBQ", "QABEQ", "P")),
    ("{AB}-[BC]-{CD}>.", ("CCA", "AAACCA"), ("CCAC", "AAACCAD")),
    ("<A>.", ("A",), ("B", "AA", "BAB", "BA", "AB")),
    ("<{A}>.", ("B", "C"), ("A", "BC", "BC", "AA")),
    ("<[AB]>.", ("A", "B"), ("C", "AB", "BA", "AA", "BB", "Q")),
    ("<{AB}>.", ("C", "Q"), ("A", "B", "AB", "AC", "BA")),
    ("<A-B>.", ("AB", ), ("AAB", "ABB", "Q", "QABQ")),
    ("<[AB]-[BC]>.", ("AB", "BB", "BC", "AC"), ("PAB", "ABQ", "QABQ")),
    ("<[AB]-[BC]-[CD]>.", ("ABC", "BBD"), ("QABC", "ABCQ", "QABCQ", "QQ")),
    ("<[AB]-[BC]-[CDEFGHIKLMN]>.", ("ABC", "ACI", "BCN"), ("QABC", "ABCQ",
                                     "QACI", "ACIQ", "QACIQ", "Q", "QQQ")),
    ("<{AB}-[BC]-[CD]>.", ("CCC", "CBC", "QCD"), ("CCCC", "QCDC", "QCDQ")),
    ("<[AB]-{BC}-[CD]>.", ("AAD", "AQC"), ("QAAD", "AADQ", "QAAD", "QQQ")),
    ("<[AB]-[BC]-{CD}>.", ("ABB", "BBB"), ("QABB", "BABQ")),
    ("<{AB}-[BC]-{CD}>.", ("CCA", "FCF"), ("QQQQ", "QCCAQ", "Q")),
    ("A-[BC>].", ("AB", "AC", "A"), ("AQ", "Q")),
    ("<A-[BC>].", ("AB", "AC", "A"), ("AA", "AQ", "Q")),
    ("A-[BC>]>.", ("AB", "AC", "A"), ("ABQ",)),
    ("[<AB]-C.", ("C", "AC", "BC", "ABC"), ("QCC", "AB")),
    ("[<AB]-[CD>].", ("", "QA", "AB", "QADQ", "QAB", "ADQ", "QQB", "QQBDQ",
                      "C", "B", "CQ"), ("Q", "QQ", "QAQ", "QC", "AQ", "QCQ")),
    ("[<ABC>].", ("", "A", "AQ", "QA", "QAQ"), ()),

    # Some ranges
    ("A(3).", ("AAA", "AAAB", "AAAAA", "BAAA"), ("AA", "AABA", "BBB", "")),
    ("x(3).", ("AAA", "BBB", "ABC"), ("AA", "B", "")),
    ("A(3)-B(3).", ("AAABBB", "BBBAAABBB"), ("ABABABAB", "AABBB")),
    ("A(2,3)-B(1,3).", ("AABB", "AAB", "AABBB", "AAAB"), ("ABBB", "QQQQQ")),
    ("A-x(1,5)-B.", ("ABB", "ABBBBBB", "ACCCCCB", "ACB"), ("ACCCCCCB", "AB")),
    ("[AB](3,4).", ("AAA", "AAAA", "ABAB", "BABA", "BBAB"), ("ABC", "CBA",
                                                             "QQQQ")),
    ("{AB}(3,5)-A(2).", ("CCCAA", "QCQAQCQAAQ"), ("ABAAA",)),
)

for pattern, good_list, bad_list in unchanged_patterns:
    p = Pattern.Prosite(pattern = pattern)
    print "Patterns:", repr(pattern), repr(p.re.pattern), \
          repr(p.grouped_re.pattern)
    assert p.tostring() == pattern, (p.tostring(), pattern)
    for x in good_list:
        if p.re.search(x) is None:
            print "\t", "should be good re, but failed:", repr(x)
    for x in bad_list:
        if p.re.search(x) is not None:
            print "\t", "should be bad re, but succeeded:", repr(x)
    for x in good_list:
        m = p.search(Seq.Seq(x))
        if m is None:
            print "\t", "should be good group re, but failed:", repr(x)
        else:
            #print "x=", repr(x)
            mapped_pattern = m.mapped_pattern()
            offset = m.start()
            for i in range(len(mapped_pattern)):
                mpi = mapped_pattern[i]
                #print "***", str(x[offset+i]), str(mpi)
                if mpi.ignore:
                    assert x[offset+i] not in mpi.letters, repr(mpi.letters)
                else:
                    if mpi.letters != "x":
                        assert x[offset+i] in mpi.letters, repr(mpi.letters)
                
    for x in bad_list:
        if p.search(Seq.Seq(x)) is not None:
            print "\t", "should be bad group re, but succeeded:", repr(x)

    
bad_patterns = (
    "A",
    "A<",
    ">A",
    "A<.",
    ">A.",
    "[AB]<.",
    ">[AB].",
    "AB.",
    "A-B<.",
    "A(B).",
    "A.B",
    "A.B.",
    "A-B",
    "[A-B]",
    "{[A-B]}.",
    "[{AB}].",
    "(A).",
    "[A>]-B.",
    "A-[<B].",
    "A(3,).",
    "A(-1).",
    #"A(5,2).",  # too complicated to check via re
    "A(2, 5).",
    "A(B,C).",
    "A-[<G].",
    "[G>]-A.",
    "1.",
    "[1].",
    "[A1B].",
    #"[AA].",   # too complicated to check via re
    "{1}.",
    "[<].",
    "[>].",
    "[A<B].",
    "[A>B].",
    )
for pat in bad_patterns:
    print pat, "-",
    try:
        Pattern.compile(pat)
    except TypeError:
        print "correctly failed"
    else:
        print "INCORRECTLY PASSED"
