import string
import math

STRONG_BONDS = ["G", "C"]
WEAK_BONDS = ["A", "T", "U"]

def Tm_staluc(s,dnac=50,saltc=50,rna=0):
    """Returns DNA/DNA tm using nearest neighbor thermodynamics. dnac is
    DNA concentration [nM] and saltc is salt concentration [mM].
    rna=0 is for DNA/DNA (default), for RNA, rna should be 1.
    Sebastian Bassi <sbassi@genesdigitales.com>"""
    
    #Credits: 
    #Main author: Sebastian Bassi <sbassi@genesdigitales.com>
    #Overcount function: Greg Singer <singerg@tcd.ie>
    #Based on the work of Nicolas Le Novere <lenov@ebi.ac.uk> 
    #Bioinformatics. 17:1226-1227(2001)

    #This function returns better results than EMBOSS DAN because it uses 
    #updated thermodinamics values and take into account inicialization 
    #parameters from SantaLucia works (1998).
    
    #Things to do:
    #+Add a function to detect complementary sequences. Change K according 
    # to result.
    #+Add support for heteroduplex (see Sugimoto et al. 1995).
    #+Correction for Mg2+. Now supports only monovalent ions.
    #+Put thermodinamics table in a external file for users to change at will
    #+Add support for danglings ends (see Le Novele. 2001) and mismatches.
    
    dh=0 #DeltaH. Enthalpy
    ds=0 #deltaS Entropy

    def tercorr(stri):
        deltah=0
        deltas=0
        if rna==0:
            #DNA/DNA
            #Allawi and SantaLucia (1997). Biochemistry 36 : 10581-10594
	    if stri[0] in STRONG_BONDS:
                deltah=deltah-0.1
                deltas=deltas+2.8
            elif stri[0] in WEAK_BONDS:
                deltah=deltah-2.3
                deltas=deltas-4.1
            if stri[0] in STRONG_BONDS:
		deltah=deltah-0.1
                deltas=deltas+2.8
            elif stri[0] in WEAK_BONDS:
                deltah=deltah-2.3
                deltas=deltas-4.1
            dhL=dh+deltah
            dsL=ds+deltas
            return dsL,dhL
        elif rna==1:
            #RNA/RNA
            if stri[0] in STRONG_BONDS:
                deltah=deltah-3.61
                deltas=deltas-1.5
            elif stri[0] in WEAK_BONDS:
                deltah=deltah-3.72
                deltas=deltas+10.5
            if stri[0] in STRONG_BONDS:
                deltah=deltah-3.61
                deltas=deltas-1.5
            elif stri[0] in WEAK_BONDS:
                deltah=deltah-3.72
                deltas=deltas+10.5
            dhL=dh+deltah
            dsL=ds+deltas
            # print "delta h=",dhL
            return dsL,dhL

    def countdinucs(s):
        """Counts dinucleotide frequencies in a sequence"""
        dinucs={}
        map(dinucs.__setitem__,[a+b for a in 'ACGT' for b in 'ACGT'],[0]*16)
        for i in range(len(s)-1):
            dn=s[i:i+2]
            dinucs[dn]+=1
        return dinucs

    sup=string.upper(s)
    R=1.987 # universal gas constant in Cal/degrees C*Mol
    vsTC,vh=tercorr(sup)
    vs=vsTC
    dinuc=countdinucs(sup)

    
    k=(dnac/4.0)*1e-8
    #With complementary check on, the 4.0 should be changed to a variable.
    
    if rna==0:
        #DNA/DNA
        #Allawi and SantaLucia (1997). Biochemistry 36 : 10581-10594
        vh=vh+dinuc["AA"]*7.9+dinuc["TT"]*7.9+dinuc["AT"]*7.2+dinuc["TA"]*7.2+\
         dinuc["CA"]*8.5+dinuc["TG"]*8.5+dinuc["GT"]*8.4+dinuc["AC"]*8.4+\
         dinuc["CT"]*7.8+dinuc["AG"]*7.8+dinuc["GA"]*8.2+dinuc["TC"]*8.2+\
         dinuc["CG"]*10.6+dinuc["GC"]*10.6+dinuc["GG"]*8+dinuc["CC"]*8
        vs=vs+dinuc["AA"]*22.2+dinuc["TT"]*22.2+dinuc["AT"]*20.4+dinuc["TA"]*21.3+\
         dinuc["CA"]*22.7+dinuc["TG"]*22.7+dinuc["GT"]*22.4+dinuc["AC"]*22.4+\
         dinuc["CT"]*21.0+dinuc["AG"]*21.0+dinuc["GA"]*22.2+dinuc["TC"]*22.2+\
         dinuc["CG"]*27.2+dinuc["GC"]*27.2+dinuc["GG"]*19.9+dinuc["CC"]*19.9
        ds=vs
        dh=vh
    else:
        #RNA/RNA hybridisation of Xia et al (1998)
        #Biochemistry 37: 14719-14735         
        vh=dinuc["AA"]*6.6+dinuc["TT"]*6.6+dinuc["AT"]*5.7+dinuc["TA"]*8.1+\
         dinuc["CA"]*10.5+dinuc["TG"]*10.5+dinuc["GT"]*10.2+dinuc["AC"]*10.2+\
         dinuc["CT"]*7.6+dinuc["AG"]*7.6+dinuc["GA"]*13.3+dinuc["TC"]*13.3+\
         dinuc["CG"]*8.0+dinuc["GC"]*14.2+dinuc["GG"]*12.2+dinuc["CC"]*12.2+\
         dinuc["AA"]*18.4+dinuc["TT"]*18.4+dinuc["AT"]*15.5+dinuc["TA"]*16.9
        vs=vs+dinuc["CA"]*27.8+dinuc["TG"]*27.8+dinuc["GT"]*26.2+dinuc["AC"]*26.2+\
         dinuc["CT"]*19.2+dinuc["AG"]*19.2+dinuc["GA"]*35.5+dinuc["TC"]*35.5+\
         dinuc["CG"]*19.4+dinuc["GC"]*34.9+dinuc["GG"]*29.7+dinuc["CC"]*29.7
        ds=vs
        dh=vh

    ds=ds-0.368*(len(s)-1)*math.log(saltc/1e3)
    tm=((1000* (-dh))/(-ds+(R * (math.log(k)))))-273.15
    # print "ds="+str(ds)
    # print "dh="+str(dh)
    return tm
