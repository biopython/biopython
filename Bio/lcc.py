# Copyright 2003 by Sebastian Bassi. sbassi@genesdigitales.com
# All rights reserved.  This code is part of the Biopython 
# distribution and governed by its license.
# Please see the LICENSE file that should have been included as part
# of this package.

import math
from string import count

crom=0
compone=[0]
lccsal=[0]

def lcc_mult(seq,wsize,start,end):
    """Return a list called lccsal, the LCC, a complexity measure 
from a sequence, called seq."""
    l2=math.log(2)
    tamseq=end-start
    global compone
    global lccsal
    compone=[0]
    lccsal=[0]
    for i in range(wsize):
        compone.append(((i+1)/float(wsize))*((math.log((i+1)/float(wsize)))/l2))
    window=seq[0:wsize]
    cant_a=count(window,'A')
    cant_c=count(window,'C')
    cant_t=count(window,'T')
    cant_g=count(window,'G')
    term_a=compone[cant_a]
    term_c=compone[cant_c]
    term_t=compone[cant_t]
    term_g=compone[cant_g]
    lccsal[0]=(-(term_a+term_c+term_t+term_g))
    tail=seq[0]
    for x in range (tamseq-wsize):
        window=seq[x+1:wsize+x+1]
        if tail==window[-1]:
            lccsal.append(lccsal[-1])
            #break
        elif tail=='A':
            cant_a=cant_a-1
            if window[-1]=='C':
                cant_c=cant_c+1
                term_a=compone[cant_a]
                term_c=compone[cant_c]
                lccsal.append(-(term_a+term_c+term_t+term_g))
            elif window[-1]=='T':
                cant_t=cant_t+1
                term_a=compone[cant_a]
                term_t=compone[cant_t]
                lccsal.append(-(term_a+term_c+term_t+term_g))
            elif window[-1]=='G':
                cant_g=cant_g+1
                term_a=compone[cant_a]
                term_g=compone[cant_g]
                lccsal.append(-(term_a+term_c+term_t+term_g))
        elif tail=='C':
            cant_c=cant_c-1
            if window[-1]=='A':
                cant_a=cant_a+1
                term_a=compone[cant_a]
                term_c=compone[cant_c]
                lccsal.append(-(term_a+term_c+term_t+term_g))
            elif window[-1]=='T':
                cant_t=cant_t+1
                term_c=compone[cant_c]
                term_t=compone[cant_t]
                lccsal.append(-(term_a+term_c+term_t+term_g))
            elif window[-1]=='G':
                cant_g=cant_g+1
                term_c=compone[cant_c]
                term_g=compone[cant_g]
                lccsal.append(-(term_a+term_c+term_t+term_g))
        elif tail=='T':
            cant_t=cant_t-1
            if window[-1]=='A':
                cant_a=cant_a+1
                term_a=compone[cant_a]
                term_t=compone[cant_t]
                lccsal.append(-(term_a+term_c+term_t+term_g))
            elif window[-1]=='C':
                cant_c=cant_c+1
                term_c=compone[cant_c]
                term_t=compone[cant_t]
                lccsal.append(-(term_a+term_c+term_t+term_g))
            elif window[-1]=='G':
                cant_g=cant_g+1
                term_t=compone[cant_t]
                term_g=compone[cant_g]
                lccsal.append(-(term_a+term_c+term_t+term_g))
        elif tail=='G':
            cant_g=cant_g-1
            if window[-1]=='A':
                cant_a=cant_a+1
                term_a=compone[cant_a]
                term_g=compone[cant_g]
                lccsal.append(-(term_a+term_c+term_t+term_g))
            elif window[-1]=='C':
                cant_c=cant_c+1
                term_c=compone[cant_c]
                term_g=compone[cant_g]
                lccsal.append(-(term_a+term_c+term_t+term_g))
            elif window[-1]=='T':
                cant_t=cant_t+1
                term_t=compone[cant_t]
                term_g=compone[cant_g]
                lccsal.append(-(term_a+term_c+term_t+term_g))
        tail=window[0]
    return lccsal

def lcc_simp(seq,start,end):
    """Return LCC, a complexity measure from a sequence (seq.)"""
    wsize=end-start
    l2=math.log(2)
    window=seq[start:end]
    if count(window,'A')==0:
        term_a=0
	# This check is usefull in order to avoid calculate log of 0.
    else:
        term_a=((count(window,'A'))/float(wsize))*((math.log((count(window,'A'))/float(wsize)))/l2)
    if count(window,'C')==0:
        term_c=0
    else:
        term_c=((count(window,'C'))/float(wsize))*((math.log((count(window,'C'))/float(wsize)))/l2)
    if count(window,'T')==0:
        term_t=0
    else:
        term_t=((count(window,'T'))/float(wsize))*((math.log((count(window,'T'))/float(wsize)))/l2)
    if count(window,'G')==0:
        term_g=0
    else:
        term_g=((count(window,'G'))/float(wsize))*((math.log((count(window,'G'))/float(wsize)))/l2)
    lccsal=-(term_a+term_c+term_t+term_g)
    return lccsal
