#!/usr/bin/env python

# ** What is this code? **
#
# This is an example script for using new biopython modules for
# generating NOE crosspeak peaklists from diagonal peaklist within
# the framework of nmrview.  
# nmrview is not required to run this script, only an installed 
# version of python.


# ** What's important? **
#
# The xpktools.py and NOEtools.py modules are what I'm trying to
# demonstrate.  They provide methods and a data class for performing
# some general analysis on NMR data taken directly from peaklist.
# The code in this script will demonstrate how they are used.


# ** Who wrote this code? **
# Robert Bussell, Jr.
# rgb2003@med.cornell.edu


# ** Running this script **
#
# To run this script on a UNIX/Linux system, make it executable and
# modify the first line of this script to point to python if necessary.
# First try running the code with the peaklist that I provide to get
# the feel of how things work, then you can use your own peaklist if you
# modify the variables under the "INITS" code block to make it work
# with your data.  The modules xpktools and NOEtools can be called from
# your own scripts when you have them in place on your computer.
# NOTE: It is very important to have an intact peaklist.  If you copy and
# paste mine into a file be prepared to remove inappropriate line breaks.


# ** Output of this script **
# 
# This script generates a human readable standard output version of the
# NOE coordinates as well as an nmrview peaklist out_example.xpk.






# ***********************************************************************

# ***** LOAD MODULES *****
import getopt
import string
import sys

# -- don't need to modify sys.path with the *tools in Biopython
# -- just need Biopython installed somewhere in the PYTHONPATH
#sys.path=[sys.path,"./"]
#sys.path=[sys.path,"/usr/people/robert/development/xpktools"]
from Bio.NMR import xpktools # Contains data classes and functions for .xpk files
from Bio.NMR import NOEtools # A module specific for generate NOE predictions 

# * * * * * * * * * * MAIN * * * * * * * * * *

# ***** INITS *****

inc=1                     # The NOE increment (n where i->i+n and i->i-n are noes)
infn="./noed.xpk"         # Input peaklist
outfn="./out_example.xpk" # Output peaklist
detectatom="H1"           # Directly detected atom
relayatom="N15"           # J-coupling from here to detected atom
fromatom="15N2"           # The other labelled nucleus 




# *-*-*  First the peaklist is read into a data class from xpktools
# *-*-*  that contains methods for easily extracting information from
# *-*-*  the peaklist file

peaklist=xpktools.Peaklist(infn) # infn is the name of the xpk file 




# *-*-* The class attribute residue_dict allows the data lines
# *-*-* to be separated from the header and returned here to
# *-*-* the dictionary <dict> as a list indexed by the assignment
# *-*-* of any of the nuclei in the file -- here, the detected atom
# *-*-* is used 

dict=peaklist.residue_dict(detectatom)




# *-*-* As well as the data, the dictionary contains two other entries,
# *-*-* corresponding to the maximum and minimum residues indexed 

MAXRES=dict["maxres"]
MINRES=dict["minres"]



#****** CALCULATE AND WRITE CROSSPEAK PEAKLIST *****



# *-*-* The peaklist class has a method for writing out the header
# *-*-* information in a format recognizable by nmrview

peaklist.write_header(outfn) # Write the header to the output file



# *-*-* Predict the i->i+inc and i->i-inc noe positions if possible
# *-*-* Write each one to the output file as they are calculated

count=0 # A counter that number the output data lines in order
res=MINRES # minimum residue number in the set 
outlist=[] # Holds the output data


while (res<=MAXRES):



# *-*-* Predicting the NOE positions based on peak assignment data
# *-*-* is done by supplying the peaklist to and specifying the label
# *-*-* of the origin and detected atom in the NOE transfer as well as
# *-*-* the residues between which the NOE transfer takes places.
 
  noe1=NOEtools.predictNOE(peaklist,"15N2","H1",res,res+inc)
  noe2=NOEtools.predictNOE(peaklist,"15N2","H1",res,res-inc)



# *-*-* The output of predictNOE is in the form of an xpk entry line
# *-*-* suitable for printing to an output file
# *-*-* Additionally, it is possible to extract information easily from
# *-*-* these output lines by using the xpktools.XpkEntry class

  entry1=xpktools.XpkEntry(noe1,peaklist.datalabels)

  if noe1!="":

  # *-*-* Here I'm using the XpkEntry class to gain access to
  # *-*-* specific fields in the that make the information
  # *-*-* more readable and suitable for creating data tables
  # *-*-* This output will be printed to the screen.
  # *-*-* The data table contains the assignment, coordinates and
  # *-*-* intensity of the resonance.

    print string.split(entry1.fields["15N2.L"],".")[0], "-->", string.split(entry1.fields["N15.L"],".")[0], "\t", entry1.fields["H1.P"], entry1.fields["N15.P"], entry1.fields["15N2.P"],entry1.fields["int"]


    noe1=noe1+"\012"
    noe1=xpktools.replace_entry(noe1,1,count)
    outlist.append(noe1)
    count += 1

    if noe2!="":
      noe2=noe2+"\012"
      noe2=xpktools.replace_entry(noe2,1,count)
      outlist.append(noe2)
      count += 1
  res += 1

# Open the output file and write the data
outfile=open(outfn,'a')
outfile.writelines(outlist) # Write the output lines to the file
outfile.close()
