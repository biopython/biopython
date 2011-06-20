# xpktools.py: A python module containing function definitions and classes
#          useful for manipulating data from nmrview .xpk peaklist files.
#
# ********** INDEX of functions and classes **********
#
#    XpkEntry class: A class suited for handling single lines of
#        non-header data from an nmrview .xpk file.  This class
#        provides methods for extracting data by the field name
#        which is listed in the last line of the peaklist header.

import sys

# * * * * * INITIALIZATIONS * * * * *
HEADERLEN=6
# * * * * * _______________ * * * * *

class XpkEntry(object):
    # Usage: XpkEntry(xpkentry,xpkheadline) where xpkentry is the line
    #        from an nmrview .xpk file and xpkheadline is the line from
    #        the header file that gives the names of the entries
    #        which is typcially the sixth line of the header (counting fm 1)
    # Variables are accessed by either their name in the header line as in
    #   self.field["H1.P"] will return the H1.P entry for example.
    #   self.field["entrynum"] returns the line number (1st field of line)

    def __init__(self,entry,headline):
       self.fields={}   # Holds all fields from input line in a dictionary
                        # keys are data labels from the .xpk header 
       datlist  = entry.split()
       headlist = headline.split()

       i=0  
       for i in range(len(datlist)-1):
         self.fields[headlist[i]]=datlist[i+1]
       i=i+1

       try:
           self.fields["entrynum"]=datlist[0]
       except IndexError, e:
           pass

class Peaklist(object):
    # This class reads in an entire xpk file and returns
    # Header file lines are available as attributes
    # The data lines are available as a list
    def __init__(self,infn):
    
        self.data=[]    # init the data line list

        infile=open(infn,'r')

        # Read in the header lines
        self.firstline=infile.readline().split("\012")[0]
        self.axislabels=infile.readline().split("\012")[0]
        self.dataset=infile.readline().split("\012")[0]
        self.sw=infile.readline().split("\012")[0]
        self.sf=infile.readline().split("\012")[0]
        self.datalabels=infile.readline().split("\012")[0]

        # Read in the data lines to a list 
        line=infile.readline()
        while line:
            self.data.append(line.split("\012")[0])
        line=infile.readline()

    def residue_dict(self,index):
        # Generate a dictionary idexed by residue number or a nucleus
        # The nucleus should be given as the input argument in the
        # same form as it appears in the xpk label line (H1, 15N for example)

        maxres=-1; minres=-1

        # Cast the data lines into the xpentry class
        self.dict={}
        for i in range(len(self.data)):
            line=self.data[i]
            ind=XpkEntry(line,self.datalabels).fields[index+".L"]
            key=ind.split(".")[0]

            res=int(key)

            if (maxres==-1):
                maxres=res
            if (minres==-1):
                minres=res

            maxres=max([maxres,res])
            minres=min([minres,res])

            if str(res) in self.dict:
                # Append additional data to list under same key
                templst=self.dict[str(res)]
                templst.append(line)
                self.dict[str(res)]=templst

            else:
                # This is a new residue, start a new list
                self.dict[str(res)]=[line]  # Use [] for list type

        self.dict["maxres"]=maxres
        self.dict["minres"]=minres

        return self.dict

    def write_header(self,outfn):
        outfile=_try_open_write(outfn)
        outfile.write(self.firstline);outfile.write("\012")
        outfile.write(self.axislabels);outfile.write("\012")
        outfile.write(self.dataset);outfile.write("\012")
        outfile.write(self.sw);outfile.write("\012")
        outfile.write(self.sf);outfile.write("\012")
        outfile.write(self.datalabels);outfile.write("\012")
        outfile.close() 

def _try_open_read(fn):
# Try to open a file for reading.  Exit on IOError
  try:
    infile=open(fn,'r')
  except IOError, e:
    print "file", fn, "could not be opened for reading - quitting."
    sys.exit(0)
  return infile

def _try_open_write(fn):
# Try to open a file for writing.  Exit on IOError
  try:
    infile=open(fn,'w')
  except IOError, e:
    print "file", fn, "could not be opened for writing - quitting."
    sys.exit(0)
  return infile


def replace_entry(line,fieldn,newentry):
        # Replace an entry in a string by the field number
        # No padding is implemented currently.  Spacing will change if
        #  the original field entry and the new field entry are of
        #  different lengths.
        # This method depends on xpktools._find_start_entry

        start=_find_start_entry(line,fieldn)
        leng=len(line[start:].split()[0])
        newline=line[:start]+str(newentry)+line[(start+leng):]
        return newline

def _find_start_entry(line,n):
        # find the starting point character for the n'th entry in
        # a space delimited line.  n is counted starting with 1
        # The n=1 field by definition begins at the first character
        # This function is used by replace_entry

        infield=0       # A flag that indicates that the counter is in a field

        if (n==1):
                return 0        # Special case

        # Count the number of fields by counting spaces
        c=1
        leng=len(line)

        # Initialize variables according to whether the first character
        #  is a space or a character
        if (line[0]==" "):
                infield=0
                field=0
        else:
                infield=1
                field=1


        while (c<leng and field<n):
                if (infield):
                        if (line[c]==" " and not (line[c-1]==" ")):
                                infield=0
                else:
                        if (not line[c]==" "):
                                infield=1
                                field=field+1

                c=c+1

        return c-1


def data_table(fn_list, datalabel, keyatom):
# Generate and generate a data table from a list of
# input xpk files <fn_list>.  The data element reported is
# <datalabel> and the index for the data table is by the 
# nucleus indicated by <keyatom>.

  outlist=[]

  [dict_list,label_line_list]=_read_dicts(fn_list,keyatom)

  # Find global max and min residue numbers
  minr=dict_list[0]["minres"]; maxr=dict_list[0]["maxres"]
 
  for dictionary in dict_list:
    if (maxr < dictionary["maxres"]):
      maxr = dictionary["maxres"]
    if (minr > dictionary["minres"]):
      minr = dictionary["minres"]

  res=minr
  while res <= maxr:        # s.t. res numbers
    count=0
    line=str(res)
    for dictionary in dict_list:      # s.t. dictionaries
      label=label_line_list[count]
      if str(res) in dictionary:
        line=line+"\t"+XpkEntry(dictionary[str(res)][0],label).fields[datalabel]
      else:
        line=line+"\t"+"*"
      count=count+1
    line=line+"\n"
    outlist.append(line)
    res=res+1

  return outlist

def _sort_keys(dictionary):
  keys=dictionary.keys()
  sorted_keys=keys.sort()
  return sorted_keys

def _read_dicts(fn_list, keyatom):
# Read multiple files into a list of residue dictionaries
  dict_list=[]; datalabel_list=[]
  for fn in fn_list:
    peaklist=Peaklist(fn); dict=peaklist.residue_dict(keyatom)
    dict_list.append(dict)
    datalabel_list.append(peaklist.datalabels)

  return [dict_list, datalabel_list]
