# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""A module containing tools to manipulate data from nmrview peaklist files.

"""
__docformat__ = "restructuredtext en"

import sys

HEADERLEN=6

class XpkEntry(object):
    """Provide dictonary access to single entry from nmrview .xpk file.

    This class is suited for handling single lines of non-header data 
    from an nmrview .xpk file. This class provides methods for extracting
    data by the field name which is listed in the last line of the 
    peaklist header.

    Parameters
    ----------
    xpkentry : str
        The line from an nmrview .xpk file.
    xpkheadline : str
        The line from the header file that gives the names of the entries.
        This is typically the sixth line of the header, 1-origin.

    Attributes
    ----------
    fields : dict
        Dictionary of fields where key is in header line, value is an entry.
        Variables are accessed by either their name in the header line as in
        self.field["H1.P"] will return the H1.P entry for example.
        self.field["entrynum"] returns the line number (1st field of line)

    """

    def __init__(self, entry, headline):
       self.fields={}   
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
    """Provide access to header lines and data from a nmrview xpk file.

    Header file lines and file data are available as attributes.

    Parameters
    ----------
    infn : str
        The input nmrview filename.

    Attributes
    ----------
    firstline  : str
        The first line in the header.
    axislabels : str
        The axis labels. 
    dataset    : str
        The label of the dataset.
    sw         : str
        The sw coordinates.
    sf         : str
        The sf coordinates.
    datalabels : str
        The labels of the entries.

    data : list
        File data after header lines.

    Examples
    --------
    >>> import Bio.NMR.xpktools
    >>> peaklist = xpktools.Peaklist('../../Doc/examples/noed.xpk')
    >>> peaklist.firstline
    'label dataset sw sf '
    >>> peaklist.dataset
    'test.nv'
    >>> peaklist.sf
    '{599.8230 } { 60.7860 } { 60.7860 }'
    >>> peaklist.datalabels
    ' H1.L  H1.P  H1.W  H1.B  H1.E  H1.J  15N2.L  15N2.P  15N2.W  15N2.B  15N2.E  15N2.J  N15.L  N15.P  N15.W  N15.B  N15.E  N15.J  vol  int  stat '


    """

    def __init__(self, infn):
        # init the data line list
        self.data = []    

        infile = open(infn,'r')

        # Read in the header lines
        self.firstline  = infile.readline().split("\012")[0]
        self.axislabels = infile.readline().split("\012")[0]
        self.dataset    = infile.readline().split("\012")[0]
        self.sw         = infile.readline().split("\012")[0]
        self.sf         = infile.readline().split("\012")[0]
        self.datalabels = infile.readline().split("\012")[0]

        # Read in the data lines to a list 
        line = infile.readline()
        while line:
            self.data.append(line.split("\012")[0])
            line = infile.readline()

    def residue_dict(self, index):
        """Return a dict of lines in `data` indexed by residue number or a nucleus.

        The nucleus should be given as the input argument in the same form as 
        it appears in the xpk label line (H1, 15N for example)

        Parameters
        ----------
        index : str
            The nucleus to index data by.

        Returns
        -------
        resdict : dict
            Mappings of index nucleus to data line.

        Examples
        --------
        >>> import Bio.NMR.xpktools
        >>> peaklist = xpktools.Peaklist('../../Doc/examples/noed.xpk')
        >>> residue_d = peaklist.residue_dict('H1')
        >>> residue_d.keys()
        ['10', 'maxres', 'minres', '3', '5', '4', '7', '6', '9', '8']
        >>> residue_d['10']
        ['8  10.hn   7.663   0.021   0.010   ++   0.000   10.n   118.341   0.324   0.010   +E   0.000   10.n   118.476   0.324   0.010   +E   0.000  0.49840 0.49840 0']

        """

        maxres, minres = -1, -1

        # Cast the data lines into the xpentry class
        self.resdict = {}
        for i in range(len(self.data)):
            line = self.data[i]
            ind = XpkEntry(line,self.datalabels).fields[index+".L"]
            key = ind.split(".")[0]

            res = int(key)

            if (maxres == -1):
                maxres = res
            if (minres == -1):
                minres = res

            maxres = max([maxres,res])
            minres = min([minres,res])

            if str(res) in self.resdict:
                # Append additional data to list under same key
                templst = self.resdict[str(res)]
                templst.append(line)
                self.resdict[str(res)] = templst

            else:
                # This is a new residue, start a new list
                self.resdict[str(res)] = [line]  # Use [] for list type

        self.resdict["maxres"] = maxres
        self.resdict["minres"] = minres

        return self.resdict

    def write_header(self,outfn):
        """Write header lines from input file to `outfn`."""
        outfile = _try_open_write(outfn)
        outfile.write(self.firstline); outfile.write("\012")
        outfile.write(self.axislabels); outfile.write("\012")
        outfile.write(self.dataset); outfile.write("\012")
        outfile.write(self.sw); outfile.write("\012")
        outfile.write(self.sf); outfile.write("\012")
        outfile.write(self.datalabels); outfile.write("\012")
        outfile.close() 

def _try_open_read(fn):
    """Try to open a file for reading, raising SystemExit if IOError."""
    try:
        infile = open(fn,'r')
    except IOError, e:
        print "file", fn, "could not be opened for reading - quitting."
        sys.exit(0)
    return infile

def _try_open_write(fn):
    """Try to open a file for writing, raising SystemExit if IOError."""
    try:
        infile = open(fn,'w')
    except IOError, e:
        print "file", fn, "could not be opened for writing - quitting."
        sys.exit(0)
    return infile


def replace_entry(line, fieldn, newentry):
    """Replace an entry in a string by the field number.

    No padding is implemented currently.  Spacing will change if the 
    original field entry and the new field entry are of different lengths.

    """

    start = _find_start_entry(line,fieldn)
    leng = len(line[start:].split()[0])
    newline = line[:start]+str(newentry)+line[(start+leng):]
    return newline

def _find_start_entry(line, n):
    """Find the starting character for entry `n` in a space delimited `line`.

    n is counted starting with 1.
    The n=1 field by definition begins at the first character.

    Returns
    -------
    starting character : str
        The starting character for entry `n`.

    """

    # A flag that indicates that the counter is in a field
    infield = 0       

    # Special case
    if (n == 1):
            return 0        

    # Count the number of fields by counting spaces
    c = 1
    leng = len(line)

    # Initialize variables according to whether the first character
    # is a space or a character
    if (line[0] == " "):
            infield = 0
            field = 0
    else:
            infield = 1
            field = 1

    while (c<leng and field<n):
            if (infield):
                    if (line[c] == " " and not (line[c-1] == " ")):
                            infield = 0
            else:
                    if (not line[c] == " "):
                            infield = 1
                            field = field+1

            c = c+1

    return c-1

def data_table(fn_list, datalabel, keyatom):
    """Generate a data table from a list of input .xpk files. 

    Parameters
    ----------
    fn_list : list
        List of .xpk file names.
    datalabel : str
        The data element reported.
    keyatom : str
        The name of the nucleus used as an index for the data table.

    Returns
    -------
    outlist : list
       List of table rows indexed by `keyatom`. 
        
    """

    outlist=[]

    [dict_list, label_line_list] = _read_dicts(fn_list, keyatom)

    # Find global max and min residue numbers
    minr = dict_list[0]["minres"]
    maxr = dict_list[0]["maxres"]
 
    for dictionary in dict_list:
        if (maxr < dictionary["maxres"]):
            maxr = dictionary["maxres"]
        if (minr > dictionary["minres"]):
             minr = dictionary["minres"]

    res = minr
    while res <= maxr:        # s.t. res numbers
        count = 0
        line = str(res)
        for dictionary in dict_list:      # s.t. dictionaries
            label=label_line_list[count]
            if str(res) in dictionary:
                line = line+"\t"+XpkEntry(dictionary[str(res)][0],label).fields[datalabel]
            else:
                line = line+"\t"+"*"
            count = count+1
        line = line+"\n"
        outlist.append(line)
        res = res + 1

    return outlist

def _sort_keys(dictionary):
    """Sort keys in dictionary and return sorted_keys."""
    keys=dictionary.keys()
    sorted_keys=keys.sort()

    """ ..warning:: This returns None, since sort() is in-place.
    """
    return sorted_keys

def _read_dicts(fn_list, keyatom):
    """Read multiple files into a list of residue dictionaries."""
    dict_list, datalabel_list = [], []
    for fn in fn_list:
        peaklist = Peaklist(fn)
        indict = peaklist.residue_dict(keyatom)
        dict_list.append(indict)
        datalabel_list.append(peaklist.datalabels)

    return [dict_list, datalabel_list]
