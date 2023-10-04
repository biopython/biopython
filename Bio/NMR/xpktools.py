# Copyright 2004 by Bob Bussell.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tools to manipulate data from nmrview .xpk peaklist files."""


HEADERLEN = 6


class XpkEntry:
    """Provide dictionary access to single entry from nmrview .xpk file.

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
        """Initialize the class."""
        # Holds all fields from input line in a dictionary
        # keys are data labels from the .xpk header
        datlist = entry.split()
        headlist = headline.split()

        self.fields = dict(zip(headlist, datlist[1:]))

        try:
            self.fields["entrynum"] = datlist[0]
        except IndexError:
            pass


class Peaklist:
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
    >>> from Bio.NMR.xpktools import Peaklist
    >>> peaklist = Peaklist('../Doc/examples/nmr/noed.xpk')
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
        """Initialize the class."""
        with open(infn) as infile:
            # Read in the header lines
            self.firstline = infile.readline().split("\012")[0]
            self.axislabels = infile.readline().split("\012")[0]
            self.dataset = infile.readline().split("\012")[0]
            self.sw = infile.readline().split("\012")[0]
            self.sf = infile.readline().split("\012")[0]
            self.datalabels = infile.readline().split("\012")[0]

            # Read in the data lines to a list
            self.data = [line.split("\012")[0] for line in infile]

    def residue_dict(self, index):
        """Return a dict of lines in 'data' indexed by residue number or a nucleus.

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
        >>> from Bio.NMR.xpktools import Peaklist
        >>> peaklist = Peaklist('../Doc/examples/nmr/noed.xpk')
        >>> residue_d = peaklist.residue_dict('H1')
        >>> sorted(residue_d.keys())
        ['10', '3', '4', '5', '6', '7', '8', '9', 'maxres', 'minres']
        >>> residue_d['10']
        ['8  10.hn   7.663   0.021   0.010   ++   0.000   10.n   118.341   0.324   0.010   +E   0.000   10.n   118.476   0.324   0.010   +E   0.000  0.49840 0.49840 0']

        """
        maxres = -1
        minres = -1

        # Cast the data lines into the xpentry class
        self.dict = {}
        for line in self.data:
            ind = XpkEntry(line, self.datalabels).fields[index + ".L"]
            key = ind.split(".")[0]

            res = int(key)

            if maxres == -1:
                maxres = res
            if minres == -1:
                minres = res

            maxres = max([maxres, res])
            minres = min([minres, res])
            res = str(res)

            try:
                # Append additional data to list under same key
                self.dict[res].append(line)
            except KeyError:
                # This is a new residue, start a new list
                self.dict[res] = [line]  # Use [] for list type

        self.dict["maxres"] = maxres
        self.dict["minres"] = minres

        return self.dict

    def write_header(self, outfn):
        """Write header lines from input file to handle ``outfn``."""
        with open(outfn, "w") as outfile:
            outfile.write(self.firstline)
            outfile.write("\012")
            outfile.write(self.axislabels)
            outfile.write("\012")
            outfile.write(self.dataset)
            outfile.write("\012")
            outfile.write(self.sw)
            outfile.write("\012")
            outfile.write(self.sf)
            outfile.write("\012")
            outfile.write(self.datalabels)
            outfile.write("\012")


def replace_entry(line, fieldn, newentry):
    """Replace an entry in a string by the field number.

    No padding is implemented currently.  Spacing will change if
    the original field entry and the new field entry are of
    different lengths.
    """
    # This method depends on xpktools._find_start_entry

    start = _find_start_entry(line, fieldn)
    leng = len(line[start:].split()[0])
    newline = line[:start] + str(newentry) + line[(start + leng) :]
    return newline


def _find_start_entry(line, n):
    """Find the starting character for entry ``n`` in a space delimited ``line`` (PRIVATE).

    n is counted starting with 1.
    The n=1 field by definition begins at the first character.

    Returns
    -------
    starting character : str
        The starting character for entry ``n``.

    """
    # This function is used by replace_entry

    if n == 1:
        return 0  # Special case

    # Count the number of fields by counting spaces
    c = 1
    leng = len(line)

    # Initialize variables according to whether the first character
    #  is a space or a character
    if line[0] == " ":
        infield = False
        field = 0
    else:
        infield = True
        field = 1

    while c < leng and field < n:
        if infield:
            if line[c] == " " and line[c - 1] != " ":
                infield = False
            else:
                if line[c] != " ":
                    infield = True
                    field += 1

        c += 1

    return c - 1


def data_table(fn_list, datalabel, keyatom):
    """Generate a data table from a list of input xpk files.

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
       List of table rows indexed by ``keyatom``.

    """
    # TODO - Clarify this docstring, add an example?
    outlist = []

    dict_list, label_line_list = _read_dicts(fn_list, keyatom)

    # Find global max and min residue numbers
    minr = dict_list[0]["minres"]
    maxr = dict_list[0]["maxres"]

    for dictionary in dict_list:
        if maxr < dictionary["maxres"]:
            maxr = dictionary["maxres"]
        if minr > dictionary["minres"]:
            minr = dictionary["minres"]

    res = minr
    while res <= maxr:  # s.t. res numbers
        count = 0
        key = str(res)
        line = key
        for dictionary in dict_list:  # s.t. dictionaries
            label = label_line_list[count]
            if key in dictionary:
                line = (
                    line + "\t" + XpkEntry(dictionary[key][0], label).fields[datalabel]
                )
            else:
                line += "\t*"
            count += 1
        line += "\n"
        outlist.append(line)
        res += 1

    return outlist


def _read_dicts(fn_list, keyatom):
    """Read multiple files into a list of residue dictionaries (PRIVATE)."""
    dict_list = []
    datalabel_list = []
    for fn in fn_list:
        peaklist = Peaklist(fn)
        dictionary = peaklist.residue_dict(keyatom)
        dict_list.append(dictionary)
        datalabel_list.append(peaklist.datalabels)

    return [dict_list, datalabel_list]


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
