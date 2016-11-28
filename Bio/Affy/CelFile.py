# Copyright 2004 by Harry Zuzan and Adam Kurkiewicz. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Classes for accessing the information in Affymetrix cel files version 3 and 4

Functions:
read      Read a cel file and store its contents in a Record

Classes:
Record    Contains the information from a cel file
"""

from __future__ import print_function
import sys
import struct
import numpy

# for debugging
# import pprint
# pp = pprint.PrettyPrinter(indent=4)


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile")


class Record(object):
    """Stores the information in a cel file

    Example usage:

    >>> from Bio.Affy import CelFile
    >>> with open('Affy/affy_v3_example.CEL', "r") as handle:
    ...     c = CelFile.read(handle, strict=False)
    ...
    >>> print(c.ncols, c.nrows)
    5 5
    >>> print(c.intensities)
    [[   234.    170.  22177.    164.  22104.]
     [   188.    188.  21871.    168.  21883.]
     [   188.    193.  21455.    198.  21300.]
     [   188.    182.  21438.    188.  20945.]
     [   193.  20370.    174.  20605.    168.]]
    >>> print(c.stdevs)
    [[   24.     34.5  2669.     19.7  3661.2]
     [   29.8    29.8  2795.9    67.9  2792.4]
     [   29.8    88.7  2976.5    62.   2914.5]
     [   29.8    76.2  2759.5    49.2  2762. ]
     [   38.8  2611.8    26.6  2810.7    24.1]]
    >>> print(c.npix)
    [[25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]]

    """
    def __init__(self):
        self.version = None
        self.GridCornerUL = None
        self.GridCornerUR = None
        self.GridCornerLR = None
        self.GridCornerLL = None
        self.DatHeader = None
        self.Algorithm = None
        self.AlgorithmParameters = None
        self.NumberCells = None
        self.intensities = None
        self.stdevs = None
        self.npix = None
        self.nrows = None
        self.ncols = None
        self.nmask = None
        self.mask = None
        self.noutliers = None
        self.outliers = None
        self.modified = None


def read(handle, strict=False):
    # If we fail to read the magic number, then it will remain None, and thus
    # we will invoke read3 (if mode is not strict), or raise IOError if mode is
    # strict.
    magicNumber = None
    # We check if the handle is a file-like object. If it isn't, and the mode
    # is strict, we raise an error. If it isn't and the mode isn't strict, we
    # continue (perhaps somebody has got a CEL-file-like iterable, which used
    # to work with previous versions of biopython and we don't want to maintain
    # backwards compatibility).
    try:
        mode = handle.mode
        # By definition an Affymetrix v4 CEL file has 64 as the first 4 bytes.
        # Note that we use little-endian irrespective of the platform, again by
        # definition.
        position = handle.tell()
        magicNumber = struct.unpack('<i', handle.read(4))[0]
    except (AttributeError, TypeError):
        if strict:
            raise IOError("You have to pass a file-like object. You can get "
                          "a file-like object by using `open`.")
    finally:
        try:
            # reset the offset, to avoid breaking either v3 or v4.
            handle.seek(position)
        except AttributeError:
            pass

    if magicNumber != 64:
        # In v4 we're always strict, as we don't have to worry about backwards
        # compatibility
        if strict and mode != "r":
            raise IOError("You're trying to open an Affymetrix v3 CEL file."
                          "You have to use a read mode, like this"
                          "`open(filename, \"r\")`.")
        return read3(handle)
    else:
        if mode != "rb":
            raise IOError("You're trying to open an Affymetrix v4 CEL file."
                          "You have to use a read binary mode, like this"
                          "`open(filename \"rb\")`.")
        return read4(handle)


# read Affymetrix files version 4.
def read4(f):
    # We follow the documentation here:
    # http://www.affymetrix.com/estore/support/developer/powertools/changelog/gcos-agcc/cel.html.affx
    record = Record()
    preHeaders = ["magic", "version", "columns", "rows", "cellNo", "headerLen"]
    preHeadersMap = dict()
    headersMap = dict()

    # load pre-headers
    for name in preHeaders:
        preHeadersMap[name] = struct.unpack('<i', f.read(4))[0]
    char = f.read(preHeadersMap["headerLen"])
    header = char.decode("ascii", "ignore")
    for header in header.split("\n"):
        if "=" in header:
            header = header.split("=")
            headersMap[header[0]] = "=".join(header[1:])

    # for debugging
    # pp.pprint("preHeadersMap")
    # pp.pprint(preHeadersMap)
    # pp.pprint("headersMap")
    # pp.pprint(headersMap)

    record.version = preHeadersMap["version"]
    if record.version != 4:
        raise(ValueError("Parse Error"))
    record.GridCornerUL = headersMap["GridCornerUL"]
    record.GridCornerUR = headersMap["GridCornerUR"]
    record.GridCornerLR = headersMap["GridCornerLR"]
    record.GridCornerLL = headersMap["GridCornerLL"]
    record.DatHeader = headersMap["DatHeader"]
    record.Algorithm = headersMap["Algorithm"]
    record.AlgorithmParameters = headersMap["AlgorithmParameters"]
    record.NumberCells = preHeadersMap["cellNo"]
    # record.intensities are set below
    # record.stdevs are set below
    # record.npix are set below
    record.nrows = int(headersMap["Rows"])
    record.ncols = int(headersMap["Cols"])

    # These cannot be reliably set in v4, because of discrepancies between real
    # data and the documented format.
    record.nmask = None
    record.mask = None
    record.noutliers = None
    record.outliers = None
    record.modified = None

    # Real data never seems to have anything but zeros here, but we don't want
    # to take chances. Raising an error is better than returning unreliable
    # data.

    if headersMap["Axis-invertX"] != "0":
        raise(ValueError("Parse Error"))

    if headersMap["AxisInvertY"] != "0":
        raise(ValueError("Parse Error"))

    if headersMap["AxisInvertY"] != "0":
        raise(ValueError("Parse Error"))

    # This is unfortunately undocumented, but it turns out that real data has
    # the `record.AlgorithmParameters` repeated in the data section, until an
    # EOF, i.e. b'\x04'.
    char = b'\x00'
    safetyValve = 10**4
    for i in range(safetyValve):
        char = f.read(1)
        # For debugging
        # print([i for i in char], end="")
        if char == b'\x04':
            break
        if i == safetyValve:
            raise(ValueError("Parse Error"))

    # After that there are precisely 15 bytes padded. Again, undocumented.
    padding = f.read(15)

    # That's how we pull out the values (triplets of the form float, float,
    # signed short).
    structa = struct.Struct("< f f h")

    # There are 10 bytes in our struct.
    structSize = 10

    # We initialise the most important: intensities, stdevs and npixs.
    record.intensities = numpy.empty(record.NumberCells, dtype=float)
    record.stdevs = numpy.empty(record.NumberCells, dtype=float)
    record.npix = numpy.empty(record.NumberCells, dtype=int)

    b = f.read(structSize * record.NumberCells)
    for i in range(record.NumberCells):
        binaryFragment = b[i * structSize: (i + 1) * structSize]
        intensity, stdevs, npix = structa.unpack(binaryFragment)
        record.intensities[i] = intensity
        record.stdevs[i] = stdevs
        record.npix[i] = npix

    # reshape without copying.
    def reshape(array):
        view = array.view()
        view.shape = (record.nrows, record.ncols)
        return view

    record.intensities = reshape(record.intensities)
    record.stdevs = reshape(record.stdevs)
    record.npix = reshape(record.npix)

    return record


def read3(handle):
    """
    Read the information in a cel file, and store it in a Record.
    """
    # Needs error handling.
    # Needs to know the chip design.
    record = Record()
    section = ""
    for line in handle:
        if not line.strip():
            continue
        # Set current section
        if line[:5] == "[CEL]":
            section = "CEL"
        elif line[:8] == "[HEADER]":
            section = "HEADER"
        elif line[:11] == "[INTENSITY]":
            section = "INTENSITY"
            record.intensities = numpy.zeros((record.nrows, record.ncols))
            record.stdevs = numpy.zeros((record.nrows, record.ncols))
            record.npix = numpy.zeros((record.nrows, record.ncols), int)
        elif line[:7] == "[MASKS]":
            section = "MASKS"
            record.mask = numpy.zeros((record.nrows, record.ncols))
        elif line[:10] == "[OUTLIERS]":
            section = "OUTLIERS"
            record.outliers = numpy.zeros((record.nrows, record.ncols))
        elif line[:10] == "[MODIFIED]":
            section = "MODIFIED"
            record.modified = numpy.zeros((record.nrows, record.ncols))
        elif line[0] == "[":
            # This would be an unknown section
            section = ""
        elif section == "CEL":
            keyword, value = line.split("=", 1)
            if keyword == 'Version':
                record.version = int(value)
        elif section == "HEADER":
            # Set record.ncols and record.nrows, remaining data goes into
            # record.header dict
            keyword, value = line.split("=", 1)
            if keyword == "Cols":
                record.ncols = int(value)
            elif keyword == "Rows":
                record.nrows = int(value)
            elif keyword == 'GridCornerUL':
                x, y = value.split()
                record.GridCornerUL = (int(x), int(y))
            elif keyword == 'GridCornerUR':
                x, y = value.split()
                record.GridCornerUR = (int(x), int(y))
            elif keyword == 'GridCornerLR':
                x, y = value.split()
                record.GridCornerLR = (int(x), int(y))
            elif keyword == 'GridCornerLL':
                x, y = value.split()
                record.GridCornerLL = (int(x), int(y))
            elif keyword == 'DatHeader':
                record.DatHeader = value.strip('\n\r')
            elif keyword == 'Algorithm':
                record.Algorithm = value.strip('\n\r')
            elif keyword == 'AlgorithmParameters':
                record.AlgorithmParameters = value.strip('\n\r')
        elif section == "INTENSITY":
            if "NumberCells" in line:
                record.NumberCells = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.intensities[x, y] = float(words[2])
                record.stdevs[x, y] = float(words[3])
                record.npix[x, y] = int(words[4])
        elif section == "MASKS":
            if "NumberCells" in line:
                record.nmask = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.mask[x, y] = int(1)
        elif section == "OUTLIERS":
            if "NumberCells" in line:
                record.noutliers = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.outliers[x, y] = int(1)
        elif section == "MODIFIED":
            if "NumberCells" in line:
                record.nmodified = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.modified[x, y] = float(words[2])
        else:
            continue
    return record

if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
