# Copyright 2004 by Bob Bussell.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""NOEtools: For predicting NOE coordinates from assignment data.

The input and output are modelled on nmrview peaklists.
This modules is suitable for directly generating an nmrview
peaklist with predicted crosspeaks directly from the
input assignment peaklist.
"""

from . import xpktools

__docformat__ = "restructuredtext en"


def predictNOE(peaklist, originNuc, detectedNuc, originResNum, toResNum):
    """Predict the i->j NOE position based on self peak (diagonal) assignments

    Parameters
    ----------
    peaklist : xprtools.Peaklist
        List of peaks from which to derive predictions
    originNuc : str
        Name of originating nucleus.
    originResNum : int
        Index of originating residue.
    detectedNuc : str
        Name of detected nucleus.

    toResNum : int
        Index of detected residue.

    Returns
    -------
    returnLine : str
        The .xpk file entry for the predicted crosspeak.

    Examples
    --------
    Using predictNOE(peaklist,"N15","H1",10,12)
    where peaklist is of the type xpktools.peaklist
    would generate a .xpk file entry for a crosspeak
    that originated on N15 of residue 10 and ended up
    as magnetization detected on the H1 nucleus of
    residue 12


    Notes
    =====
    The initial peaklist is assumed to be diagonal (self peaks only)
    and currently there is no checking done to insure that this
    assumption holds true.  Check your peaklist for errors and
    off diagonal peaks before attempting to use predictNOE.
    """

    returnLine = ""  # The modified line to be returned to the caller

    datamap = _data_map(peaklist.datalabels)

    # Construct labels for keying into dictionary
    originAssCol = datamap[originNuc + ".L"] + 1
    originPPMCol = datamap[originNuc + ".P"] + 1
    detectedPPMCol = datamap[detectedNuc + ".P"] + 1

    # Make a list of the data lines involving the detected
    if str(toResNum) in peaklist.residue_dict(detectedNuc) \
    and str(originResNum) in peaklist.residue_dict(detectedNuc):
        detectedList = peaklist.residue_dict(detectedNuc)[str(toResNum)]
        originList = peaklist.residue_dict(detectedNuc)[str(originResNum)]
        returnLine = detectedList[0]

        for line in detectedList:
            aveDetectedPPM = _col_ave(detectedList, detectedPPMCol)
            aveOriginPPM = _col_ave(originList, originPPMCol)
            originAss = originList[0].split()[originAssCol]

        returnLine = xpktools.replace_entry(returnLine, originAssCol + 1, originAss)
        returnLine = xpktools.replace_entry(returnLine, originPPMCol + 1, aveOriginPPM)

    return returnLine


def _data_map(labelline):
    # Generate a map between datalabels and column number
    #   based on a labelline
    i = 0  # A counter
    datamap = {}  # The data map dictionary
    labelList = labelline.split()  # Get the label line

    # Get the column number for each label
    for i in range(len(labelList)):
        datamap[labelList[i]] = i

    return datamap


def _col_ave(list, col):
    # Compute average values from a particular column in a string list
    total = 0.0
    n = 0
    for element in list:
        total += float(element.split()[col])
        n += 1
    return total / n
