# Copyright 2011 by Wibowo Arindrarto (w.arindrarto@gmail.com)
# Revisions copyright 2011-2016 by Peter Cock.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO parser for the ABI format.

ABI is the format used by Applied Biosystem's sequencing machines to store
sequencing results.

For more details on the format specification, visit:
http://www6.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf

"""
import datetime
import struct
import sys

from os.path import basename

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import SequenceIterator


# dictionary for determining which tags goes into SeqRecord annotation
# each key is tag_name + tag_number
# if a tag entry needs to be added, just add its key and its key
# for the annotations dictionary as the value
# dictionary for tags that require preprocessing before use in creating
# seqrecords
_EXTRACT = {
    "TUBE1": "sample_well",
    "DySN1": "dye",
    "GTyp1": "polymer",
    "MODL1": "machine_model",
}


# Complete data structure representing 98% of the API. The general section
# represents the part of the API that's common to ALL instruments, whereas the
# instrument specific sections are labelled as they are in the ABIF spec
#
# Keys don't seem to clash from machine to machine, so when we parse, we look
# for ANY key, and store that in the raw ABIF data structure attached to the
# annotations, with the assumption that anyone parsing the data can look up the
# spec themself
#
# Key definitions are retained in case end users want "nice" labels pre-made
# for them for all of the available fields.
_INSTRUMENT_SPECIFIC_TAGS = {}

# fmt: off
_INSTRUMENT_SPECIFIC_TAGS["general"] = {
    "APFN2": "Sequencing Analysis parameters file name",
    "APXV1": "Analysis Protocol XML schema version",
    "APrN1": "Analysis Protocol settings name",
    "APrV1": "Analysis Protocol settings version",
    "APrX1": "Analysis Protocol XML string",
    "CMNT1": "Sample Comment",
    "CTID1": "Container Identifier, a.k.a. plate barcode",
    "CTNM1": "Container name, usually identical to CTID, but not necessarily so",
    "CTTL1": "Comment Title",
    "CpEP1": "Capillary type electrophoresis. 1 for a capillary based machine. 0 for a slab gel based machine.",
    "DATA1": "Channel 1 raw data",
    "DATA2": "Channel 2 raw data",
    "DATA3": "Channel 3 raw data",
    "DATA4": "Channel 4 raw data",
    "DATA5": "Short Array holding measured volts/10 (EP voltage) during run",
    "DATA6": "Short Array holding measured milliAmps trace (EP current) during run",
    "DATA7": "Short Array holding measured milliWatts trace (Laser EP Power) during run",
    "DATA8": "Short Array holding measured oven Temperature (polymer temperature) trace during run",
    "DATA9": "Channel 9 processed data",
    "DATA10": "Channel 10 processed data",
    "DATA11": "Channel 11 processed data",
    "DATA12": "Channel 12 processed data",
    # Prism 3100/3100-Avant may provide DATA105
    #          3130/3130-XL may provide DATA105
    # 3530/3530-XL may provide DATA105-199, 9-12, 205-299
    "DSam1": "Downsampling factor",
    "DySN1": "Dye set name",
    "Dye#1": "Number of dyes",
    "DyeN1": "Dye 1 name",
    "DyeN2": "Dye 2 name",
    "DyeN3": "Dye 3 name",
    "DyeN4": "Dye 4 name",
    "DyeW1": "Dye 1 wavelength",
    "DyeW2": "Dye 2 wavelength",
    "DyeW3": "Dye 3 wavelength",
    "DyeW4": "Dye 4 wavelength",
    # 'DyeN5-N': 'Dye 5-N Name',
    # 'DyeW5-N': 'Dye 5-N Wavelength',
    "EPVt1": "Electrophoresis voltage setting (volts)",
    "EVNT1": "Start Run event",
    "EVNT2": "Stop Run event",
    "EVNT3": "Start Collection event",
    "EVNT4": "Stop Collection event",
    "FWO_1": 'Base Order. Sequencing Analysis Filter wheel order. Fixed for 3500 at "GATC"',
    "GTyp1": "Gel or polymer Type",
    "InSc1": "Injection time (seconds)",
    "InVt1": "Injection voltage (volts)",
    "LANE1": "Lane/Capillary",
    "LIMS1": "Sample tracking ID",
    "LNTD1": "Length to detector",
    "LsrP1": "Laser Power setting (micro Watts)",
    "MCHN1": "Instrument name and serial number",
    "MODF1": "Data collection module file",
    "MODL1": "Model number",
    "NAVG1": "Pixels averaged per lane",
    "NLNE1": "Number of capillaries",
    "OfSc1": "List of scans that are marked off scale in Collection. (optional)",
    # OvrI and OrvV are listed as "1-N", and "One for each dye (unanalyzed
    # and/or analyzed data)"
    "OvrI1": "List of scan number indexes that have values greater than 32767 but did not "
             "saturate the camera. In Genemapper samples, this can have indexes with "
             "values greater than 32000. In sequencing samples, this cannot have "
             "indexes with values greater than 32000.",
    "OvrI2": "List of scan number indexes that have values greater than 32767 but did not "
             "saturate the camera. In Genemapper samples, this can have indexes with "
             "values greater than 32000. In sequencing samples, this cannot have "
             "indexes with values greater than 32000.",
    "OvrI3": "List of scan number indexes that have values greater than 32767 but did not "
             "saturate the camera. In Genemapper samples, this can have indexes with "
             "values greater than 32000. In sequencing samples, this cannot have "
             "indexes with values greater than 32000.",
    "OvrI4": "List of scan number indexes that have values greater than 32767 but did not "
             "saturate the camera. In Genemapper samples, this can have indexes with "
             "values greater than 32000. In sequencing samples, this cannot have "
             "indexes with values greater than 32000.",
    "OvrV1": "List of color data values found at the locations listed in the OvrI tag. "
             "There must be exactly as many numbers in this array as in the OvrI array.",
    "OvrV2": "List of color data values found at the locations listed in the OvrI tag. "
             "There must be exactly as many numbers in this array as in the OvrI array.",
    "OvrV3": "List of color data values found at the locations listed in the OvrI tag. "
             "There must be exactly as many numbers in this array as in the OvrI array.",
    "OvrV4": "List of color data values found at the locations listed in the OvrI tag. "
             "There must be exactly as many numbers in this array as in the OvrI array.",
    "PDMF1": "Sequencing Analysis Mobility file name chosen in collection",
    "RMXV1": "Run Module XML schema version",
    "RMdN1": "Run Module name (same as MODF)",
    "RMdX1": "Run Module XML string",
    "RPrN1": "Run Protocol name",
    "RPrV1": "Run Protocol version",
    "RUND1": "Run Started Date",
    "RUND2": "Run Stopped Date",
    "RUND3": "Data Collection Started Date",
    "RUND4": "Data Collection Stopped date",
    "RUNT1": "Run Started Time",
    "RUNT2": "Run Stopped Time",
    "RUNT3": "Data Collection Started Time",
    "RUNT4": "Data Collection Stopped Time",
    "Rate1": "Scanning Rate. Milliseconds per frame.",
    "RunN1": "Run Name",
    "SCAN1": "Number of scans",
    "SMED1": "Polymer lot expiration date",
    "SMLt1": "Polymer lot number",
    "SMPL1": "Sample name",
    "SVER1": "Data collection software version",
    "SVER3": "Data collection firmware version",
    "Satd1": "Array of longs representing the scan numbers of data points, which are flagged as saturated by data collection (optional)",
    "Scal1": "Rescaling divisor for color data",
    "Scan1": "Number of scans (legacy - use SCAN)",
    "TUBE1": "Well ID",
    "Tmpr1": "Run temperature setting",
    "User1": "Name of user who created the plate (optional)",
}

#  No instrument specific tags
# _INSTRUMENT_SPECIFIC_TAGS['abi_prism_3100/3100-Avant'] = {
# }

_INSTRUMENT_SPECIFIC_TAGS["abi_3130/3130xl"] = {
    "CTOw1": "Container owner",
    "HCFG1": "Instrument Class",
    "HCFG2": "Instrument Family",
    "HCFG3": "Official Instrument Name",
    "HCFG4": "Instrument Parameters",
    "RMdVa1": "Run Module version",
}

_INSTRUMENT_SPECIFIC_TAGS["abi_3530/3530xl"] = {
    "AAct1": "Primary Analysis Audit Active indication. True if system auditing was enabled during the last write of this file, "
             "false if system auditing was disabled.",
    "ABED1": "Anode buffer expiration date using ISO 8601 format using the patterns YYYY-MM-DDTHH:MM:SS.ss+/-HH:MM. Hundredths of a second are optional.",
    "ABID1": "Anode buffer tray first installed date",
    "ABLt1": "Anode buffer lot number",
    "ABRn1": "Number of runs (injections) processed with the current Anode Buffer (runs allowed - runs remaining)",
    "ABTp1": "Anode buffer type",
    "AEPt1": "Analysis Ending scan number for basecalling on initial analysis",
    "AEPt2": "Analysis Ending scan number for basecalling on last analysis",
    "APCN1": "Amplicon name",
    "ARTN1": "Analysis Return code. Produced only by 5 Prime basecaller 1.0b3",
    "ASPF1": "Flag to indicate whether adaptive processing worked or not",
    "ASPt1": "Analysis Starting scan number for first analysis",
    "ASPt2": "Analysis Starting scan number for last analysis",
    "AUDT2": "Audit log used across 3500 software (optional)",
    "AVld1": "Assay validation flag (true or false)",
    "AmbT1": "Record of ambient temperature readings",
    "AsyC1": "The assay contents (xml format)",
    "AsyN1": "The assay name",
    "AsyV1": "The assay version",
    "B1Pt1": "Reference scan number for mobility and spacing curves for first analysis",
    "B1Pt2": "Reference scan number for mobility and spacing curves for last analysis",
    "BCTS1": "Basecaller timestamp. Time of completion of most recent analysis",
    "BcRn1": "Basecalling qc code",
    "BcRs1": "Basecalling warnings, a concatenated comma separated string",
    "BcRs2": "Basecalling errors, a concatenated comma separated string",
    "CAED1": "Capillary array expiration",
    "CALt1": "Capillary array lot number",
    "CARn1": "Number of injections processed (including the one of which this sample was a part) through the capillary array",
    "CASN1": "Capillary array serial number",
    "CBED1": "Cathode buffer expiration date",
    "CBID1": "Cathode buffer tray first installed date",
    "CBLt1": "Cathode buffer lot number",
    "CBRn1": "Number of runs (injections) processed with the current Cathode Buffer (runs allowed - runs remaining)",
    "CBTp1": "Cathode buffer type",
    "CLRG1": "Start of the clear range (inclusive).",
    "CLRG2": "Clear range length",
    "CRLn1": "Contiguous read length",
    "CRLn2": 'One of "Pass", "Fail", or "Check"',
    "CTOw1": "The name entered as the Owner of a plate, in the plate editor",
    "CkSm1": "File checksum",
    "DCEv1": "A list of door-close events, separated by semicolon. Door open events are generally paired with door close events.",
    "DCHT1": "Reserved for backward compatibility. The detection cell heater temperature setting from the Run Module. Not used for 3500.",
    "DOEv1": "A list of door-open events, separated by semicolon. Door close events are generally paired with door open events.",
    "ESig2": "Electronic signature record used across 3500 software",
    "FTab1": "Feature table. Can be created by Nibbler for Clear Range.",
    "FVoc1": "Feature table vocabulary. Can be created by Nibbler for Clear Range.",
    "Feat1": "Features. Can be created by Nibbler for Clear Range.",
    "HCFG1": "The Instrument Class. All upper case, no spaces. Initial valid value: CE",
    "HCFG2": "The Instrument Family. All upper case, no spaces. Valid values: 31XX or 37XX for UDC, 35XX (for 3500)",
    "HCFG3": "The official instrument name. Mixed case, minus any special formatting. Initial valid values: 3130, 3130xl, 3730, 3730xl, 3500, 3500xl.",
    "HCFG4": "Instrument parameters. Contains key-value pairs of instrument configuration information, separated by semicolons. "
             "Four parameters are included initially: UnitID=<UNITD number>, CPUBoard=<board type>, "
             "ArraySize=<# of capillaries>, SerialNumber=<Instrument Serial#>.",
    "InjN1": "Injection name",
    "LAST1": "Parameter settings information",
    "NOIS1": "The estimate of rms baseline noise (S/N ratio) for each dye for a successfully analyzed sample. "
             "Corresponds in order to the raw data in tags DATA 1-4. KB basecaller only.",
    "P1AM1": "Amplitude of primary peak, which is not necessarily equal to corresponding signal strength at that position",
    "P1RL1": "Deviation of primary peak position from (PLoc,2), times 100, rounded to integer",
    "P1WD1": "Full-width Half-max of primary peak, times 100, rounded to integer. "
             "Corresponding signal intensity is not necessarily equal to one half of primary peak amplitude",
    "P2AM1": "Amplitude of secondary peak, which is not necessarily equal to corresponding signal strength at that position",
    "P2BA1": "Base of secondary peak",
    "P2RL1": "Deviation of secondary peak position from (PLoc,2), times 100, rounded to integer",
    "PBAS1": "Array of sequence characters edited by user",
    "PBAS2": "Array of sequence characters as called by Basecaller",
    "PCON1": "Array of quality Values (0-255) as edited by user",
    "PCON2": "Array of quality values (0-255) as called by Basecaller",
    "PDMF2": "Mobility file name chosen in most recent analysis (identical to PDMF1)",
    "PLOC1": "Array of peak locations edited by user",
    "PLOC2": "Array of peak locations as called by Basecaller",
    "PRJT1": "SeqScape 2.0 project template name",
    "PROJ4": "SeqScape 2.0 project name",
    "PSZE1": "Plate size. The number of sample positions in the container. Current allowed values: 96, 384.",
    "PTYP1": "Plate type. Current allowed values: 96-Well, 384-Well.",
    "PuSc1": "Median pupscore",
    "QV201": "QV20+ value",
    "QV202": 'One of "Pass", "Fail", or "Check"',
    "QcPa1": "QC parameters",
    "QcRn1": "Trimming and QC code",
    "QcRs1": "QC warnings, a concatenated comma separated string",
    "QcRs2": "QC errors, a concatenated comma separated string",
    "RGOw1": "The name entered as the Owner of a Results Group, in the Results Group Editor. Implemented as the user name from the results group.",
    "RInj1": "Reinjection number. The reinjection number that this sample belongs to. Not present if there was no reinjection.",
    "RNmF1": "Raman normalization factor",
    "RevC1": "for whether the sequence has been complemented",
    "RunN1": "Run name (which, for 3500, is different from injection name)",
    "S/N%1": "Signal strength for each dye",
    "SMID1": "Polymer first installed date",
    "SMRn1": "Number of runs (injections) processed with the current polymer (runs allowed - runs remaining)",
    "SPAC1": "Average peak spacing used in last analysis",
    "SPAC2": "Basecaller name - corresponds to name of bcp file.",
    "SPAC3": "Average peak spacing last calculated by the Basecaller.",
    "SPEC1": "Sequencing Analysis Specimen Name",
    "SVER2": "Basecaller version number",
    "SVER4": "Sample File Format Version String",
    "ScPa1": "The parameter string of size caller",
    "ScSt1": "Raw data start point. Set to 0 for 3500 data collection.",
    "SpeN1": "Active spectral calibration name",
    "TrPa1": "Trimming parameters",
    "TrSc1": "Trace score.",
    "TrSc2": 'One of "Pass", "Fail", or "Check"',
    "phAR1": "Trace peak aria ratio",
    "phCH1": 'Chemistry type ("term", "prim", "unknown"), based on DYE_1 information',
    "phDY1": 'Dye ("big", "d-rhod", "unknown"), based on mob file information',
    "phQL1": "Maximum Quality Value",
    "phTR1": "Set Trim region",
    "phTR2": "Trim probability",
}

_INSTRUMENT_SPECIFIC_TAGS["abi_3730/3730xl"] = {
    "BufT1": "Buffer tray heater temperature (degrees C)",
}
# fmt: on

# dictionary for data unpacking format
_BYTEFMT = {
    1: "b",  # byte
    2: "s",  # char
    3: "H",  # word
    4: "h",  # short
    5: "i",  # long
    6: "2i",  # rational, legacy unsupported
    7: "f",  # float
    8: "d",  # double
    10: "h2B",  # date
    11: "4B",  # time
    12: "2i2b",  # thumb
    13: "B",  # bool
    14: "2h",  # point, legacy unsupported
    15: "4h",  # rect, legacy unsupported
    16: "2i",  # vPoint, legacy unsupported
    17: "4i",  # vRect, legacy unsupported
    18: "s",  # pString
    19: "s",  # cString
    20: "2i",  # tag, legacy unsupported
}
# header data structure (excluding 4 byte ABIF marker)
_HEADFMT = ">H4sI2H3I"
# directory data structure
_DIRFMT = ">4sI2H4I"

__global_tag_listing = []
for tag in _INSTRUMENT_SPECIFIC_TAGS.values():
    __global_tag_listing += tag.keys()


def _get_string_tag(opt_bytes_value, default=None):
    """Return the string value of the given an optional raw bytes tag value.

    If the bytes value is None, return the given default value.

    """
    if opt_bytes_value is None:
        return default
    try:
        return opt_bytes_value.decode()
    except UnicodeDecodeError:
        return opt_bytes_value.decode(encoding=sys.getdefaultencoding())


class AbiIterator(SequenceIterator):
    """Parser for Abi files."""

    def __init__(self, source, trim=False):
        """Return an iterator for the Abi file format."""
        self.trim = trim
        super().__init__(source, mode="b", fmt="ABI")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        # check if input file is a valid Abi file
        marker = handle.read(4)
        if not marker:
            # handle empty file gracefully
            raise ValueError("Empty file.")

        if marker != b"ABIF":
            raise OSError(f"File should start ABIF, not {marker!r}")
        records = self.iterate(handle)
        return records

    def iterate(self, handle):
        """Parse the file and generate SeqRecord objects."""
        # dirty hack for handling time information
        times = {"RUND1": "", "RUND2": "", "RUNT1": "", "RUNT2": ""}

        # initialize annotations
        annot = dict(zip(_EXTRACT.values(), [None] * len(_EXTRACT)))

        # parse header and extract data from directories
        header = struct.unpack(_HEADFMT, handle.read(struct.calcsize(_HEADFMT)))

        # Set default sample ID value, which we expect to be present in most
        # cases in the SMPL1 tag, but may be missing.
        sample_id = "<unknown id>"

        raw = {}
        seq = qual = None
        for tag_name, tag_number, tag_data in _abi_parse_header(header, handle):
            key = tag_name + str(tag_number)

            raw[key] = tag_data

            # PBAS2 is base-called sequence, only available in 3530
            if key == "PBAS2":
                seq = tag_data.decode()
            # PCON2 is quality values of base-called sequence
            elif key == "PCON2":
                qual = [ord(val) for val in tag_data.decode()]
            # SMPL1 is sample id entered before sequencing run, it must be
            # a string.
            elif key == "SMPL1":
                sample_id = _get_string_tag(tag_data)
            elif key in times:
                times[key] = tag_data
            else:
                if key in _EXTRACT:
                    annot[_EXTRACT[key]] = tag_data

        # set time annotations
        annot["run_start"] = f"{times['RUND1']} {times['RUNT1']}"
        annot["run_finish"] = f"{times['RUND2']} {times['RUNT2']}"

        # raw data (for advanced end users benefit)
        annot["abif_raw"] = raw

        # fsa check
        is_fsa_file = all(tn not in raw for tn in ("PBAS1", "PBAS2"))

        if is_fsa_file:
            try:
                file_name = basename(handle.name).replace(".fsa", "")
            except AttributeError:
                file_name = ""

            sample_id = _get_string_tag(raw.get("LIMS1"), sample_id)
            description = _get_string_tag(raw.get("CTID1"), "<unknown description>")
            record = SeqRecord(
                Seq(""),
                id=sample_id,
                name=file_name,
                description=description,
                annotations=annot,
            )

        else:
            # use the file name as SeqRecord.name if available
            try:
                file_name = basename(handle.name).replace(".ab1", "")
            except AttributeError:
                file_name = ""
            record = SeqRecord(
                Seq(seq),
                id=sample_id,
                name=file_name,
                description="",
                annotations=annot,
            )
        if qual:
            # Expect this to be missing for FSA files.
            record.letter_annotations["phred_quality"] = qual
        elif not is_fsa_file and not qual and self.trim:
            raise ValueError(
                "The 'abi-trim' format can not be used for files without"
                " quality values."
            )

        if self.trim and not is_fsa_file:
            record = _abi_trim(record)

        record.annotations["molecule_type"] = "DNA"
        yield record


def _AbiTrimIterator(handle):
    """Return an iterator for the Abi file format that yields trimmed SeqRecord objects (PRIVATE)."""
    return AbiIterator(handle, trim=True)


def _abi_parse_header(header, handle):
    """Return directory contents (PRIVATE)."""
    # header structure (after ABIF marker):
    # file version, tag name, tag number,
    # element type code, element size, number of elements
    # data size, data offset, handle (not file handle)
    head_elem_size = header[4]
    head_elem_num = header[5]
    head_offset = header[7]
    index = 0

    while index < head_elem_num:
        start = head_offset + index * head_elem_size
        # add directory offset to tuple
        # to handle directories with data size <= 4 bytes
        handle.seek(start)
        dir_entry = struct.unpack(_DIRFMT, handle.read(struct.calcsize(_DIRFMT))) + (
            start,
        )
        index += 1
        # only parse desired dirs
        key = dir_entry[0].decode()
        key += str(dir_entry[1])

        tag_name = dir_entry[0].decode()
        tag_number = dir_entry[1]
        elem_code = dir_entry[2]
        elem_num = dir_entry[4]
        data_size = dir_entry[5]
        data_offset = dir_entry[6]
        tag_offset = dir_entry[8]
        # if data size <= 4 bytes, data is stored inside tag
        # so offset needs to be changed
        if data_size <= 4:
            data_offset = tag_offset + 20
        handle.seek(data_offset)
        data = handle.read(data_size)
        yield tag_name, tag_number, _parse_tag_data(elem_code, elem_num, data)


def _abi_trim(seq_record):
    """Trims the sequence using Richard Mott's modified trimming algorithm (PRIVATE).

    Arguments:
        - seq_record - SeqRecord object to be trimmed.

    Trimmed bases are determined from their segment score, which is a
    cumulative sum of each base's score. Base scores are calculated from
    their quality values.

    More about the trimming algorithm:
    http://www.phrap.org/phredphrap/phred.html
    http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Quality_trimming.html
    """
    start = False  # flag for starting position of trimmed sequence
    segment = 20  # minimum sequence length
    trim_start = 0  # init start index
    cutoff = 0.05  # default cutoff value for calculating base score

    if len(seq_record) <= segment:
        return seq_record
    else:
        # calculate base score
        score_list = [
            cutoff - (10 ** (qual / -10.0))
            for qual in seq_record.letter_annotations["phred_quality"]
        ]

        # calculate cumulative score
        # if cumulative value < 0, set it to 0
        # first value is set to 0, because of the assumption that
        # the first base will always be trimmed out
        cummul_score = [0]
        for i in range(1, len(score_list)):
            score = cummul_score[-1] + score_list[i]
            if score < 0:
                cummul_score.append(0)
            else:
                cummul_score.append(score)
                if not start:
                    # trim_start = value when cumulative score is first > 0
                    trim_start = i
                    start = True

        # trim_finish = index of highest cumulative score,
        # marking the end of sequence segment with highest cumulative score
        trim_finish = cummul_score.index(max(cummul_score))

        return seq_record[trim_start:trim_finish]


def _parse_tag_data(elem_code, elem_num, raw_data):
    """Return single data value (PRIVATE).

    Arguments:
     - elem_code - What kind of data
     - elem_num - How many data points
     - raw_data - abi file object from which the tags would be unpacked

    """
    if elem_code in _BYTEFMT:
        # because '>1s' unpack differently from '>s'
        if elem_num == 1:
            num = ""
        else:
            num = str(elem_num)
        fmt = ">" + num + _BYTEFMT[elem_code]

        assert len(raw_data) == struct.calcsize(fmt)
        data = struct.unpack(fmt, raw_data)

        # no need to use tuple if len(data) == 1
        # also if data is date / time
        if elem_code not in [10, 11] and len(data) == 1:
            data = data[0]

        # account for different data types
        if elem_code == 2:
            return data
        elif elem_code == 10:
            return str(datetime.date(*data))
        elif elem_code == 11:
            return str(datetime.time(*data[:3]))
        elif elem_code == 13:
            return bool(data)
        elif elem_code == 18:
            return data[1:]
        elif elem_code == 19:
            return data[:-1]
        else:
            return data
    else:
        return None


if __name__ == "__main__":
    pass
