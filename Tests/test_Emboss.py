# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Runs a few EMBOSS tools to check our wrappers and parsers."""

import os
import sys
import unittest
import subprocess
from io import StringIO

from Bio.Emboss.Applications import WaterCommandline, NeedleCommandline
from Bio.Emboss.Applications import SeqretCommandline, SeqmatchallCommandline
from Bio import SeqIO
from Bio import AlignIO
from Bio import MissingExternalDependencyError
from Bio.Application import _escape_filename
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord


# ###############################################################

# Try to avoid problems when the OS is in another language
os.environ["LANG"] = "C"

exes_wanted = ["water", "needle", "seqret", "transeq", "seqmatchall", "embossversion"]
exes = {}  # Dictionary mapping from names to exe locations

if "EMBOSS_ROOT" in os.environ:
    # Windows default installation path is C:\mEMBOSS which contains the exes.
    # EMBOSS also sets an environment variable which we will check for.
    path = os.environ["EMBOSS_ROOT"]
    if os.path.isdir(path):
        for name in exes_wanted:
            if os.path.isfile(os.path.join(path, name + ".exe")):
                exes[name] = os.path.join(path, name + ".exe")
        del name
    else:
        raise MissingExternalDependencyError(
            "$EMBOSS_ROOT=%r which does not exist!" % path
        )
    del path
if sys.platform != "win32":
    from subprocess import getoutput

    for name in exes_wanted:
        # This will "just work" if installed on the path as normal on Unix
        output = getoutput("%s -help" % name)
        if "not found" not in output and "not recognized" not in output:
            exes[name] = name
        del output
    del name

if len(exes) < len(exes_wanted):
    raise MissingExternalDependencyError(
        "Install EMBOSS if you want to use Bio.Emboss."
    )


def get_emboss_version():
    """Return a tuple of three ints, e.g. (6,1,0)."""
    # Windows and Unix versions of EMBOSS seem to differ in
    # which lines go to stdout and stderr - so merge them.
    child = subprocess.Popen(
        _escape_filename(exes["embossversion"]),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        shell=(sys.platform != "win32"),
    )
    stdout, stderr = child.communicate()
    child.stdout.close()  # This is both stdout and stderr
    del child
    assert stderr is None  # Send to stdout instead
    for line in stdout.split("\n"):
        if line.strip() == "Report the current EMBOSS version number":
            # e.g.
            # $ embossversion
            # Report the current EMBOSS version number
            # 6.5.7.0
            pass
        elif line.strip() == "Reports the current EMBOSS version number":
            # e.g.
            # $ embossversion
            # Reports the current EMBOSS version number
            # 6.3.1
            pass
        elif line.startswith("Writes the current EMBOSS version number"):
            pass
        elif line.count(".") == 2:
            return tuple(int(v) for v in line.strip().split("."))
        elif line.count(".") == 3:
            # e.g. I installed mEMBOSS-6.2.0.1-setup.exe
            # which reports 6.2.0.1 - for this return (6,2,0)
            return tuple(int(v) for v in line.strip().split("."))[:3]
        else:
            # Either we can't understand the output, or this is really
            # an error message not caught earlier (e.g. not in English)
            raise MissingExternalDependencyError(
                "Install EMBOSS if you want to use Bio.Emboss (%s)." % line
            )
    # In case there was no output at all...
    raise MissingExternalDependencyError("Could not get EMBOSS version")


# To avoid confusing known errors from old versions of EMBOSS ...
emboss_version = get_emboss_version()
if emboss_version < (6, 1, 0):
    raise MissingExternalDependencyError("Test requires EMBOSS 6.1.0 patch 3 or later.")

#################################################################


# Top level function as this makes it easier to use for debugging:
def emboss_piped_SeqIO_convert(records, old_format, new_format):
    """Run seqret, returns records (as a generator)."""
    # Setup, this assumes for all the format names used
    # Biopython and EMBOSS names are consistent!
    cline = SeqretCommandline(
        exes["seqret"],
        sformat=old_format,
        osformat=new_format,
        auto=True,  # no prompting
        filter=True,
    )
    # Run the tool,
    child = subprocess.Popen(
        str(cline),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        shell=(sys.platform != "win32"),
    )
    SeqIO.write(records, child.stdin, old_format)
    child.stdin.close()
    child.stderr.close()
    # TODO - Is there a nice way to return an iterator AND
    # automatically close the handle?
    records = list(SeqIO.parse(child.stdout, new_format))
    child.stdout.close()
    return records


# Top level function as this makes it easier to use for debugging:
def emboss_piped_AlignIO_convert(alignments, old_format, new_format):
    """Run seqret, returns alignments (as a generator)."""
    # Setup, this assumes for all the format names used
    # Biopython and EMBOSS names are consistent!
    cline = SeqretCommandline(
        exes["seqret"],
        sformat=old_format,
        osformat=new_format,
        auto=True,  # no prompting
        filter=True,
    )
    # Run the tool,
    with subprocess.Popen(
        str(cline),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        shell=(sys.platform != "win32"),
    ) as child:
        AlignIO.write(alignments, child.stdin, old_format)
        child.stdin.close()
        aligns = list(AlignIO.parse(child.stdout, new_format))
    return aligns


class SeqRetTests(unittest.TestCase):
    """Base class providing SeqRecord comparison method."""

    def compare_records(self, old_records, new_records, msg=None):
        self.assertEqual(len(old_records), len(new_records), msg)
        for old, new in zip(old_records, new_records):
            # Note the name matching is a bit fuzzy, e.g. truncation and
            # no spaces in PHYLIP files.
            self.assertTrue(
                (old.id in new.id)
                or (new.id in old.id)
                or (old.id.replace(" ", "_") == new.id.replace(" ", "_"))
                or (old.name == new.name),
                msg,
            )
            self.assertEqual(len(old.seq), len(new.seq), msg)
            if str(old.seq).upper() != str(new.seq).upper():
                if str(old.seq).replace("X", "N") == str(new.seq):
                    self.fail("%s: X -> N (protein forced into nucleotide?)" % msg)
                else:
                    self.assertEqual(old.seq, new.seq, msg)
            if old.features and new.features:
                self.assertEqual(len(old.features), len(new.features), msg)
            # TODO - check annotation


class SeqRetSeqIOTests(SeqRetTests):
    """Check EMBOSS seqret against Bio.SeqIO for converting files."""

    def tearDown(self):
        clean_up()

    def check_SeqIO_to_EMBOSS(self, in_filename, in_format, skip_formats=()):
        """Check SeqIO writes files seqret can read back."""
        records = list(SeqIO.parse(in_filename, in_format))
        for temp_format in ["genbank", "embl", "fasta"]:
            if temp_format in skip_formats:
                continue
            new_records = list(
                emboss_piped_SeqIO_convert(records, temp_format, "fasta")
            )
            msg = "converting %s from %s to %s" % (in_filename, in_format, temp_format)
            self.compare_records(records, new_records, msg)

    def check_EMBOSS_to_SeqIO(self, filename, old_format, skip_formats=()):
        """Check SeqIO can read read seqret's conversion output."""
        # TODO: Why can't we read EMBOSS's swiss output?
        self.assertTrue(os.path.isfile(filename))
        old_records = list(SeqIO.parse(filename, old_format))
        for new_format in ["genbank", "fasta", "pir", "embl", "ig"]:
            if new_format in skip_formats:
                continue
            cline = SeqretCommandline(
                exes["seqret"],
                sequence=filename,
                sformat=old_format,
                osformat=new_format,
                auto=True,  # no prompting
                stdout=True,
            )
            # Run the tool,
            with subprocess.Popen(
                str(cline),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                shell=(sys.platform != "win32"),
            ) as child:
                child.stdin.close()
                new_records = list(SeqIO.parse(child.stdout, new_format))
            msg = "converting %s from %s to %s" % (filename, old_format, new_format)
            self.compare_records(old_records, new_records, msg)

    def check_SeqIO_with_EMBOSS(self, filename, old_format, skip_formats=()):
        # Check EMBOSS can read Bio.SeqIO output...
        self.check_SeqIO_to_EMBOSS(filename, old_format, skip_formats)
        # Check Bio.SeqIO can read EMBOSS seqret output...
        self.check_EMBOSS_to_SeqIO(filename, old_format, skip_formats)

    def test_abi(self):
        """Check SeqIO agrees with EMBOSS' Abi to FASTQ conversion."""
        # This lets use check the id, sequence, and quality scores
        for filename in ["Abi/3730.ab1", "Abi/empty.ab1"]:
            old = SeqIO.read(filename, "abi")
            cline = SeqretCommandline(
                exes["seqret"],
                sequence=filename,
                sformat="abi",
                osformat="fastq-sanger",
                auto=True,  # no prompting
                stdout=True,
            )
            # Run the tool,
            with subprocess.Popen(
                str(cline),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                shell=(sys.platform != "win32"),
            ) as child:
                child.stdin.close()
                new = SeqIO.read(child.stdout, "fastq-sanger")
            if emboss_version == (6, 4, 0) and new.id == "EMBOSS_001":
                # Avoid bug in EMBOSS 6.4.0 (patch forthcoming)
                pass
            else:
                self.assertEqual(old.id, new.id)
            self.assertEqual(str(old.seq), str(new.seq))
            if emboss_version < (6, 3, 0) and new.letter_annotations[
                "phred_quality"
            ] == [1] * len(old):
                # Apparent bug in EMBOSS 6.2.0.1 on Windows
                pass
            else:
                self.assertEqual(old.letter_annotations, new.letter_annotations)

    def test_genbank(self):
        """Check SeqIO & EMBOSS reading each other's conversions of a GenBank file."""
        self.check_SeqIO_with_EMBOSS("GenBank/cor6_6.gb", "genbank")

    def test_genbank2(self):
        """Check SeqIO & EMBOSS reading each other's conversions of another GenBank file."""
        self.check_SeqIO_with_EMBOSS("GenBank/NC_000932.gb", "genbank")

    def test_embl(self):
        """Check SeqIO & EMBOSS reading each other's conversions of an EMBL file."""
        self.check_SeqIO_with_EMBOSS("EMBL/U87107.embl", "embl")

    def test_ig(self):
        """Check SeqIO & EMBOSS reading each other's conversions of an ig file."""
        # NOTE - EMBOSS considers "genbank" to be for nucleotides only,
        # and will turn "X" into "N" for GenBank output.
        self.check_SeqIO_to_EMBOSS(
            "IntelliGenetics/VIF_mase-pro.txt", "ig", skip_formats=["genbank", "embl"],
        )
        # TODO - What does a % in an ig sequence mean?
        # e.g. "IntelliGenetics/vpu_nucaligned.txt"
        # and  "IntelliGenetics/TAT_mase_nuc.txt"
        # EMBOSS seems to ignore them.

    def test_pir(self):
        """Check SeqIO & EMBOSS reading each other's conversions of a PIR file."""
        # Skip genbank here, EMBOSS mangles the LOCUS line:
        self.check_SeqIO_with_EMBOSS(
            "NBRF/clustalw.pir", "pir", skip_formats=["genbank"]
        )
        # Skip EMBL here, EMBOSS mangles the ID line
        # Skip GenBank, EMBOSS 6.0.1 on Windows won't output proteins as GenBank
        self.check_SeqIO_with_EMBOSS(
            "NBRF/DMB_prot.pir", "pir", skip_formats=["embl", "genbank"]
        )

    def test_clustalw(self):
        """Check SeqIO & EMBOSS reading each other's conversions of a Clustalw file."""
        self.check_SeqIO_with_EMBOSS(
            "Clustalw/hedgehog.aln", "clustal", skip_formats=["embl", "genbank"]
        )
        self.check_SeqIO_with_EMBOSS(
            "Clustalw/opuntia.aln", "clustal", skip_formats=["embl", "genbank"]
        )


class SeqRetAlignIOTests(SeqRetTests):
    """Check EMBOSS seqret against Bio.AlignIO for converting files."""

    def tearDown(self):
        clean_up()

    def compare_alignments(self, old_list, new_list, msg=None):
        self.assertEqual(len(old_list), len(new_list), msg)
        for old, new in zip(old_list, new_list):
            self.compare_records(old, new, msg)

    def check_EMBOSS_to_AlignIO(self, filename, old_format, skip_formats=()):
        """Check AlignIO can read seqret's conversion of the file."""
        self.assertTrue(os.path.isfile(filename), filename)
        old_aligns = list(AlignIO.parse(filename, old_format))
        formats = ["clustal", "phylip", "ig", "msf"]
        if len(old_aligns) == 1:
            formats.extend(["fasta", "nexus"])
        for new_format in formats:
            if new_format in skip_formats:
                continue
            cline = SeqretCommandline(
                exes["seqret"],
                sequence=filename,
                sformat=old_format,
                osformat=new_format,
                auto=True,  # no prompting
                stdout=True,
            )
            # Run the tool,
            with subprocess.Popen(
                str(cline),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                shell=(sys.platform != "win32"),
            ) as child:
                child.stdin.close()
                new_aligns = list(AlignIO.parse(child.stdout, new_format))
            msg = "converting %s from %s to %s" % (filename, old_format, new_format)
            self.compare_alignments(old_aligns, new_aligns, msg)

    def check_AlignIO_to_EMBOSS(self, in_filename, in_format, skip_formats=()):
        """Check Bio.AlignIO can write files seqret can read."""
        old_aligns = list(AlignIO.parse(in_filename, in_format))

        formats = ["clustal", "phylip"]
        if len(old_aligns) == 1:
            formats.extend(["fasta", "nexus"])
        for temp_format in formats:
            if temp_format in skip_formats:
                continue
            # PHYLIP is a simple format which explicitly supports
            # multiple alignments (unlike FASTA).
            try:
                new_aligns = list(
                    emboss_piped_AlignIO_convert(old_aligns, temp_format, "phylip")
                )
            except ValueError as e:
                self.assertIn(
                    str(e),
                    (
                        "Need the molecule type to be defined",
                        "Repeated name 'AT3G20900.' (originally 'AT3G20900.1-SEQ'), possibly due to truncation",
                        "Repeated name 'gi|1377497' (originally 'gi|13774975|gb|AAK39115.1|AF35'), possibly due to truncation",
                        "Repeated name 'gi_1393639' (originally 'gi_13936397_dbj_BAB47195.'), possibly due to truncation",
                    ),
                )
                continue
            msg = "converting %s from %s to %s" % (in_filename, in_format, temp_format)
            self.compare_alignments(old_aligns, new_aligns, msg)

    def check_AlignIO_with_EMBOSS(self, filename, old_format, skip_formats=()):
        # Check EMBOSS can read Bio.AlignIO output...
        self.check_AlignIO_to_EMBOSS(filename, old_format, skip_formats)
        # Check Bio.AlignIO can read EMBOSS seqret output...
        self.check_EMBOSS_to_AlignIO(filename, old_format, skip_formats)

    def test_align_clustalw(self):
        """Check AlignIO & EMBOSS reading each other's conversions of a ClustalW file."""
        self.check_AlignIO_with_EMBOSS("Clustalw/hedgehog.aln", "clustal")
        self.check_AlignIO_with_EMBOSS("Clustalw/opuntia.aln", "clustal")
        self.check_AlignIO_with_EMBOSS(
            "Clustalw/odd_consensus.aln", "clustal", skip_formats=["nexus"]
        )  # TODO - why not nexus?
        self.check_AlignIO_with_EMBOSS("Clustalw/protein.aln", "clustal")
        self.check_AlignIO_with_EMBOSS("Clustalw/promals3d.aln", "clustal")

    def test_clustalw(self):
        """Check AlignIO & EMBOSS reading each other's conversions of a PHYLIP file."""
        self.check_AlignIO_with_EMBOSS("Phylip/horses.phy", "phylip")
        self.check_AlignIO_with_EMBOSS("Phylip/hennigian.phy", "phylip")
        self.check_AlignIO_with_EMBOSS("Phylip/reference_dna.phy", "phylip")
        self.check_AlignIO_with_EMBOSS("Phylip/reference_dna2.phy", "phylip")
        self.check_AlignIO_with_EMBOSS("Phylip/interlaced.phy", "phylip")
        self.check_AlignIO_with_EMBOSS("Phylip/interlaced2.phy", "phylip")
        self.check_AlignIO_with_EMBOSS("Phylip/random.phy", "phylip")


class PairwiseAlignmentTests(unittest.TestCase):
    """Run pairwise alignments with water and needle, and parse them."""

    def tearDown(self):
        clean_up()

    def pairwise_alignment_check(self, query_seq, targets, alignments, local=True):
        """Check pairwise alignment data is sane."""
        # The datasets should be small, so making iterators into lists is OK
        targets = list(targets)
        alignments = list(alignments)
        self.assertEqual(len(targets), len(alignments))
        for target, alignment in zip(targets, alignments):
            self.assertEqual(len(alignment), 2)
            # self.assertEqual(target.id, alignment[1].id) #too strict
            if alignment[1].id not in target.id and alignment[1].id not in target.name:
                raise AssertionError(
                    "%s vs %s or %s" % (alignment[1].id, target.id, target.name)
                )
            if local:
                # Local alignment
                self.assertIn(str(alignment[0].seq).replace("-", ""), query_seq)
                self.assertIn(
                    str(alignment[1].seq).replace("-", "").upper(),
                    str(target.seq).upper(),
                )
            else:
                # Global alignment
                self.assertEqual(str(query_seq), str(alignment[0].seq).replace("-", ""))
                self.assertEqual(
                    str(target.seq).upper(),
                    str(alignment[1].seq).replace("-", "").upper(),
                )
        return True

    def run_water(self, cline):
        # Run the tool,
        stdout, stderr = cline()
        self.assertTrue(
            stderr.strip().startswith("Smith-Waterman local alignment"), stderr
        )
        if cline.outfile:
            self.assertEqual(stdout.strip(), "")
            self.assertTrue(
                os.path.isfile(cline.outfile),
                "Missing output file %r from:\n%s" % (cline.outfile, cline),
            )
        else:
            # Don't use this yet... could return stdout handle instead?
            return stdout

    def test_water_file(self):
        """Run water with the asis trick, output to a file."""
        # Setup, try a mixture of keyword arguments and later additions:
        cline = WaterCommandline(cmd=exes["water"], gapopen="10", gapextend="0.5")
        # Try using both human readable names, and the literal ones:
        cline.set_parameter("asequence", "asis:ACCCGGGCGCGGT")
        cline.set_parameter("-bsequence", "asis:ACCCGAGCGCGGT")
        # Try using a property set here:
        cline.outfile = "Emboss/temp with space.water"
        self.assertEqual(str(eval(repr(cline))), str(cline))
        # Run the tool,
        self.run_water(cline)
        # Check we can parse the output...
        align = AlignIO.read(cline.outfile, "emboss")
        self.assertEqual(len(align), 2)
        self.assertEqual(str(align[0].seq), "ACCCGGGCGCGGT")
        self.assertEqual(str(align[1].seq), "ACCCGAGCGCGGT")
        # Clean up,
        os.remove(cline.outfile)

    def test_water_piped(self):
        """Run water with asis trick, output piped to stdout."""
        cline = WaterCommandline(
            cmd=exes["water"],
            asequence="asis:ACCCGGGCGCGGT",
            bsequence="asis:ACCCGAGCGCGGT",
            gapopen=10,
            gapextend=0.5,
            auto=True,
            filter=True,
        )
        self.assertEqual(
            str(cline),
            exes["water"]
            + " -auto -filter"
            + " -asequence=asis:ACCCGGGCGCGGT"
            + " -bsequence=asis:ACCCGAGCGCGGT"
            + " -gapopen=10 -gapextend=0.5",
        )
        # Run the tool,
        child = subprocess.Popen(
            str(cline),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            shell=(sys.platform != "win32"),
        )
        child.stdin.close()
        # Check we could read it's output
        align = AlignIO.read(child.stdout, "emboss")
        self.assertEqual(len(align), 2)
        self.assertEqual(str(align[0].seq), "ACCCGGGCGCGGT")
        self.assertEqual(str(align[1].seq), "ACCCGAGCGCGGT")
        # Check no error output:
        self.assertEqual(child.stderr.read(), "")
        self.assertEqual(0, child.wait())
        child.stdout.close()
        child.stderr.close()

    def test_needle_file(self):
        """Run needle with the asis trick, output to a file."""
        # Setup,
        cline = NeedleCommandline(cmd=exes["needle"])
        cline.set_parameter("-asequence", "asis:ACCCGGGCGCGGT")
        cline.set_parameter("-bsequence", "asis:ACCCGAGCGCGGT")
        cline.set_parameter("-gapopen", "10")
        cline.set_parameter("-gapextend", "0.5")
        # EMBOSS would guess this, but let's be explicit:
        cline.set_parameter("-snucleotide", "True")
        cline.set_parameter("-outfile", "Emboss/temp with space.needle")
        self.assertEqual(str(eval(repr(cline))), str(cline))
        # Run the tool,
        stdout, stderr = cline()
        # Check it worked,
        self.assertTrue(
            stderr.strip().startswith("Needleman-Wunsch global alignment"), stderr
        )
        self.assertEqual(stdout.strip(), "")
        filename = cline.outfile
        self.assertTrue(
            os.path.isfile(filename),
            "Missing output file %r from:\n%s" % (filename, cline),
        )
        # Check we can parse the output...
        align = AlignIO.read(filename, "emboss")
        self.assertEqual(len(align), 2)
        self.assertEqual(str(align[0].seq), "ACCCGGGCGCGGT")
        self.assertEqual(str(align[1].seq), "ACCCGAGCGCGGT")
        # Clean up,
        os.remove(filename)

    def test_needle_piped(self):
        """Run needle with asis trick, output piped to stdout."""
        cline = NeedleCommandline(
            cmd=exes["needle"],
            asequence="asis:ACCCGGGCGCGGT",
            bsequence="asis:ACCCGAGCGCGGT",
            gapopen=10,
            gapextend=0.5,
            auto=True,
            filter=True,
        )
        self.assertEqual(
            str(cline),
            exes["needle"]
            + " -auto -filter"
            + " -asequence=asis:ACCCGGGCGCGGT"
            + " -bsequence=asis:ACCCGAGCGCGGT"
            + " -gapopen=10 -gapextend=0.5",
        )
        # Run the tool,
        child = subprocess.Popen(
            str(cline),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            shell=(sys.platform != "win32"),
        )
        child.stdin.close()
        # Check we could read it's output
        align = AlignIO.read(child.stdout, "emboss")
        self.assertEqual(len(align), 2)
        self.assertEqual(str(align[0].seq), "ACCCGGGCGCGGT")
        self.assertEqual(str(align[1].seq), "ACCCGAGCGCGGT")
        # Check no error output:
        self.assertEqual(child.stderr.read(), "")
        self.assertEqual(0, child.wait())
        child.stdout.close()
        child.stderr.close()

    def test_water_file2(self):
        """Run water with the asis trick and nucleotide FASTA file, output to a file."""
        # Setup,
        query = "ACACACTCACACACACTTGGTCAGAGATGCTGTGCTTCTTGGAAGCAAGGNCTCAAAGGCAAGGTGCACGCAGAGGGACGTTTGAGTCTGGGATGAAGCATGTNCGTATTATTTATATGATGGAATTTCACGTTTTTATG"
        out_file = "Emboss/temp_test2.water"
        in_file = "Fasta/f002"
        self.assertTrue(os.path.isfile(in_file))
        if os.path.isfile(out_file):
            os.remove(out_file)
        cline = WaterCommandline(cmd=exes["water"])
        cline.set_parameter("-asequence", "asis:%s" % query)
        cline.set_parameter("-bsequence", in_file)
        cline.set_parameter("-gapopen", "10")
        cline.set_parameter("-gapextend", "0.5")
        cline.set_parameter("-outfile", out_file)
        self.assertEqual(str(eval(repr(cline))), str(cline))
        # Run the tool,
        self.run_water(cline)
        # Check we can parse the output and it is sensible...
        self.pairwise_alignment_check(
            query,
            SeqIO.parse(in_file, "fasta"),
            AlignIO.parse(out_file, "emboss"),
            local=True,
        )
        # Clean up,
        os.remove(out_file)

    def test_water_file3(self):
        """Run water with the asis trick and GenBank file, output to a file."""
        # Setup,
        query = "TGTTGTAATGTTTTAATGTTTCTTCTCCCTTTAGATGTACTACGTTTGGA"
        out_file = "Emboss/temp_test3.water"
        in_file = "GenBank/cor6_6.gb"
        self.assertTrue(os.path.isfile(in_file))
        if os.path.isfile(out_file):
            os.remove(out_file)
        cline = WaterCommandline(cmd=exes["water"])
        cline.set_parameter("asequence", "asis:%s" % query)
        cline.set_parameter("bsequence", in_file)
        # TODO - Tell water this is a GenBank file!
        cline.set_parameter("gapopen", "1")
        cline.set_parameter("gapextend", "0.5")
        cline.set_parameter("outfile", out_file)
        self.assertEqual(str(eval(repr(cline))), str(cline))
        # Run the tool,
        self.run_water(cline)
        # Check we can parse the output and it is sensible...
        self.pairwise_alignment_check(
            query,
            SeqIO.parse(in_file, "genbank"),
            AlignIO.parse(out_file, "emboss"),
            local=True,
        )
        # Clean up,
        os.remove(out_file)

    def test_water_file4(self):
        """Run water with the asis trick and SwissProt file, output to a file."""
        # Setup,
        query = "DVCTGKALCDPVTQNIKTYPVKIENLRVMI"
        out_file = "Emboss/temp_test4.water"
        in_file = "SwissProt/P0A186.txt"
        self.assertTrue(os.path.isfile(in_file))
        if os.path.isfile(out_file):
            os.remove(out_file)
        cline = WaterCommandline(cmd=exes["water"])
        cline.set_parameter("-asequence", "asis:%s" % query)
        cline.set_parameter("-bsequence", in_file)
        # EMBOSS should work this out, but let's be explicit:
        cline.set_parameter("-sprotein", True)
        # TODO - Tell water this is a SwissProt file!
        cline.set_parameter("-gapopen", "20")
        cline.set_parameter("-gapextend", "5")
        cline.set_parameter("-outfile", out_file)
        self.assertEqual(str(eval(repr(cline))), str(cline))
        # Run the tool,
        self.run_water(cline)
        # Check we can parse the output and it is sensible...
        self.pairwise_alignment_check(
            query,
            SeqIO.parse(in_file, "swiss"),
            AlignIO.parse(out_file, "emboss"),
            local=True,
        )
        # Clean up,
        os.remove(out_file)

    def test_needle_piped2(self):
        """Run needle with asis trick, and nucleotide FASTA file, output piped to stdout."""
        # TODO - Support needle in Bio.Emboss.Applications
        # (ideally with the -auto and -filter arguments)
        # Setup,
        query = "ACACACTCACACACACTTGGTCAGAGATGCTGTGCTTCTTGGAA"
        cline = exes["needle"]
        cline += " -asequence asis:" + query
        cline += " -bsequence Fasta/f002"
        cline += " -auto"  # no prompting
        cline += " -filter"  # use stdout
        # Run the tool,
        child = subprocess.Popen(
            str(cline),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            shell=(sys.platform != "win32"),
        )
        child.stdin.close()
        # Check we can parse the output and it is sensible...
        self.pairwise_alignment_check(
            query,
            SeqIO.parse("Fasta/f002", "fasta"),
            AlignIO.parse(child.stdout, "emboss"),
            local=False,
        )
        # Check no error output:
        self.assertEqual(child.stderr.read(), "")
        self.assertEqual(0, child.wait())
        child.stdout.close()
        child.stderr.close()

    def test_water_needs_output(self):
        """Run water without output file or stdout/filter should give error."""
        cline = WaterCommandline(
            cmd=exes["water"],
            asequence="asis:ACCCGGGCGCGGT",
            bsequence="asis:ACCCGAGCGCGGT",
            gapopen=10,
            gapextend=0.5,
            auto=True,
        )
        self.assertTrue(cline.auto)
        self.assertTrue(not cline.stdout)
        self.assertTrue(not cline.filter)
        self.assertIsNone(cline.outfile)
        self.assertRaises(ValueError, str, cline)

    def test_needle_needs_output(self):
        """Run needle without output file or stdout/filter should give error."""
        cline = NeedleCommandline(
            cmd=exes["needle"],
            asequence="asis:ACCCGGGCGCGGT",
            bsequence="asis:ACCCGAGCGCGGT",
            gapopen=10,
            gapextend=0.5,
            auto=True,
        )
        self.assertTrue(cline.auto)
        self.assertTrue(not cline.stdout)
        self.assertTrue(not cline.filter)
        self.assertIsNone(cline.outfile)
        self.assertRaises(ValueError, str, cline)

    def test_seqtmatchall_piped(self):
        """Run seqmatchall with pair output piped to stdout."""
        cline = SeqmatchallCommandline(
            cmd=exes["seqmatchall"],
            sequence="Fasta/f002",
            aformat="pair",
            wordsize=9,
            auto=True,
            stdout=True,
        )
        self.assertEqual(
            str(cline),
            exes["seqmatchall"]
            + " -auto -stdout"
            + " -sequence=Fasta/f002"
            + " -wordsize=9 -aformat=pair",
        )
        # Run the tool,
        child = subprocess.Popen(
            str(cline),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            shell=(sys.platform != "win32"),
        )
        child.stdin.close()
        # Check we could read it's output
        for align in AlignIO.parse(child.stdout, "emboss"):
            self.assertEqual(len(align), 2)
            self.assertEqual(align.get_alignment_length(), 9)
        # Check no error output:
        self.assertEqual(child.stderr.read(), "")
        self.assertEqual(0, child.wait())
        child.stdout.close()
        child.stderr.close()


# Top level function as this makes it easier to use for debugging:
def emboss_translate(sequence, table=None, frame=None):
    """Call transeq, returns protein sequence as string."""
    # TODO - Support transeq in Bio.Emboss.Applications?
    # (doesn't seem worthwhile as Biopython can do translations)

    if not sequence:
        raise ValueError(sequence)

    # Setup,
    cline = exes["transeq"]

    if len(sequence) < 100:
        filename = None
        cline += " -sequence asis:%s" % sequence
    else:
        # There are limits on command line string lengths...
        # use a temp file instead.
        filename = "Emboss/temp_transeq.txt"
        SeqIO.write(SeqRecord(sequence, id="Test"), filename, "fasta")
        cline += " -sequence %s" % filename

    cline += " -auto"  # no prompting
    cline += " -filter"  # use stdout
    if table is not None:
        cline += " -table %s" % str(table)
    if frame is not None:
        cline += " -frame %s" % str(frame)
    # Run the tool,
    child = subprocess.Popen(
        str(cline),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        shell=(sys.platform != "win32"),
    )
    out, err = child.communicate()
    # Check no error output:
    if err != "":
        raise ValueError(str(cline) + "\n" + err)

    # Check we could read it's output
    record = SeqIO.read(StringIO(out), "fasta")

    if 0 != child.wait():
        raise ValueError(str(cline))

    if filename:
        os.remove(filename)
        if not record.id.startswith("Test"):
            raise ValueError(str(cline))
    else:
        if not record.id.startswith("asis"):
            raise ValueError(str(cline))
    return str(record.seq)


# Top level function as this makes it easier to use for debugging:
def check_translation(sequence, translation, table=None):
    if table is None:
        t = 1
    else:
        t = table
    if (
        translation != str(sequence.translate(t))
        or translation != str(translate(sequence, t))
        or translation != translate(str(sequence), t)
    ):
        # More details...
        for i, amino in enumerate(translation):
            codon = sequence[i * 3 : i * 3 + 3]
            if amino != str(codon.translate(t)):
                raise ValueError(
                    "%s -> %s not %s (table %s)" % (codon, amino, codon.translate(t), t)
                )
        # Shouldn't reach this line:
        raise ValueError("%s -> %s (table %s)" % (sequence, translation, t))
    return True


class TranslationTests(unittest.TestCase):
    """Run pairwise alignments with water and needle, and parse them."""

    def tearDown(self):
        clean_up()

    def test_simple(self):
        """Run transeq vs Bio.Seq for simple translations (including alt tables)."""
        examples = [
            Seq("ACGTGACTGACGTAGCATGCCACTAGG"),
            # Unamibguous TA? codons:
            Seq("TAATACTATTAG"),
            # Most of the ambiguous TA? codons:
            Seq("TANTARTAYTAMTAKTAHTABTADTAV"),
            # Problem cases,
            #
            # Seq("TAW"),
            # W = A or T, but EMBOSS does TAW -> X
            # TAA -> Y, TAT ->Y, so in Biopython TAW -> Y
            #
            # Seq("TAS"),
            # S = C or G, but EMBOSS does TAS -> Y
            # TAG -> *, TAC ->Y, so in Biopython TAS -> X (Y or *)
            #
            # Seq("AAS"),
            # On table 9, EMBOSS gives N, we give X.
            # S = C or G, so according to my reading of
            # table 9 on the NCBI page, AAC=N, AAG=K
            # suggesting this is a bug in EMBOSS.
            #
            Seq("ACGGGGGGGGTAAGTGGTGTGTGTGTAGT"),
        ]

        for sequence in examples:
            # EMBOSS treats spare residues differently... avoid this issue
            if len(sequence) % 3 != 0:
                sequence = sequence[: -(len(sequence) % 3)]
            self.assertEqual(len(sequence) % 3, 0)
            self.assertGreater(len(sequence), 0)
            self.check(sequence)

    def check(self, sequence):
        """Compare our translation to EMBOSS's using all tables.

        Takes a Seq object (and a filename containing it).
        """
        translation = emboss_translate(sequence)
        self.assertTrue(check_translation(sequence, translation))

        for table in [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23]:
            translation = emboss_translate(sequence, table)
            self.assertTrue(check_translation(sequence, translation, table))
        return True

    def translate_all_codons(self, letters):
        sequence = Seq(
            "".join(c1 + c3 + c3 for c1 in letters for c2 in letters for c3 in letters)
        )
        self.check(sequence)

    # def test_all_ambig_dna_codons(self):
    #    """transeq vs Bio.Seq on ambiguous DNA codons (inc. alt tables)."""
    #    self.translate_all_codons(ambiguous_dna_letters)

    def test_all_unambig_dna_codons(self):
        """Run transeq vs Bio.Seq on unambiguous DNA codons (inc. alt tables)."""
        self.translate_all_codons("ATCGatcg")

    def test_all_unambig_rna_codons(self):
        """Run transeq vs Bio.Seq on unambiguous RNA codons (inc. alt tables)."""
        self.translate_all_codons("AUCGaucg")

    def test_mixed_unambig_rna_codons(self):
        """Run transeq vs Bio.Seq on unambiguous DNA/RNA codons (inc. alt tables)."""
        self.translate_all_codons("ATUCGatucg")


def clean_up():
    """Fallback clean up method to remove temp files."""
    for filename in os.listdir("Emboss"):
        if filename.startswith("temp_"):
            try:
                os.remove(filename)
            except Exception:  # TODO - Which exceptions?
                pass


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
    clean_up()
