# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests of the transcription and translation methods of Seq objects."""

import unittest

from Bio import Seq
from Bio.Data.CodonTable import TranslationError


class TestTranscriptionTranslation(unittest.TestCase):
    def test_transcription(self):
        s = "ATA"
        dna = Seq.Seq(s)
        rna = dna.transcribe()
        self.assertEqual(rna, "AUA")
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
        dna = Seq.Seq(s)
        rna = dna.transcribe()
        self.assertEqual(
            rna,
            "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU",
        )
        s = "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
        rna = Seq.Seq(s)
        dna = rna.back_transcribe()
        self.assertEqual(
            dna,
            "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT",
        )

    def test_translation(self):
        s = ""
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        self.assertEqual(protein, "")
        s = "TAA"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        self.assertEqual(protein, "")
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCA"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        self.assertEqual(protein, "ENSFSLDFL")
        s = "GAA"
        dna = Seq.Seq(s)
        protein = dna.translate(15, to_stop=True)
        self.assertEqual(protein, "E")
        s = "ATA"
        dna = Seq.Seq(s)
        protein = dna.translate("Vertebrate Mitochondrial", to_stop=True)
        self.assertEqual(protein, "M")
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATAT"
        dna = Seq.Seq(s)
        protein = dna.translate("SGC8", to_stop=True)
        self.assertEqual(protein, "ENSFSLDFLWNPSPSNDAWDSSY")

    def test_dna_rna_translation(self):
        s = "TCAAAAAGGTGCATCTAGATG"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        self.assertEqual(protein, "SKRCI")
        gapped_protein = dna.translate()
        self.assertEqual(gapped_protein, "SKRCI*M")
        # The table used here has "AGG" as a stop codon:
        p2 = dna.translate(table=2, to_stop=True)
        self.assertEqual(p2, "SK")
        p2 = dna.translate(table=2)
        self.assertEqual(p2, "SK*CI*M")
        p2 = dna.translate(table=2, stop_symbol="+")
        self.assertEqual(p2, "SK+CI+M")
        r = s.replace("T", "U")
        rna = Seq.Seq(r)
        protein = rna.translate(to_stop=True)
        self.assertEqual(protein, "SKRCI")
        gapped_protein = rna.translate()
        self.assertEqual(gapped_protein, "SKRCI*M")

    def test_ambiguous(self):
        s = "RATGATTARAATYTA"
        dna = Seq.Seq(s)
        protein = dna.translate("Vertebrate Mitochondrial")
        self.assertEqual(protein, "BD*NL")
        stop_protein = dna.translate("SGC1", to_stop=True)
        self.assertEqual(stop_protein, "BD")

    # ------------------------------------------------------------------
    # Gap translation tests (added by Biopython PR #5186).
    #
    # Context: Multiple-sequence alignment files (FASTA output from mafft,
    # muscle, etc.) contain gap characters ('-') used to maintain reading-
    # frame alignment across insertions and deletions.  Without this PR,
    # Biopython raises TranslationError on ALL partial-gap codons (any
    # codon containing '-' except the full-gap '---'), even when gap='-'
    # is explicitly provided.  This forces downstream tools to implement
    # per-codon workarounds — splitting 3822-character protein-coding
    # sequences into individual triplets, substituting '-' with 'N', and
    # calling translate() on each one separately.
    #
    # The PR resolves this by treating gap characters as 'N' (unknown
    # nucleotide) within partial-gap codons when gap='-' is provided,
    # leveraging Biopython's existing IUPAC ambiguity resolution tables.
    #
    # Side-by-side comparison of all 60 partial-gap codon patterns:
    #
    #   Without PR: ALL 60 partial-gap codons raise TranslationError
    #   With PR:    8 resolve to a specific amino acid (four-fold
    #               degenerate families), 52 resolve to 'X' (ambiguous),
    #               0 raise TranslationError
    #
    # Real-world use case: mutation_scatter_plot pipeline
    # (https://github.com/host-patho-evo/mutation_scatter_plot)
    # processes 64 GB GISAID SARS-CoV-2 spike protein alignments.
    # The per-codon workaround (splitting + caching) was the only way
    # to use Biopython's translate() with alternative genetic code tables
    # on alignment-padded data.  See also:
    #   - Biopython PR #4992 (original proposal, declined)
    #   - Biopython PR #4994 (gap default alignment, merged)
    #   - Biopython issue #5036
    # ------------------------------------------------------------------

    def test_gap_full_codon(self):
        """Full-gap codon '---' translates to '-' (pre-existing behavior)."""
        self.assertEqual(Seq.Seq("---").translate(gap="-"), "-")

    def test_gap_full_codon_in_sequence(self):
        """Full-gap codon embedded in a sequence."""
        self.assertEqual(Seq.Seq("TCA---TCG").translate(gap="-"), "S-S")

    # -- Four-fold degenerate families: third-position gap resolves to AA --

    def test_gap_TC_dash_serine(self):
        """TC- → TCN → S: all four TC[ACGT] encode Serine."""
        self.assertEqual(Seq.Seq("TC-").translate(gap="-"), "S")

    def test_gap_AC_dash_threonine(self):
        """AC- → ACN → T: all four AC[ACGT] encode Threonine."""
        self.assertEqual(Seq.Seq("AC-").translate(gap="-"), "T")

    def test_gap_CC_dash_proline(self):
        """CC- → CCN → P: all four CC[ACGT] encode Proline."""
        self.assertEqual(Seq.Seq("CC-").translate(gap="-"), "P")

    def test_gap_CG_dash_arginine(self):
        """CG- → CGN → R: all four CG[ACGT] encode Arginine."""
        self.assertEqual(Seq.Seq("CG-").translate(gap="-"), "R")

    def test_gap_CT_dash_leucine(self):
        """CT- → CTN → L: all four CT[ACGT] encode Leucine."""
        self.assertEqual(Seq.Seq("CT-").translate(gap="-"), "L")

    def test_gap_GC_dash_alanine(self):
        """GC- → GCN → A: all four GC[ACGT] encode Alanine."""
        self.assertEqual(Seq.Seq("GC-").translate(gap="-"), "A")

    def test_gap_GG_dash_glycine(self):
        """GG- → GGN → G: all four GG[ACGT] encode Glycine."""
        self.assertEqual(Seq.Seq("GG-").translate(gap="-"), "G")

    def test_gap_GT_dash_valine(self):
        """GT- → GTN → V: all four GT[ACGT] encode Valine."""
        self.assertEqual(Seq.Seq("GT-").translate(gap="-"), "V")

    # -- Ambiguous partial-gap codons: resolve to 'X' --

    def test_gap_AT_dash_ambiguous(self):
        """AT- → ATN → X: ATA/ATC/ATT=Ile vs ATG=Met → ambiguous."""
        self.assertEqual(Seq.Seq("AT-").translate(gap="-"), "X")

    def test_gap_AA_dash_ambiguous(self):
        """AA- → AAN → X: AAA=Lys vs AAC/AAT=Asn vs AAG=Lys → ambiguous."""
        self.assertEqual(Seq.Seq("AA-").translate(gap="-"), "X")

    def test_gap_first_position(self):
        """-TC → NTC → X: first position unknown → ambiguous."""
        self.assertEqual(Seq.Seq("-TC").translate(gap="-"), "X")

    def test_gap_second_position(self):
        """T-C → TNC → X: second position unknown → ambiguous."""
        self.assertEqual(Seq.Seq("T-C").translate(gap="-"), "X")

    def test_gap_two_dashes(self):
        """A-- → ANN → X: two unknowns → ambiguous."""
        self.assertEqual(Seq.Seq("A--").translate(gap="-"), "X")

    def test_gap_two_dashes_leading(self):
        """--A → NNA → X: two unknowns → ambiguous."""
        self.assertEqual(Seq.Seq("--A").translate(gap="-"), "X")

    # -- Without gap='-', partial-gap codons must still raise --

    def test_no_gap_kwarg_still_raises(self):
        """Without gap='-', string-level translate('TC-') must still raise.

        Note: Seq('TC-').translate() defaults to gap='-' since PR #4994,
        so this test uses the string-level function instead.
        """
        with self.assertRaises(TranslationError):
            Seq.translate("TC-")

    def test_no_gap_kwarg_full_gap_still_raises(self):
        """Without gap='-', string-level translate('---') must still raise."""
        with self.assertRaises(TranslationError):
            Seq.translate("---")

    # -- Out-of-phase gap triplets (the reason global substitution fails) --

    def test_out_of_phase_gap_triplet(self):
        """AA---A: gaps span codon boundary → AA-|--A → XX, not a gap codon.

        This demonstrates why .replace('---', '???').replace('-', 'N')
        .replace('???', '---') is incorrect: the '---' substring is not
        in phase with the codon frame.
        """
        self.assertEqual(Seq.Seq("AA---A").translate(gap="-"), "XX")

    # -- Mixed sequences: real-world alignment patterns --

    def test_mixed_alignment_sequence(self):
        """Realistic alignment: standard + full-gap + partial-gap codons.

        ATC=I, ---=-, GCA=A, AT-=X(ambiguous), TC-=S
        """
        seq = "ATC---GCAAT-TC-"
        self.assertEqual(Seq.Seq(seq).translate(gap="-"), "I-AXS")

    def test_spike_protein_deletion_pattern(self):
        """Pattern from SARS-CoV-2 BA.2.86 spike alignment.

        ATG=M, TCA=S, ---=-(deletion), GGC=G, TC-=S(four-fold)
        """
        seq = "ATGTCA---GGCTC-"
        self.assertEqual(Seq.Seq(seq).translate(gap="-"), "MS-GS")


class TestTranslationPerformance(unittest.TestCase):
    """Performance benchmarks: per-codon workaround vs. native translate().

    These benchmarks demonstrate the real-world performance cost of the
    per-codon workaround that downstream tools (e.g. mutation_scatter_plot)
    are forced to use because Biopython's translate() crashes on partial-gap
    codons.

    The workaround — adapted from mutation_scatter_plot.alt_translate() —
    splits sequences into individual triplets and calls translate() on
    each one, optionally with functools.lru_cache.

    Existing software tools that handle gapped alignment translation correctly:
    - SeaView (http://doua.prabi.fr/software/seaview) — "View as proteins"
    - MEGA (https://www.megasoftware.net/) — "Align by Codon"
    - Jalview (https://www.jalview.org/) — linked cDNA/protein views
    - MACSE (https://www.agap-ge2pop.org/macse/) — codon-aware alignment
    - BioEdit (https://bioedit.software.informer.com/) — manual editor
    - SnapGene (https://www.snapgene.com/) — commercial sequence editor
    - ApE (https://jorgensen.biology.utah.edu/wayned/ape/) — plasmid editor
    """

    @staticmethod
    def _alt_translate_no_cache(seq, table=1):
        """Per-codon workaround WITHOUT lru_cache (worst-case scenario)."""
        result = []
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            if not codon or len(codon) < 3:
                continue
            if "-" in codon and codon != "---":
                codon = codon.replace("-", "N")
            result.append(Seq.translate(codon, gap="-"))
        return "".join(result)

    @staticmethod
    def _make_cached_translator():
        """Per-codon workaround WITH lru_cache (production scenario)."""
        import functools

        @functools.lru_cache(maxsize=128)
        def _cached_translate_codon(codon, table=1):
            if "-" in codon and codon != "---":
                codon = codon.replace("-", "N")
            return Seq.translate(codon, gap="-")

        def _alt_translate_cached(seq, table=1):
            result = []
            for i in range(0, len(seq), 3):
                codon = seq[i : i + 3]
                if not codon or len(codon) < 3:
                    continue
                result.append(_cached_translate_codon(codon, table))
            return "".join(result)

        return _alt_translate_cached

    @staticmethod
    def _build_test_sequence(n_codons=1000, gap_fraction=0.05):
        """Build a realistic gapped alignment sequence.

        ~95% standard ACGT codons, ~5% partial-gap codons (typical of
        SARS-CoV-2 spike protein alignments from GISAID).
        """
        import random

        random.seed(42)
        bases = "ACGT"
        seq = []
        for _ in range(n_codons):
            if random.random() < gap_fraction:
                # Partial-gap codon (1 or 2 dashes)
                codon = [random.choice(bases) for _ in range(3)]
                n_dashes = random.choice([1, 2])
                positions = random.sample(range(3), n_dashes)
                for p in positions:
                    codon[p] = "-"
                seq.append("".join(codon))
            else:
                seq.append("".join(random.choice(bases) for _ in range(3)))
        return "".join(seq)

    def test_workaround_correctness(self):
        """Verify the per-codon workaround produces same results as native."""
        seq = self._build_test_sequence(n_codons=200)
        native_result = str(Seq.Seq(seq).translate(gap="-"))
        workaround_result = self._alt_translate_no_cache(seq)
        self.assertEqual(native_result, workaround_result)

    def test_cached_workaround_correctness(self):
        """Verify cached per-codon workaround produces same results."""
        seq = self._build_test_sequence(n_codons=200)
        native_result = str(Seq.Seq(seq).translate(gap="-"))
        cached_fn = self._make_cached_translator()
        cached_result = cached_fn(seq)
        self.assertEqual(native_result, cached_result)

    def test_performance_comparison(self):
        """Benchmark: native vs per-codon vs per-codon+cache.

        This test measures wall-clock time for translating a realistic
        1274-codon spike protein sequence 100 times (simulating batch
        processing of many records).

        Expected results (approximate):
          native translate():        1× (baseline, this PR)
          per-codon + lru_cache:    ~5-10× slower
          per-codon without cache: ~50-100× slower
        """
        import timeit

        # Build a spike-protein-length sequence (1274 codons = 3822 nt)
        seq = self._build_test_sequence(n_codons=1274, gap_fraction=0.03)
        n_repeats = 100

        # 1. Native translate() (this PR)
        def native():
            return str(Seq.Seq(seq).translate(gap="-"))

        t_native = timeit.timeit(native, number=n_repeats)

        # 2. Per-codon workaround WITHOUT cache
        def workaround_no_cache():
            return self._alt_translate_no_cache(seq)

        t_no_cache = timeit.timeit(workaround_no_cache, number=n_repeats)

        # 3. Per-codon workaround WITH lru_cache
        cached_fn = self._make_cached_translator()
        # Warm up the cache
        cached_fn(seq)

        def workaround_cached():
            return cached_fn(seq)

        t_cached = timeit.timeit(workaround_cached, number=n_repeats)

        # Report results
        print(
            f"\n  Performance comparison ({n_repeats} iterations, "
            f"{len(seq)} nt sequence):"
        )
        print(f"    Native translate():       {t_native:.3f}s (1.0×)")
        if t_native > 0:
            print(
                f"    Per-codon (no cache):     {t_no_cache:.3f}s "
                f"({t_no_cache / t_native:.1f}×)"
            )
            print(
                f"    Per-codon (lru_cache):    {t_cached:.3f}s "
                f"({t_cached / t_native:.1f}×)"
            )

        # Native should be faster than uncached workaround
        self.assertLess(
            t_native,
            t_no_cache,
            "Native translate() should be faster than "
            "per-codon workaround without cache",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
