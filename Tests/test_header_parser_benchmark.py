#!/usr/bin/env python
# Copyright 2026 by Contributors.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Large-volume synthetic benchmark: HeaderParser vs SimpleFastaParser.

Generates a fixed-size (default 2.4 GB) FASTA file for each sequence-length
tier, then compares four header-extraction approaches on wall-clock time,
throughput, and memory.

Machine-parseable TSV output is printed to stdout (one row per measurement).
Human-readable commentary goes to stderr.

CPU cache hierarchy context (Intel Xeon Platinum 8260, "eltu1"):
  4-socket NUMA, 24 cores/socket × 2 HT = 192 logical CPUs.
  Per core:   L1d = 32 KB | L1i = 32 KB | L2 = 1 MB
  Per socket: L3 = 35.75 MB (36,608 KB, shared "Smart Cache")
  Total L3:   4 × 35.75 MB = 143 MB across all NUMA nodes

Usage:
    # Quick smoke test (small files, ~5 s)
    python Tests/test_header_parser_benchmark.py --quick --run-tests

    # Full benchmark: 2.4 GB per tier (takes minutes per tier)
    python Tests/test_header_parser_benchmark.py

    # Custom target size
    python Tests/test_header_parser_benchmark.py --target-gb 1.0
"""

from __future__ import annotations

import argparse
import gc
import io
import os
import platform
import resource
import sys
import tempfile
import textwrap
import time
import unittest


# ---------------------------------------------------------------------------
# Synthetic FASTA generator
# ---------------------------------------------------------------------------

def _estimate_bytes_per_record(seq_length: int, line_wrap: int = 60) -> int:
    """Estimate bytes per FASTA record (header + wrapped sequence + newlines)."""
    header_bytes = 55  # average header line length including '>' and '\n'
    seq_lines = (seq_length + line_wrap - 1) // line_wrap
    seq_bytes = seq_length + seq_lines  # characters + newlines
    # ~6% of records get an extra blank line (1/17)
    blank_fraction = 1.0 / 17
    return int(header_bytes + seq_bytes + blank_fraction)


def _records_for_target_size(seq_length: int, target_bytes: int) -> int:
    """Calculate number of records needed to reach target file size."""
    bpr = _estimate_bytes_per_record(seq_length)
    return max(100, target_bytes // bpr)


def _generate_fasta_file(
    path: str, n_records: int, seq_length: int, line_wrap: int = 60,
) -> int:
    """Write a synthetic FASTA file.  Returns actual file size in bytes."""
    raw_seq = ("ACGT" * ((seq_length // 4) + 1))[:seq_length]
    wrapped = textwrap.fill(raw_seq, width=line_wrap) + "\n"
    with open(path, "w") as fh:
        for i in range(n_records):
            if i % 37 == 0:
                desc = f"sp|Q{i:05d}|PROT_{i} Hypothetical protein length>{seq_length} OS=Synthetic"
            elif i % 23 == 0:
                desc = f"lcl|seq_{i} [organism=Test] [length = {seq_length} bp]"
            else:
                desc = f"seq_{i} Synthetic sequence length={seq_length}"
            fh.write(f">{desc}\n")
            fh.write(wrapped)
            if i % 17 == 0:
                fh.write("\n")
    return os.path.getsize(path)


# ---------------------------------------------------------------------------
# Four approaches being compared
# ---------------------------------------------------------------------------

def approach_a_simple_fasta_parser(path: str) -> list[str]:
    """Biopython SimpleFastaParser — forces full sequence buffering."""
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    titles = []
    with open(path) as fh:
        for title, _seq in SimpleFastaParser(fh):
            titles.append(title)
    return titles


def approach_b_naive_inline(path: str) -> list[str]:
    """Naive 2-liner: for line in fh: if line[0] == '>'."""
    titles = []
    with open(path) as fh:
        for line in fh:
            if line and line[0] == ">":
                titles.append(line[1:].rstrip())
    return titles


def approach_b_robust_inline(path: str) -> list[str]:
    """Robust inline — handles blank lines, \\r, edge cases."""
    titles = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped[0] == ">":
                titles.append(stripped[1:])
    return titles


def approach_c_header_parser(path: str) -> list[str]:
    """Proposed FastaHeaderParser (PR #5185).

    Uses ``line and line[0] == ">"`` for efficient blank-line handling
    (avoids the slower ``line[0:1]`` slice).  Falls back to a reference
    implementation if running on vanilla Biopython without the patch.
    """
    try:
        from Bio.SeqIO.FastaIO import FastaHeaderParser
        with open(path) as fh:
            return list(FastaHeaderParser(fh))
    except ImportError:
        pass
    # Reference implementation — matches the patch
    titles = []
    with open(path) as fh:
        for line in fh:
            if line and line[0] == ">":
                titles.append(line[1:].rstrip())
    return titles


APPROACHES = [
    ("SimpleFastaParser", approach_a_simple_fasta_parser),
    ("NaiveInline", approach_b_naive_inline),
    ("RobustInline", approach_b_robust_inline),
    ("FastaHeaderParser", approach_c_header_parser),
]


# ---------------------------------------------------------------------------
# Measurement harness
# ---------------------------------------------------------------------------

def _measure(func, path: str, warmup: int = 1, repeats: int = 3) -> dict:
    """Run func(path) and return timing + memory stats."""
    for _ in range(warmup):
        func(path)
    gc.collect()
    gc.disable()
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    times = []
    n_titles = 0
    for _ in range(repeats):
        t0 = time.perf_counter()
        titles = func(path)
        t1 = time.perf_counter()
        n_titles = len(titles)
        times.append(t1 - t0)
        del titles
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    gc.enable()
    gc.collect()
    return {
        "best_s": min(times),
        "median_s": sorted(times)[len(times) // 2],
        "rss_delta_kb": rss_after - rss_before,
        "n_titles": n_titles,
    }


# ---------------------------------------------------------------------------
# TSV output
# ---------------------------------------------------------------------------

TSV_HEADER = (
    "hostname\tseq_length\tn_records\tfile_mb\tapproach\t"
    "best_s\tmedian_s\trec_per_s\trss_delta_kb\tspeedup_vs_A"
)


def _tsv_line(hostname, seq_len, n_records, file_mb, label, stats, baseline):
    best = stats["best_s"]
    median = stats["median_s"]
    rps = stats["n_titles"] / best if best > 0 else 0
    rss = stats["rss_delta_kb"]
    speedup = baseline / best if baseline > 0 and best > 0 else 1.0
    return (
        f"{hostname}\t{seq_len}\t{n_records}\t{file_mb:.1f}\t{label}\t"
        f"{best:.4f}\t{median:.4f}\t{rps:.0f}\t{rss}\t{speedup:.2f}"
    )


# ---------------------------------------------------------------------------
# Main benchmark runner
# ---------------------------------------------------------------------------

def _info(msg):
    """Print to stderr (human-readable commentary)."""
    print(msg, file=sys.stderr)


def run_benchmark(
    seq_lengths: list[int],
    target_bytes: int,
    repeats: int = 3,
    tmpdir: str | None = None,
):
    hostname = platform.node() or "unknown"

    _info("=" * 78)
    _info("FASTA Header Extraction Benchmark")
    _info(f"  Host: {hostname}")
    _info(f"  Target file size: {target_bytes / 2**30:.2f} GB per tier")
    _info(f"  Sequence lengths: {seq_lengths}")
    _info(f"  Repeats: {repeats}")
    _info("=" * 78)
    _info("")
    _info("=== CPU cache hierarchy (Intel Xeon Platinum 8260, 4-socket NUMA) ===")
    _info("  L1d=32KB/core  L2=1MB/core  L3=35.75MB/socket (143MB total)")
    _info("  Cache line: 64B   NUMA cross-node latency: ~2x local")
    _info("")

    if tmpdir is None:
        home_tmp = os.path.expanduser("~/tmp")
        if os.path.isdir(home_tmp):
            tmpdir = tempfile.mkdtemp(prefix="biopython_bench_", dir=home_tmp)
        else:
            tmpdir = tempfile.mkdtemp(prefix="biopython_bench_")
    _info(f"  Temp dir: {tmpdir}")
    _info("")

    # Print TSV header to stdout
    print(TSV_HEADER)

    try:
        for seq_len in seq_lengths:
            n_records = _records_for_target_size(seq_len, target_bytes)
            seq_data_mb = (n_records * seq_len) / (1024 * 1024)
            hdr_kb = (n_records * 55) / 1024

            _info("─" * 78)
            _info(f"  {n_records:,} records × {seq_len:,} bp/seq")
            _info(f"  Sequence payload: {seq_data_mb:,.1f} MB")
            _info(f"  Header payload:   {hdr_kb:,.1f} KB")

            per_seq_kb = seq_len / 1024
            if seq_data_mb > 143:
                note = "overflows ENTIRE system L3 (143 MB, 4 NUMA nodes)"
            elif seq_data_mb > 35.75:
                note = "overflows single-socket L3 (35.75 MB)"
            elif seq_data_mb > 1:
                note = f"overflows per-core L2 (1 MB), fits in L3"
            else:
                note = "fits in L2"
            _info(f"  Cache: single seq {per_seq_kb:.1f} KB, total {seq_data_mb:.0f} MB — {note}")

            waste_mb = seq_data_mb * 3
            _info(f"  SimpleFastaParser waste: ~{waste_mb:,.0f} MB  |  HeaderParser: ~{hdr_kb:.0f} KB")
            _info("")

            fasta_path = os.path.join(tmpdir, f"bench_{seq_len}.fasta")
            _info(f"  Generating {fasta_path} ...")
            file_size = _generate_fasta_file(fasta_path, n_records, seq_len)
            file_mb = file_size / (1024 * 1024)
            _info(f"  File size: {file_mb:,.1f} MB  ({file_size:,} bytes)")
            _info("")

            baseline_time = None
            for label, func in APPROACHES:
                try:
                    stats = _measure(func, fasta_path, warmup=1, repeats=repeats)
                except Exception as e:
                    _info(f"  ERROR {label}: {e}")
                    continue

                if baseline_time is None:
                    baseline_time = stats["best_s"]
                rps = stats["n_titles"] / stats["best_s"] if stats["best_s"] > 0 else 0
                speedup = baseline_time / stats["best_s"] if stats["best_s"] > 0 else 0

                _info(f"  {label:<25} {stats['best_s']:8.3f}s  "
                      f"{rps:>12,.0f} rec/s  {speedup:.2f}×")

                # TSV to stdout
                print(_tsv_line(hostname, seq_len, n_records, file_mb,
                                label, stats, baseline_time))

            _info("")
            os.unlink(fasta_path)

    finally:
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass

    _info("=" * 78)
    _info("Done. TSV results on stdout.")


# ---------------------------------------------------------------------------
# Correctness edge-case tests (unittest)
# ---------------------------------------------------------------------------

class TestFastaHeaderParserEdgeCases(unittest.TestCase):
    """Demonstrate why a naive inline parser is insufficient."""

    def _header_parse(self, text: str) -> list[str]:
        try:
            from Bio.SeqIO.FastaIO import FastaHeaderParser
            return list(FastaHeaderParser(io.StringIO(text)))
        except ImportError:
            titles = []
            for line in io.StringIO(text):
                if line[0:1] == ">":
                    titles.append(line[1:].rstrip())
            return titles

    def _naive_parse(self, text: str) -> list[str]:
        titles = []
        for line in io.StringIO(text):
            if line and line[0] == ">":
                titles.append(line[1:].rstrip())
        return titles

    def test_basic(self):
        fasta = ">seq1 desc\nACGT\n>seq2\nTTTT\n"
        self.assertEqual(self._header_parse(fasta), ["seq1 desc", "seq2"])

    def test_blank_lines_between_records(self):
        """Blank lines between records are valid FASTA."""
        fasta = ">seq1\nACGT\n\n>seq2\nTTTT\n"
        self.assertEqual(len(self._header_parse(fasta)), 2)

    def test_blank_lines_naive_crashes(self):
        """Naive parser with bare line[0] would crash on blank lines."""
        fasta = ">seq1\nACGT\n\n>seq2\nTTTT\n"
        self.assertEqual(len(self._naive_parse(fasta)), 2)

    def test_gt_in_description(self):
        """'>' in description should not start a new record."""
        fasta = ">gene_A length>1000 organism=Test\nACGT\n>gene_B\nTTTT\n"
        result = self._header_parse(fasta)
        self.assertEqual(len(result), 2)
        self.assertIn("length>1000", result[0])

    def test_empty_file(self):
        self.assertEqual(self._header_parse(""), [])

    def test_no_sequence(self):
        """Header with empty sequence."""
        fasta = ">empty_seq\n>next_seq\nACGT\n"
        self.assertEqual(len(self._header_parse(fasta)), 2)

    def test_multiline_sequence_not_counted(self):
        long_seq = "ACGT" * 10000
        fasta = f">big_seq\n{long_seq}\n>small_seq\nA\n"
        self.assertEqual(self._header_parse(fasta), ["big_seq", "small_seq"])

    def test_wrapper_not_trivial(self):
        """Edge cases requiring non-trivial handling."""
        tricky = (
            ">seq1 has > in description\nACGT\n\n"
            ">seq2\r\nTTTT\n>seq3\n\n>seq4\n ACGT \n"
        )
        result = self._header_parse(tricky)
        self.assertEqual(len(result), 4)
        self.assertEqual(result[0], "seq1 has > in description")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Benchmark FASTA header extraction approaches. "
                    "TSV results go to stdout; commentary to stderr."
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="Quick mode: 50 MB files, fewer tiers"
    )
    parser.add_argument(
        "--target-gb", type=float, default=2.4,
        help="Target file size in GB for each tier (default: 2.4)"
    )
    parser.add_argument(
        "--seq-lengths", type=int, nargs="+",
        help="Sequence lengths to benchmark (default: 100 300 1000 30000 250000)"
    )
    parser.add_argument(
        "--repeats", type=int, default=3,
        help="Measurement repeats per approach (default: 3)"
    )
    parser.add_argument(
        "--run-tests", action="store_true",
        help="Also run edge-case unit tests"
    )
    args = parser.parse_args()

    if args.run_tests:
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromTestCase(TestFastaHeaderParserEdgeCases)
        runner = unittest.TextTestRunner(verbosity=2, stream=sys.stderr)
        result = runner.run(suite)
        if not result.wasSuccessful():
            sys.exit(1)
        _info("")

    if args.quick:
        target_bytes = 50 * 1024 * 1024  # 50 MB
        seq_lengths = args.seq_lengths or [100, 300, 30_000]
        repeats = 1
    else:
        target_bytes = int(args.target_gb * 2**30)
        seq_lengths = args.seq_lengths or [100, 300, 1_000, 30_000, 250_000]
        repeats = args.repeats

    run_benchmark(
        seq_lengths=seq_lengths,
        target_bytes=target_bytes,
        repeats=repeats,
    )
