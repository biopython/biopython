#!/usr/bin/env python3
"""Benchmark SeqIO FASTA index construction: original vs optimised.

This script measures the wall-clock time to build an in-memory index
(scanning all record boundaries and extracting IDs) for a synthetic
FASTA file of configurable size.

Two approaches are compared:

  Original  – compiled regex ``marker_re.match(line)`` on every line
               plus ``handle.tell()`` before every ``readline()``.
               This is the code path in upstream Biopython's
               ``SequentialSeqFileRandomAccess.__iter__``.

  Optimised – first-byte short-circuit before regex (avoids the regex
               engine on ~99 % of lines) plus arithmetic offset
               tracking for non-BGZF files (eliminates ``tell()``
               syscalls on sequence lines entirely).

Usage::

    python test_index_benchmark.py --target-mb 500 --seq-length 100
    python test_index_benchmark.py --target-mb 2000 --seq-length 1000

Output is tab-separated for easy pasting into spreadsheets or PRs.
"""

import argparse
import os
import platform
import re
import socket
import sys
import tempfile
import time
import unittest

# ── synthetic FASTA generation ────────────────────────────────────────────────


def _generate_fasta(path, n_records, seq_length):
    """Write a synthetic FASTA file with *n_records* records.

    Each record has a header ``>seq_XXXXXXXX synthetic record`` followed
    by the nucleotide sequence wrapped at 70 characters per line.
    """
    seq_line = "A" * 70 + "\n"
    full_lines, remainder = divmod(seq_length, 70)
    with open(path, "w") as fh:
        for i in range(n_records):
            fh.write(f">seq_{i:08d} synthetic record\n")
            for _ in range(full_lines):
                fh.write(seq_line)
            if remainder:
                fh.write("A" * remainder + "\n")


def _estimate_records(target_mb, seq_length):
    """Estimate how many records are needed for *target_mb* of FASTA."""
    header_len = len(">seq_00000000 synthetic record\n")
    full_lines, remainder = divmod(seq_length, 70)
    seq_bytes = full_lines * 71 + (remainder + 1 if remainder else 0)
    record_bytes = header_len + seq_bytes
    return int(target_mb * 1024 * 1024 / record_bytes)


# ── original (baseline) indexer ───────────────────────────────────────────────


def _index_original(path):
    """Faithful reproduction of the original __iter__ code path.

    Uses ``marker_re.match(line)`` on EVERY line and
    ``handle.tell()`` before every ``readline()`` — the exact
    upstream code path before this optimisation.
    """
    marker_re = re.compile(b"^>")
    marker_offset = 1
    count = 0
    with open(path, "rb") as handle:
        # Skip header
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        while marker_re.match(line):
            id_bytes = line[marker_offset:].strip().split(None, 1)[0]
            length = len(line)
            while True:
                end_offset = handle.tell()
                line = handle.readline()
                if marker_re.match(line) or not line:
                    count += 1
                    start_offset = end_offset
                    break
                else:
                    length += len(line)
    return count


# ── optimised indexer ─────────────────────────────────────────────────────────


def _index_optimised(path):
    """Use the optimised SequentialSeqFileRandomAccess.__iter__.

    Uses first-byte short-circuit before regex plus arithmetic
    offset tracking (no ``tell()`` on sequence lines for non-BGZF
    handles).  The handle is opened via ``_open_for_random_access``
    exactly as ``SeqIO.index()`` would.
    """
    from Bio.SeqIO._index import SequentialSeqFileRandomAccess

    idx = SequentialSeqFileRandomAccess(path, "fasta")
    count = 0
    for _ in idx:
        count += 1
    idx._handle.close()
    return count


# ── benchmark runner ──────────────────────────────────────────────────────────


def _run_benchmark(target_mb, seq_length, repeats=3):
    """Generate test file, run both indexers, report timing."""
    n_records = _estimate_records(target_mb, seq_length)

    tmpdir = tempfile.mkdtemp(prefix="bench_idx_")
    fasta_path = os.path.join(tmpdir, "bench.fasta")

    print(
        f"Generating {n_records:,} records × {seq_length} bp ...",
        file=sys.stderr,
        flush=True,
    )
    _generate_fasta(fasta_path, n_records, seq_length)
    file_mb = os.path.getsize(fasta_path) / (1024 * 1024)
    print(f"File: {file_mb:.1f} MB", file=sys.stderr, flush=True)

    hostname = socket.gethostname()

    # Warm up filesystem cache
    with open(fasta_path, "rb") as fh:
        while fh.read(1 << 20):
            pass

    results = {}
    for label, func in [
        ("Original", _index_original),
        ("Optimised", _index_optimised),
    ]:
        times = []
        n_indexed = 0
        for _ in range(repeats):
            t0 = time.perf_counter()
            n_indexed = func(fasta_path)
            t1 = time.perf_counter()
            times.append(t1 - t0)
        best = min(times)
        median = sorted(times)[len(times) // 2]
        results[label] = {
            "best": best,
            "median": median,
            "n_indexed": n_indexed,
            "rec_per_s": n_indexed / best if best > 0 else 0,
        }

    # Cleanup
    os.unlink(fasta_path)
    os.rmdir(tmpdir)

    speedup = (
        results["Original"]["best"] / results["Optimised"]["best"]
        if results["Optimised"]["best"] > 0
        else 0
    )

    # TSV output
    header = (
        "hostname\tseq_length\tn_records\tfile_mb\t"
        "approach\tbest_s\tmedian_s\trec_per_s\tspeedup_vs_original"
    )
    print(header)
    for label in ["Original", "Optimised"]:
        r = results[label]
        sp = 1.0 if label == "Original" else speedup
        print(
            f"{hostname}\t{seq_length}\t{r['n_indexed']}\t{file_mb:.1f}\t"
            f"{label}\t{r['best']:.4f}\t{r['median']:.4f}\t"
            f"{int(r['rec_per_s'])}\t{sp:.2f}"
        )

    return results


# ── unit tests ────────────────────────────────────────────────────────────────


class TestIndexOptimisation(unittest.TestCase):
    """Verify the optimised indexer produces identical results."""

    def _make_original_results(self, path):
        """Run the original indexer and return (id, offset, length) list."""
        marker_re = re.compile(b"^>")
        results = []
        with open(path, "rb") as handle:
            while True:
                start_offset = handle.tell()
                line = handle.readline()
                if marker_re.match(line) or not line:
                    break
            while marker_re.match(line):
                id_bytes = line[1:].strip().split(None, 1)[0]
                length = len(line)
                while True:
                    end_offset = handle.tell()
                    line = handle.readline()
                    if marker_re.match(line) or not line:
                        results.append((id_bytes.decode(), start_offset, length))
                        start_offset = end_offset
                        break
                    else:
                        length += len(line)
        return results

    def _check(self, n_records, seq_length):
        from Bio.SeqIO._index import SequentialSeqFileRandomAccess

        tmpdir = tempfile.mkdtemp(prefix="test_idx_")
        path = os.path.join(tmpdir, "test.fasta")
        _generate_fasta(path, n_records, seq_length)
        try:
            ref = self._make_original_results(path)
            idx = SequentialSeqFileRandomAccess(path, "fasta")
            opt = list(idx)
            idx._handle.close()
            self.assertEqual(
                len(ref),
                len(opt),
                f"n_records={n_records}, seq_length={seq_length}",
            )
            for i, (r, o) in enumerate(zip(ref, opt)):
                self.assertEqual(r, o, f"Record {i}")
        finally:
            os.unlink(path)
            os.rmdir(tmpdir)

    def test_single_record(self):
        """Single record, short sequence."""
        self._check(1, 10)

    def test_100_records_short(self):
        """100 records × 50 bp."""
        self._check(100, 50)

    def test_100_records_1kbp(self):
        """100 records × 1000 bp (multi-line sequence)."""
        self._check(100, 1000)

    def test_1000_records_300bp(self):
        """1000 records × 300 bp."""
        self._check(1000, 300)

    def test_empty_file(self):
        """Empty FASTA file."""
        from Bio.SeqIO._index import SequentialSeqFileRandomAccess

        tmpdir = tempfile.mkdtemp(prefix="test_idx_")
        path = os.path.join(tmpdir, "empty.fasta")
        with open(path, "w") as fh:
            pass
        try:
            ref = self._make_original_results(path)
            idx = SequentialSeqFileRandomAccess(path, "fasta")
            opt = list(idx)
            idx._handle.close()
            self.assertEqual(ref, [])
            self.assertEqual(opt, [])
        finally:
            os.unlink(path)
            os.rmdir(tmpdir)

    def test_blank_lines(self):
        """FASTA with blank lines between records."""
        from Bio.SeqIO._index import SequentialSeqFileRandomAccess

        tmpdir = tempfile.mkdtemp(prefix="test_idx_")
        path = os.path.join(tmpdir, "blanks.fasta")
        with open(path, "w") as fh:
            fh.write(">rec1\nACGT\n\n>rec2\nTGCA\n\n")
        try:
            ref = self._make_original_results(path)
            idx = SequentialSeqFileRandomAccess(path, "fasta")
            opt = list(idx)
            idx._handle.close()
            self.assertEqual(len(ref), len(opt))
            for i, (r, o) in enumerate(zip(ref, opt)):
                self.assertEqual(r, o, f"Record {i}")
        finally:
            os.unlink(path)
            os.rmdir(tmpdir)

    def test_get_raw_consistency(self):
        """get_raw returns correct bytes for each record."""
        from Bio.SeqIO._index import SequentialSeqFileRandomAccess

        tmpdir = tempfile.mkdtemp(prefix="test_idx_")
        path = os.path.join(tmpdir, "raw.fasta")
        _generate_fasta(path, 50, 200)
        try:
            idx = SequentialSeqFileRandomAccess(path, "fasta")
            records = list(idx)
            for rec_id, offset, length in records:
                raw = idx.get_raw(offset)
                self.assertTrue(
                    raw.startswith(b">"),
                    f"Record {rec_id} raw doesn't start with >",
                )
                self.assertEqual(len(raw), length, f"Record {rec_id} length mismatch")
            idx._handle.close()
        finally:
            os.unlink(path)
            os.rmdir(tmpdir)

    def test_tell_calls_reduced(self):
        """Verify tell() calls are reduced to ~1 for non-BGZF."""
        from Bio.SeqIO._index import SequentialSeqFileRandomAccess

        tmpdir = tempfile.mkdtemp(prefix="test_idx_")
        path = os.path.join(tmpdir, "tell_count.fasta")
        _generate_fasta(path, 1000, 100)
        try:
            idx = SequentialSeqFileRandomAccess(path, "fasta")

            class CountingHandle:
                """Wrapper to count tell() calls."""

                def __init__(self, h):
                    self._h = h
                    self.tell_count = 0

                def seek(self, *a):
                    return self._h.seek(*a)

                def readline(self):
                    return self._h.readline()

                def read(self, n=-1):
                    return self._h.read(n)

                def tell(self):
                    self.tell_count += 1
                    return self._h.tell()

            ch = CountingHandle(idx._handle)
            idx._handle = ch
            records = list(idx)
            self.assertEqual(len(records), 1000)
            # Original would call tell() ~3000 times (3 lines per record)
            # Optimised should call tell() only ~1 time (initial seek)
            self.assertLess(
                ch.tell_count,
                10,
                f"Expected very few tell() calls, got {ch.tell_count}",
            )
            idx._handle = ch._h
            idx._handle.close()
        finally:
            os.unlink(path)
            os.rmdir(tmpdir)


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--target-mb",
        type=float,
        default=500,
        help="Target FASTA file size in MB [default: 500]",
    )
    parser.add_argument(
        "--seq-length",
        type=int,
        default=100,
        help="Nucleotide sequence length per record [default: 100]",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Timing repetitions [default: 3]",
    )
    parser.add_argument(
        "--run-tests",
        action="store_true",
        help="Run unit tests instead of benchmark",
    )
    args = parser.parse_args()

    if args.run_tests:
        unittest.main(argv=[""], exit=True, verbosity=2)
    else:
        _run_benchmark(args.target_mb, args.seq_length, args.repeats)
