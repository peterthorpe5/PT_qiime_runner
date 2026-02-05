#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QIIME 2 amplicon runner (DADA2 or Deblur).

Overview
--------
End-to-end QIIME 2 workflow for paired-end amplicon data with a switchable
denoiser: DADA2 (paired) or Deblur (single). Keeps QIIME artefacts/provenance
intact, emits TSV summaries, and prepares outputs for downstream analysis
(including a phyloseq-ready folder named 'phyloseq_output/').

Design choices
--------------
- Import as paired-end using PairedEndFastqManifestPhred33V2.
- Strongly encourage primer removal using q2-cutadapt trim-paired.
- If DADA2 is chosen, do not pre-merge; use denoise-paired.
- If Deblur is chosen, either join pairs inside QIIME (vsearch) or use
  forward-only reads; then denoise at a fixed trim length.

Notes
-----
- All tables written by this script are tab-separated (TSV) when applicable.
- Named arguments are required; no positional arguments are used.
- UK English spelling is used throughout documentation strings.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional
from tempfile import NamedTemporaryFile
import stat
import os
import errno
import logging
import csv
import time
import pandas as pd
import gzip
import math
from typing import Iterable, Optional, Dict
try:
    import psutil
except Exception:
    psutil = None


# Wall-clock start for runtime/elapsed logging
_SCRIPT_START_TIME = time.time()


V34_FWD = "CCTACGGGNGGCWGCAG"
V34_REV = "GACTACHVGGGTATCTAATCC"
V34_FWD_OVH = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" + V34_FWD
V34_REV_OVH = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" + V34_REV

class Paths:
    """Container for key filesystem paths used in the run.

    Attributes
    ----------
    root : Path
        Root directory for all results.
    artifacts : Path
        Directory for QIIME 2 artefacts (.qza).
    visuals : Path
        Directory for QIIME 2 visualisations (.qzv).
    taxonomy : Path
        Directory for taxonomy artefacts.
    phylogeny : Path
        Directory for phylogeny artefacts.
    exports : Path
        Directory for optional exported files.
    logs : Path
        Directory for log files.
    phyloseq_output : Path
        Directory for phyloseq outputs (artefacts + metadata TSV).
    report_tsv : Path
        Tab-separated summary of key QC checkpoints.
    qza, qzv : Path
        Back-compat aliases to artifacts/visuals.
    """

    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.artifacts = self.root / "artifacts_qza"
        self.visuals = self.root / "visuals_qzv"
        self.taxonomy = self.root / "taxonomy"
        self.phylogeny = self.root / "phylogeny"
        self.exports = self.root / "exports"
        self.logs = self.root / "logs"
        self.phyloseq_output = self.root / "phyloseq_output"
        self.report_tsv = self.root / "run_report.tsv"
        # Aliases for older code paths.qza / paths.qzv
        self.qza = self.artifacts
        self.qzv = self.visuals

    def mkdirs(self) -> None:
        """Create all output directories if they do not already exist."""
        for p in (
            self.artifacts,
            self.visuals,
            self.taxonomy,
            self.phylogeny,
            self.exports,
            self.logs,
            self.phyloseq_output,
        ):
            p.mkdir(parents=True, exist_ok=True)


# ----------------------------- logging -------------------------- #

def setup_logging(*, out_dir: Path, run_label: str) -> logging.Logger:
    """
    Configure structured logging to both stderr (human) and a file (machine).

    The file log captures DEBUG+ with timestamps; the stderr stream shows
    INFO+ with compact formatting. Intended for batch runs on HPC/GPFS.

    Parameters
    ----------
    out_dir : pathlib.Path
        The run's output directory (e.g., results/<RUN>).
    run_label : str
        Identifier for the run; used in initial metadata lines.

    Returns
    -------
    logging.Logger
        Configured logger instance ('q2_amplicon_runner').
    """
    log_dir = Path(out_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "run_debug.log"

    logger = logging.getLogger("q2_amplicon_runner")
    logger.setLevel(logging.DEBUG)
    # Avoid duplicate handlers if reinitialised
    logger.handlers.clear()
    logger.propagate = False

    # ---- stderr (human) ----
    stream_handler = logging.StreamHandler(stream=sys.stderr)
    stream_handler.setLevel(logging.INFO)
    stream_fmt = logging.Formatter("%(levelname)s: %(message)s")
    stream_handler.setFormatter(stream_fmt)

    # ---- file (machine) ----
    file_handler = logging.FileHandler(filename=log_file, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
    )
    file_handler.setFormatter(file_fmt)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    # Context banner
    logger.info("Run label: %s", run_label)
    logger.info("Output directory: %s", Path(out_dir).resolve())
    logger.debug("Python version: %s", " ".join(map(str, sys.version_info)))
    logger.debug("Command line: %s", " ".join(sys.argv))

    return logger


def validate_metadata_content(
    *,
    metadata_tsv: Path,
    min_cols: int = 2,
) -> None:
    """
    Validate that a QIIME-style metadata TSV looks sane.

    Checks:
      - No obvious non-TSV payloads (e.g., lines starting with a Python shebang).
      - First row exists and has at least `min_cols` columns.
      - All non-comment rows have the same number of columns as the header.

    Parameters
    ----------
    metadata_tsv : pathlib.Path
        Candidate metadata TSV.
    min_cols : int, default 2
        Minimum number of columns required in the header.

    Raises
    ------
    ValueError
        If any check fails.
    """
    if not metadata_tsv.exists():
        raise ValueError(f"Metadata file not found: {metadata_tsv}")

    with metadata_tsv.open("r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError(f"Metadata is empty: {metadata_tsv}")

    header = lines[0].rstrip("\n\r")
    n_cols = len(header.split("\t")) if header else 0
    if n_cols < min_cols:
        raise ValueError(
            f"Metadata header has too few columns ({n_cols} < {min_cols}) in {metadata_tsv}"
        )

    # Obvious non-TSV payloads
    bad_markers = ("#!/usr/bin/env python", "from __future__ import", "import argparse")
    for i, raw in enumerate(lines[:200], start=1):  # scan first 200 lines
        if any(raw.startswith(m) for m in bad_markers):
            raise ValueError(
                f"Detected non-TSV content at line {i} of {metadata_tsv}: {raw[:60]!r}"
            )

    # Column count consistency on data lines (skip comments)
    for i, raw in enumerate(lines[1:], start=2):
        if not raw.strip() or raw.startswith("#"):
            continue
        nf = len(raw.rstrip("\n\r").split("\t"))
        if nf != n_cols:
            raise ValueError(
                f"Column count mismatch at line {i} ({nf} vs header {n_cols}) in {metadata_tsv}"
            )


def _rmtree_onerror(func, path, exc_info):
    """
    Handle errors from shutil.rmtree by making paths writable and retrying.

    This is useful on shared filesystems where permissions and transient locks
    can cause intermittent deletion failures.
    """
    try:
        os.chmod(path, 0o700)
    except OSError:
        pass
    try:
        func(path)
    except OSError:
        # Let the caller decide whether to retry.
        raise


def safe_rmtree(path: Path, retries: int = 6, sleep_s: float = 0.5) -> bool:
    """
    Remove a directory tree robustly on shared filesystems.

    Returns True if the directory is removed (or did not exist), False otherwise.
    """
    if not path.exists():
        return True

    for attempt in range(1, retries + 1):
        try:
            shutil.rmtree(path, onerror=_rmtree_onerror)
            return True
        except OSError as exc:
            # Errno 39: directory not empty (common on NFS/parallel writes)
            if exc.errno in {errno.ENOTEMPTY, errno.EBUSY, errno.EACCES}:
                time.sleep(sleep_s * attempt)
                continue
            time.sleep(sleep_s * attempt)

    return False


def safe_rotate_then_rmtree(path: Path) -> Optional[Path]:
    """
    Rename a directory to a unique name, then attempt deletion.

    Renaming first reduces races where other processes are still writing into
    the original directory name.
    """
    if not path.exists():
        return None

    rotated = path.with_name(f"{path.name}.old.{int(time.time())}")
    try:
        path.rename(rotated)
    except OSError:
        # If rename fails, fall back to deleting in place
        rotated = path

    _ = safe_rmtree(rotated)
    return rotated



def make_qiime_metadata_copy(
    *,
    src_tsv: Path,
    out_dir: Path,
    id_col_candidates: Iterable[str] = ("sample-id", "SampleID", "sampleID", "Sample_ID", "sample_id", "Sample-ID", "Sample ID", "Sample.ID"),
    force_id_header: str = "sample-id",
    add_q2_types: bool = False,
) -> Path:
    """
    Write a persistent, normalised metadata TSV for QIIME.

    Normalises the first column name to `force_id_header` and keeps the file
    on disk for the whole run (no ephemeral NamedTemporaryFile).

    Parameters
    ----------
    src_tsv : pathlib.Path
        Source metadata TSV.
    out_dir : pathlib.Path
        Folder to write the normalised copy into.
    id_col_candidates : Iterable[str], default (...)
        Acceptable first-column names to treat as the ID column.
    force_id_header : str, default 'sample-id'
        Header to use for the first column in the output.
    add_q2_types : bool, default False
        If True, add a '#q2:types' second line with 'categorical' for all columns.

    Returns
    -------
    pathlib.Path
        Path to the persistent, QIIME-friendly metadata TSV.

    Raises
    ------
    ValueError
        If the source file fails validation or no ID column is found.
    """
    validate_metadata_content(metadata_tsv=src_tsv)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "metadata.for_qiime.tsv"

    df = pd.read_csv(src_tsv, sep="\t", dtype=str).fillna("")

    # Identify the ID column
    first_col = df.columns[0]
    if first_col not in id_col_candidates:
        raise ValueError(
            f"First column must be one of {tuple(id_col_candidates)}; found '{first_col}' in {src_tsv}"
        )

    # Rename and clean ID column
    if first_col != force_id_header:
        df = df.rename(columns={first_col: force_id_header})
    df[force_id_header] = df[force_id_header].astype(str).str.strip()

    # Ensure ID column is first
    cols = [force_id_header] + [c for c in df.columns if c != force_id_header]
    df = df[cols]

    if add_q2_types:
        header = "\t".join(df.columns)
        types_line = "#q2:types\t" + "\t".join("categorical" for _ in df.columns[1:])
        with out_path.open("w", encoding="utf-8") as handle:
            handle.write(header + "\n")
            handle.write(types_line + "\n")
        df.to_csv(out_path, sep="\t", index=False, mode="a", encoding="utf-8")
    else:
        df.to_csv(out_path, sep="\t", index=False, encoding="utf-8")

    # Re-validate the written copy (defensive)
    validate_metadata_content(metadata_tsv=out_path)
    return out_path




def log_section(*, logger: logging.Logger, title: str) -> None:
    """
    Emit a visible section divider in logs.

    Parameters
    ----------
    logger : logging.Logger
        Logger to write to.
    title : str
        Short section title.
    """
    sep = "=" * max(10, min(80, len(title) + 8))
    logger.info("%s", sep)
    logger.info("== %s ==", title)
    logger.info("%s", sep)


def log_memory_usage(
    logger: logging.Logger,
    prefix: str = "",
    extra_msg: str | None = None,
) -> None:
    """
    Log the current and peak memory usage (resident set size), plus elapsed time.

    Parameters
    ----------
    logger : logging.Logger
        Logger instance to emit the message.
    prefix : str
        Optional prefix (e.g., 'START', 'END', or a pipeline stage label).
    extra_msg : str | None
        Optional extra text appended to the log message.
    """
    cur_gb = None
    peak_gb = None

    # Try psutil for current RSS; fall back to /proc/self/status on Linux
    if psutil is not None:
        try:
            cur_gb = psutil.Process(os.getpid()).memory_info().rss / (1024 ** 3)
        except Exception:
            cur_gb = None

    if cur_gb is None:
        # very light fallback: try /proc/self/status VmRSS
        try:
            with open("/proc/self/status", "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if line.startswith("VmRSS:"):
                        kb = float(line.split()[1])
                        cur_gb = kb / (1024 ** 2)
                        break
        except Exception:
            cur_gb = None

    # Peak RSS via resource (kilobytes on Linux)
    try:
        import resource  # local import to avoid platform issues
        peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # Linux reports KB, macOS reports bytes
        if os.uname().sysname == "Linux":
            peak_gb = peak_kb / (1024 ** 2)
        else:
            peak_gb = peak_kb / (1024 ** 3)
    except Exception:
        peak_gb = None

    # Elapsed wall time
    elapsed_s = max(0.0, time.time() - _SCRIPT_START_TIME)
    elapsed_min = elapsed_s / 60.0

    parts = []
    if prefix:
        parts.append(prefix.strip())
    if cur_gb is not None:
        parts.append(f"RAM: {cur_gb:.2f} GB")
    if peak_gb is not None:
        parts.append(f"Peak: {peak_gb:.2f} GB")
    parts.append(f"Elapsed: {elapsed_min:.1f} min")
    if extra_msg:
        parts.append(extra_msg)

    logger.info(" | ".join(parts))



def _open_maybe_gzip(*, path: Path):
    """Open a text FASTQ file that may be gzipped."""
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open(mode="r", encoding="utf-8", errors="replace")


def _sample_fastq_read_lengths(
    *,
    fastq_paths: Iterable[Path],
    max_reads: int = 20000,
) -> list[int]:
    """Sample read lengths from one or more FASTQ files.

    Reads up to `max_reads` reads total across the provided FASTQs.

    Parameters
    ----------
    fastq_paths : Iterable[pathlib.Path]
        FASTQ or FASTQ.GZ files.
    max_reads : int, default 20000
        Maximum number of reads to sample across all files.

    Returns
    -------
    list[int]
        Observed read lengths.
    """
    lengths: list[int] = []
    reads_seen = 0

    for fp in fastq_paths:
        if reads_seen >= max_reads:
            break
        if not fp.exists():
            continue

        with _open_maybe_gzip(path=fp) as handle:
            while reads_seen < max_reads:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline()
                plus = handle.readline()
                qual = handle.readline()
                if not qual:
                    break
                lengths.append(len(seq.strip()))
                reads_seen += 1

    return lengths


def _percentile_int(*, values: list[int], q: float) -> int:
    """Return an integer percentile using a simple sorted-index method."""
    if not values:
        raise ValueError("No values provided for percentile calculation.")
    if q <= 0:
        return min(values)
    if q >= 100:
        return max(values)

    vals = sorted(values)
    # 0-indexed rank
    rank = int(math.floor((q / 100.0) * (len(vals) - 1)))
    return int(vals[rank])


def infer_deblur_trim_length(
    *,
    fastq_paths: Iterable[Path],
    percentile: float = 10.0,
    min_len: int = 150,
    max_len: Optional[int] = None,
    warn_below: int = 200,
    logger: Optional[logging.Logger] = None,
) -> int:
    """Infer a Deblur trim length from observed read lengths.

    Uses a low percentile so most reads are long enough after trimming.

    Returns
    -------
    int
        Recommended trim length.
    """
    lengths = _sample_fastq_read_lengths(fastq_paths=fastq_paths, max_reads=20000)
    if not lengths:
        if logger:
            logger.warning(
                "Deblur trim-length auto-inference failed (no reads sampled); "
                "falling back to 200. Consider setting --deblur_trim_length explicitly."
            )
        return 200

    p = _percentile_int(values=lengths, q=percentile)
    rec = max(min_len, int(p))
    if max_len is not None:
        rec = min(rec, int(max_len))

    if logger:
        logger.info(
            "Deblur trim-length auto-inference: sampled=%d reads; "
            "min=%d, p%.0f=%d, median=%d, max=%d; chosen=%d",
            len(lengths),
            min(lengths),
            percentile,
            p,
            _percentile_int(values=lengths, q=50.0),
            max(lengths),
            rec,
        )
        if rec < warn_below:
            logger.warning(
                "Auto-inferred Deblur trim length (%d) is low. This will force all ASVs to %d bp "
                "and can reduce taxonomic resolution. Check your joining/trimming settings or set "
                "--deblur_trim_length explicitly.",
                rec,
                rec,
            )

    return rec



# ----------------------------- helpers ----------------------------- #
def run_cmd(*, cmd: list[str], log_file: Path, logger: Optional[logging.Logger] = None) -> None:
    """
    Run a shell command with all stdout/stderr tee'd to a step log.

    Parameters
    ----------
    cmd : list of str
        Command tokens (no shell=True).
    log_file : pathlib.Path
        Path to the step-specific log file.
    logger : Optional[logging.Logger]
        If provided, logs the command start and destination log file.

    Raises
    ------
    subprocess.CalledProcessError
        If the command exits non-zero.
    """
    from datetime import datetime

    log_file.parent.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if logger is not None:
        logger.info("Running command (%s): %s", timestamp, " ".join(cmd))
        logger.debug("Step log: %s", log_file)

    with log_file.open("a", encoding="utf-8") as lf:
        lf.write(f"\n===== BEGIN COMMAND @ {timestamp} =====\n")
        lf.write("$ " + " ".join(cmd) + "\n\n")
        lf.flush()

        subprocess.run(cmd, stdout=lf, stderr=lf, check=True)

        timestamp_end = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        lf.write(f"\n===== END COMMAND @ {timestamp_end} =====\n")



def _min_sample_depth_from_biom_summary(summary_txt: Path) -> Optional[int]:
    """Parse 'Minimum per-sample count' from BIOM summarize-table output."""
    if not summary_txt.exists():
        return None
    text = summary_txt.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"Minimum per-sample count:\s+(\d+)", text)
    return int(m.group(1)) if m else None


def run_qc_bundle(
    *, paths: "Paths", metadata_tsv: Path, denoiser: str,
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Build a standard QC bundle of QIIME 2 visualisations and text summaries.

    The function consumes core artefacts (feature table, representative
    sequences, rooted tree, optional taxonomy) and writes .qzv outputs and a
    few lightweight TSV/TXT summaries under ``<paths.root>/qc``. It is safe to
    call regardless of which artefacts exist; steps are skipped if inputs are
    missing.

    New behaviour:
        - Derives per-sample depths from the exported feature table and saves
          ``sample-frequency.tsv`` (SampleID → Frequency).
        - Chooses an automatic rarefaction max depth as the 15th percentile of
          per-sample depths (retain ≥ 85% of samples), capped at 10,000 and
          with a minimum of 10. Writes a short explanation to
          ``depth_recommendation.txt`` and uses this value for
          ``alpha-rarefaction``.
        - Adds ``asv_length_summary.tsv`` (min/median/mean/max/N) alongside
          ``asv_lengths.tsv``.
        - Expands ``counts_summary.tsv`` to include min/median/max depth,
          the chosen auto max depth, and estimated retained percentage.

    Created outputs (when inputs are available)
    -------------------------------------------
      - feature-table_summarize.qzv
      - tabulate-seqs.qzv
      - dada2_stats.qzv (if denoiser == "dada2_paired")
      - deblur_stats.qzv (if denoiser == "deblur_single")
      - taxa-barplot.qzv (if taxonomy present)
      - alpha-rarefaction.qzv (if rooted tree present; auto max depth)
      - sample_depths.txt      (BIOM per-sample depth summary)
      - sample-frequency.tsv   (SampleID → Frequency, derived)
      - counts_summary.tsv     (overall counts + depth stats)
      - asv_lengths.tsv        (FeatureID → ASV length, from rep-seqs FASTA)
      - asv_length_summary.tsv (summary of ASV lengths)
      - export_table/, export_repseqs/ (raw exports used by summaries)

    Parameters
    ----------
    paths : Paths
        A container with standardised output/input directories:
        ``paths.root``, ``paths.qza``, ``paths.qzv``, ``paths.logs``,
        ``paths.taxonomy``, ``paths.phylogeny``.
    metadata_tsv : pathlib.Path
        Path to the sample metadata table (TSV).
    denoiser : str
        Either ``"dada2_paired"`` or ``"deblur_single"``; used to choose
        which stats visualisation to render if present.
    logger : logging.Logger, optional
        Logger for progress and warnings.

    Returns
    -------
    None
    """
    qza_dir = paths.qza
    logs = paths.logs

    out_qc = paths.root / "qc"
    out_qc.mkdir(parents=True, exist_ok=True)

    # Prefer filtered artefacts when present (avoids ID mismatches in QC).
    table_filtered = qza_dir / "table.lengthfilter.qza"
    rep_filtered = qza_dir / "rep-seqs.lengthfilter.qza"

    table_qza = table_filtered if table_filtered.exists() else qza_dir / "table.qza"
    repseqs_qza = rep_filtered if rep_filtered.exists() else qza_dir / "rep-seqs.qza"

    rooted_tree_qza = paths.phylogeny / "rooted-tree.qza"
    stats_qza = qza_dir / "denoise-stats.qza"        # DADA2
    deblur_stats_qza = qza_dir / "deblur_stats.qza"  # Deblur (if present)
    taxonomy_qza = paths.taxonomy / "taxonomy.qza"

    # 1) Feature-table summarise (+ metadata so sample groups show)
    if table_qza.exists():
        run_cmd(
            cmd=[
                "qiime", "feature-table", "summarize",
                "--i-table", str(table_qza),
                "--m-sample-metadata-file", str(metadata_tsv),
                "--o-visualization", str(out_qc / "feature-table_summarize.qzv"),
            ],
            log_file=logs / "qc_feature_table_summarize.log",
            logger=logger,
        )

    # 2) Tabulate representative sequences
    if repseqs_qza.exists():
        run_cmd(
            cmd=[
                "qiime", "feature-table", "tabulate-seqs",
                "--i-data", str(repseqs_qza),
                "--o-visualization", str(out_qc / "tabulate-seqs.qzv"),
            ],
            log_file=logs / "qc_tabulate_seqs.log",
            logger=logger,
        )

    # 3) Denoiser stats visual
    if denoiser == "dada2_paired" and stats_qza.exists():
        run_cmd(
            cmd=[
                "qiime", "metadata", "tabulate",
                "--m-input-file", str(stats_qza),
                "--o-visualization", str(out_qc / "dada2_stats.qzv"),
            ],
            log_file=logs / "qc_dada2_stats.log",
            logger=logger,
        )
    elif denoiser == "deblur_single" and deblur_stats_qza.exists():
        run_cmd(
            cmd=[
                "qiime", "deblur", "visualize-stats",
                "--i-deblur-stats", str(deblur_stats_qza),
                "--o-visualization", str(out_qc / "deblur_stats.qzv"),
            ],
            log_file=logs / "qc_deblur_stats.log",
            logger=logger,
        )

    # 4) Taxa barplots (only if taxonomy + table exist)
    if taxonomy_qza.exists() and table_qza.exists():
        run_cmd(
            cmd=[
                "qiime", "taxa", "barplot",
                "--i-table", str(table_qza),
                "--i-taxonomy", str(taxonomy_qza),
                "--m-metadata-file", str(metadata_tsv),
                "--o-visualization", str(out_qc / "taxa-barplot.qzv"),
            ],
            log_file=logs / "qc_taxa_barplot.log",
            logger=logger,
        )

    # 5) Plain-text quick summaries + per-sample frequencies
    export_dir = out_qc / "export_table"
    export_dir.mkdir(parents=True, exist_ok=True)

    biom_fp = export_dir / "feature-table.biom"
    table_tsv = export_dir / "feature-table.tsv"
    depth_txt = out_qc / "sample_depths.txt"
    freq_tsv = out_qc / "sample-frequency.tsv"
    depth_note = out_qc / "depth_recommendation.txt"

    if table_qza.exists():
        run_cmd(
            cmd=[
                "qiime", "tools", "export",
                "--input-path", str(table_qza),
                "--output-path", str(export_dir),
            ],
            log_file=logs / "qc_export_table.log",
            logger=logger,
        )

    # BIOM summary (overall counts) and TSV conversion for deterministic parsing.
    if biom_fp.exists():
        run_cmd(
            cmd=[
                "biom", "summarize-table",
                "-i", str(biom_fp),
                "-o", str(depth_txt),
            ],
            log_file=logs / "qc_biom_summarize.log",
            logger=logger,
        )
        run_cmd(
            cmd=[
                "biom", "convert",
                "--input-fp", str(biom_fp),
                "--output-fp", str(table_tsv),
                "--to-tsv",
            ],
            log_file=logs / "qc_biom_convert_tsv.log",
            logger=logger,
        )

    # Parse TSV matrix → per-sample depths and write sample-frequency.tsv
    sample_freqs: Dict[str, int] = {}
    try:
        if table_tsv.exists():
            with table_tsv.open("r", encoding="utf-8", errors="replace") as fh:
                header: Optional[list[str]] = None
                for line in fh:
                    if line.startswith("#") and "OTU ID" in line:
                        # Header line: "#OTU ID\tS1\tS2\t..."
                        header = line.strip().lstrip("#").split("\t")
                        sample_freqs = {sid: 0 for sid in header[1:]}
                        continue
                    if line.startswith("#") or not line.strip():
                        continue
                    if header is None:
                        continue
                    fields = line.rstrip("\n\r").split("\t")
                    values = fields[1:]
                    for sid, val in zip(header[1:], values):
                        try:
                            sample_freqs[sid] += int(float(val))
                        except Exception:  # noqa: BLE001
                            # Treat non-numeric as zero
                            pass
            with freq_tsv.open("w", encoding="utf-8") as out:
                out.write("sample-id\tfrequency\n")
                for sid, freq in sample_freqs.items():
                    out.write(f"{sid}\t{freq}\n")
    except Exception as err:  # noqa: BLE001
        if logger:
            logger.warning("Failed to compute per-sample frequencies: %s", err)
        sample_freqs = {}

    # 6) Export rep-seqs FASTA → ASV lengths and summary
    rep_export = out_qc / "export_repseqs"
    rep_export.mkdir(parents=True, exist_ok=True)

    if repseqs_qza.exists():
        run_cmd(
            cmd=[
                "qiime", "tools", "export",
                "--input-path", str(repseqs_qza),
                "--output-path", str(rep_export),
            ],
            log_file=logs / "qc_export_repseqs.log",
            logger=logger,
        )

    fasta = rep_export / "dna-sequences.fasta"
    if fasta.exists():
        tsv = out_qc / "asv_lengths.tsv"
        tsv_summary = out_qc / "asv_length_summary.tsv"
        try:
            lengths: Dict[str, int] = {}
            seq_id: Optional[str] = None
            seq_buf: list[str] = []

            with fasta.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    line = line.rstrip("\n\r")
                    if line.startswith(">"):
                        if seq_id is not None:
                            lengths[seq_id] = sum(len(x) for x in seq_buf)
                        seq_id = line[1:].split()[0]
                        seq_buf = []
                    else:
                        seq_buf.append(line)
                if seq_id is not None:
                    lengths[seq_id] = sum(len(x) for x in seq_buf)

            # Per-ASV lengths
            with tsv.open("w", encoding="utf-8") as out:
                out.write("FeatureID\tLength\n")
                for k, v in lengths.items():
                    out.write(f"{k}\t{v}\n")

            # Tiny summary (min/median/mean/max/N)
            if lengths:
                vals = sorted(lengths.values())
                n = len(vals)
                if n % 2 == 1:
                    median = vals[n // 2]
                else:
                    median = int((vals[n // 2 - 1] + vals[n // 2]) / 2)
                mean = int(round(sum(vals) / n))
                with tsv_summary.open("w", encoding="utf-8") as out:
                    out.write("n\tmin\tmedian\tmean\tmax\n")
                    out.write(f"{n}\t{min(vals)}\t{median}\t{mean}\t{max(vals)}\n")
        except Exception as err:  # noqa: BLE001
            if logger:
                logger.warning("ASV length TSVs failed: %s", err)

    # 7) Tiny overall counts file using computed per-sample depths (preferred)
    counts_tsv = out_qc / "counts_summary.tsv"
    auto_max_depth = 10_000
    try:
        # Auto-select rarefaction depth from per-sample frequencies.
        if sample_freqs:
            depths = sorted(sample_freqs.values())
            n = len(depths)
            idx = max(0, min(n - 1, int(n * 0.15)))  # 15th percentile
            candidate = depths[idx]
            candidate = max(10, min(10_000, int(candidate)))
            auto_max_depth = candidate

        # Write a note describing the decision.
        with (out_qc / "depth_recommendation.txt").open("w", encoding="utf-8") as out:
            if sample_freqs:
                out.write(
                    "Auto rarefaction depth selection\n"
                    f"- Samples: {len(sample_freqs)}\n"
                    "- Strategy: 15th percentile of per-sample depths (retain ≥ 85%)\n"
                    f"- Chosen max depth: {auto_max_depth}\n"
                    "- Caps: min=10, max=10,000\n"
                )
            else:
                out.write(
                    "Auto rarefaction depth selection\n"
                    "- Could not derive per-sample depths; using fallback 10,000.\n"
                )
    except Exception as err:  # noqa: BLE001
        if logger:
            logger.warning("Auto depth selection failed, using 10,000: %s", err)
        auto_max_depth = 10_000

    try:
        # Fill counts summary using derived depths; fallback to BIOM summary.
        n_samples = -1
        n_taxa = -1
        total_reads = -1
        min_depth = -1
        median_depth = -1
        max_depth = -1
        retained_pc = -1

        if sample_freqs:
            depths = sorted(sample_freqs.values())
            n_samples = len(depths)
            total_reads = sum(depths)
            min_depth = depths[0]
            max_depth = depths[-1]
            if n_samples % 2 == 1:
                median_depth = depths[n_samples // 2]
            else:
                median_depth = int((depths[n_samples // 2 - 1] + depths[n_samples // 2]) / 2)

            # Try to extract taxa count from BIOM summary (if available).
            if depth_txt.exists():
                text = depth_txt.read_text(encoding="utf-8", errors="ignore")
                m = re.search(r"Number of observations:\s+(\d+)", text)
                n_taxa = int(m.group(1)) if m else -1

            # Estimated retention at chosen depth.
            try:
                retained = sum(1 for d in depths if d >= auto_max_depth)
                retained_pc = int(round(100 * retained / n_samples)) if n_samples > 0 else -1
            except Exception:  # noqa: BLE001
                retained_pc = -1

        else:
            # Fallback entirely to BIOM summary if present.
            if depth_txt.exists():
                text = depth_txt.read_text(encoding="utf-8", errors="ignore")
                m = re.search(r"Number of samples:\s+(\d+)", text)
                n_samples = int(m.group(1)) if m else -1
                m = re.search(r"Number of observations:\s+(\d+)", text)
                n_taxa = int(m.group(1)) if m else -1
                m = re.search(r"Total count:\s+(\d+)", text)
                total_reads = int(m.group(1)) if m else -1

        with counts_tsv.open("w", encoding="utf-8") as out:
            out.write(
                "samples\ttaxa\ttotal_reads\tmin_depth\tmedian_depth\tmax_depth\t"
                "auto_max_depth\test_retained_percent\n"
            )
            out.write(
                f"{n_samples}\t{n_taxa}\t{total_reads}\t{min_depth}\t{median_depth}\t"
                f"{max_depth}\t{auto_max_depth}\t{retained_pc}\n"
            )
    except Exception as err:  # noqa: BLE001
        if logger:
            logger.warning("Counts summary failed: %s", err)

    # 8) Alpha-rarefaction (requires tree + table) with auto max depth
    # --- Compute min_depth and a safe steps value ---
    # Derive min_depth from per-sample frequencies if available; else default to 1
    if sample_freqs:
        depths_sorted = sorted(sample_freqs.values())
        min_depth_observed = max(1, depths_sorted[0])
    else:
        min_depth_observed = 1

    # If auto_max_depth would be below min_depth, bump it up to min_depth
    if auto_max_depth < min_depth_observed:
        auto_max_depth = min_depth_observed
        if logger:
            logger.warning("Adjusted max depth to %d to be ≥ min depth.", auto_max_depth)

    # Steps cannot exceed the number of integer depths in [min_depth, max_depth]
    steps_possible = max(1, auto_max_depth - min_depth_observed + 1)
    steps_safe = min(10, steps_possible)  # 10 is the q2 default; clamp if needed

    # Optionally skip if the window is degenerate
    if steps_safe < 1 or auto_max_depth < 1:
        if logger:
            logger.warning("Skipping alpha-rarefaction (non-positive depth window).")
    else:
        alpha_qzv = out_qc / "alpha-rarefaction.qzv"
        alpha_qzv.unlink(missing_ok=True)  # avoid overwrite errors

        def _alpha_rarefy(with_tree: bool) -> None:
            cmd = [
                "qiime", "diversity", "alpha-rarefaction",
                "--i-table", str(table_qza),
                "--m-metadata-file", str(metadata_tsv),
                "--p-min-depth", str(min_depth_observed),
                "--p-max-depth", str(auto_max_depth),
                "--p-steps", str(steps_safe),
                "--o-visualization", str(alpha_qzv),
            ]
            log_path = logs / ("qc_alpha_rarefaction.log" if with_tree else "qc_alpha_rarefaction_no_tree.log")
            if with_tree:
                cmd[4:4] = ["--i-phylogeny", str(rooted_tree_qza)]
            run_cmd(cmd=cmd, log_file=log_path, logger=logger)

        # Prefer phylogenetic run, but fall back gracefully if it fails
        try:
            if rooted_tree_qza.exists():
                _alpha_rarefy(with_tree=True)
            else:
                _alpha_rarefy(with_tree=False)
        except subprocess.CalledProcessError:
            if logger:
                logger.warning(
                    "Phylogenetic alpha-rarefaction failed; retrying without a tree. "
                    "See logs/qc_alpha_rarefaction.log."
                )
            alpha_qzv.unlink(missing_ok=True)
            _alpha_rarefy(with_tree=False)




def _deduplicate_tsv_headers(in_path: Path, out_path: Path) -> Path:
    """Rewrite a TSV so its header names are unique (A, A.2, A.3, ...)."""
    with in_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty TSV: {in_path}")
    header, body = rows[0], rows[1:]
    seen = {}
    new_header = []
    for h in header:
        count = seen.get(h, 0) + 1
        seen[h] = count
        new_header.append(h if count == 1 else f"{h}.{count}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(new_header)
        w.writerows(body)
    return out_path

def ensure_legacy_metadata_location(
    *,
    metadata_src: Path,
    legacy_relpath: Path = Path("metadata/metadata16S.tsv"),
    prefer_symlink: bool = True,
) -> Path:
    """
    Ensure a metadata file exists at a legacy, hard-coded location.

    This function creates 'metadata/metadata16S.tsv' so any downstream
    hard-coded consumers (R scripts, older glue code) keep working.

    Parameters
    ----------
    metadata_src : pathlib.Path
        The canonical metadata TSV provided by the user (absolute or relative).
    legacy_relpath : pathlib.Path, optional
        Relative path to the legacy location. Default is 'metadata/metadata16S.tsv'.
    prefer_symlink : bool, optional
        If True, attempt to create a symlink. If the filesystem forbids
        symlinks (or creation fails), a regular file copy is made instead.

    Returns
    -------
    pathlib.Path
        The absolute path of the legacy file (symlink or copy).

    Raises
    ------
    FileNotFoundError
        If `metadata_src` does not exist.
    """
    src = Path(metadata_src).expanduser().resolve()
    if not src.exists():
        raise FileNotFoundError(f"Metadata not found: {src}")

    legacy_abs = Path.cwd() / legacy_relpath
    legacy_abs.parent.mkdir(parents=True, exist_ok=True)

    # If something already exists there, replace if it points elsewhere.
    if legacy_abs.exists() or legacy_abs.is_symlink():
        try:
            current_target = legacy_abs.resolve()
        except FileNotFoundError:
            current_target = None
        if current_target != src:
            legacy_abs.unlink(missing_ok=True)
        else:
            return legacy_abs

    if prefer_symlink:
        try:
            legacy_abs.symlink_to(src)
            return legacy_abs
        except (OSError, NotImplementedError):
            pass  # fall back to copy on restricted filesystems

    shutil.copy2(src, legacy_abs)
    return legacy_abs


def write_report_row(*, report_path: Path, fields: Dict[str, str]) -> None:
    """Append a single row to the run report TSV, creating header if needed.

    Parameters
    ----------
    report_path : pathlib.Path
        Path to the report TSV.
    fields : dict
        Mapping from column name to value for this row.
    """
    report_path.parent.mkdir(parents=True, exist_ok=True)
    is_new = not report_path.exists()
    with report_path.open("a", encoding="utf-8") as fh:
        if is_new:
            fh.write("\t".join(fields.keys()) + "\n")
        fh.write("\t".join(str(v) for v in fields.values()) + "\n")


def ensure_dirs(*, out_dir: Path) -> Paths:
    """Create the run directory tree and return a populated Paths object.

    Args:
        out_dir: Root output directory for this run.

    Returns:
        A ``Paths`` instance with subdirectories created.
    """
    out_dir = Path(out_dir)
    paths = Paths(root=out_dir)   # <-- important: use 'root=', not 'out_dir='
    paths.mkdirs()
    return paths



def _canonicalise_metadata_first_col(*, metadata_tsv: Path, logger: Optional[logging.Logger] = None) -> Path:
    """
    Ensure the first column header is exactly 'sample-id'.

    Accepts common variants like 'sampleID', 'SampleID', 'sample_id', 'Sample_ID'.
    If a variant is detected, writes a temporary corrected copy and returns its path.
    Otherwise returns the original path unchanged.

    Args:
        metadata_tsv: Path to the provided metadata TSV.

    Returns:
        Path to a metadata TSV whose first column header is exactly 'sample-id'.
    """
    first, second = _read_tsv_first_two_lines(metadata_tsv)

    # Extract the very first header token (left of first tab, or entire line)
    raw_first_col = (first.split("\t")[0]).strip()

    # Normalise variants: lower-case and drop hyphens/underscores/spaces
    norm = raw_first_col.lower().replace("-", "").replace("_", "").replace(" ", "")

    if norm == "sampleid":
        # Rewrite a corrected copy with 'sample-id' as the first column header
        with metadata_tsv.open("r", encoding="utf-8", errors="replace") as src, \
            NamedTemporaryFile("w", encoding="utf-8", delete=False) as tmp:
            fixed_path = Path(tmp.name)
            header = src.readline().rstrip("\n\r")
            rest = src.read()
            # Replace only the very first column header token
            if "\t" in header:
                header_tokens = header.split("\t")
                header_tokens[0] = "sample-id"
                header_fixed = "\t".join(header_tokens)
            else:
                header_fixed = "sample-id"
            tmp.write(header_fixed + "\n")
            tmp.write(rest)
        # Warn if missing '#q2:types' (not fatal)
        if not (second or "").startswith("#q2:types"):
            logger.warning("No '#q2:types' on second line of metadata.")
        return fixed_path

    # If it’s already correct, just warn about missing types line (not fatal)
    if raw_first_col != "sample-id":
        # If it’s some *other* unexpected header, fail fast with a clear message
        raise ValueError(
            f"Metadata first column must be 'sample-id' (or a common variant). "
            f"Found: '{raw_first_col}'"
        )
    if not (second or "").startswith("#q2:types"):
        logger.warning("No '#q2:types' on second line of metadata.")
        sys.stderr.write("[warn] No '#q2:types' on second line of metadata.\n")

    return metadata_tsv



# ------------------------ metadata/manifest checks ------------------------ #
def _read_tsv_first_two_lines(path: Path) -> tuple[str, str]:
    """Read the first two lines of a TSV file.

    Parameters
    ----------
    path : pathlib.Path
        Path to the TSV file.

    Returns
    -------
    tuple of str
        First and second lines with trailing newlines removed.
    """
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline().rstrip("\n\r")
        second = fh.readline().rstrip("\n\r")
    return first, second


def _collect_sample_ids(metadata_tsv: Path) -> list[str]:
    """Collect sample IDs from a QIIME-style metadata TSV.

    Skips the header row and any subsequent comment lines that begin with '#'.
    Does not assume a '#q2:types' line is present.
    """
    ids: list[str] = []
    with metadata_tsv.open("r", encoding="utf-8", errors="replace") as fh:
        # consume header
        header = fh.readline()
        # read remaining lines
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):  # only skip if it really starts with '#'
                continue
            parts = line.rstrip("\n\r").split("\t")
            ids.append(parts[0].strip())
    return ids


def validate_manifest_and_metadata(*, manifest: Path, metadata_tsv: Path) -> None:
    """
    Validate that all sample IDs in the manifest and metadata overlap.

    This version tolerates underscores and dots as equivalent, so that
    `ES_P202` and `ES.P202` are treated as matching IDs.  It still reports
    the original identifiers in any mismatch message.

    Args:
        manifest:
            Path to the QIIME 2 manifest (tab-delimited; must include 'sample-id').
        metadata_tsv:
            Path to the metadata file (tab-delimited; must include a header row).

    Raises:
        ValueError:
            If the manifest is empty or if sample IDs do not overlap after
            normalised comparison.
    """
    # ---------------------------- Read manifest ---------------------------- #
    with manifest.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "sample-id" not in reader.fieldnames:
            raise ValueError("Manifest must include a 'sample-id' column.")
        mani_ids = [row["sample-id"].strip() for row in reader if row.get("sample-id")]
    if not mani_ids:
        raise ValueError("No rows detected in manifest (after header).")

    # ---------------------------- Read metadata ---------------------------- #
    with metadata_tsv.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Metadata file has no header row.")
        sample_col = reader.fieldnames[0]
        meta_ids = [row[sample_col].strip() for row in reader if row.get(sample_col)]

    # ----------------------- Normalise for comparison ---------------------- #
    def _cmp_norm(s: str) -> str:
        """Return a normalised version of an ID for set comparison."""
        return (s or "").strip().replace("_", ".").lower()

    set_meta_raw = set(meta_ids)
    set_mani_raw = set(mani_ids)

    set_meta_cmp = {_cmp_norm(x) for x in set_meta_raw}
    set_mani_cmp = {_cmp_norm(x) for x in set_mani_raw}

    # ----------------------------- Evaluate -------------------------------- #
    if set_meta_cmp != set_mani_cmp:
        only_in_meta = sorted(x for x in set_meta_raw if _cmp_norm(x) not in set_mani_cmp)
        only_in_mani = sorted(x for x in set_mani_raw if _cmp_norm(x) not in set_meta_cmp)

        msg_parts = ["Mismatch between metadata and manifest sample IDs."]
        if only_in_meta:
            msg_parts.append("Only in metadata: " + ", ".join(only_in_meta))
        if only_in_mani:
            msg_parts.append("Only in manifest: " + ", ".join(only_in_mani))
        raise ValueError(" \n".join(msg_parts))

    # ----------------------------- Success --------------------------------- #
    print(f"[INFO] Manifest–metadata validation OK ({len(set_mani_raw)} samples).",
          flush=True)




# ---------------------------- Flash steps ---------------------------- #

def parse_paired_manifest(*, manifest: Path) -> list[tuple[str, Path, Path]]:
    """
    Parse a PairedEndFastqManifestPhred33V2 into (sample_id, R1, R2) rows.

    Parameters
    ----------
    manifest : pathlib.Path
        Path to the paired-end manifest TSV.

    Returns
    -------
    list of (sample_id, r1_path, r2_path)
        Parsed rows with absolute Paths.

    Raises
    ------
    FileNotFoundError
        If the manifest or any referenced FASTQ file does not exist.
    ValueError
        If manifest header is unexpected.
    """
    rows: list[tuple[str, Path, Path]] = []
    with manifest.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        expected = ["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"]
        if header[:3] != expected:
            raise ValueError(f"Unexpected paired manifest header: {header[:3]} != {expected}")
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n\r").split("\t")
            sid = parts[0].strip()
            r1 = Path(parts[1].strip()).expanduser().resolve()
            r2 = Path(parts[2].strip()).expanduser().resolve()
            if not r1.exists():
                raise FileNotFoundError(f"R1 not found for {sid}: {r1}")
            if not r2.exists():
                raise FileNotFoundError(f"R2 not found for {sid}: {r2}")
            rows.append((sid, r1, r2))
    return rows




def resolve_flash_binary(*, user_arg: Optional[str]) -> str:
    """
    Resolve the FLASH/FLASH2 executable to use.

    If `user_arg` is provided, it is returned as-is (trust the user).
    Otherwise, search the current PATH for a suitable binary, preferring
    'flash2' then 'flash'. Raise a clear error if neither is found.

    Parameters
    ----------
    user_arg : Optional[str]
        Executable name or path supplied by the user (may be None).

    Returns
    -------
    str
        Executable name or absolute path to run.

    Raises
    ------
    FileNotFoundError
        If no suitable binary is found in PATH.
    """
    if user_arg:
        return user_arg

    for candidate in ("flash2", "flash"):
        found = shutil.which(candidate)
        if found:
            return found

    raise FileNotFoundError(
        "Could not locate 'flash2' or 'flash' on PATH. "
        "Install it or provide --flash_binary /full/path/to/flash."
    )


def run_flash_join(
    *,
    flash_bin: str,
    sample_id: str,
    r1_fastq: Path,
    r2_fastq: Path,
    out_dir: Path,
    min_overlap: int,
    max_overlap: int,
    mismatch_ratio: float,
    logs_dir: Path,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """
    Run FLASH/FLASH2 to merge paired reads for a single sample.

    Parameters
    ----------
    flash_bin : str
        Executable name/path, e.g. 'flash' or 'flash2'.
    sample_id : str
        Sample identifier used to name outputs.
    r1_fastq : pathlib.Path
        Path to forward FASTQ.
    r2_fastq : pathlib.Path
        Path to reverse FASTQ.
    out_dir : pathlib.Path
        Output folder for joined FASTQ.
    min_overlap : int
        Minimum overlap (-m).
    max_overlap : int
        Maximum overlap (-M).
    mismatch_ratio : float
        Maximum mismatch density in the overlap (-x).
    logs_dir : pathlib.Path
        Directory for log files.

    Returns
    -------
    pathlib.Path
        Path to the joined FASTQ (extended fragments).

    Raises
    ------
    subprocess.CalledProcessError
        If FLASH exits non-zero.
    FileNotFoundError
        If the expected joined file is not produced.
    """
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = sample_id
    log_file = logs_dir / f"04a_flash_{sample_id}.log"

    cmd = [
        flash_bin,
        "-m", str(min_overlap),
        "-M", str(max_overlap),
        "-x", str(mismatch_ratio),
        "-d", str(out_dir),          # ensure outputs go in out_dir
        "-o", prefix,                # file prefix within out_dir
        str(r1_fastq),
        str(r2_fastq),
    ]
    run_cmd(cmd=cmd, log_file=log_file, logger=logger)

    # Accept either plain or gzipped output
    joined_plain = out_dir / f"{prefix}.extendedFrags.fastq"
    joined_gz    = out_dir / f"{prefix}.extendedFrags.fastq.gz"
    if joined_plain.exists():
        return joined_plain
    if joined_gz.exists():
        return joined_gz

    raise FileNotFoundError(
        f"FLASH did not produce {joined_plain} or {joined_gz}. See log: {log_file}"
    )



def write_single_end_manifest(
    *,
    rows: list[tuple[str, Path]],
    out_manifest: Path,
) -> None:
    """
    Write a SingleEndFastqManifestPhred33V2 file from (sample_id, fastq) pairs.

    Parameters
    ----------
    rows : list of (sample_id, fastq_path)
        Joined FASTQ files per sample.
    out_manifest : pathlib.Path
        Where to write the manifest. Tab-separated, UTF-8.
    """
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with out_manifest.open("w", encoding="utf-8") as fh:
        fh.write("sample-id\tabsolute-filepath\n")
        for sid, fq in rows:
            fh.write(f"{sid}\t{fq}\n")


def qiime_import_single_end(*, manifest: Path, out_artifact: Path, logs: Path, logger: Optional[logging.Logger] = None) -> None:
    """
    Import single-end reads using a single-end manifest.
    """
    cmd = [
        "qiime", "tools", "import",
        "--type", "SampleData[SequencesWithQuality]",
        "--input-path", str(manifest),
        "--output-path", str(out_artifact),
        "--input-format", "SingleEndFastqManifestPhred33V2",
    ]
    run_cmd(cmd=cmd, log_file=logs / "01_import_single_end.log", logger=logger)


# ----------------------------- QIIME steps ----------------------------- #
def qiime_import_paired(*, manifest: Path, out_artifact: Path, logs: Path, logger: Optional[logging.Logger] = None) -> None:
    """Import paired-end reads using a paired manifest."""
    cmd = [
        "qiime", "tools", "import",
        "--type", "SampleData[PairedEndSequencesWithQuality]",
        "--input-path", str(manifest),
        "--output-path", str(out_artifact),
        "--input-format", "PairedEndFastqManifestPhred33V2",
    ]
    run_cmd(cmd=cmd, log_file=logs / "01_import.log", logger=logger)


def qiime_demux_summary(*, demux_qza: Path, out_qzv: Path, logs: Path, logger: Optional[logging.Logger] = None) -> None:
    """Create a demultiplexing summary visualisation."""
    cmd = [
        "qiime", "demux", "summarize",
        "--i-data", str(demux_qza),
        "--o-visualization", str(out_qzv),
    ]
    run_cmd(cmd=cmd, log_file=logs / "02_demux_summarize.log", logger=logger)


def ensure_qiime2_tmp(*, tmp_root: Path) -> Path:
    """
    Ensure QIIME 2 temp cache directory exists with sticky-world-writable perms (01777).

    Parameters
    ----------
    tmp_root : pathlib.Path
        The run's tmp directory (e.g., paths.root / "tmp").

    Returns
    -------
    pathlib.Path
        Path to the qiime2 temp dir.

    Raises
    ------
    PermissionError
        If we cannot set the required mode and the filesystem rejects it.
    """
    qdir = Path(tmp_root) / "qiime2"
    qdir.mkdir(parents=True, exist_ok=True)
    wanted = 0o41777  # sticky + rwx for all
    try:
        os.chmod(qdir, wanted)
    except PermissionError as e:
        # Try to recover by removing and recreating with correct mode
        try:
            qdir.rmdir()
        except OSError:
            pass
        qdir.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(qdir, wanted)
        except PermissionError:
            raise PermissionError(
                f"Cannot set sticky perms on {qdir}. Try running `chmod 1777 {qdir}` manually "
                "or choose a different --out_dir on a filesystem that supports the sticky bit."
            ) from e
    return qdir


def qiime_cutadapt_trim_paired(
    *, in_qza: Path, out_qza: Path,
    front_f: Optional[list[str]], front_r: Optional[list[str]],
    discard_untrimmed: bool, logs: Path, logger: Optional[logging.Logger] = None
) -> None:
    """
    Run q2-cutadapt trim-paired with support for multiple forward and reverse primers.

    Parameters
    ----------
    in_qza : Path
        Input demuxed reads.
    out_qza : Path
        Output trimmed sequences.
    front_f : list of str or None
        List of forward primers.
    front_r : list of str or None
        List of reverse primers.
    discard_untrimmed : bool
        Whether to discard reads without primer matches.
    logs : Path
        Directory for logs.
    """
    if front_f and front_r:
        cmd = [
            "qiime", "cutadapt", "trim-paired",
            "--i-demultiplexed-sequences", str(in_qza),
            "--p-discard-untrimmed", str(discard_untrimmed).lower(),
            "--o-trimmed-sequences", str(out_qza),
        ]

        # add repeated flags for each primer
        for f in front_f:
            cmd.extend(["--p-front-f", f])
        for r in front_r:
            cmd.extend(["--p-front-r", r])

        run_cmd(cmd=cmd, log_file=logs / "03_cutadapt_trim_paired.log", logger=logger)

    else:
        shutil.copy2(in_qza, out_qza)



def filter_asvs_by_length(*, repseqs_qza: Path, table_qza: Path,
                          out_repseqs_qza: Path, out_table_qza: Path,
                          length_range: tuple[int, int], work_dir: Path,
                          logs: Path, logger=None):
    """Keep only features whose sequence length is within [min_len, max_len]."""
    keep_ids = []
    tmp = work_dir / "lenfilter"
    tmp.mkdir(parents=True, exist_ok=True)

    # 1) export rep-seqs to fasta
    run_cmd(cmd=["qiime", "tools", "export",
                 "--input-path", str(repseqs_qza),
                 "--output-path", str(tmp)],
            log_file=logs / "08_lenfilter_export_repseqs.log", logger=logger)

    fasta = tmp / "dna-sequences.fasta"

    # 2) parse fasta to collect IDs within window
    min_len, max_len = length_range
    with fasta.open() as fh:
        seq_id = None
        seq = []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_id is not None:
                    L = sum(len(x) for x in seq)
                    if min_len <= L <= max_len:
                        keep_ids.append(seq_id)
                seq_id = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if seq_id is not None:
            L = sum(len(x) for x in seq)
            if min_len <= L <= max_len:
                keep_ids.append(seq_id)

    if not keep_ids:
        raise RuntimeError(f"No ASVs within length range {min_len}:{max_len}")

    # 3) write metadata-style ID file
    ids_tsv = tmp / "ids_to_keep.tsv"
    with ids_tsv.open("w") as out:
        out.write("Feature ID\n")
        for k in keep_ids:
            out.write(f"{k}\n")

    # 4) filter rep-seqs
    run_cmd(cmd=["qiime", "feature-table", "filter-seqs",
                 "--i-data", str(repseqs_qza),
                 "--m-metadata-file", str(ids_tsv),
                 "--o-filtered-data", str(out_repseqs_qza)],
            log_file=logs / "08_lenfilter_filter_seqs.log", logger=logger)

    # 5) filter table (features)
    run_cmd(cmd=["qiime", "feature-table", "filter-features",
                 "--i-table", str(table_qza),
                 "--m-metadata-file", str(ids_tsv),
                 "--o-filtered-table", str(out_table_qza)],
            log_file=logs / "08_lenfilter_filter_table.log", logger=logger)


def qiime_dada2_denoise_paired(
    *, in_qza: Path, table_qza: Path, repseqs_qza: Path, stats_qza: Path,
    trim_left_f: int, trim_left_r: int, trunc_len_f: int, trunc_len_r: int,
    max_ee_f: float, max_ee_r: float, threads: int, logs: Path,
    logger: Optional[logging.Logger] = None) -> None:
    """Run DADA2 paired-end denoising in QIIME 2."""
    cmd = [
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", str(in_qza),
        "--p-trim-left-f", str(trim_left_f),
        "--p-trim-left-r", str(trim_left_r),
        "--p-trunc-len-f", str(trunc_len_f),
        "--p-trunc-len-r", str(trunc_len_r),
        "--p-max-ee-f", str(max_ee_f),
        "--p-max-ee-r", str(max_ee_r),
        "--p-n-threads", str(threads),
        "--o-table", str(table_qza),
        "--o-representative-sequences", str(repseqs_qza),
        "--o-denoising-stats", str(stats_qza),
    ]
    run_cmd(cmd=cmd, log_file=logs / "04_dada2_paired.log", logger=logger)


def qiime_vsearch_join_pairs(*, in_qza: Path, out_qza: Path, logs: Path, logger: Optional[logging.Logger] = None) -> None:
    """Join paired-end reads inside QIIME 2 using vsearch."""
    cmd = [
        "qiime", "vsearch", "join-pairs",
        "--i-demultiplexed-seqs", str(in_qza),
        "--o-joined-sequences", str(out_qza),
    ]
    run_cmd(cmd=cmd, log_file=logs / "04a_vsearch_join_pairs.log", logger=logger)


def qiime_deblur_denoise_single(
    *, in_qza: Path, table_qza: Path, repseqs_qza: Path,
    trim_length: int, threads: int, logs: Path, 
    logger: Optional[logging.Logger] = None) -> None:
    """Run Deblur denoising on single-end (joined or forward-only) reads."""
    if trim_length <= 0:
        raise ValueError("Deblur requires --deblur_trim_length > 0.")
    cmd = [
        "qiime", "deblur", "denoise-16S",
        "--i-demultiplexed-seqs", str(in_qza),
        "--p-trim-length", str(trim_length),
        "--p-jobs-to-start", str(threads),
        "--o-table", str(table_qza),
        "--o-representative-sequences", str(repseqs_qza),
        "--o-stats", str(in_qza.parent / "deblur_stats.qza"),
    ]
    run_cmd(cmd=cmd, log_file=logs / "05_deblur_single.log", logger=logger)


def qiime_phylogeny_mafft_fasttree(
    *, repseqs_qza: Path, aligned_qza: Path, masked_qza: Path,
    unrooted_qza: Path, rooted_qza: Path, threads: int, logs: Path,
    logger: Optional[logging.Logger] = None) -> None:
    """Build a phylogenetic tree using MAFFT alignment and FastTree."""
    cmd = [
        "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
        "--i-sequences", str(repseqs_qza),
        "--p-n-threads", str(threads),
        "--o-alignment", str(aligned_qza),
        "--o-masked-alignment", str(masked_qza),
        "--o-tree", str(unrooted_qza),
        "--o-rooted-tree", str(rooted_qza),
    ]
    run_cmd(cmd=cmd, log_file=logs / "06_phylogeny_mafft_fasttree.log", logger=logger)


def qiime_taxonomy_sklearn(
    *, repseqs_qza: Path, classifier_qza: Path, taxonomy_qza: Path, logs: Path,
    logger: Optional[logging.Logger] = None) -> None:
    """Assign taxonomy using a pre-trained sklearn classifier artefact."""
    cmd = [
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-classifier", str(classifier_qza),
        "--i-reads", str(repseqs_qza),
        "--o-classification", str(taxonomy_qza),
    ]
    run_cmd(cmd=cmd, log_file=logs / "07_taxonomy_sklearn.log", logger=logger)


# ----------------------------- main orchestration ----------------------------- #
def build_arg_parser() -> argparse.ArgumentParser:
    """Construct the command-line interface for the runner.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser with named-only arguments.
    """
    p = argparse.ArgumentParser(
        description=("QIIME 2 amplicon runner with DADA2 or Deblur. Named arguments only."),
        allow_abbrev=False,
    )
    p.add_argument("--run_label", required=True, type=str, help="Run label.")
    p.add_argument("--manifest", required=True, type=Path, help="Paired manifest TSV.")
    p.add_argument("--metadata_tsv", required=True, type=Path, help="QIIME metadata TSV.")
    p.add_argument("--denoiser", required=True, choices=["dada2_paired", "deblur_single"], help="Choose denoiser.")
    p.add_argument("--pair_strategy", choices=["join_in_qiime", "join_with_flash", 
                                               "forward_only"], default="join_with_flash",
                   help="For Deblur: join pairs in QIIME or use forward-only reads.")
    p.add_argument("--cutadapt_f", default=None, type=str, help="Forward primer.")
    p.add_argument("--cutadapt_r", default=None, type=str, help="Reverse primer.")
    p.add_argument(
        "--cutadapt_discard_untrimmed",
        action="store_true",
        help="Discard reads lacking primer matches."
    )

    p.add_argument(
        "--no-cutadapt_discard_untrimmed",
        dest="cutadapt_discard_untrimmed",
        action="store_false",
        help="Do not discard reads lacking primer matches."
    )

    p.set_defaults(cutadapt_discard_untrimmed=False)

    # DADA2 
    p.add_argument("--trim_left_f", default=0, type=int, help="DADA2 trim-left F.")
    p.add_argument("--trim_left_r", default=0, type=int, help="DADA2 trim-left R.")
    p.add_argument("--trunc_len_f", default=0, type=int, help="DADA2 trunc-len F.")
    p.add_argument("--trunc_len_r", default=0, type=int, help="DADA2 trunc-len R.")
    # max_ee of 2.0 was too strict in our tests, leading to very few reads passing. 
    # 5.0 is more reasonable for typical Illumina data, but users can adjust as needed.
    p.add_argument("--max_ee_f", default=5.0, type=float, help="DADA2 maxEE F.")
    p.add_argument("--max_ee_r", default=5.0, type=float, help="DADA2 maxEE R.")

    p.add_argument("--asv_length_range", default="100:970", type=str,
               help="Length window for ASVs, e.g. '100:970'. If omitted, no length filter is applied.")

    # trimming
    p.add_argument("--primer_preset",
            choices=["auto","v34_locus_only","v34_overhangs"],
            default="v34_overhangs",
            help=("Primer preset. fallback to 'v34' (locus-only):  'auto' uses any explicit --cutadapt_front_* if provided, "
                "else falls back to 'v34' (locus-only). Use 'v34_overhangs' when Nextera "
                "overhangs are still present."))

    # Deblur knob
    p.add_argument("--deblur_trim_length", default=0, type=int, help="Deblur trim length.")
    # flash merge options
    p.add_argument("--flash_min_overlap", type=int, default=20,
                help="FLASH: minimum required overlap (-m).")
    p.add_argument("--flash_max_overlap", type=int, default=300,
                help="FLASH: maximum overlap (-M).")
    p.add_argument("--flash_mismatch_ratio", type=float, default=0.10,
                help="FLASH: maximum mismatch density in overlap (-x).")
    p.add_argument(
        "--flash_binary",
        type=str,
        default=None,  # auto-detect when omitted
        help="FLASH/FLASH2 executable. If omitted, auto-detect from PATH.",
    )
    # Common
    p.add_argument("--threads", default=4, type=int, help="Threads for heavy steps.")
    p.add_argument("--classifier_qza", default=None, type=Path,
                   help="Pre-trained sklearn classifier (.qza) matched to region.")
    p.add_argument("--out_dir", required=True, type=Path, help="Output directory for the run.")
    return p


def main() -> None:
    """Entry point for the QIIME 2 amplicon runner.

    This function parses arguments, prepares folders, validates metadata
    and manifest, runs the chosen QIIME 2 steps, and assembles outputs for
    downstream use including a phyloseq-ready folder.
    """
    args = build_arg_parser().parse_args()
    logger = setup_logging(out_dir=args.out_dir, run_label=args.run_label)

    # Log explicit start time and initial RAM
    logger.info("Start time: %s", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(_SCRIPT_START_TIME)))
    log_memory_usage(logger=logger, prefix="START")

    args.out_dir = Path(args.out_dir).expanduser().resolve()
    # Always wipe logs before starting
    logs_dir = args.out_dir / "logs"

    # Instead of shutil.rmtree(logs_dir) - this fails on our filesystem when files are 
    # in use (e.g. if the user tries to re-run without closing previous logs). 
    # Instead, we try to rename it out of the way first, then remove the old one asynchronously. 
    # This way we avoid blocking on file locks and ensure a clean logs directory for the new run.
    rotated = safe_rotate_then_rmtree(Path(logs_dir))
    if rotated is not None and Path(rotated).exists():
        logger.warning("Could not fully remove logs directory: %s", rotated)


    logger = setup_logging(out_dir=args.out_dir, run_label=args.run_label)
    args.manifest = Path(args.manifest).expanduser().resolve()
    args.metadata_tsv = Path(args.metadata_tsv).expanduser().resolve()

    logger.info("CWD=%s", Path.cwd())
    logger.info("Denoiser=%s | Pair strategy=%s", args.denoiser, args.pair_strategy)

    paths = Paths(Path(args.out_dir))
    paths.mkdirs()  # ensure folders exist before any command
    logger.info("out_dir=%s denoiser=%s", paths.root, args.denoiser)


    tmp = paths.root / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(tmp)
    os.environ["TEMP"] = str(tmp)
    os.environ["TMP"] = str(tmp)
    os.environ["QIIMETMPDIR"] = str(tmp)
    os.environ["XDG_CACHE_HOME"] = str(paths.root / ".cache")
    q2tmp = ensure_qiime2_tmp(tmp_root=tmp)
    logger.info("QIIME2_TMP=%s", q2tmp)

    meta_for_qiime = make_qiime_metadata_copy(
    src_tsv=args.metadata_tsv,
    out_dir=paths.root / "tmp",
    add_q2_types=False,  # set True if you want the types line
            )
    # Make legacy path (symlink/copy) point to the canonical copy
    legacy_meta = ensure_legacy_metadata_location(
        metadata_src=meta_for_qiime,
        legacy_relpath=Path("metadata/metadata16S.tsv"),
        prefer_symlink=True,
    )

    # Internal validation (manifest/metadata) — validate against the canonical copy
    validate_manifest_and_metadata(manifest=args.manifest, 
                                   metadata_tsv=meta_for_qiime)

    # loads of faff here tyriing to sort which primers
    # Resolve primer settings

    # ---------------- Primer selection (Option 1: soft anchoring) ---------------- #
    # Rules:
    #   - If the user supplies --cutadapt_f / --cutadapt_r, trust them fully.
    #   - If primer_preset is v34 or v34_overhangs, supply *unanchored* primers.
    #   - If primer_preset is auto and user gave no primers, fall back to unanchored V3–V4.
    #   - No caret (^) anchoring is ever inserted automatically.

    # ---------------- Primer selection (Option 1: soft anchoring) ---------------- #

    # Final primer strings to pass to cutadapt
    primer_f = None
    primer_r = None

    if args.cutadapt_f or args.cutadapt_r:
        # User explicitly supplied primers
        if not (args.cutadapt_f and args.cutadapt_r):
            logger.error("Provide BOTH --cutadapt_f and --cutadapt_r, or neither.")
            sys.exit(2)
        primer_f = [args.cutadapt_f]
        primer_r = [args.cutadapt_r]
        logger.info("Using user-supplied cutadapt primers.")
    else:

        # New full adapter sets
        ALL_FORWARD = [
            "AGTCAGTCAGCCGGACTACHVGGGTWTCTAAT",  # 515F
            "GTGCCAGCMGCCGCGGTAA",               # 515FC
            "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG",  # 241F
            "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",                    # FWoverhang
        ]

        ALL_REVERSE = [
            "GGACTACHVGGGTWTCTAAT",              # 806R
            "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC",  # 785R
            "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",                        # RVoverhang
        ]

        # override final primers
        primer_f = ALL_FORWARD
        primer_r = ALL_REVERSE
        logger.info("Using multiple cutadapt primers (F=%d, R=%d).", len(primer_f), len(primer_r))


    # Log the final primers actually used
    logger.info("Final primers: F=%s", primer_f)
    logger.info("Final primers: R=%s", primer_r)
    logger.info("Cutadapt will search %d forward and %d reverse primers.", len(primer_f), len(primer_r))




    # sanity: both or neither
    if bool(args.cutadapt_f) ^ bool(args.cutadapt_r):
        logger.error("Provide BOTH --cutadapt_f and --cutadapt_r (or use --primer_preset).")
        sys.exit(2)


    paths = ensure_dirs(out_dir=args.out_dir)
    logger.info("out_dir=%s denoiser=%s", paths.root, args.denoiser)

    # Optional: external checker script (non-interactive)
    checker = shutil.which("check_map.sh")
    if checker:
        try:
            run_cmd(
                cmd=[checker, "-m", str(args.metadata_tsv), "-M", str(args.manifest), "--no-fix"],
                log_file=paths.logs / "00_check_map.log", logger=logger,
            )
        except subprocess.CalledProcessError:
            logger.error("Metadata check failed (check_map.sh). See log for details.")
            sys.exit(2)

    # 1) Import
    try:
        run_cmd(cmd=["qiime", "info"], log_file=paths.logs / "00_qiime_info.log", logger=logger)
    except Exception as e:
        logger.debug("qiime info failed: %s", e)

    demux_paired = paths.qza / "paired-demux.qza"
    qiime_import_paired(manifest=args.manifest, out_artifact=demux_paired, logs=paths.logs, 
                        logger=logger)

    # 2) Demux summary (pre-trim)
    demux_qzv_pre = paths.qzv / "demux_pre_trim.qzv"
    qiime_demux_summary(demux_qza=demux_paired, out_qzv=demux_qzv_pre, logs=paths.logs, 
                        logger=logger)

    # 3) Cutadapt (optional)
    demux_trimmed = paths.qza / "paired-demux-trimmed.qza"
    qiime_cutadapt_trim_paired(
        in_qza=demux_paired,
        out_qza=demux_trimmed,
        front_f=primer_f,
        front_r=primer_r,
        discard_untrimmed=args.cutadapt_discard_untrimmed,
        logs=paths.logs,
    )

    # 4) Demux summary (post-trim)
    demux_qzv_post = paths.qzv / "demux_post_trim.qzv"
    qiime_demux_summary(demux_qza=demux_trimmed, out_qzv=demux_qzv_post, logs=paths.logs,
                        logger=logger)

    # 5) Denoise branch
    table_qza = paths.qza / "table.qza"
    repseqs_qza = paths.qza / "rep-seqs.qza"
    stats_qza = paths.qza / "denoise-stats.qza"

    if args.denoiser == "dada2_paired":
        qiime_dada2_denoise_paired(
            in_qza=demux_trimmed,
            table_qza=table_qza,
            repseqs_qza=repseqs_qza,
            stats_qza=stats_qza,
            trim_left_f=args.trim_left_f,
            trim_left_r=args.trim_left_r,
            trunc_len_f=args.trunc_len_f,
            trunc_len_r=args.trunc_len_r,
            max_ee_f=args.max_ee_f,
            max_ee_r=args.max_ee_r,
            threads=args.threads,
            logs=paths.logs,
        )
    else:
        # Deblur path: join pairs in QIIME or use forward-only (join is preferred)
        if args.pair_strategy == "join_in_qiime":
            # existing vsearch join
            joined_qza = paths.qza / "joined-single-end.qza"
            qiime_vsearch_join_pairs(
                in_qza=demux_trimmed,
                out_qza=joined_qza,
                logs=paths.logs,
            )
            deblur_in = joined_qza

        elif args.pair_strategy == "forward_only":
            # use the trimmed paired artefact but Deblur will ignore R2
            deblur_in = demux_trimmed

        elif args.pair_strategy == "join_with_flash":
            # Pre-join on raw FASTQs listed in the paired manifest, then import single-end
            flash_exe = resolve_flash_binary(user_arg=args.flash_binary)
            logger.info("FLASH_BINARY=%s", flash_exe)

            flash_dir = paths.root / "flash_joined"
            flash_dir.mkdir(parents=True, exist_ok=True)
            logger.info("FLASH_OUT_DIR=%s", flash_dir)

            triples = parse_paired_manifest(manifest=args.manifest)
            logger.info("FLASH_SAMPLES=%s", ",".join([sid for sid, _, _ in triples]))

            joined_rows: list[tuple[str, Path]] = []
            for sid, r1, r2 in triples:
                joined = run_flash_join(
                    flash_bin=flash_exe,
                    sample_id=sid,
                    r1_fastq=r1,
                    r2_fastq=r2,
                    out_dir=flash_dir,
                    min_overlap=args.flash_min_overlap,
                    max_overlap=args.flash_max_overlap,
                    mismatch_ratio=args.flash_mismatch_ratio,
                    logs_dir=paths.logs,
                    logger=logger,
                )
                joined_rows.append((sid, joined))
            if not joined.exists() or joined.stat().st_size == 0:
                raise RuntimeError(f"FLASH produced no merged reads for sample {sid}. "
                                "Try adjusting overlap/mismatch, or use '--pair_strategy forward_only'.")
            
            se_manifest = flash_dir / "single_end_manifest.tsv"
            write_single_end_manifest(rows=joined_rows, out_manifest=se_manifest)
            logger.info("FLASH_SINGLE_END_MANIFEST=%s", se_manifest)

            # Import single-end artefact
            demux_single = paths.qza / "single-demux.qza"
            qiime_import_single_end(manifest=se_manifest, out_artifact=demux_single, logs=paths.logs, 
                                    logger=logger)

            # Optional quick summary of single-end demux
            demux_qzv_single = paths.qzv / "demux_single.qzv"
            qiime_demux_summary(demux_qza=demux_single, out_qzv=demux_qzv_single, logs=paths.logs,
                                logger=logger)

            deblur_in = demux_single

        else:
            raise ValueError(f"Unknown pair_strategy: {args.pair_strategy}")

        # Deblur as before
        # Determine Deblur trim length (explicit or auto)
        deblur_trim = int(args.deblur_trim_length)

        if deblur_trim == 0:
            if args.pair_strategy == "join_with_flash":
                # Infer from the actual merged reads that will be imported
                fastqs = [fq for _, fq in joined_rows]
                deblur_trim = infer_deblur_trim_length(
                    fastq_paths=fastqs,
                    percentile=10.0,
                    min_len=150,
                    max_len=None,
                    warn_below=200,
                    logger=logger,
                )
            elif args.pair_strategy == "forward_only":
                # Infer from forward reads in the paired manifest
                triples = parse_paired_manifest(manifest=args.manifest)
                fastqs = [r1 for _, r1, _ in triples]
                deblur_trim = infer_deblur_trim_length(
                    fastq_paths=fastqs,
                    percentile=10.0,
                    min_len=150,
                    max_len=None,
                    warn_below=200,
                    logger=logger,
                )
            else:
                # join_in_qiime: you *can* export joined_qza and sample it, but fallback is OK
                deblur_trim = 200
                logger.warning(
                    "Auto Deblur trim-length with pair_strategy=join_in_qiime is not implemented; "
                    "using fallback %d. Set --deblur_trim_length explicitly if needed.",
                    deblur_trim,
                )

        logger.info("Deblur trim length in use: %d", deblur_trim)

        qiime_deblur_denoise_single(
            in_qza=deblur_in,
            table_qza=table_qza,
            repseqs_qza=repseqs_qza,
            trim_length=deblur_trim,
            threads=args.threads,
            logs=paths.logs,
        )

########


    if args.asv_length_range:
        try:
            Lmin, Lmax = map(int, args.asv_length_range.split(":"))
        except Exception:
            raise ValueError("--asv_length_range must look like 'MIN:MAX', e.g. 440:470")
        filt_rep = paths.qza / "rep-seqs.lengthfilter.qza"
        filt_tab = paths.qza / "table.lengthfilter.qza"
        filter_asvs_by_length(
            repseqs_qza=repseqs_qza, table_qza=table_qza,
            out_repseqs_qza=filt_rep, out_table_qza=filt_tab,
            length_range=(Lmin, Lmax), work_dir=paths.root / "tmp",
            logs=paths.logs, logger=logger
        )
        # continue pipeline using filtered artefacts:
        repseqs_qza = filt_rep
        table_qza   = filt_tab
        logger.info("Applied ASV length filter: %d:%d bp", Lmin, Lmax)


    # 6) Phylogeny
    aligned = paths.phylogeny / "aligned-rep-seqs.qza"
    masked = paths.phylogeny / "masked-aligned-rep-seqs.qza"
    unrooted = paths.phylogeny / "unrooted-tree.qza"
    rooted = paths.phylogeny / "rooted-tree.qza"
    qiime_phylogeny_mafft_fasttree(
        repseqs_qza=repseqs_qza,
        aligned_qza=aligned,
        masked_qza=masked,
        unrooted_qza=unrooted,
        rooted_qza=rooted,
        threads=args.threads,
        logs=paths.logs,
    )

    # 7) Taxonomy (sklearn only in this version)
    taxonomy_qza: Optional[Path] = None
    if args.classifier_qza is not None:
        taxonomy_qza = paths.taxonomy / "taxonomy.qza"
        qiime_taxonomy_sklearn(
            repseqs_qza=repseqs_qza,
            classifier_qza=args.classifier_qza,
            taxonomy_qza=taxonomy_qza,
            logs=paths.logs,
        )
    else:
        logger.warning("No classifier provided; skipping taxonomy assignment.")

    # 8) Prepare phyloseq_output folder
    shutil.copy2(meta_for_qiime, paths.phyloseq_output / "metadata.tsv")
    shutil.copy2(table_qza, paths.phyloseq_output / "table.qza")
    shutil.copy2(repseqs_qza, paths.phyloseq_output / "rep-seqs.qza")
    shutil.copy2(rooted, paths.phyloseq_output / "rooted-tree.qza")
    if taxonomy_qza is not None:
        shutil.copy2(taxonomy_qza, paths.phyloseq_output / "taxonomy.qza")

    # 9) Minimal run report
    write_report_row(
        report_path=paths.report_tsv,
        fields={
            "run_label": args.run_label,
            "denoiser": args.denoiser,
            "threads": str(args.threads),
            "used_cutadapt": str(bool(args.cutadapt_f and args.cutadapt_r)),
            "pair_strategy": "join_in_qiime" if args.denoiser == "deblur_single" else "n/a",
            "deblur_trim_length": str(args.deblur_trim_length if args.denoiser == "deblur_single" else 0),
            "classifier": str(args.classifier_qza) if args.classifier_qza else "none",
        },
    )
    # QC bundle (lives in results/<RUN>/qc)
    run_qc_bundle(paths=paths, metadata_tsv=meta_for_qiime, denoiser=args.denoiser, logger=logger)
    logger.info("QC bundle written to: %s", paths.root / "qc")
    end_ts = time.time()
    logger.info("End time: %s", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_ts)))
    log_memory_usage(logger=logger, prefix="END", extra_msg="Pipeline complete")




if __name__ == "__main__":
    main()
