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
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional
from tempfile import NamedTemporaryFile
import stat
import os
import logging
import csv
import time
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
    log_file.parent.mkdir(parents=True, exist_ok=True)
    if logger is not None:
        logger.info("▶ %s", " ".join(cmd))
        logger.debug("Step log: %s", log_file)

    with log_file.open("a", encoding="utf-8") as lf:
        lf.write("$ " + " ".join(cmd) + "\n")
        lf.flush()
        subprocess.run(cmd, stdout=lf, stderr=lf, check=True)


def run_qc_bundle(
    *, paths: Paths, metadata_tsv: Path, denoiser: str,
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Build a standard QC bundle of QIIME 2 visualisations and text summaries.

    The function consumes core artefacts (feature table, representative
    sequences, rooted tree, optional taxonomy) and writes .qzv outputs and a
    few lightweight TSV/TXT summaries under ``<paths.root>/qc``. It is safe to
    call regardless of which artefacts exist; steps are skipped if inputs are
    missing.

    Created outputs (when inputs are available):
      - feature-table_summarize.qzv
      - tabulate-seqs.qzv
      - dada2_stats.qzv (if denoiser == "dada2_paired")
      - deblur_stats.qzv (if denoiser == "deblur_single")
      - taxa-barplot.qzv (if taxonomy present)
      - alpha-rarefaction.qzv (if rooted tree present; max depth = 10,000)
      - sample_depths.txt  (BIOM per-sample depth summary)
      - asv_lengths.tsv    (FeatureID → ASV length, from rep-seqs FASTA)
      - export_table/, export_repseqs/ (raw exports used by summaries)
    """
    qza_dir = paths.qza
    qzv_dir = paths.qzv
    logs = paths.logs

    out_qc = paths.root / "qc"
    out_qc.mkdir(parents=True, exist_ok=True)

    # Prefer filtered artefacts when present (avoids ID mismatches in QC).
    table_filtered = qza_dir / "table.lengthfilter.qza"
    rep_filtered = qza_dir / "rep-seqs.lengthfilter.qza"

    table_qza = table_filtered if table_filtered.exists() else qza_dir / "table.qza"
    repseqs_qza = (
        rep_filtered if rep_filtered.exists() else qza_dir / "rep-seqs.qza"
    )

    rooted_tree_qza = paths.phylogeny / "rooted-tree.qza"
    stats_qza = qza_dir / "denoise-stats.qza"          # DADA2
    deblur_stats_qza = qza_dir / "deblur_stats.qza"    # Deblur (if present)
    taxonomy_qza = paths.taxonomy / "taxonomy.qza"

    # 1) Feature-table summarize (+ metadata so sample groups show)
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

    # 5) Alpha-rarefaction (requires tree + table)
    if rooted_tree_qza.exists() and table_qza.exists():
        run_cmd(
            cmd=[
                "qiime", "diversity", "alpha-rarefaction",
                "--i-table", str(table_qza),
                "--i-phylogeny", str(rooted_tree_qza),
                "--m-metadata-file", str(metadata_tsv),
                "--p-max-depth", "10000",
                "--o-visualization", str(out_qc / "alpha-rarefaction.qzv"),
            ],
            log_file=logs / "qc_alpha_rarefaction.log",
            logger=logger,
        )

    # 6) Plain-text quick summaries
    # 6a) per-sample depths via BIOM
    export_dir = out_qc / "export_table"
    export_dir.mkdir(parents=True, exist_ok=True)

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

    biom_fp = export_dir / "feature-table.biom"
    if biom_fp.exists():
        run_cmd(
            cmd=[
                "biom", "summarize-table",
                "-i", str(biom_fp),
                "-o", str(out_qc / "sample_depths.txt"),
            ],
            log_file=logs / "qc_biom_summarize.log",
            logger=logger,
        )

    # 6b) ASV lengths: export rep-seqs FASTA, compute lengths
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

            with tsv.open("w", encoding="utf-8") as out:
                out.write("FeatureID\tLength\n")
                for k, v in lengths.items():
                    out.write(f"{k}\t{v}\n")
        except Exception as err:  # noqa: BLE001
            if logger:
                logger.warning("ASV length TSV failed: %s", err)

    # 6c) tiny overall counts file (parsed from BIOM summary)
    counts_tsv = out_qc / "counts_summary.tsv"
    try:
        summary_txt = out_qc / "sample_depths.txt"
        if summary_txt.exists():
            text = summary_txt.read_text(encoding="utf-8", errors="ignore")
            m = re.search(r"Number of samples:\s+(\d+)", text)
            n_samples = int(m.group(1)) if m else -1
            m = re.search(r"Number of observations:\s+(\d+)", text)
            n_taxa = int(m.group(1)) if m else -1
            m = re.search(r"Total count:\s+(\d+)", text)
            total = int(m.group(1)) if m else -1

            with counts_tsv.open("w", encoding="utf-8") as out:
                out.write("samples\ttaxa\ttotal_reads\n")
                out.write(f"{n_samples}\t{n_taxa}\t{total}\n")
    except Exception:
        # Best-effort only; fine to skip on failure.
        pass



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
    front_f: Optional[str], front_r: Optional[str],
    discard_untrimmed: bool, logs: Path, logger: Optional[logging.Logger] = None
) -> None:
    """Run q2-cutadapt trim-paired if primers are supplied; else copy artefact."""
    if front_f and front_r:
        cmd = [
            "qiime", "cutadapt", "trim-paired",
            "--i-demultiplexed-sequences", str(in_qza),
            "--p-front-f", front_f,
            "--p-front-r", front_r,
            "--p-discard-untrimmed", str(discard_untrimmed).lower(),
            "--o-trimmed-sequences", str(out_qza),
        ]
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
    p.add_argument("--cutadapt_discard_untrimmed", default=True,
                   type=lambda x: str(x).lower() in {"1", "true", "yes"},
                   help="Discard reads lacking primer matches.")
    # DADA2 
    p.add_argument("--trim_left_f", default=0, type=int, help="DADA2 trim-left F.")
    p.add_argument("--trim_left_r", default=0, type=int, help="DADA2 trim-left R.")
    p.add_argument("--trunc_len_f", default=0, type=int, help="DADA2 trunc-len F.")
    p.add_argument("--trunc_len_r", default=0, type=int, help="DADA2 trunc-len R.")
    p.add_argument("--max_ee_f", default=2.0, type=float, help="DADA2 maxEE F.")
    p.add_argument("--max_ee_r", default=2.0, type=float, help="DADA2 maxEE R.")

    p.add_argument("--asv_length_range", default="440:470", type=str,
               help="Length window for ASVs, e.g. '440:470'. If omitted, no length filter is applied.")

    # trimming
    p.add_argument("--primer_preset",
            choices=["auto","v34","v34_overhangs"],
            default="auto",
            help=("Primer preset. 'auto' uses any explicit --cutadapt_front_* if provided, "
                "else falls back to 'v34' (locus-only). Use 'v34_overhangs' when Nextera "
                "overhangs are still present."))

    # Deblur knob
    p.add_argument("--deblur_trim_length", default=0, type=int, help="Deblur trim length.")
    # flash merge options
    p.add_argument("--flash_min_overlap", type=int, default=20,
                help="FLASH: minimum required overlap (-m).")
    p.add_argument("--flash_max_overlap", type=int, default=300,
                help="FLASH: maximum overlap (-M).")
    p.add_argument("--flash_mismatch_ratio", type=float, default=0.25,
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

    # loads of faff here tyriing to sort which primers
    # Resolve primer settings
    used_default = False
    if args.primer_preset == "v34_overhangs" and not (args.cutadapt_f or args.cutadapt_r):
        args.cutadapt_f = f"^{V34_FWD_OVH}"
        args.cutadapt_r = f"^{V34_REV_OVH}"
        logger.warning("Using V3–V4 primers WITH Nextera overhangs (anchored).")
    elif not (args.cutadapt_f or args.cutadapt_r):
        # auto -> default to locus-only
        args.cutadapt_f = f"^{V34_FWD}"
        args.cutadapt_r = f"^{V34_REV}"
        used_default = True

    # sanity: both or neither
    if bool(args.cutadapt_f) ^ bool(args.cutadapt_r):
        logger.error("Provide BOTH --cutadapt_f and --cutadapt_r (or use --primer_preset).")
        sys.exit(2)

    if used_default:
        logger.warning(
            "Using default V3–V4 locus-only primers (anchored): F=%s | R=%s. "
            "If reads still include Nextera overhangs, rerun with --primer_preset v34_overhangs.",
            args.cutadapt_f, args.cutadapt_r
        )



    # Make legacy metadata path available for any hard-coded consumers
    legacy_meta = ensure_legacy_metadata_location(
        metadata_src=args.metadata_tsv,
        legacy_relpath=Path("metadata/metadata16S.tsv"),
        prefer_symlink=True,
    )

    paths = ensure_dirs(out_dir=args.out_dir)
    logger.info("out_dir=%s denoiser=%s", paths.root, args.denoiser)

    # Internal validation (manifest/metadata)
    fixed_meta = _canonicalise_metadata_first_col(metadata_tsv=args.metadata_tsv, logger=logger)
    validate_manifest_and_metadata(manifest=args.manifest, metadata_tsv=fixed_meta)


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
        front_f=args.cutadapt_f,
        front_r=args.cutadapt_r,
        discard_untrimmed=bool(args.cutadapt_discard_untrimmed),
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
        qiime_deblur_denoise_single(
            in_qza=deblur_in,
            table_qza=table_qza,
            repseqs_qza=repseqs_qza,
            trim_length=args.deblur_trim_length,
            threads=args.threads,
            logs=paths.logs,
        )


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
    shutil.copy2(fixed_meta, paths.phyloseq_output / "metadata.tsv") 
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
    run_qc_bundle(paths=paths, metadata_tsv=fixed_meta, denoiser=args.denoiser, logger=logger)
    logger.info("QC bundle written to: %s", paths.root / "qc")
    end_ts = time.time()
    logger.info("End time: %s", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_ts)))
    log_memory_usage(logger=logger, prefix="END", extra_msg="Pipeline complete")




if __name__ == "__main__":
    main()
