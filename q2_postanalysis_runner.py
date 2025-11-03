#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
q2_postanalysis_runner
======================

Standard post-analysis bundle for QIIME 2 outputs.

Given a directory containing:
  - table.qza
  - rep-seqs.qza
  - rooted-tree.qza
  - metadata.tsv
  - taxonomy.qza (optional)

This script will:

1) Pick or compute a rarefaction depth (or use --sampling_depth).
2) Run alpha/beta diversity, PCoA + Emperor plots.
3) Run alpha/beta group significance tests (if --group_column given).
4) Make taxonomic barplots at multiple levels (if taxonomy is present).
5) Optionally run ANCOM on a chosen metadata column.
6) (Optional) Export all .qzv visualizations to HTML and try PDF rendering.
7) Write a compact README.md and a summary.tsv with high-level QC stats.

Requires a QIIME 2 environment on PATH.

Example
-------
python q2_postanalysis_runner.py \
  --input_dir results/JH102/phyloseq_output \
  --out_dir results/JH102/post_analysis \
  --group_column Treatment \
  --sampling_depth 12000 \
  --run_ancom \
  --export_visuals \
  --threads 24
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import subprocess
from pathlib import Path
from statistics import median_low
import sys
from typing import Dict, Iterable, List, Optional, Tuple
from html import escape
from datetime import datetime
import pandas as pd 


# ----------------------------- utilities ----------------------------- #


def setup_logging(out_dir: Path) -> logging.Logger:
    """Configure structured logging to stderr and to a log file.

    Args:
        out_dir: Output directory where `postanalysis.log` will be written.

    Returns:
        A configured logger instance (`q2_postanalysis`).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = out_dir / "postanalysis.log"

    logger = logging.getLogger("q2_postanalysis")
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)

    stream = logging.StreamHandler(stream=sys.stderr)
    stream.setLevel(logging.INFO)
    stream.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))

    fileh = logging.FileHandler(str(log_file), mode="w", encoding="utf-8")
    fileh.setLevel(logging.DEBUG)
    fileh.setFormatter(
        logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    )

    logger.addHandler(stream)
    logger.addHandler(fileh)
    logger.info("Logging to %s", log_file)
    return logger


def run_optional(
    *,
    cmd: list[str],
    log: logging.Logger,
    logfile: Optional[Path] = None,
    warn: str = "Step failed; skipping.",
) -> bool:
    """
    Run a subprocess command; on failure, log a warning and continue.

    Parameters
    ----------
    cmd : list[str]
        Command tokens to execute.
    log : logging.Logger
        Logger used to emit messages.
    logfile : pathlib.Path, optional
        If provided, stdout/stderr are appended to this file.
    warn : str
        Warning message to emit on failure.

    Returns
    -------
    bool
        True if the command succeeded, False if it failed.
    """
    log.info("▶ %s", " ".join(cmd))
    try:
        if logfile:
            logfile.parent.mkdir(parents=True, exist_ok=True)
            with logfile.open("a", encoding="utf-8") as lf:
                lf.write("$ " + " ".join(cmd) + "\n")
                lf.flush()
                subprocess.run(cmd, stdout=lf, stderr=lf, check=True)
        else:
            subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as exc:
        log.warning("%s (exit=%s)", warn, exc.returncode)
        return False


def run(cmd: List[str], log: logging.Logger, logfile: Optional[Path] = None) -> None:
    """Run a subprocess command and optionally tee to a step log.

    Args:
        cmd: Command tokens to execute (no shell=True).
        log: Logger used to emit the starting message.
        logfile: If provided, stdout/stderr are appended to this file.

    Raises:
        subprocess.CalledProcessError: If the command exits non-zero.
    """
    log.info("▶ %s", " ".join(cmd))
    if logfile:
        logfile.parent.mkdir(parents=True, exist_ok=True)
        with logfile.open("a", encoding="utf-8") as lf:
            lf.write("$ " + " ".join(cmd) + "\n")
            lf.flush()
            subprocess.run(cmd, stdout=lf, stderr=lf, check=True)
    else:
        subprocess.run(cmd, check=True)



from typing import Optional, Sequence
from pathlib import Path
import subprocess
import logging


def _find_sample_frequency_file(
    *,
    export_dir: Path,
    log: logging.Logger,
) -> Optional[Path]:
    """
    Locate the per-sample frequency table inside an exported QIIME visualisation.

    QIIME often places the table under ``<export_dir>/data/sample-frequency.tsv``.
    Some templates omit TSVs entirely; in that case this function returns ``None``.

    Parameters
    ----------
    export_dir : pathlib.Path
        Directory produced by ``qiime tools export`` for the feature-table summary.
    log : logging.Logger
        Logger used for debug/progress messages.

    Returns
    -------
    pathlib.Path or None
        Path to a plausible frequency table if found, otherwise ``None``.
    """
    candidates: Sequence[Path] = (
        export_dir / "data" / "sample-frequency.tsv",
        export_dir / "sample-frequency.tsv",
        export_dir / "data" / "sample-frequencies.tsv",
        export_dir / "data" / "sample-frequency-detail.tsv",
        export_dir / "data" / "sample-frequency-detail.csv",
    )
    for path in candidates:
        if path.exists():
            return path

    data_dir = export_dir / "data"
    if data_dir.exists():
        for path in sorted(list(data_dir.glob("*.tsv")) + list(data_dir.glob("*.csv"))):
            try:

                pd.read_csv(
                    path,
                    sep="\t" if path.suffix.lower() == ".tsv" else ",",
                    nrows=5,
                    dtype=str,
                )
                log.debug("Considering %s as a sample-frequency candidate.", path)
                return path
            except Exception:  # nosec B902 - broad on purpose, we just skip candidates
                continue
    return None


def _load_depths_from_export(*, sample_freq_tsv: Path) -> pd.DataFrame:
    """Load per-sample read depths from QIIME's sample-frequency.tsv.

    The export produced by `qiime feature-table summarize` typically includes a
    TSV with two columns: Sample ID and Frequency. Column names can vary
    slightly across QIIME versions, so this function matches them defensively.

    Args:
        sample_freq_tsv: Path to the exported `sample-frequency.tsv`.

    Returns:
        DataFrame with columns: 'sample-id' (str) and 'depth' (int).

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If required columns cannot be identified.
    """

    if not sample_freq_tsv.exists():
        raise FileNotFoundError(f"Missing {sample_freq_tsv}")

    df = pd.read_csv(sample_freq_tsv, sep="\t", dtype=str)
    canon = {c.strip().lower(): c for c in df.columns}

    sid_col = next((canon[k] for k in canon if "sample" in k and "id" in k), None)
    freq_col = next((canon[k] for k in canon if "frequency" in k), None)
    if sid_col is None or freq_col is None:
        raise ValueError("Could not find 'Sample ID' and 'Frequency' columns in sample-frequency.tsv")

    out = df[[sid_col, freq_col]].copy()
    out.columns = ["sample-id", "depth"]
    out["depth"] = pd.to_numeric(out["depth"], errors="coerce").fillna(0).astype(int)
    return out


def _load_metadata_for_depth(*, metadata_tsv: Path) -> pd.DataFrame:
    """Load metadata TSV and ensure the first column is 'sample-id'.

    Args:
        metadata_tsv: Path to the analysis metadata TSV.

    Returns:
        DataFrame containing at least 'sample-id'.
    """

    md = pd.read_csv(metadata_tsv, sep="\t", dtype=str, comment="#")
    first = md.columns[0]
    if first.lower().replace("_", "-") != "sample-id":
        md = md.rename(columns={first: "sample-id"})
    return md


def _retention_curve(*, depths: "pd.Series", candidates: "List[int]") -> "pd.DataFrame":
    """Compute overall retention at each candidate depth.

    Args:
        depths: Per-sample read depths (ints).
        candidates: Candidate rarefaction depths.

    Returns:
        DataFrame with columns: depth, n_kept, frac_kept.
    """
    n_total = int((depths > 0).sum())
    rows = []
    for d in candidates:
        kept = int((depths >= d).sum())
        frac = (kept / n_total) if n_total else 0.0
        rows.append({"depth": int(d), "n_kept": kept, "frac_kept": frac})
    return pd.DataFrame(rows)


def _retention_by_group(
    *, df: "pd.DataFrame", group_column: str, candidates: "List[int]"
) -> "pd.DataFrame":
    """Compute per-group retention at each candidate depth.

    Args:
        df: DataFrame with columns: sample-id, depth, and the group column.
        group_column: Metadata column to balance (e.g., 'Genotype').
        candidates: Candidate depths.

    Returns:
        DataFrame with columns: depth, group, n_group, n_kept, frac_kept.
    """
    rows = []
    for g, sub in df.groupby(group_column, dropna=False):
        sub_total = int((sub["depth"] > 0).sum())
        for d in candidates:
            kept = int((sub["depth"] >= d).sum())
            frac = (kept / sub_total) if sub_total else 0.0
            rows.append(
                {
                    "depth": int(d),
                    "group": str(g),
                    "n_group": sub_total,
                    "n_kept": kept,
                    "frac_kept": frac,
                }
            )
    return pd.DataFrame(rows)


def _choose_depth_balanced(
    *,
    overall: "pd.DataFrame",
    by_group: "pd.DataFrame",
    min_overall: float,
    min_each_group: float,
) -> int:
    """Select the highest depth meeting overall AND per-group retention thresholds.

    Falls back to highest depth meeting overall only; if none, returns the
    minimum candidate depth.

    Args:
        overall: Overall retention table (depth, frac_kept).
        by_group: Per-group retention table (depth, group, frac_kept).
        min_overall: Minimum overall retention fraction (e.g., 0.85).
        min_each_group: Minimum per-group retention fraction (e.g., 0.75).

    Returns:
        Recommended integer depth.
    """
    valid = []
    for d, sub_o in overall.groupby("depth"):
        ok_overall = float(sub_o["frac_kept"].iloc[0]) >= min_overall
        sub_g = by_group[by_group["depth"] == d]
        ok_groups = sub_g["frac_kept"].min() >= min_each_group if not sub_g.empty else False
        if ok_overall and ok_groups:
            valid.append(int(d))

    if valid:
        return max(valid)

    # Fallback: overall only
    candidates = overall[overall["frac_kept"] >= min_overall]["depth"].astype(int).tolist()
    if candidates:
        return max(candidates)

    # Final fallback: smallest candidate
    return int(overall["depth"].min())


def recommend_rarefaction_depth(
    *,
    summary_export_dir: Path,
    metadata_tsv: Path,
    group_column: Optional[str],
    out_dir: Path,
    min_overall: float,
    min_each_group: float,
    log: logging.Logger,
) -> Optional[int]:
    """Recommend a rarefaction depth balancing overall and per-group retention.

    This function reads the exported `sample-frequency.tsv`, merges with the
    metadata, evaluates retention across a sensible candidate grid, writes
    supporting TSVs, and returns the recommended depth.

    If any step is unavailable (e.g., missing files or no group column), it
    returns None and the caller may fall back to a simpler heuristic.

    Args:
        summary_export_dir: Directory created from exporting feature-table summary.
        metadata_tsv: Path to run metadata TSV (first column is sample-id).
        group_column: Column name to balance. If None or not present, returns None.
        out_dir: Folder to write `sample_depths_with_groups.tsv`, `retention_overall.tsv`,
            `retention_by_group.tsv`, and `recommendation.tsv`.
        min_overall: Minimum overall retention fraction required.
        min_each_group: Minimum per-group retention fraction required.
        log: Logger for progress messages.

    Returns:
        Integer depth recommendation, or None if not derivable.
    """
    try:
        # 1) Prefer a QZV-exported frequency file if present
        freq_file = _find_sample_frequency_file(export_dir=summary_export_dir, log=log)
        if freq_file:
            depths = _load_depths_from_export(sample_freq_tsv=freq_file)

            log.info("Loaded per-sample depths from %s", freq_file)
        else:
            # 2) Fallback: compute depths from feature table directly (BIOM → TSV)
            log.info("No sample-frequency.tsv in QZV export; deriving depths from table.qza instead.")
            # table.qza is next to metadata in the input dir; we pass it via caller in main()
            # We cannot see it here, so we’ll re-derive its path in main() and pass via out closure.
            # Easiest: infer `table.qza` location from metadata path (both in the same folder).
            table_qza = metadata_tsv.parent / "table.qza"
            depths = _load_depths_fallback_from_table(table_qza=table_qza, out_dir=out_dir.parent, log=log)
            if depths is None:
                log.warning("Depth recommendation: could not derive depths from table.qza fallback.")
                return None

        meta = _load_metadata_for_depth(metadata_tsv=metadata_tsv)

        depths["sample-id"] = (
                    depths["sample-id"]
                    .astype(str)
                    .str.strip()
                    )

        meta["sample-id"] = (
                    meta["sample-id"]
                    .astype(str)
                    .str.strip()
                )



        if not group_column or group_column not in meta.columns:
            log.info(
                "Depth recommendation: group column not provided/found; "
                "falling back to percentile-based inference."
            )
            return None

        df = depths.merge(meta, on="sample-id", how="inner")
        n_depth = int(depths["sample-id"].nunique())
        n_meta = int(meta["sample-id"].nunique())
        n_overlap = int(df["sample-id"].nunique())
        log.info(
            "Depths×Metadata sample IDs: %d in table, %d in metadata, %d overlapping.",
            n_depth, n_meta, n_overlap,
        )
        if n_overlap == 0:
            # Optional: show a few examples to debug naming/whitespace issues
            ex_depth = list(depths["sample-id"].astype(str).head(5))
            ex_meta = list(meta["sample-id"].astype(str).head(5))
            log.warning(
                "No overlap between depths and metadata sample IDs. "
                "Examples (depths): %s | (metadata): %s",
                ex_depth, ex_meta,
            )

        if df.empty:
            log.warning("Depth recommendation: no overlapping sample IDs between depths and metadata.")
            return None


        # Candidate grid: unique depths or quantile grid if very large.
        uniq = sorted(set(df["depth"].tolist()))
        if len(uniq) > 200:
            qs = list(range(5, 100, 5))
            qs_vals = [int(pd.Series(df["depth"]).quantile(q / 100.0)) for q in qs]
            candidates = sorted(set([x for x in qs_vals if x > 0] + [min(uniq), max(uniq)]))
        else:
            candidates = [int(x) for x in uniq if x > 0]

        overall = _retention_curve(depths=df["depth"], candidates=candidates)
        by_group = _retention_by_group(df=df, group_column=group_column, candidates=candidates)

        out_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_dir / "sample_depths_with_groups.tsv", sep="\t", index=False)
        overall.to_csv(out_dir / "retention_overall.tsv", sep="\t", index=False)
        by_group.to_csv(out_dir / "retention_by_group.tsv", sep="\t", index=False)

        rec = _choose_depth_balanced(
            overall=overall,
            by_group=by_group,
            min_overall=min_overall,
            min_each_group=min_each_group,
        )
        pd.DataFrame(
            [
                {
                    "recommended_depth": rec,
                    "min_overall": min_overall,
                    "min_each_group": min_each_group,
                }
            ]
        ).to_csv(out_dir / "recommendation.tsv", sep="\t", index=False)

        log.info(
            "Depth recommendation: %d (overall ≥ %.2f, each %s ≥ %.2f).",
            rec,
            min_overall,
            group_column,
            min_each_group,
        )
        return int(rec)

    except Exception as e:
        log.warning("Depth recommendation failed (%s); falling back to percentile heuristic.", e)
        return None



def write_html_index(out_dir: Path) -> None:
    """Create a minimal landing page linking to QZVs and exported HTML/PDF.

    The page lists:
      - All QZV visualizations in ``visuals/`` with direct links.
      - All exported QZV folders in ``exports/`` (if present), linking to their
        ``index.html`` and an accompanying PDF when available.

    Args:
        out_dir: Root post-analysis directory (contains ``visuals/`` and optionally
            ``exports/``).
    """
    visuals = out_dir / "visuals"
    exports = out_dir / "exports"
    items_qzv = sorted(visuals.glob("*.qzv")) if visuals.exists() else []
    export_dirs = sorted([p for p in exports.glob("*") if p.is_dir()]) if exports.exists() else []

    def rel(p: Path) -> str:
        return escape(str(p.relative_to(out_dir)).replace("\\", "/"))

    # Simple, readable HTML with tiny bit of CSS for clarity.
    lines = [
        "<!doctype html>",
        "<html lang='en'>",
        "<head>",
        "  <meta charset='utf-8' />",
        f"  <title>QIIME 2 Post-analysis — {escape(out_dir.name)}</title>",
        "  <meta name='viewport' content='width=device-width,initial-scale=1' />",
        "  <style>",
        "    body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Ubuntu,Cantarell,'Helvetica Neue',Arial,sans-serif;line-height:1.45;margin:2rem;}",
        "    h1{font-size:1.6rem;margin:0 0 0.5rem}",
        "    h2{font-size:1.2rem;margin:2rem 0 0.5rem}",
        "    .note{color:#555;margin:0.25rem 0 1rem}",
        "    ul{padding-left:1.2rem}",
        "    li{margin:0.25rem 0}",
        "    code{background:#f6f8fa;padding:0.1rem 0.3rem;border-radius:4px}",
        "    .muted{color:#888}",
        "  </style>",
        "</head>",
        "<body>",
        f"  <h1>QIIME 2 post-analysis: <code>{escape(out_dir.name)}</code></h1>",
        "  <p class='note'>This index links to the primary visualisations and their exported HTML/PDF (when available).</p>",
    ]

    # QZV visualisations
    lines.append("  <h2>Visualizations (.qzv)</h2>")
    if items_qzv:
        lines.append("  <ul>")
        for q in items_qzv:
            lines.append(f"    <li><a href='{rel(q)}'>{escape(q.name)}</a></li>")
        lines.append("  </ul>")
    else:
        lines.append("  <p class='muted'>No .qzv files found in <code>visuals/</code>.</p>")

    # Exports with HTML/PDF
    lines.append("  <h2>Exports (HTML / PDF)</h2>")
    if export_dirs:
        lines.append("  <ul>")
        for d in export_dirs:
            idx = d / "index.html"
            pdf = d.with_suffix(".pdf")
            parts = [escape(d.name)]
            if idx.exists():
                parts.append(f"<a href='{rel(idx)}'>HTML</a>")
            else:
                parts.append("<span class='muted'>HTML n/a</span>")
            if pdf.exists():
                parts.append(f"<a href='{rel(pdf)}'>PDF</a>")
            else:
                parts.append("<span class='muted'>PDF n/a</span>")
            lines.append(f"    <li>{' — '.join(parts)}</li>")
        lines.append("  </ul>")
    else:
        lines.append("  <p class='muted'>No exports found in <code>exports/</code>.</p>")

    lines += [
        "  <hr />",
        "  <p class='muted'>Generated automatically by q2_postanalysis_runner.py</p>",
        "</body>",
        "</html>",
    ]

    (out_dir / "index.html").write_text("\n".join(lines), encoding="utf-8")



def ensure_files(base: Path, names: Iterable[str]) -> None:
    """Ensure a set of files exists in a directory.

    Args:
        base: Parent directory to check.
        names: Filenames expected to exist under `base`.

    Raises:
        FileNotFoundError: If any of the required files are missing.
    """
    missing = [n for n in names if not (base / n).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required files: {', '.join(missing)}")


def export_qzv(qzv: Path, outdir: Path, log: logging.Logger) -> Path:
    """Export a QIIME 2 visualization (.qzv) to an HTML folder.

    Args:
        qzv: Path to the .qzv file.
        outdir: Directory to create and place exported assets in.
        log: Logger for progress messages.

    Returns:
        The path to the created export directory.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    run(
        ["qiime", "tools", "export", "--input-path", str(qzv), "--output-path", str(outdir)],
        log,
    )
    return outdir


def try_html_to_pdf(export_dir: Path, log: logging.Logger, engine: str = "auto") -> None:
    """Render an exported QZV (index.html) to PDF with the requested engine.

    Args:
        export_dir: Directory produced by `qiime tools export` (has index.html).
        log: Logger for progress messages.
        engine: One of {'auto','chrome','wkhtmltopdf','none'}.
    """
    index_html = export_dir / "index.html"
    if not index_html.exists():
        log.debug("No index.html in %s; skipping PDF render.", export_dir)
        return

    pdf_path = export_dir.with_suffix(".pdf")

    def try_chrome() -> bool:
        for browser in ("chromium", "chromium-browser", "google-chrome", "google-chrome-stable", "chrome"):
            exe = shutil.which(browser)
            if not exe:
                continue
            log.info("Rendering PDF with Chrome/Chromium: %s", browser)
            cmd = [
                exe,
                "--headless",
                "--disable-gpu",
                f"--print-to-pdf={pdf_path}",
                str(index_html),
            ]
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if pdf_path.exists():
                    log.info("PDF written: %s", pdf_path)
                    return True
            except subprocess.CalledProcessError as e:
                log.debug("Chrome render failed (%s): %s", browser, e)
                continue
        return False

    def try_wkhtmltopdf() -> bool:
        exe = shutil.which("wkhtmltopdf")
        if not exe:
            return False
        log.info("Rendering PDF with wkhtmltopdf")
        # Note: wkhtmltopdf has limited JS support; Vega/Emperor may be incomplete.
        cmd = [
            exe,
            "--quiet",
            "--enable-local-file-access",
            "--print-media-type",
            "--page-size", "A4",
            str(index_html),
            str(pdf_path),
        ]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if pdf_path.exists():
                log.info("PDF written: %s", pdf_path)
                return True
        except subprocess.CalledProcessError as e:
            log.debug("wkhtmltopdf render failed: %s", e)
        return False

    if engine == "none":
        log.info("PDF rendering disabled (--pdf_engine none).")
        return
    elif engine == "chrome":
        ok = try_chrome()
        if not ok:
            log.info("Chrome/Chromium not available or failed; skipping PDF for %s", export_dir)
        return
    elif engine == "wkhtmltopdf":
        ok = try_wkhtmltopdf()
        if not ok:
            log.info("wkhtmltopdf not available or failed; skipping PDF for %s", export_dir)
        return
    else:  # auto
        if try_chrome():
            return
        if try_wkhtmltopdf():
            return
        log.info("No PDF engine available; skipping PDF for %s", export_dir)


def infer_sampling_depth_from_summary(summary_export_dir: Path) -> Optional[int]:
    """Infer a conservative rarefaction depth from a feature-table summary export.

    The function parses `sample-frequency.tsv` and returns the ~15th percentile
    read depth (rounded down), which is a reasonable default for many cohorts.

    Args:
        summary_export_dir: Directory created by exporting `feature-table summarize`.

    Returns:
        An integer rarefaction depth if derivable, otherwise `None`.
    """
    freq_tsv = summary_export_dir / "sample-frequency.tsv"
    if not freq_tsv.exists():
        return None

    depths: List[int] = []
    with freq_tsv.open("r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        _ = next(reader, None)  # header
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            try:
                depths.append(int(row[1]))
            except Exception:
                continue

    if not depths:
        return None

    depths.sort()
    idx = max(0, int(len(depths) * 0.15) - 1)
    return depths[idx]


def summarize_depths(summary_export_dir: Path) -> Dict[str, int]:
    """Summarize basic depth statistics from `sample-frequency.tsv`.

    Args:
        summary_export_dir: Directory created by exporting `feature-table summarize`.

    Returns:
        A dictionary with keys: n_samples, min_depth, median_depth, max_depth.
        If unavailable, returns zeros for all fields.
    """
    freq_tsv = summary_export_dir / "sample-frequency.tsv"
    if not freq_tsv.exists():
        return {"n_samples": 0, "min_depth": 0, "median_depth": 0, "max_depth": 0}

    depths: List[int] = []
    with freq_tsv.open("r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        _ = next(reader, None)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            try:
                depths.append(int(row[1]))
            except Exception:
                continue

    if not depths:
        return {"n_samples": 0, "min_depth": 0, "median_depth": 0, "max_depth": 0}

    depths.sort()
    return {
        "n_samples": len(depths),
        "min_depth": depths[0],
        "median_depth": int(median_low(depths)),
        "max_depth": depths[-1],
    }


def _export_table_biom(*, table_qza: Path, export_dir: Path, log: logging.Logger) -> Path:
    """Export a QIIME 2 feature table artefact to BIOM format.

    Args:
        table_qza: Path to table.qza.
        export_dir: Directory to write the exported BIOM.
        log: Logger for progress messages.

    Returns:
        Path to the exported BIOM file (feature-table.biom).
    """
    export_dir.mkdir(parents=True, exist_ok=True)
    biom_path = export_dir / "feature-table.biom"
    if not biom_path.exists():
        run(
            ["qiime", "tools", "export", "--input-path", str(table_qza), "--output-path", str(export_dir)],
            log,
        )
    return biom_path


def _collapsed_table_stats(*, table_qza: Path, work_dir: Path, log: logging.Logger) -> dict:
    """
    Export a collapsed table.qza -> BIOM -> TSV and return basic stats:
    - n_features (rows)
    - n_nonzero_features (rows with any nonzero counts)
    - total_frequency (sum of all counts)
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    # Export to BIOM (reuses your helper)
    biom_path = _export_table_biom(table_qza=table_qza, export_dir=work_dir, log=log)
    tsv_path = work_dir / "feature-table.tsv"

    # Convert BIOM -> TSV
    subprocess.run(
        ["biom", "convert", "--to-tsv", "--input-fp", str(biom_path), "--output-fp", str(tsv_path)],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )

    # Detect header line (“#OTU ID …” often sits on line 2)
    with tsv_path.open("r", encoding="utf-8") as fh:
        _ = fh.readline()
        maybe_header = fh.readline()

    header_line_index = 1 if (maybe_header.startswith("#OTU ID") or maybe_header.startswith("#OTUID")) else 0
    df = pd.read_csv(tsv_path, sep="\t", header=header_line_index, dtype=str, engine="python")

    if df.shape[1] < 2:
        return {"n_features": 0, "n_nonzero_features": 0, "total_frequency": 0}

    # samples are columns 1..N
    numeric = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0)
    row_sums = numeric.sum(axis=1)
    n_features = int(df.shape[0])
    n_nonzero = int((row_sums > 0).sum())
    total_freq = int(numeric.values.sum())
    return {"n_features": n_features, "n_nonzero_features": n_nonzero, "total_frequency": total_freq}


def write_taxa_summary_tsv(*, out_dir: Path, art_dir: Path, tax_levels: List[int], log: logging.Logger) -> None:
    """
    For each collapsed level table_L{L}.qza in `art_dir`, compute stats and write
    out_dir/taxa_summary.tsv with columns: level, n_features, n_nonzero_features, total_frequency.
    """
    rows = []
    for level in sorted(tax_levels):
        qza = art_dir / f"table_L{level}.qza"
        if not qza.exists():
            log.info("No collapsed table for level %s; skipping in taxa_summary.", level)
            continue
        stats = _collapsed_table_stats(table_qza=qza, work_dir=art_dir / f"table_L{level}_export", log=log)
        rows.append({
            "level": level,
            "n_features": stats["n_features"],
            "n_nonzero_features": stats["n_nonzero_features"],
            "total_frequency": stats["total_frequency"],
        })

    out_fp = out_dir / "taxa_summary.tsv"
    with out_fp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["level", "n_features", "n_nonzero_features", "total_frequency"])
        for r in rows:
            w.writerow([r["level"], r["n_features"], r["n_nonzero_features"], r["total_frequency"]])
    log.info("Wrote %s", out_fp)



def _group_counts_after_rarefy(
    *,
    rarefied_table: Path,
    metadata_tsv: Path,
    group_column: str,
) -> dict[str, int]:
    """
    Estimate per-group sample counts among samples retained after rarefaction.

    Returns a dict mapping group -> count. If detection fails, returns {}.
    """
    try:
        # Export BIOM
        tmp = rarefied_table.parent / "tmp_counts"
        tmp.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["qiime", "tools", "export", "--input-path", str(rarefied_table), "--output-path", str(tmp)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        biom_fp = tmp / "feature-table.biom"
        tsv_fp = tmp / "feature-table.tsv"
        subprocess.run(
            ["biom", "convert", "--to-tsv", "--input-fp", str(biom_fp), "--output-fp", str(tsv_fp)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        # Second line header lists sample IDs
        with tsv_fp.open("r", encoding="utf-8") as fh:
            _ = fh.readline()
            header = fh.readline().strip().split("\t")[1:]

        md = pd.read_csv(metadata_tsv, sep="\t", dtype=str, comment="#")
        md.rename(columns={md.columns[0]: "sample-id"}, inplace=True)
        md["sample-id"] = md["sample-id"].astype(str).str.strip()
        kept = md[md["sample-id"].isin(header)]
        return kept[group_column].value_counts(dropna=False).to_dict()
    except Exception:
        return {}


def _biom_to_depths_via_tsv(
    *,
    biom_path: Path,
    work_dir: Path,
) -> "pd.DataFrame":
    """
    Convert a BIOM matrix to TSV and compute per-sample read depths.

    This preserves the BIOM TSV header (which often starts with '#OTU ID')
    so that sample IDs are parsed correctly.

    Parameters
    ----------
    biom_path : pathlib.Path
        Path to ``feature-table.biom``.
    work_dir : pathlib.Path
        Working directory to store the temporary TSV output.

    Returns
    -------
    pandas.DataFrame
        Data frame with columns ``sample-id`` (str) and ``depth`` (int).

    Raises
    ------
    ValueError
        If the converted TSV has an unexpected layout.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = work_dir / "feature-table.tsv"

    cmd = [
        "biom",
        "convert",
        "--to-tsv",
        "--input-fp",
        str(biom_path),
        "--output-fp",
        str(tsv_path),
    ]
    subprocess.run(  # nosec B603
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Detect whether the second line is the header ('#OTU ID ...')
    with tsv_path.open("r", encoding="utf-8") as fh:
        first = fh.readline()
        second = fh.readline()

    header_line_index = 0
    if second.startswith("#OTU ID") or second.startswith("#OTUID"):
        # Use the second line as header; pandas will strip the leading '#'
        header_line_index = 1

    df = pd.read_csv(
        tsv_path,
        sep="\t",
        header=header_line_index,
        dtype=str,
        engine="python",
    )

    # Normalise first column name to something like 'Feature ID'
    first_col = df.columns[0]
    sample_cols = [c for c in df.columns[1:]]
    if not sample_cols:
        raise ValueError(f"Unexpected BIOM TSV layout in {tsv_path}")

    # Convert counts to numeric and sum per sample
    numeric = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    depths = numeric.sum(axis=0).astype(int).reset_index()
    depths.columns = ["sample-id", "depth"]
    return depths



def _load_depths_fallback_from_table(
    *, table_qza: Path, out_dir: Path, log: logging.Logger
) -> Optional["pd.DataFrame"]:
    """Fallback path to obtain per-sample depths from table.qza when QZV export lacks TSVs.

    Args:
        table_qza: Path to table.qza.
        out_dir: Output directory where temporary/exported files will live.
        log: Logger for progress messages.

    Returns:
        DataFrame with columns 'sample-id' and 'depth', or None on failure.
    """
    try:
        biom_dir = out_dir / "exports" / "table_export"
        biom_path = _export_table_biom(table_qza=table_qza, export_dir=biom_dir, log=log)
        work_dir = out_dir / "exports" / "table_export" / "tmp"
        depths = _biom_to_depths_via_tsv(biom_path=biom_path, work_dir=work_dir)
        log.info("Derived per-sample depths from feature table (BIOM → TSV).")
        return depths
    except Exception as e:
        log.warning("Fallback depth derivation from table.qza failed: %s", e)
        return None


def _alpha_vectors_for_metric(
    *,
    metric: str,
    core_dir: Path,
    art_dir: Path,
    rarefied_table: Path,
    table_raw: Path,
    log: logging.Logger,
    logs_dir: Path,
) -> Path:
    """
    Resolve (or compute) the alpha-diversity vector for a given metric.

    For metrics produced by core-metrics-phylogenetic we prefer the rarefied
    outputs. Only if a metric is not present do we compute it; in that case,
    we compute from the rarefied table to keep inputs consistent.

    Parameters
    ----------
    metric : str
        Alpha-diversity metric name (e.g., 'observed_features', 'shannon',
        'faith_pd', 'pielou_e').
    core_dir : pathlib.Path
        Directory where core-metrics outputs for this depth were written.
    art_dir : pathlib.Path
        Directory for ad-hoc artefacts.
    rarefied_table : pathlib.Path
        Rarefied table produced by core-metrics.
    table_raw : pathlib.Path
        Original (unrarefied) table (kept for reference).
    log : logging.Logger
        Logger for progress messages.
    logs_dir : pathlib.Path
        Directory for per-step logs.

    Returns
    -------
    pathlib.Path
        Path to the alpha-diversity vector .qza for this metric.
    """
    # Mapping from requested metric to core-metrics vector file name
    metric_to_core = {
        "observed_features": "observed_features_vector.qza",
        "shannon": "shannon_vector.qza",
        "faith_pd": "faith_pd_vector.qza",
        "pielou_e": "evenness_vector.qza",  # Pielou == evenness from core-metrics
    }

    core_name = metric_to_core.get(metric)
    if core_name:
        core_vec = core_dir / core_name
        if core_vec.exists():
            return core_vec

    # Fallback: compute the metric from the rarefied table for consistency
    out_vec = art_dir / f"{metric}_vector.qza"
    run(
        [
            "qiime", "diversity", "alpha",
            "--i-table", str(rarefied_table),
            "--p-metric", metric,
            "--o-alpha-diversity", str(out_vec),
        ],
        log,
        logs_dir / f"alpha_{metric}.log",
    )
    return out_vec



def _is_phylogenetic_beta(metric: str) -> bool:
    """Return True if a beta metric is phylogenetic (needs a tree).

    Args:
        metric: Beta-diversity metric name.

    Returns:
        True for UniFrac variants; False otherwise.
    """
    m = metric.strip().lower()
    return m in {"unweighted_unifrac", "weighted_unifrac",
                 "generalised_unifrac", "generalized_unifrac"}



def write_readme(
    out_dir: Path,
    sampling_depth: int,
    group_column: Optional[str],
    alpha_metrics: List[str],
    beta_metrics: List[str],
    tax_levels: List[int],
    taxonomy_present: bool,
) -> None:
    """Write a concise README.md describing the post-analysis bundle.

    Args:
        out_dir: Root output directory for post-analysis.
        sampling_depth: Rarefaction depth used for diversity analyses.
        group_column: Metadata column used for group significance tests (if any).
        alpha_metrics: Alpha metrics computed.
        beta_metrics: Beta metrics computed.
        tax_levels: Taxonomic levels collapsed for barplots.
        taxonomy_present: Whether taxonomy.qza was present.
    """
    lines = [
        "# QIIME 2 Post-analysis Bundle",
        "",
        f"- **Sampling depth:** `{sampling_depth}`",
        f"- **Group column:** `{group_column or 'NA'}`",
        f"- **Alpha metrics:** {', '.join(alpha_metrics)}",
        f"- **Beta metrics:** {', '.join(beta_metrics)}",
        f"- **Taxonomy present:** {taxonomy_present}",
        f"- **Taxa barplot levels:** {', '.join(map(str, tax_levels)) if taxonomy_present else 'n/a'}",
        "",
        "## Key folders",
        "- `artifacts/` – intermediate QIIME 2 artifacts (distance matrices, PCoAs, collapsed tables, etc.)",
        "- `visuals/` – QZV visualizations (and exported HTML/PDF if `--export_visuals` was used)",
        "- `logs/` – per-step QIIME command logs",
        "- `exports/` – exports of select QZVs for inspection/archiving (optional)",
        "",
        "## What was run",
        "1. Feature table and representative sequences summarization",
        "2. Core phylogenetic diversity metrics (alpha/beta) with PCoA + Emperor",
        "3. Alpha and beta group significance (if group column provided)",
        "4. Taxonomy collapse and barplots (if taxonomy present)",
        "5. ANCOM differential abundance (if requested)",
        "",
        "This bundle is designed for quick QC and first-pass exploration.",
        "For custom analyses, consider downstream R (phyloseq, vegan) or Python workflows.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def write_summary_tsv(
    out_dir: Path,
    depth_stats: Dict[str, int],
    sampling_depth: int,
    taxonomy_present: bool,
    alpha_metrics: List[str],
    beta_metrics: List[str],
    tax_levels: List[int],
    group_column: Optional[str],
) -> None:
    """Write a summary TSV with high-level QC stats and configuration.

    Args:
        out_dir: Root output directory where `summary.tsv` will be written.
        depth_stats: Dictionary with n_samples, min_depth, median_depth, max_depth.
        sampling_depth: Rarefaction depth used.
        taxonomy_present: Whether taxonomy.qza was present.
        alpha_metrics: Alpha metrics computed.
        beta_metrics: Beta metrics computed.
        tax_levels: Taxonomy levels used for barplots (if any).
        group_column: Grouping column name used for significance tests (if any).
    """
    rows = [
        ("n_samples", depth_stats.get("n_samples", 0)),
        ("min_depth", depth_stats.get("min_depth", 0)),
        ("median_depth", depth_stats.get("median_depth", 0)),
        ("max_depth", depth_stats.get("max_depth", 0)),
        ("sampling_depth", sampling_depth),
        ("taxonomy_present", int(bool(taxonomy_present))),
        ("alpha_metrics", ",".join(alpha_metrics)),
        ("beta_metrics", ",".join(beta_metrics)),
        ("tax_levels", ",".join(map(str, tax_levels)) if taxonomy_present else "NA"),
        ("group_column", group_column or "NA"),
    ]
    with (out_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["key", "value"])
        writer.writerows(rows)


# ----------------------------- main ----------------------------- #


def main() -> None:
    """CLI entrypoint for running the QIIME 2 post-analysis bundle."""
    parser = argparse.ArgumentParser(
        description="Run a standard QIIME 2 post-analysis bundle on denoised outputs."
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="Folder containing table.qza, rep-seqs.qza, rooted-tree.qza, metadata.tsv (and optional taxonomy.qza).",
    )
    parser.add_argument(
        "--out_dir",
        required=False,
        type=Path,
        help="Output directory for post-analysis (default: <input_dir>/post_analysis).",
    )
    parser.add_argument(
        "--sampling_depth",
        type=int,
        default=None,
        help="Rarefaction depth. If not provided, it is estimated from feature-table summarize.",
    )
    parser.add_argument(
        "--group_column",
        type=str,
        default=None,
        help="Metadata column for group significance tests (PERMANOVA, alpha group significance, ANCOM).",
    )
    parser.add_argument(
        "--alpha_metrics",
        type=str,
        nargs="+",
        default=["observed_features", "shannon", "faith_pd", "pielou_e"],
        help="Alpha diversity metrics to compute.",
    )
    parser.add_argument(
        "--beta_metrics",
        type=str,
        nargs="+",
        default=["braycurtis", "jaccard", "unweighted_unifrac", "weighted_unifrac"],
        help="Beta diversity metrics to compute.",
    )
    parser.add_argument(
        "--tax_levels",
        type=int,
        nargs="+",
        default=[2, 3, 4, 5, 6],
        help="Taxonomy levels for collapsed barplots.",
    )

    parser.add_argument(
        "--no_ancom",
        dest="run_ancom",
        action="store_false",
        help="Skip ANCOM analysis (it now runs by default).",
    )
    parser.set_defaults(run_ancom=True)

    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads for heavy steps.",
    )

    parser.add_argument(
        "--min_depth_floor",
        type=int,
        default=5,
        help="Minimum rarefaction depth floor; if the recommended depth is lower, "
            "this floor will be used instead (default: 5).",
    )

    parser.add_argument(
        "--min_overall",
        type=float,
        default=0.85,
        help="Minimum overall retention fraction required for depth recommendation (default: 0.85).",
    )
    parser.add_argument(
        "--min_each_group",
        type=float,
        default=0.75,
        help="Minimum per-group retention fraction required for depth recommendation (default: 0.75).",
    )
    parser.add_argument(
        "--depth_report_dir",
        type=Path,
        default=None,
        help="Optional subfolder to write depth recommendation TSVs (default: <out_dir>/depth_recommendation).",
    )


    parser.add_argument(
            "--pdf_engine",
            type=str,
            choices=["auto", "chrome", "wkhtmltopdf", "none"],
            default="auto",
            help=(
                "How to render PDFs from exported QZVs. "
                "'chrome' uses headless Chromium/Chrome if found; "
                "'wkhtmltopdf' uses wkhtmltopdf; "
                "'auto' tries Chrome first then wkhtmltopdf; "
                "'none' disables PDF rendering."
            ),
        )

    parser.add_argument(
        "--export_visuals",
        action="store_true",
        help="Export .qzv visuals to HTML folders; try to render PDFs if Chrome is present.",
    )

    args = parser.parse_args()
    base = args.input_dir.resolve()
    out = (args.out_dir or (base / "post_analysis")).resolve()

    viz = out / "visuals"
    art = out / "artifacts"
    logs = out / "logs"
    for d in (viz, art, logs):
        d.mkdir(parents=True, exist_ok=True)

    log = setup_logging(out)
    log.info("INPUT_DIR=%s", base)
    log.info("OUT_DIR=%s", out)

    # Required inputs
    ensure_files(base, ["table.qza", "rep-seqs.qza", "rooted-tree.qza", "metadata.tsv"])
    table = base / "table.qza"
    repseqs = base / "rep-seqs.qza"
    tree = base / "rooted-tree.qza"
    metadata = base / "metadata.tsv"
    taxonomy_present = (base / "taxonomy.qza").exists()
    taxonomy = base / "taxonomy.qza" if taxonomy_present else None

    # 0) Basic summaries first (also helps choose depth)
    ft_summary_qzv = viz / "feature_table_summary.qzv"
    run(
        [
            "qiime",
            "feature-table",
            "summarize",
            "--i-table",
            str(table),
            "--o-visualization",
            str(ft_summary_qzv),
            "--m-sample-metadata-file",
            str(metadata),
        ],
        log,
        logs / "feature_table_summarize.log",
    )

    repseqs_tab_qzv = viz / "repseqs_tabulate_seqs.qzv"
    run(
        [
            "qiime",
            "feature-table",
            "tabulate-seqs",
            "--i-data",
            str(repseqs),
            "--o-visualization",
            str(repseqs_tab_qzv),
        ],
        log,
        logs / "tabulate_seqs.log",
    )

    # Export summary and compute depth recommendation (balanced by group if possible)
    exp_dir = out / "exports" / "feature_table_summary"
    export_qzv(ft_summary_qzv, exp_dir, log)
    depth_stats = summarize_depths(exp_dir)

    # Respect explicit user depth if provided
    depth = args.sampling_depth

    # Where to write the depth report (TSVs)
    depth_report_dir = args.depth_report_dir or (out / "depth_recommendation")

    if depth is None:
        # Try the balanced recommender first (requires valid group column)
        depth = recommend_rarefaction_depth(
            summary_export_dir=exp_dir,
            metadata_tsv=metadata,
            group_column=args.group_column,
            out_dir=depth_report_dir,
            min_overall=args.min_overall,
            min_each_group=args.min_each_group,
            log=log,
        )

    if depth is None:
        # Fall back to robust ~15th percentile heuristic
        depth = infer_sampling_depth_from_summary(exp_dir)
        if depth is not None:
            log.info("Depth (percentile heuristic): %d", depth)

    if depth is None:
        # Final safety default
        log.warning("Could not infer sampling depth; defaulting to 10000.")
        depth = 10000
    
    # If a floor is set, clamp small depths up to the floor (with a warning)
    if depth < args.min_depth_floor:
        log.warning(
            "Recommended depth %d is below the floor %d; using the floor instead. "
            "Consider dropping ultra-low samples or relaxing retention thresholds.",
            depth, args.min_depth_floor,
        )
        depth = args.min_depth_floor


    log.info("Using sampling depth: %d", depth)


    # 1) Core metrics (phylogenetic)

    core_dir = art / f"core_metrics_d{depth}"
    if core_dir.exists():
        ts = datetime.now().strftime("%Y%m%d-%H%M%S")
        alt = art / f"core_metrics_d{depth}_{ts}"
        log.warning(
            "Core metrics output dir %s already exists; using %s instead.",
            core_dir, alt
        )
        core_dir = alt


    run(
        [
            "qiime",
            "diversity",
            "core-metrics-phylogenetic",
            "--i-phylogeny",
            str(tree),
            "--i-table",
            str(table),
            "--p-sampling-depth",
            str(depth),
            "--m-metadata-file",
            str(metadata),
            "--p-n-jobs-or-threads",
            str(args.threads),
            "--output-dir",
            str(core_dir),
        ],
        log,
        logs / "core_metrics_phylogenetic.log",
    )


    # 2) Alpha diversity vectors + group significance (consistent with rarefaction)
    skip_alpha = False
    rarefied = core_dir / "rarefied_table.qza"
    if args.group_column:
        gc = _group_counts_after_rarefy(
            rarefied_table=rarefied,
            metadata_tsv=metadata,
            group_column=args.group_column,
        )
        if gc:
            too_small = {k: v for k, v in gc.items() if int(v) < 2}
            if too_small:
                log.warning(
                    "Skipping alpha-group-significance: too few samples per group after rarefaction: %s",
                    too_small,
                )
                skip_alpha = True
        keep_groups = [g for g, n in gc.items() if int(n) >= 2]
        skip_beta_group = False
        if args.group_column:
            if not keep_groups:
                log.warning("All groups are singletons after rarefaction; skipping beta-group-significance.")
                skip_beta_group = True
            else:
                # Build a safe SQL-ish WHERE for QIIME (single quotes, escaped)
                esc = lambda s: str(s).replace("'", "''")
                where_groups = ", ".join(f"'{esc(g)}'" for g in keep_groups)
                where_clause = f"{args.group_column} IN ({where_groups})"
                log.info("PERMANOVA will be run on groups with ≥2 samples: %s", keep_groups)


    if not skip_alpha:
        
        for metric in args.alpha_metrics:
            vec = _alpha_vectors_for_metric(
                metric=metric,
                core_dir=core_dir,
                art_dir=art,
                rarefied_table=rarefied,
                table_raw=table,
                log=log,
                logs_dir=logs,
            )
            if args.group_column:
                viz_alpha = viz / f"alpha_{metric}_by_{args.group_column}.qzv"
                run_optional(
                    cmd=[
                        "qiime", "diversity", "alpha-group-significance",
                        "--i-alpha-diversity", str(vec),
                        "--m-metadata-file", str(metadata),
                        "--p-ignore-missing-samples",  # <<< makes it tolerant
                        "--o-visualization", str(viz_alpha),
                    ],
                    log=log,
                    logfile=logs / f"alpha_group_{metric}.log",
                    warn=f"alpha-group-significance failed for {metric}; skipping.",
                )


    # 3) Beta diversity: distances, PCoA + Emperor + PERMANOVA if group column

    # after core-metrics runs:
    rarefied = core_dir / "rarefied_table.qza"

    core_map = {
        "braycurtis": (
            "bray_curtis_distance_matrix.qza",
            "bray_curtis_pcoa_results.qza",
            "bray_curtis_emperor.qzv",
        ),
        "jaccard": (
            "jaccard_distance_matrix.qza",
            "jaccard_pcoa_results.qza",
            "jaccard_emperor.qzv",
        ),
        "unweighted_unifrac": (
            "unweighted_unifrac_distance_matrix.qza",
            "unweighted_unifrac_pcoa_results.qza",
            "unweighted_unifrac_emperor.qzv",
        ),
        "weighted_unifrac": (
            "weighted_unifrac_distance_matrix.qza",
            "weighted_unifrac_pcoa_results.qza",
            "weighted_unifrac_emperor.qzv",
        ),
    }

    for metric in args.beta_metrics:
        # 1) Try to reuse artefacts made by core-metrics-phylogenetic
        dist = pcoa = emp = None
        if metric in core_map:
            d_name, p_name, e_name = core_map[metric]
            d_path, p_path, e_path = core_dir / d_name, core_dir / p_name, core_dir / e_name
            if d_path.exists() and p_path.exists() and e_path.exists():
                dist, pcoa, emp = d_path, p_path, e_path

        # 2) If missing, compute from the *rarefied* table (consistent with alpha)
        if dist is None:
            dist = art / f"{metric}_distance.qza"
            if _is_phylogenetic_beta(metric):
                run(
                    [
                        "qiime", "diversity", "beta-phylogenetic",
                        "--i-table", str(rarefied),
                        "--i-phylogeny", str(tree),
                        "--p-metric", metric,
                        "--p-n-jobs-or-threads", str(args.threads),
                        "--o-distance-matrix", str(dist),
                    ],
                    log, logs / f"beta_{metric}.log",
                )
            else:
                run(
                    [
                        "qiime", "diversity", "beta",
                        "--i-table", str(rarefied),
                        "--p-metric", metric,
                        "--o-distance-matrix", str(dist),
                    ],
                    log, logs / f"beta_{metric}.log",
                )

        if pcoa is None:
            pcoa = art / f"{metric}_pcoa.qza"
            run(
                ["qiime", "diversity", "pcoa",
                "--i-distance-matrix", str(dist),
                "--o-pcoa", str(pcoa)],
                log, logs / f"pcoa_{metric}.log",
            )

        if emp is None:
            emp = viz / f"{metric}_emperor.qzv"
            run(
                [
                    "qiime", "emperor", "plot",
                    "--i-pcoa", str(pcoa),
                    "--m-metadata-file", str(metadata),
                    "--o-visualization", str(emp),
                ],
                log, logs / f"emperor_{metric}.log",
            )

        if args.group_column and not skip_beta_group:
            # Filter the distance matrix to eligible groups
            dist_filt = art / f"{metric}_distance_filtered.qza"
            run(
                [
                    "qiime", "diversity", "filter-distance-matrix",
                    "--i-distance-matrix", str(dist),
                    "--m-metadata-file", str(metadata),
                    "--p-where", where_clause,
                    "--o-filtered-distance-matrix", str(dist_filt),
                ],
                log, logs / f"filter_dm_{metric}.log",
            )

            bg = viz / f"{metric}_PERMANOVA_{args.group_column}.qzv"
            run(
                [
                    "qiime", "diversity", "beta-group-significance",
                    "--i-distance-matrix", str(dist_filt),
                    "--m-metadata-file", str(metadata),
                    "--m-metadata-column", args.group_column,
                    "--o-visualization", str(bg),
                    "--p-pairwise",
                ],
                log, logs / f"beta_group_{metric}.log",
            )


    # 4) Taxonomy barplots at specified levels (if taxonomy available)
    if taxonomy_present and taxonomy:
        for level in args.tax_levels:
            collapsed = art / f"table_L{level}.qza"
            collapsed = art / f"table_L{level}.qza"
            run(
                [
                    "qiime", "taxa", "collapse",
                    "--i-table", str(table),
                    "--i-taxonomy", str(taxonomy),
                    "--p-level", str(level),
                    "--o-collapsed-table", str(collapsed),
                ],
                log, logs / f"taxa_collapse_L{level}.log",
            )
            write_taxa_summary_tsv(out_dir=out, art_dir=art, tax_levels=args.tax_levels, log=log)

    else:
        log.info("taxonomy.qza not found; skipping taxonomy barplots.")

    # 5) ANCOM (optional; requires group column)
    if args.run_ancom and args.group_column:
        comp = art / "composition_table.qza"
        run(
            [
                "qiime",
                "composition",
                "add-pseudocount",
                "--i-table",
                str(table),
                "--o-composition-table",
                str(comp),
            ],
            log,
            logs / "composition_add_pseudocount.log",
        )
        ancom_viz = viz / f"ancom_{args.group_column}.qzv"
        run(
            [
                "qiime",
                "composition",
                "ancom",
                "--i-table",
                str(comp),
                "--m-metadata-file",
                str(metadata),
                "--m-metadata-column",
                args.group_column,
                "--o-visualization",
                str(ancom_viz),
            ],
            log,
            logs / "ancom.log",
        )
    elif args.run_ancom and not args.group_column:
        log.warning("--run_ancom was set but no --group_column provided; skipping ANCOM.")

    # 6) Optional: export visuals (+ attempt PDF rendering)
    if args.export_visuals:
        export_root = out / "exports"
        export_root.mkdir(parents=True, exist_ok=True)
        for qzv in viz.glob("*.qzv"):
            dest = export_root / qzv.stem
            export_qzv(qzv, dest, log)
            try_html_to_pdf(dest, log, args.pdf_engine)


    # 7) README + summary TSV
    write_readme(
        out_dir=out,
        sampling_depth=depth,
        group_column=args.group_column,
        alpha_metrics=args.alpha_metrics,
        beta_metrics=args.beta_metrics,
        tax_levels=args.tax_levels,
        taxonomy_present=taxonomy_present,
    )
    write_summary_tsv(
        out_dir=out,
        depth_stats=depth_stats,
        sampling_depth=depth,
        taxonomy_present=taxonomy_present,
        alpha_metrics=args.alpha_metrics,
        beta_metrics=args.beta_metrics,
        tax_levels=args.tax_levels,
        group_column=args.group_column,
    )

    # 8) HTML landing page
    write_html_index(out)
    log.info("Wrote index.html")

    log.info("Done. Outputs in: %s", out)


if __name__ == "__main__":
    main()
