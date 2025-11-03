#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract QIIME 2 alpha-group significance and PERMANOVA results from .qzv files.

Overview
--------
This script scans a folder of QIIME 2 visualisations (.qzv), unzips each one
in memory, and looks for plausible results tables (TSV) that contain p-values
from:
  • alpha-group significance (typically Kruskal–Wallis)
  • beta-group significance (PERMANOVA; global and optional pairwise)

It is designed to be conservative and robust:
  • It never assumes fixed internal filenames; instead it inspects every TSV
    found inside a .qzv and looks for columns that look like p-values.
  • It falls back to extracting `metric` and `group_column` from the .qzv
    filename patterns commonly produced by the post-analysis runner, e.g.:
        alpha_shannon_by_Treatment.qzv
        braycurtis_PERMANOVA_Treatment.qzv
  • If a .qzv contains no usable stats, it is simply skipped.

Outputs
-------
A single tab-separated file with columns:
    kind            alpha | beta
    test            e.g. Kruskal-Wallis | PERMANOVA
    metric          e.g. shannon | braycurtis
    group_column    e.g. Genotype
    group1          for pairwise rows (else NA)
    group2          for pairwise rows (else NA)
    statistic       stat value if available (else NA)
    p_value
    q_value         adjusted p if available (else NA)
    n_groups        if derivable (else NA)
    qzv_file        source visualisation filename

Requirements
------------
Python 3.9+, pandas. No QIIME environment is required (parsing is zip-based).

Usage
-----
python extract_qzv_stats.py \
    --visuals_dir results/JH102_DADA2/post_analysis/visuals \
    --out_tsv results/JH102_DADA2/post_analysis/stat_summaries/diversity_stats_summary.tsv
"""

from __future__ import annotations

import argparse
import io
import logging
import re
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


def setup_logging(*, verbose: bool) -> logging.Logger:
    """
    Configure logging for console output.

    Parameters
    ----------
    verbose : bool
        If True, sets level to DEBUG; otherwise INFO.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger("qzv_stats")
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    h.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(h)
    return logger


def read_tsv_from_zip(*, zf: zipfile.ZipFile, member: str) -> Optional[pd.DataFrame]:
    """
    Read a TSV file stored within a zip member into a DataFrame.

    Parameters
    ----------
    zf : zipfile.ZipFile
        Opened zip archive.
    member : str
        Member path inside the archive.

    Returns
    -------
    pandas.DataFrame or None
        DataFrame if parsed successfully; otherwise None.
    """
    try:
        with zf.open(member) as fh:
            raw = fh.read()
        df = pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str, engine="python")
        # normalise column names (strip, lower)
        df.columns = [c.strip() for c in df.columns]
        return df
    except Exception:
        return None


def infer_context_from_filename(*, filename: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Infer (kind, metric, group_column) from a typical .qzv filename.

    Parameters
    ----------
    filename : str
        Base name of the .qzv file.

    Returns
    -------
    tuple
        kind ('alpha' or 'beta' or None), metric (or None), group_column (or None).
    """
    base = Path(filename).stem
    # alpha_<metric>_by_<group>.qzv
    m = re.match(r"alpha_(?P<metric>[^_]+)_by_(?P<group>.+)$", base)
    if m:
        return "alpha", m.group("metric"), m.group("group")

    # <metric>_PERMANOVA_<group>.qzv
    m = re.match(r"(?P<metric>[^_]+)_PERMANOVA_(?P<group>.+)$", base, flags=re.IGNORECASE)
    if m:
        return "beta", m.group("metric"), m.group("group")

    # <metric>_emperor.qzv etc. (no group info)
    m = re.match(r"(?P<metric>[^_]+)_(emperor|pcoa)$", base, flags=re.IGNORECASE)
    if m:
        return "beta", m.group("metric"), None

    return None, None, None


def pick_first_present(*, df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    """
    Pick the first column present in a DataFrame from a list of candidates.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to inspect.
    candidates : Iterable[str]
        Candidate column names to try (case-sensitive).

    Returns
    -------
    str or None
        The first matching column name, or None if none are present.
    """
    for c in candidates:
        if c in df.columns:
            return c
        # try case-insensitive match
        for dc in df.columns:
            if dc.lower() == c.lower():
                return dc
    return None


def parse_alpha_stats(*, df: pd.DataFrame) -> Optional[Dict[str, str]]:
    """
    Parse alpha-group significance summary row from a DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Candidate table extracted from an alpha .qzv.

    Returns
    -------
    dict or None
        Mapping with keys 'test', 'statistic', 'p_value', 'q_value', 'n_groups' if found.
    """
    # Common QIIME table tends to have a single KW summary row.
    # Look for p-value & statistic-like fields.
    pcol = pick_first_present(df=df, candidates=["p-value", "p value", "p_value", "p"])
    statcol = pick_first_present(df=df, candidates=["H", "H (Kruskal-Wallis)", "Statistic", "statistic"])
    qcol = pick_first_present(df=df, candidates=["q-value", "FDR p-value", "p-value (FDR-corrected)", "p (FDR)", "q"])
    ncol = pick_first_present(df=df, candidates=["Number of groups", "n_groups", "k"])

    if pcol is None and statcol is None:
        return None

    # Use the first row if there are multiple
    row = df.iloc[0].to_dict()
    return {
        "test": "Kruskal-Wallis",
        "statistic": str(row.get(statcol, "NA")) if statcol else "NA",
        "p_value": str(row.get(pcol, "NA")) if pcol else "NA",
        "q_value": str(row.get(qcol, "NA")) if qcol else "NA",
        "n_groups": str(row.get(ncol, "NA")) if ncol else "NA",
    }


def parse_permanova_stats(*, df: pd.DataFrame) -> List[Dict[str, str]]:
    """
    Parse PERMANOVA results (global and/or pairwise) from a DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Candidate table extracted from a beta-group .qzv.

    Returns
    -------
    list of dict
        Each dict has keys: 'test','statistic','p_value','q_value','group1','group2'.
    """
    # Heuristics for columns seen in QIIME exports
    pcol = pick_first_present(df=df, candidates=["p-value", "p value", "p_value", "p"])
    statcol = pick_first_present(df=df, candidates=["pseudo-F", "F", "F-value", "statistic"])
    g1 = pick_first_present(df=df, candidates=["Group 1", "group 1", "group1", "group_A", "A"])
    g2 = pick_first_present(df=df, candidates=["Group 2", "group 2", "group2", "group_B", "B"])
    qcol = pick_first_present(df=df, candidates=["q-value", "FDR p-value", "p-value (FDR-corrected)", "p (FDR)", "q"])

    # If there's a single-row global table (no groups), still capture it.
    out: List[Dict[str, str]] = []
    if pcol is None and statcol is None:
        return out

    # If group columns are present, treat as pairwise; else as global.
    if g1 and g2:
        for _, r in df.iterrows():
            out.append({
                "test": "PERMANOVA",
                "statistic": str(r.get(statcol, "NA")) if statcol else "NA",
                "p_value": str(r.get(pcol, "NA")) if pcol else "NA",
                "q_value": str(r.get(qcol, "NA")) if qcol else "NA",
                "group1": str(r.get(g1, "NA")),
                "group2": str(r.get(g2, "NA")),
            })
    else:
        # Use first row as the global test summary
        r0 = df.iloc[0]
        out.append({
            "test": "PERMANOVA",
            "statistic": str(r0.get(statcol, "NA")) if statcol else "NA",
            "p_value": str(r0.get(pcol, "NA")) if pcol else "NA",
            "q_value": str(r0.get(qcol, "NA")) if qcol else "NA",
            "group1": "NA",
            "group2": "NA",
        })
    return out


def harvest_from_qzv(*, qzv_path: Path, logger: logging.Logger) -> List[Dict[str, str]]:
    """
    Harvest alpha-group and PERMANOVA stats from a single .qzv file.

    Parameters
    ----------
    qzv_path : pathlib.Path
        Path to the .qzv file.
    logger : logging.Logger
        Logger for progress messages.

    Returns
    -------
    list of dict
        Zero or more parsed result rows.
    """
    kind, metric, group_col = infer_context_from_filename(filename=qzv_path.name)
    rows: List[Dict[str, str]] = []

    with zipfile.ZipFile(qzv_path, "r") as zf:
        # Consider any TSV under data/ or root as a candidate
        tsv_members = [m for m in zf.namelist() if m.lower().endswith(".tsv")]
        for member in tsv_members:
            df = read_tsv_from_zip(zf=zf, member=member)
            if df is None or df.empty:
                continue

            # Try alpha first if filename suggests so, else both
            parsed_alpha = parse_alpha_stats(df=df)
            if parsed_alpha and (kind == "alpha" or kind is None):
                rows.append({
                    "kind": "alpha",
                    "test": parsed_alpha["test"],
                    "metric": metric or "unknown",
                    "group_column": group_col or "unknown",
                    "group1": "NA",
                    "group2": "NA",
                    "statistic": parsed_alpha["statistic"],
                    "p_value": parsed_alpha["p_value"],
                    "q_value": parsed_alpha["q_value"],
                    "n_groups": parsed_alpha["n_groups"],
                    "qzv_file": qzv_path.name,
                })
                continue  # avoid double-counting same file as beta

            # Try PERMANOVA
            parsed_beta = parse_permanova_stats(df=df)
            if parsed_beta:
                for r in parsed_beta:
                    rows.append({
                        "kind": "beta",
                        "test": r["test"],
                        "metric": metric or "unknown",
                        "group_column": group_col or "unknown",
                        "group1": r["group1"],
                        "group2": r["group2"],
                        "statistic": r["statistic"],
                        "p_value": r["p_value"],
                        "q_value": r["q_value"],
                        "n_groups": "NA",
                        "qzv_file": qzv_path.name,
                    })

    if not rows:
        logger.debug("No stats harvested from %s", qzv_path.name)
    return rows


def main() -> None:
    """
    Command-line entrypoint. Scans a visuals directory of .qzv files and
    writes a single TSV summary with extracted statistics.
    """
    ap = argparse.ArgumentParser(
        description="Extract alpha-group significance and PERMANOVA stats from QIIME 2 .qzv files."
    )
    ap.add_argument(
        "--visuals_dir",
        required=True,
        type=Path,
        help="Directory containing .qzv outputs (e.g., post_analysis/visuals).",
    )
    ap.add_argument(
        "--out_tsv",
        required=True,
        type=Path,
        help="Path to write the combined TSV summary.",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose debugging logs.",
    )
    args = ap.parse_args()

    logger = setup_logging(verbose=args.verbose)
    visuals_dir = args.visuals_dir
    out_tsv = args.out_tsv
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    qzvs = sorted(visuals_dir.glob("*.qzv"))
    if not qzvs:
        logger.warning("No .qzv files found in %s", visuals_dir)
        # still write an empty table with header for reproducibility
        pd.DataFrame(
            columns=[
                "kind", "test", "metric", "group_column", "group1", "group2",
                "statistic", "p_value", "q_value", "n_groups", "qzv_file"
            ]
        ).to_csv(out_tsv, sep="\t", index=False)
        return

    all_rows: List[Dict[str, str]] = []
    for q in qzvs:
        try:
            rows = harvest_from_qzv(qzv_path=q, logger=logger)
            all_rows.extend(rows)
        except zipfile.BadZipFile:
            logger.warning("Skipping corrupt .qzv: %s", q.name)
        except Exception as e:
            logger.warning("Skipping %s due to error: %s", q.name, e)

    if not all_rows:
        logger.info("No statistics discovered; writing an empty summary.")
        cols = ["kind", "test", "metric", "group_column", "group1", "group2",
                "statistic", "p_value", "q_value", "n_groups", "qzv_file"]
        pd.DataFrame(columns=cols).to_csv(out_tsv, sep="\t", index=False)
        return

    df = pd.DataFrame(all_rows)
    # Light normalisation
    for col in ("p_value", "q_value", "statistic"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    df.to_csv(out_tsv, sep="\t", index=False)
    logger.info("Wrote %d rows to %s", len(df), out_tsv)


if __name__ == "__main__":
    main()
