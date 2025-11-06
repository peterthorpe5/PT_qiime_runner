#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a QIIME 2 PairedEndFastqManifestPhred33V2 from FASTQs, optionally using metadata.

This utility supports three modes:

A) Metadata mapping by key:
   Provide --metadata-tsv, --filename-key-column and --sample-id-column.
   The value in --filename-key-column (e.g. 'DNA.code') must appear in FASTQ
   filenames and will be mapped to the desired 'sample-id' from metadata.

B) Metadata present, no key column:
   Provide --metadata-tsv and --sample-id-column (omit --filename-key-column).
   Assumes filenames already contain the sample-id as the stem before the R1/R2
   marker; the metadata is used to limit/validate sample IDs and to name symlinks.

C) No metadata:
   Omit --metadata-tsv. The sample-id is inferred directly from the filename
   stem before the R1/R2 marker.

Optionally, tidy symlinks are created as <sample-id>_R1.fastq.gz / _R2.fastq.gz.
All tabular files are TSV. All arguments are named (no positional arguments).
UK English spelling is used in docstrings.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


print("[DEBUG] make_paired_manifest.py loaded – boundary-aware matcher active", flush=True)


def normalise_sample_id(text: str) -> str:
    """Return a normalised sample identifier safe for filenames and downstream tools.

    The normalisation performs:
    - Strip leading/trailing whitespace.
    - Replace spaces and underscores with '.'.
    - Replace bracketed numeric '(123)' with '.123'.
    - Coerce any remaining disallowed characters to '.' (allowed: [A-Za-z0-9.-]).

    Args:
        text: The original sample identifier.

    Returns:
        A normalised sample identifier.
    """
    s = text.strip()
    s = s.replace(" ", ".").replace("_", ".")
    s = re.sub(r"\(([0-9]+)\)", r".\1", s)
    s = re.sub(r"[^A-Za-z0-9.\-]+", ".", s)
    return s


def _prefix_boundary_match(stem: str, sid: str) -> bool:
    """
    Return True if 'sid' is a prefix of 'stem' with a word-like boundary.

    A boundary is end-of-string or one of [._-] immediately after the prefix.
    Tries four variants to bridge '.' and '_' differences:

      1) stem startswith sid
      2) (stem with '_'→'.') startswith sid
      3) stem startswith (sid with '.'→'_')
      4) (stem with '_'→'.') startswith (sid with '.'→'_')

    Examples
    --------
    Ema.10  ↔  Ema_10_S153    → match
    ES.B57RepA  ↔  ES_B57RepA_S106 → match
    Ema.1   ↔  Ema_10_S153    → no match (protected by boundary)
    """
    def _match(s: str, prefix: str) -> bool:
        # anchor at start; require boundary after the prefix
        pat = r'^' + re.escape(prefix) + r'(?:$|[._-])'
        return re.match(pat, s) is not None

    sid_dot2us = sid.replace(".", "_")
    stem_us2dot = stem.replace("_", ".")
    return (
        _match(stem, sid)
        or _match(stem_us2dot, sid)
        or _match(stem, sid_dot2us)
        or _match(stem_us2dot, sid_dot2us)
    )


def read_metadata(
    *,
    metadata_tsv: Path,
    sample_id_column: str,
    filename_key_column: Optional[str],
) -> Tuple[Optional[Dict[str, str]], Optional[List[str]]]:
    """Read metadata and return mappings for Mode A or a sample-id list for Mode B.

    If ``filename_key_column`` is provided (Mode A), a mapping from that column
    to normalised sample IDs is returned. Otherwise (Mode B), a list of normalised
    sample IDs is returned.

    Args:
        metadata_tsv: Path to the QIIME-style metadata TSV.
        sample_id_column: Column providing the desired final 'sample-id' values.
        filename_key_column: Column expected to occur in filenames (e.g. 'DNA.code').
            If ``None``, mapping-by-key is disabled and Mode B is selected.

    Returns:
        A tuple ``(key_to_sample, sample_ids)`` where exactly one is not ``None``:
        - ``key_to_sample``: dict mapping filename key -> normalised sample-id (Mode A).
        - ``sample_ids``: list of normalised sample-ids from metadata (Mode B).

    Raises:
        ValueError: If required columns are missing or duplicates are detected.
    """
    with metadata_tsv.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        if sample_id_column not in cols:
            raise ValueError(
                f"Column '{sample_id_column}' not found in metadata. "
                f"Columns present: {', '.join(cols)}"
            )

        if filename_key_column:
            if filename_key_column not in cols:
                raise ValueError(
                    f"Column '{filename_key_column}' not found in metadata. "
                    f"Columns present: {', '.join(cols)}"
                )
            key_to_sample: Dict[str, str] = {}
            seen_samples: set[str] = set()
            seen_keys: set[str] = set()
            for row in reader:
                sample = normalise_sample_id(row.get(sample_id_column, ""))
                key = (row.get(filename_key_column, "") or "").strip()
                if not sample or not key:
                    continue
                if sample in seen_samples:
                    raise ValueError(f"Duplicate sample-id in metadata: '{sample}'")
                if key in seen_keys:
                    raise ValueError(
                        f"Duplicate filename-key in metadata: '{key}' "
                        f"(column '{filename_key_column}')"
                    )
                seen_samples.add(sample)
                seen_keys.add(key)
                key_to_sample[key] = sample
            return key_to_sample, None

        sample_ids: List[str] = []
        seen_samples = set()
        for row in reader:
            sample = normalise_sample_id(row.get(sample_id_column, ""))
            if not sample:
                continue
            if sample in seen_samples:
                raise ValueError(f"Duplicate sample-id in metadata: '{sample}'")
            seen_samples.add(sample)
            sample_ids.append(sample)
        return None, sample_ids


def find_fastqs(
    *,
    reads_dir: Path,
    exts: Tuple[str, ...] = (".fastq.gz", ".fq.gz", ".fastq", ".fq"),
) -> List[Path]:
    """Recursively discover FASTQ files under a directory.

    Args:
        reads_dir: Root directory to scan.
        exts: Allowed filename extensions.

    Returns:
        A list of absolute ``Path`` objects to FASTQ files.
    """
    paths: List[Path] = []
    for p in reads_dir.rglob("*"):
        if p.is_file() and any(str(p).endswith(e) for e in exts):
            paths.append(p.resolve())
    return paths


def split_r1_r2(
    *,
    files: Iterable[Path],
    r1_pattern: str,
    r2_pattern: str,
) -> Tuple[List[Path], List[Path]]:
    """Split an iterable of paths into forward (R1) and reverse (R2) lists.

    Args:
        files: FASTQ paths to examine.
        r1_pattern: Substring designating forward reads in filenames.
        r2_pattern: Substring designating reverse reads in filenames.

    Returns:
        A pair ``(r1_list, r2_list)`` of ``Path`` lists.
    """
    r1, r2 = [], []
    for f in files:
        name = f.name
        if r1_pattern in name:
            r1.append(f)
        elif r2_pattern in name:
            r2.append(f)
    return r1, r2


def infer_stem(name: str, stem_regex: re.Pattern[str]) -> Optional[str]:
    """Extract the sample stem from a filename using a compiled regex.

    The regex must contain exactly one capture group which represents the stem.

    Args:
        name: Filename to parse.
        stem_regex: Compiled regex with one capture group for the stem.

    Returns:
        The extracted stem, or ``None`` if no match is found.
    """
    m = stem_regex.search(name)
    return m.group(1) if m else None


def group_by_stem(
    *,
    files: Iterable[Path],
    stem_regex: re.Pattern[str],
    r1_pattern: str,
    r2_pattern: str,
) -> Dict[str, Dict[str, List[Path]]]:
    """Group FASTQs by inferred stem and separate into R1/R2 lists.

    Args:
        files: FASTQ paths to group.
        stem_regex: Compiled regex with one capture group for the stem.
        r1_pattern: Substring designating forward reads.
        r2_pattern: Substring designating reverse reads.

    Returns:
        A dictionary mapping ``stem -> {'R1': [...], 'R2': [...]}``.
    """
    groups: Dict[str, Dict[str, List[Path]]] = {}
    for f in files:
        stem = infer_stem(f.name, stem_regex)
        if not stem:
            continue
        bucket = groups.setdefault(stem, {"R1": [], "R2": []})
        if r1_pattern in f.name:
            bucket["R1"].append(f)
        elif r2_pattern in f.name:
            bucket["R2"].append(f)
    return groups


def write_manifest(*, rows: List[Tuple[str, Path, Path]], out_path: Path) -> None:
    """Write a PairedEndFastqManifestPhred33V2 TSV.

    Args:
        rows: List of tuples ``(sample_id, r1_path, r2_path)``.
        out_path: Destination TSV path.

    Returns:
        None
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for sid, r1, r2 in rows:
            fh.write(f"{sid}\t{r1.resolve()}\t{r2.resolve()}\n")


def safe_symlink(src: Path, dest: Path) -> None:
    """Create or replace a symlink from ``dest`` to ``src``.

    Args:
        src: Existing source path.
        dest: Symlink path to create (will be replaced if present).

    Returns:
        None

    Raises:
        OSError: If creation fails.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    os.symlink(src, dest)


def build_rows_mode_a(
    *,
    groups: Dict[str, Dict[str, List[Path]]],
    key_to_sample: Dict[str, str],
    symlink_dir: Optional[Path],
) -> List[Tuple[str, Path, Path]]:
    """Build manifest rows by matching stems using metadata filename keys (Mode A).

    Each metadata key is matched as a substring within the stem. Multiple lanes
    are supported; pairs are ordered and, when symlinking, suffixed as _L###.

    Args:
        groups: Mapping of stem to R1/R2 lists.
        key_to_sample: Mapping from filename key -> normalised sample-id.
        symlink_dir: Directory for tidy symlinks, or ``None`` to disable.

    Returns:
        A list of ``(sample_id, r1_path, r2_path)`` tuples.

    Raises:
        RuntimeError: If any metadata keys cannot be matched to stems.
    """
    rows: List[Tuple[str, Path, Path]] = []
    missing: List[str] = []
    for key, sample in key_to_sample.items():
        matched = [(stem, rr) for stem, rr in groups.items() if key in stem]
        if not matched:
            missing.append(f"{key} -> {sample}")
            continue
        for _, rr in matched:
            r1s = sorted(rr["R1"])
            r2s = sorted(rr["R2"])
            n = min(len(r1s), len(r2s))
            for i in range(n):
                r1, r2 = r1s[i], r2s[i]
                rows.append((sample, r1, r2))
                if symlink_dir:
                    suffix = "" if n == 1 else f"_L{i+1:03d}"
                    safe_symlink(r1, symlink_dir / f"{sample}_R1{suffix}.fastq.gz")
                    safe_symlink(r2, symlink_dir / f"{sample}_R2{suffix}.fastq.gz")
    if missing:
        raise RuntimeError(
            "Unmatched metadata keys (no stems contained these substrings):\n  - "
            + "\n  - ".join(missing)
        )
    return rows


def build_rows_mode_b(
    *,
    groups: dict[str, dict[str, list[Path]]],
    sample_ids: list[str],
    symlink_dir: Path | None,
) -> list[tuple[str, Path, Path]]:
    """
    Build manifest rows when we have metadata but no explicit filename key.

    Logic
    -----
    - Treat each metadata sample-id as a prefix of the FASTQ stem.
    - If no direct prefix match, retry with '.' <-> '_' substitutions.
    - Enforce a boundary after the prefix to avoid 'Ema.1' matching 'Ema_10...'.
    - Emit all R1/R2 pairs for any matched stems (multi-lane friendly).
    """
    print(f"[DEBUG] build_rows_mode_b: {len(sample_ids)} sample IDs from metadata", flush=True)
    print(f"[DEBUG] build_rows_mode_b: {len(groups)} FASTQ stems detected", flush=True)

    for sid in sample_ids[:10]:  # show first 10 only
        matched = [(stem, rr) for stem, rr in groups.items() if _prefix_boundary_match(stem, sid)]
        if matched:
            print(f"[DEBUG] {sid}: matched {len(matched)} stems (e.g. {list(dict(matched).keys())[:3]})", flush=True)
        else:
            print(f"[DEBUG] {sid}: no matches", flush=True)


    rows: list[tuple[str, Path, Path]] = []
    missing: list[str] = []

    for sid in sample_ids:
        matches = [(stem, rr) for stem, rr in groups.items()
                   if _prefix_boundary_match(stem, sid)]
        if not matches:
            missing.append(sid)
            continue

        for _, rr in matches:
            r1s = sorted(rr["R1"])
            r2s = sorted(rr["R2"])
            n = min(len(r1s), len(r2s))
            for i in range(n):
                r1, r2 = r1s[i], r2s[i]
                rows.append((sid, r1, r2))
                if symlink_dir:
                    suffix = "" if n == 1 else f"_L{i+1:03d}"
                    safe_symlink(r1, symlink_dir / f"{sid}_R1{suffix}.fastq.gz")
                    safe_symlink(r2, symlink_dir / f"{sid}_R2{suffix}.fastq.gz")

    if missing:
        raise RuntimeError(
            "Samples from metadata not found as stems in filenames:\n  - "
            + "\n  - ".join(missing)
        )
    return rows



def build_rows_mode_c(
    *,
    groups: Dict[str, Dict[str, List[Path]]],
    symlink_dir: Optional[Path],
) -> List[Tuple[str, Path, Path]]:
    """Build manifest rows by inferring sample-id directly from stems (Mode C).

    Args:
        groups: Mapping of stem to R1/R2 lists.
        symlink_dir: Directory for tidy symlinks, or ``None`` to disable.

    Returns:
        A list of ``(sample_id, r1_path, r2_path)`` tuples.
    """
    rows: List[Tuple[str, Path, Path]] = []
    for stem, rr in sorted(groups.items()):
        sid = normalise_sample_id(stem)
        r1s = sorted(rr["R1"])
        r2s = sorted(rr["R2"])
        n = min(len(r1s), len(r2s))
        for i in range(n):
            r1, r2 = r1s[i], r2s[i]
            rows.append((sid, r1, r2))
            if symlink_dir:
                suffix = "" if n == 1 else f"_L{i+1:03d}"
                safe_symlink(r1, symlink_dir / f"{sid}_R1{suffix}.fastq.gz")
                safe_symlink(r2, symlink_dir / f"{sid}_R2{suffix}.fastq.gz")
    return rows


def build_manifest(
    *,
    reads_dir: Path,
    manifest_out: Path,
    r1_pattern: str,
    r2_pattern: str,
    stem_regex: str,
    metadata_tsv: Optional[Path],
    sample_id_column: Optional[str],
    filename_key_column: Optional[str],
    symlink_dir: Optional[Path],
    dry_run: bool,
) -> None:
    """Create the paired-end manifest (and optional symlinks) across Modes A/B/C.

    Args:
        reads_dir: Root directory scanned for FASTQs (recursively).
        manifest_out: Destination TSV path for the manifest.
        r1_pattern: Substring marking forward reads in filenames.
        r2_pattern: Substring marking reverse reads in filenames.
        stem_regex: Regex string with one capture group to extract the stem.
        metadata_tsv: Optional metadata TSV path.
        sample_id_column: Column providing 'sample-id' when metadata is used.
        filename_key_column: Column that appears in filenames (Mode A only).
        symlink_dir: If set, create symlinks for each pair at this path.
        dry_run: If ``True``, do not write files; print intended actions.

    Returns:
        None

    Raises:
        RuntimeError: If FASTQs are missing or metadata entries are unmatched.
        ValueError: If metadata arguments are inconsistent.
    """
    fastqs = find_fastqs(reads_dir=reads_dir)
    if not fastqs:
        raise RuntimeError(f"No FASTQs found under: {reads_dir}")

    stem_re = re.compile(stem_regex)
    groups = group_by_stem(
        files=fastqs, stem_regex=stem_re, r1_pattern=r1_pattern, r2_pattern=r2_pattern
    )

    # Warn if stems lack R1/R2 after grouping (common when patterns are too strict)
    empty_pairs = [s for s, rr in groups.items() if not rr["R1"] or not rr["R2"]]
    if empty_pairs:
        print(f"[WARN] {len(empty_pairs)} stems lack complete R1/R2 pairs. "
            f"Consider using --r1-pattern _R1 and --r2-pattern _R2. "
            f"Example: {empty_pairs[:3]}", flush=True)


    if metadata_tsv:
        if not sample_id_column:
            raise ValueError("--sample-id-column is required when --metadata-tsv is provided.")
        key_to_sample, sample_ids = read_metadata(
            metadata_tsv=metadata_tsv,
            sample_id_column=sample_id_column,
            filename_key_column=filename_key_column,
        )
        rows = (
            build_rows_mode_a(groups=groups, key_to_sample=key_to_sample or {}, symlink_dir=symlink_dir)
            if key_to_sample is not None
            else build_rows_mode_b(groups=groups, sample_ids=sample_ids or [], symlink_dir=symlink_dir)
        )
    else:
        rows = build_rows_mode_c(groups=groups, symlink_dir=symlink_dir)

    if dry_run:
        print(f"[dry-run] would write manifest: {manifest_out}")
        for sid, r1, r2 in rows:
            print(f"[dry-run] {sid}\t{r1}\t{r2}")
    else:
        write_manifest(rows=rows, out_path=manifest_out)


def build_parser() -> argparse.ArgumentParser:
    """Return the command-line argument parser for this utility.

    Returns:
        An ``argparse.ArgumentParser`` configured for named-only arguments.
    """
    p = argparse.ArgumentParser(
        description=(
            "Create a QIIME 2 PairedEndFastqManifestPhred33V2 from FASTQs, "
            "with optional metadata mapping."
        ),
        allow_abbrev=False,
    )
    p.add_argument(
        "--reads-dir",
        required=True,
        type=Path,
        help="Directory containing FASTQs (scanned recursively).",
    )
    p.add_argument(
        "--manifest-out",
        required=True,
        type=Path,
        help="Output path for the paired manifest TSV.",
    )
    p.add_argument(
        "--r1-pattern",
        default="_R1_",
        type=str,
        help="Substring marking forward reads (default: _R1_).",
    )
    p.add_argument(
        "--r2-pattern",
        default="_R2_",
        type=str,
        help="Substring marking reverse reads (default: _R2_).",
    )
    p.add_argument(
        "--stem-regex",
        default=r"^(.+?)_R[12]",
        type=str,
        help=r"Regex with one capture group to extract sample stem (default: '^(.+?)_R[12]').",
    )
    p.add_argument(
        "--metadata-tsv",
        type=Path,
        default=None,
        help="Metadata TSV (optional).",
    )
    p.add_argument(
        "--sample-id-column",
        type=str,
        default=None,
        help="Metadata column providing desired sample IDs (required if --metadata-tsv is set).",
    )
    p.add_argument(
        "--filename-key-column",
        type=str,
        default=None,
        help=(
            "Metadata column expected in filenames (optional). "
            "If omitted while using --metadata-tsv, Mode B is used."
        ),
    )
    p.add_argument(
        "--symlink-dir",
        type=Path,
        default=None,
        help=(
            "If set, create symlinks <sample-id>_R1.fastq.gz / _R2.fastq.gz here "
            "(multi-lane pairs receive _L###)."
        ),
    )
    p.add_argument(
        "--dry-run",
        default=False,
        type=lambda x: str(x).lower() in {"1", "true", "yes"},
        help="Print intended actions without writing files (default: false).",
    )
    return p


def main() -> None:
    """Entry point: parse arguments and build the manifest (and symlinks)."""
    args = build_parser().parse_args()
    try:
        build_manifest(
            reads_dir=args.reads_dir,
            manifest_out=args.manifest_out,
            r1_pattern=args.r1_pattern,
            r2_pattern=args.r2_pattern,
            stem_regex=args.stem_regex,
            metadata_tsv=args.metadata_tsv,
            sample_id_column=args.sample_id_column,
            filename_key_column=args.filename_key_column,
            symlink_dir=args.symlink_dir,
            dry_run=bool(args.dry_run),
        )
    except (ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
