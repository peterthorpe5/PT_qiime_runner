#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(qiime2R)
  library(phyloseq)
  library(stringr)
  library(Biostrings)
  library(ggplot2)
  library(ape)
})

# -------------------------------------------------------------------
# Helpers (PEP 8–style docstrings in comments)
# -------------------------------------------------------------------

# """
# normalise_ids(x: character) -> character
#
# Trim whitespace and collapse internal spaces to underscores.
# Safe normalisation for sample ID comparisons.
# """
normalise_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  gsub("\\s+", "_", x)
}

# """
# safe_read_metadata(file_path: str, id_col: str|NA) -> phyloseq::sample_data
#
# Read QIIME 2 metadata TSV, handling '#q2:types' and similar comment lines.
# If id_col is provided, use that column for rownames; otherwise use the first column.
# Applies simple ID normalisation. Returns sample_data().
# """
safe_read_metadata <- function(file_path, id_col = NA) {
  md <- read.delim(
    file = file_path, header = TRUE, sep = "\t",
    check.names = FALSE, quote = "", comment.char = "#",
    stringsAsFactors = FALSE
  )
  if (ncol(md) == 0) stop("metadata.tsv appears to have no columns.")

  if (!is.na(id_col)) {
    if (!id_col %in% colnames(md)) {
      stop(sprintf("Requested --metadata_id_col '%s' not found in metadata.tsv.", id_col))
    }
    sid_col <- id_col
  } else {
    sid_col <- colnames(md)[1]
  }

  md[[sid_col]] <- normalise_ids(md[[sid_col]])
  rownames(md) <- md[[sid_col]]

  if (sid_col == colnames(md)[1]) {
    colnames(md)[1] <- "sample-id"
  }
  sample_data(md)
}

# """
# orient_feature_table_by_overlap(x, meta_ids) -> list(mat, taxa_are_rows)
#
# Decide orientation by which axis overlaps the metadata sample IDs.
# Returns a numeric matrix and a logical flag for taxa_are_rows.
# """
orient_feature_table_by_overlap <- function(x, meta_ids) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x))   x <- as.matrix(x)
  mode(x) <- "numeric"

  col_overlap <- sum(normalise_ids(colnames(x)) %in% meta_ids)
  row_overlap <- sum(normalise_ids(rownames(x)) %in% meta_ids)

  if (col_overlap > 0 && row_overlap == 0) {
    return(list(mat = x, taxa_are_rows = TRUE))   # taxa rows, samples columns
  } else if (row_overlap > 0 && col_overlap == 0) {
    return(list(mat = t(x), taxa_are_rows = TRUE)) # transpose to taxa rows
  } else if (col_overlap > 0 && row_overlap > 0) {
    stop("Ambiguous orientation: metadata IDs appear in BOTH rows and columns of the table.")
  } else {
    stop("Cannot find metadata sample IDs in either rows or columns of the table.")
  }
}

# """
# asv_lengths_from_qza(seqs_qza) -> data.frame
#
# Extract ASV IDs and lengths from a rep-seqs.qza object (qiime2R::read_qza()).
# Returns a data.frame with columns: FeatureID, Length.
# """
asv_lengths_from_qza <- function(seqs_qza) {
  obj <- seqs_qza$data
  ids <- names(obj)
  if (is.null(ids) || any(!nzchar(ids))) {
    stop("rep-seqs.qza has no valid ASV IDs in names().")
  }
  lens <- if (inherits(obj, "DNAStringSet")) {
    Biostrings::width(obj)
  } else {
    nchar(as.character(obj))
  }
  data.frame(FeatureID = ids, Length = as.integer(lens), stringsAsFactors = FALSE)
}

# """
# apply_length_window(keep_ids_source: data.frame, min_len: int, max_len: int) -> character[]
#
# Given a data.frame with FeatureID and Length, return the subset of FeatureID
# within the inclusive [min_len, max_len] window.
# """
apply_length_window <- function(keep_ids_source, min_len, max_len) {
  keep <- keep_ids_source$FeatureID[
    keep_ids_source$Length >= min_len & keep_ids_source$Length <= max_len
  ]
  unique(as.character(keep))
}

# """
# align_taxonomy_to(ps: phyloseq, tax_qza: list) -> phyloseq
#
# Align and merge taxonomy with an existing phyloseq object.
# Accepts a QIIME taxonomy with either 'Taxon' or rank-wise columns.
# Drops taxa not present in ps; reorders to match ps taxa_names().
# """
align_taxonomy_to <- function(ps, tax_qza) {
  tax_df <- tax_qza$data

  if ("Feature.ID" %in% colnames(tax_df)) {
    rownames(tax_df) <- tax_df$Feature.ID
  } else if (is.null(rownames(tax_df)) || any(!nzchar(rownames(tax_df)))) {
    stop("taxonomy.qza lacks 'Feature.ID' and rownames are empty; cannot align taxonomy.")
  }

  if ("Taxon" %in% colnames(tax_df)) {
    ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    split_list <- stringr::str_split(tax_df$Taxon, ";\\s*")
    max_len <- min(max(vapply(split_list, length, integer(1))), length(ranks))
    split_mat <- do.call(rbind, lapply(split_list, function(x) {
      x <- x[seq_len(min(length(x), max_len))]
      length(x) <- max_len
      x
    }))
    colnames(split_mat) <- ranks[seq_len(max_len)]
    rownames(split_mat) <- rownames(tax_df)
    tax_tab <- as.matrix(split_mat)
  } else {
    drop_cols <- intersect(colnames(tax_df), c("Feature.ID", "Confidence"))
    tax_tab <- as.matrix(tax_df[, setdiff(colnames(tax_df), drop_cols), drop = FALSE])
  }

  common <- intersect(taxa_names(ps), rownames(tax_tab))
  if (length(common) == 0) {
    warning("No overlapping ASV IDs between feature table and taxonomy; returning ps without taxonomy.")
    return(ps)
  }
  if (length(common) < ntaxa(ps)) {
    dropped <- setdiff(taxa_names(ps), common)
    warning(sprintf("Dropping %d taxa without taxonomy (e.g., %s)",
                    length(dropped), paste(utils::head(dropped, 5), collapse = ", ")))
    ps <- prune_taxa(common, ps)
  }
  tax_tab <- tax_tab[taxa_names(ps), , drop = FALSE]
  merge_phyloseq(ps, tax_table(tax_tab))
}

# """
# align_tree_to(ps: phyloseq, tree_qza: list) -> phyloseq
#
# Prune and reorder a rooted tree so that its tip labels exactly match taxa_names(ps).
# Uses ape::keep.tip and then reorders tip.label to ps order before merging.
# """
align_tree_to <- function(ps, tree_qza) {
  tr <- tree_qza$data
  otus <- taxa_names(ps)
  tips <- tr$tip.label

  common <- intersect(otus, tips)
  if (length(common) == 0) {
    warning("No overlapping ASV IDs between feature table and tree; returning ps without tree.")
    return(ps)
  }
  if (length(common) < length(otus)) {
    dropped <- setdiff(otus, common)
    warning(sprintf("Dropping %d taxa absent from tree (e.g., %s)",
                    length(dropped), paste(utils::head(dropped, 5), collapse = ", ")))
    ps <- prune_taxa(common, ps)
  }
  if (length(common) < length(tips)) {
    tr <- ape::keep.tip(tr, taxa_names(ps))
  }
  tr$tip.label <- taxa_names(ps)
  merge_phyloseq(ps, phy_tree(tr))
}

# """
# plot_asv_length_hist(df: data.frame, out_pdf: str, title: str) -> None
#
# Plot a histogram of ASV reference sequence lengths (column 'Length') and save as PDF.
# """
plot_asv_length_hist <- function(df, out_pdf, title) {
  stopifnot(all(c("Length") %in% colnames(df)))
  p <- ggplot(df, aes(x = Length)) +
    geom_histogram(bins = 40) +
    labs(x = "ASV reference sequence length (bp)", y = "Count", title = title) +
    theme_bw()
  ggsave(filename = out_pdf, plot = p, width = 6, height = 4, units = "in")
}

# -------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------

option_list <- list(
  make_option(c("--input_dir"), type = "character",
              help = "Directory with table.qza, rep-seqs.qza, (rooted-tree.qza), (taxonomy.qza), metadata.tsv"),
  make_option(c("--output_rds"), type = "character",
              help = "Path to save the phyloseq object (.rds)"),
  make_option(c("--length_range"), type = "character", default = NA,
              help = "Optional ASV length window, e.g. 240:260. Applied BEFORE construction."),
  make_option(c("--metadata_id_col"), type = "character", default = NA,
              help = "Optional column in metadata.tsv to use as sample IDs. If unset, the first column is used."),
  make_option(c("--allow_partial_samples"), type = "character", default = "false",
              help = "Proceed with intersection if some samples are missing. Choices: true|false (default false).")
)

opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = FALSE)
stopifnot(!is.null(opt$input_dir), !is.null(opt$output_rds))

inp <- normalizePath(opt$input_dir, mustWork = TRUE)
tax_path  <- file.path(inp, "taxonomy.qza")
tree_path <- file.path(inp, "rooted-tree.qza")
req       <- c("table.qza", "rep-seqs.qza", "metadata.tsv")
miss      <- req[!file.exists(file.path(inp, req))]
if (length(miss) > 0) {
  stop(paste0("Missing required files in --input_dir: ", paste(miss, collapse = ", ")))
}

allow_partial <- tolower(opt$allow_partial_samples) %in% c("true","t","1","yes","y")

cat("Reading QIIME 2 artefacts ...\n")
ft_qza   <- read_qza(file.path(inp, "table.qza"))
seqs_qza <- read_qza(file.path(inp, "rep-seqs.qza"))
meta     <- safe_read_metadata(file.path(inp, "metadata.tsv"), id_col = opt$metadata_id_col)

have_tax  <- file.exists(tax_path)
have_tree <- file.exists(tree_path)
if (have_tax)  tax_qza  <- read_qza(tax_path)
if (have_tree) tree_qza <- read_qza(tree_path)

# Feature table orientation by overlap with metadata
meta_ids <- normalise_ids(rownames(data.frame(meta)))
oriented <- orient_feature_table_by_overlap(ft_qza$data, meta_ids)
ft       <- oriented$mat
taxa_rows_flag <- oriented$taxa_are_rows  # TRUE by construction

# Diagnostics: write raw sample IDs from table and metadata
utils::write.table(
  data.frame(sample_id = normalise_ids(colnames(ft))),
  file = file.path(inp, "diagnostic_table_samples.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  data.frame(sample_id = meta_ids),
  file = file.path(inp, "diagnostic_metadata_samples.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Prune to shared samples (metadata ↔ table columns)
table_samples <- normalise_ids(colnames(ft))
common_samples <- intersect(table_samples, meta_ids)
if (length(common_samples) == 0) {
  stop(paste0(
    "No overlapping samples between feature table and metadata after normalisation.\n",
    "Wrote 'diagnostic_table_samples.tsv' and 'diagnostic_metadata_samples.tsv' in ", inp, "."
  ))
}
if (!allow_partial && length(common_samples) < length(table_samples)) {
  stop(sprintf(
    "Found only %d overlapping samples out of %d. Re-run with --allow_partial_samples true if intersection is acceptable.",
    length(common_samples), length(table_samples)
  ))
}
ft   <- ft[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]
samp <- sample_data(meta)

# ASV lengths (BEFORE filtering) and plot
len_before <- asv_lengths_from_qza(seqs_qza)
utils::write.table(len_before, file = file.path(inp, "asv_lengths_before.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
plot_asv_length_hist(
  df = len_before,
  out_pdf = file.path(inp, "asv_length_hist_before.pdf"),
  title = "ASV length distribution (before filtering)"
)

# Optional length window (APPLIED BEFORE construction)
if (!is.na(opt$length_range)) {
  parts <- stringr::str_split(opt$length_range, ":", n = 2)[[1]]
  ok <- length(parts) == 2 && all(suppressWarnings(!is.na(as.integer(parts))))
  if (!ok) {
    warning("Ignoring malformed --length_range; expected 'min:max' with integers.")
  } else {
    min_len <- as.integer(parts[1]); max_len <- as.integer(parts[2])
    keep_by_len <- apply_length_window(len_before, min_len, max_len)
    keep_ids <- intersect(keep_by_len, rownames(ft))
    if (length(keep_ids) == 0) {
      warning("Length filter removed all taxa; continuing without applying the filter.")
    } else {
      ft <- ft[keep_ids, , drop = FALSE]
      cat(sprintf("Applied length filter early: %d:%d (retained %d ASVs)\n",
                  min_len, max_len, length(keep_ids)))
    }
  }
}

# Drop empty rows/cols defensively
ft <- ft[rowSums(ft) > 0, , drop = FALSE]
ft <- ft[, colSums(ft) > 0, drop = FALSE]

# Construct base phyloseq
otu <- otu_table(ft, taxa_are_rows = taxa_rows_flag)
ps  <- phyloseq(otu, samp)

# Taxonomy (optional, aligned and merged)
if (have_tax) {
  ps <- align_taxonomy_to(ps, tax_qza)
} else {
  warning("taxonomy.qza not found; building phyloseq without taxonomy.")
}

# Tree (optional, aligned and merged)
if (have_tree) {
  ps <- align_tree_to(ps, tree_qza)
} else {
  warning("rooted-tree.qza not found; building phyloseq without a tree.")
}

# Final tidy pruning
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

# ASV lengths AFTER filtering (restrict to taxa_names(ps)) and plot
len_after <- len_before[len_before$FeatureID %in% taxa_names(ps), , drop = FALSE]
utils::write.table(len_after, file = file.path(inp, "asv_lengths_after.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
plot_asv_length_hist(
  df = len_after,
  out_pdf = file.path(inp, "asv_length_hist_after.pdf"),
  title = "ASV length distribution (after filtering)"
)

# Diagnostics
utils::write.table(
  data.frame(sample = sample_names(ps), reads = sample_sums(ps)),
  file = file.path(inp, "sample_sums.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  data.frame(FeatureID = taxa_names(ps), total = taxa_sums(ps)),
  file = file.path(inp, "taxa_sums.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(sprintf(
  "Summary: %d samples | %d taxa | min/median/max reads: %s/%s/%s\n",
  nsamples(ps), ntaxa(ps),
  ifelse(nsamples(ps)>0, min(sample_sums(ps)), 0),
  ifelse(nsamples(ps)>0, median(sample_sums(ps)), 0),
  ifelse(nsamples(ps)>0, max(sample_sums(ps)), 0)
))

# Save
saveRDS(ps, file = opt$output_rds)
cat("Saved phyloseq object to: ", opt$output_rds, "\n", sep = "")
cat("Done.\n")
