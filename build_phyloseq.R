#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(qiime2R)
  library(phyloseq)
  library(stringr)
  library(Biostrings)
})

# Build and save a phyloseq object from QIIME 2 artefacts.
# Expects in --input_dir (from q2_amplicon_runner.py):
#   - table.qza
#   - rep-seqs.qza
#   - rooted-tree.qza   (optional)
#   - taxonomy.qza      (optional)
#   - metadata.tsv
#
# All args are named; no positionals.

option_list <- list(
  make_option(c("--input_dir"), type = "character",
              help = "Directory with table.qza, rep-seqs.qza, (rooted-tree.qza), (taxonomy.qza), metadata.tsv"),
  make_option(c("--output_rds"), type = "character",
              help = "Path to save the phyloseq object (.rds)"),
  make_option(c("--length_range"), type = "character", default = NA,
              help = "Optional ASV length window (e.g., 240:260). Applied after construction.")
)
opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = FALSE)

stopifnot(!is.null(opt$input_dir), !is.null(opt$output_rds))
inp <- normalizePath(opt$input_dir, mustWork = TRUE)

req_files <- c("table.qza", "rep-seqs.qza", "metadata.tsv")
missing <- req_files[!file.exists(file.path(inp, req_files))]
if (length(missing) > 0) {
  stop(paste0("Missing required files in --input_dir: ", paste(missing, collapse = ", ")))
}

cat("Reading QIIME 2 artefacts ...\n")

# --- read artefacts ---
ft_qza   <- read_qza(file.path(inp, "table.qza"))        # feature table
seqs_qza <- read_qza(file.path(inp, "rep-seqs.qza"))     # representative sequences

tree_path <- file.path(inp, "rooted-tree.qza")
have_tree <- file.exists(tree_path)
if (have_tree) {
  tree_qza <- read_qza(tree_path)                        # rooted tree
}

# Metadata: treat leading '#' as comments so '#q2:types' is skipped
meta <- read.delim(
  file = file.path(inp, "metadata.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "#",
  stringsAsFactors = FALSE
)
if (ncol(meta) == 0) stop("metadata.tsv appears to have no columns.")
colnames(meta)[1] <- "sample-id"
rownames(meta) <- meta[["sample-id"]]
samp <- sample_data(meta)

# --- ensure table orientation is taxa-as-rows ---
ft <- ft_qza$data
if (is.matrix(ft) || inherits(ft, "Matrix")) {
  # Heuristic: taxa-as-rows is typical. If fewer rows than cols, flip.
  if (nrow(ft) < ncol(ft)) {
    ft <- t(ft)
  }
} else if (is.data.frame(ft)) {
  ft <- as.matrix(ft)
  if (nrow(ft) < ncol(ft)) ft <- t(ft)
} else {
  # Try a best-effort coerce
  ft <- as.matrix(ft)
  if (nrow(ft) < ncol(ft)) ft <- t(ft)
}
otu <- otu_table(ft, taxa_are_rows = TRUE)

# --- construct base phyloseq (tree optional) ---
if (have_tree) {
  ps <- phyloseq(otu, samp, phy_tree(tree_qza$data))
} else {
  warning("rooted-tree.qza not found; building phyloseq without a tree.")
  ps <- phyloseq(otu, samp)
}

# --- advisory: samples present in table but missing in metadata ---
otu_samples  <- sample_names(ps)
meta_samples <- rownames(meta)
if (!all(otu_samples %in% meta_samples)) {
  missing_in_meta <- setdiff(otu_samples, meta_samples)
  warning(sprintf(
    "Samples present in table but missing in metadata: %s",
    paste(missing_in_meta, collapse = ", ")
  ))
}
# prune to shared samples (keeps ps consistent)
common_samples <- intersect(sample_names(ps), rownames(meta))
if (length(common_samples) < nsamples(ps)) {
  dropped <- setdiff(sample_names(ps), common_samples)
  warning(sprintf("Dropping %d samples absent from metadata (e.g., %s)",
                  length(dropped), paste(utils::head(dropped, 5), collapse = ", ")))
  ps <- prune_samples(common_samples, ps)
}

# --- taxonomy (optional) ---
tax_path <- file.path(inp, "taxonomy.qza")
if (file.exists(tax_path)) {
  tax_df <- read_qza(tax_path)$data

  # Row key for features
  if ("Feature.ID" %in% colnames(tax_df)) {
    rownames(tax_df) <- tax_df$Feature.ID
  } else if (is.null(rownames(tax_df)) || any(!nzchar(rownames(tax_df)))) {
    stop("taxonomy.qza lacks 'Feature.ID' and row names are empty; cannot align taxonomy.")
  }

  # Expand "Taxon" to rank columns if present
  if ("Taxon" %in% colnames(tax_df)) {
    ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    split_list <- str_split(tax_df$Taxon, ";\\s*")
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
    # assume rank-wise columns already
    tax_tab <- as.matrix(tax_df[, setdiff(colnames(tax_df), c("Feature.ID","Confidence")), drop = FALSE])
  }

  # align taxa
  common_taxa <- intersect(taxa_names(ps), rownames(tax_tab))
  if (length(common_taxa) == 0) {
    stop("No overlapping ASV IDs between table and taxonomy.")
  }
  if (length(common_taxa) < ntaxa(ps)) {
    dropped <- setdiff(taxa_names(ps), common_taxa)
    warning(sprintf("Dropping %d taxa without taxonomy (first few: %s)",
                    length(dropped), paste(utils::head(dropped, 5), collapse = ", ")))
  }
  ps <- prune_taxa(common_taxa, ps)
  ps <- merge_phyloseq(ps, tax_table(tax_tab[common_taxa, , drop = FALSE]))
} else {
  warning("taxonomy.qza not found; building phyloseq without taxonomy.")
}

# --- optional ASV length filter (based on rep-seqs.qza) ---
if (!is.na(opt$length_range)) {
  parts <- str_split(opt$length_range, ":", n = 2)[[1]]
  ok <- length(parts) == 2 && all(suppressWarnings(!is.na(as.integer(parts))))
  if (ok) {
    min_len <- as.integer(parts[1]); max_len <- as.integer(parts[2])

    seq_obj <- seqs_qza$data
    if (inherits(seq_obj, "DNAStringSet")) {
      asv_ids <- names(seq_obj)
      asv_len <- Biostrings::width(seq_obj)
    } else {
      asv_ids <- names(seq_obj)
      if (is.null(asv_ids) || any(!nzchar(asv_ids))) {
        stop("rep-seqs.qza has an unexpected structure; cannot determine ASV IDs.")
      }
      asv_len <- nchar(as.character(seq_obj))
    }

    keep_ids <- asv_ids[asv_len >= min_len & asv_len <= max_len]
    keep_ids <- intersect(keep_ids, taxa_names(ps))
    ps <- prune_taxa(keep_ids, ps)

    cat(sprintf("Applied length filter: %d:%d (retained %d ASVs)\n",
                min_len, max_len, length(keep_ids)))
  } else {
    warning("Ignoring malformed --length_range; expected 'min:max' with integers.")
  }
}

# --- tidy up zeroes after any pruning ---
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

# --- quick summary + diagnostics ---
cat(sprintf(
  "Summary: %d samples | %d taxa | min/median/max reads: %s/%s/%s\n",
  nsamples(ps), ntaxa(ps),
  ifelse(nsamples(ps)>0, min(sample_sums(ps)), 0),
  ifelse(nsamples(ps)>0, median(sample_sums(ps)), 0),
  ifelse(nsamples(ps)>0, max(sample_sums(ps)), 0)
))

# sample read totals
utils::write.table(
  data.frame(sample = sample_names(ps), reads = sample_sums(ps)),
  file = file.path(inp, "sample_sums.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# taxa totals
utils::write.table(
  data.frame(FeatureID = taxa_names(ps), total = taxa_sums(ps)),
  file = file.path(inp, "taxa_sums.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ASV lengths (from rep-seqs)
seq_obj <- seqs_qza$data
if (inherits(seq_obj, "DNAStringSet")) {
  asv_ids <- names(seq_obj)
  asv_len <- Biostrings::width(seq_obj)
} else {
  asv_ids <- names(seq_obj)
  asv_len <- nchar(as.character(seq_obj))
}
len_df  <- data.frame(FeatureID = asv_ids, Length = as.integer(asv_len), stringsAsFactors = FALSE)
utils::write.table(len_df, file = file.path(inp, "asv_lengths.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# --- save ---
saveRDS(ps, file = opt$output_rds)
cat("Saved phyloseq object to: ", opt$output_rds, "\n", sep = "")
cat("Done.\n")
