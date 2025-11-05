library(qiime2R)
library(dplyr)
library(ggplot2)
library(readr)

plot_asv_length_hist <- function(rep_seqs_qza, out_pdf) {
  #' Plot a histogram of ASV reference sequence lengths from a QIIME 2 rep-seqs artefact.
  #'
  #' @param rep_seqs_qza Path to rep-seqs.qza.
  #' @param out_pdf Output PDF file path.
  #' @return Invisible NULL; writes a PDF.
  x <- read_qza(rep_seqs_qza)$data
  # x is a DNAStringSet-like; convert to lengths
  lens <- Biostrings::width(x)  # requires Biostrings loaded by qiime2R
  df <- tibble::tibble(Length = as.integer(lens))

  p <- ggplot(df, aes(x = Length)) +
    geom_histogram(bins = 40) +
    labs(x = "ASV reference sequence length (bp)", y = "Count",
         title = "ASV length distribution") +
    theme_bw()

  ggsave(filename = out_pdf, plot = p, width = 6, height = 4, units = "in")
  invisible(NULL)
}

# Example:
 plot_asv_length_hist("results/JH102/phyloseq_output/rep-seqs.qza",
                      "results/JH102/qc/dada2_asv_length_hist.pdf")



 plot_asv_length_hist("results/JH102_DEBLUR/phyloseq_output/rep-seqs.qza",
                      "results/JH102/qc/deblur_asv_length_hist.pdf")
