QIIME 2 Amplicon Analysis Pipeline

This pipeline processes paired-end amplicon sequencing data (e.g. 16S rRNA) using QIIME 2 with a choice of DADA2 or Deblur denoisers. It produces standard QIIME artefacts and a phyloseq_output/ directory ready for downstream analysis in R.


Install


mamba env create -n qiime2-amplicon-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

mamba install -n qiime2-amplicon-2024.10 -c conda-forge r-qiime2r

conda activate qiime2-amplicon-2024.10

mamba install flash

 mamba install -c conda-forge r-qiime2r r-stringi r-tidyverse


mamba install -c conda-forge wkhtmltopdf

R -q -e 'install.packages("remotes", repos="https://cloud.r-project.org")'

R -q -e 'remotes::install_github("jbisanz/qiime2R")'

download the classifier from:  https://library.qiime2.org/data-resources#naive-bayes-classifiers 







Contents

q2_amplicon_runner.py – main pipeline (Python)

build_phyloseq.R – converts artefacts to a phyloseq object

config.yaml – example configuration (for documentation and reproducibility)

check_map.sh – optional metadata checker

phyloseq_output/ – created automatically with artefacts and metadata

1. Prepare your inputs
Required files

FASTQs: paired-end .fastq.gz files.

Manifest: a PairedEndFastqManifestPhred33V2 TSV describing sample IDs and FASTQ locations.

Metadata: a QIIME-formatted metadata TSV (first column sample-id, second line #q2:types).

Optional

Pre-trained classifier: .qza trained for your primer pair and region (e.g. SILVA 138 V3-V4).

Primer sequences: if you want to trim primers with Cutadapt.

2. Validate metadata

Run the metadata checker before starting:

./check_map.sh -m metadata/metadata16S.tsv -M manifests/paired_manifest.tsv --no-fix


If any errors are reported, fix them before continuing.

3. Run the pipeline

Example: DADA2 paired-end analysis with primer trimming.

python q2_amplicon_runner.py \
  --run_label V3V4_STUDY1 \
  --manifest manifests/paired_manifest.tsv \
  --metadata_tsv metadata/metadata16S.tsv \
  --denoiser dada2_paired \
  --cutadapt_front_f CCTACGGGNGGCWGCAG \
  --cutadapt_front_r GACTACHVGGGTATCTAATCC \
  --trim_left_f 0 \
  --trim_left_r 0 \
  --trunc_len_f 240 \
  --trunc_len_r 200 \
  --max_ee_f 2.0 \
  --max_ee_r 2.0 \
  --threads 24 \
  --classifier_qza taxonomy/silva-138.1-99-515-806-nb-classifier.qza \
  --out_dir results/V3V4_STUDY1


Example: Deblur analysis with read-joining.

python q2_amplicon_runner.py \
  --run_label V3V4_STUDY1_DEBLUR \
  --manifest manifests/paired_manifest.tsv \
  --metadata_tsv metadata/metadata16S.tsv \
  --denoiser deblur_single \
  --pair_strategy join_in_qiime \
  --deblur_trim_length 250 \
  --threads 24 \
  --classifier_qza taxonomy/silva-138.1-99-515-806-nb-classifier.qza \
  --out_dir results_deblur/V3V4_STUDY1_DEBLUR

4. Outputs

Inside the --out_dir folder you’ll find:

Folder	Contents
artifacts_qza/	QIIME artefacts (.qza) for tables, sequences, stats
visuals_qzv/	QIIME visualisations (.qzv)
taxonomy/	Taxonomic assignment results
phylogeny/	Aligned sequences and phylogenetic trees
phyloseq_output/	Artefacts + metadata for R/phyloseq
logs/	Logs of each QIIME command
run_report.tsv	Summary of key parameters used
5. Build a phyloseq object

Once the run is complete, convert the artefacts into a ready-to-use phyloseq object:

Rscript build_phyloseq.R \
  --input_dir results/V3V4_STUDY1/phyloseq_output \
  --output_rds results/V3V4_STUDY1/phyloseq_output/phyloseq.rds \
  --length_range 240:260


This saves phyloseq.rds, which you can load directly in R:

ps <- readRDS("results/V3V4_STUDY1/phyloseq_output/phyloseq.rds")

6. Next steps

Once you have a phyloseq object, you can proceed with:

Alpha/beta diversity analyses

Differential abundance testing

Taxonomic summaries and visualisations

Notes & tips

Primer trimming (cutadapt) is highly recommended for improved denoising accuracy.

Choose truncation lengths based on the demux quality plots to ensure enough overlap for merging.

If DADA2 becomes computationally heavy for very large datasets, consider Deblur for scalability.

Always ensure metadata sample IDs match those in the manifest and FASTQs.
