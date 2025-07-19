# AMCA
AMCA: GUI-Based QIIME2 Microbiome Analysis Workflow

This project provides a user-friendly, Zenity GUI-driven pipeline for performing comprehensive microbiome analysis using QIIME2. Designed to simplify the analysis process, the pipeline guides users through host read depletion, sequence import, denoising, diversity and taxonomic analyses, and biomarker discovery using LEfSe—all with step-by-step dialog prompts.

Features

Optional host read depletion using Bowtie2 (Human GRCh38 or custom index)
Zenity GUI for user interaction at every step
Manifest generation and FASTQ import to QIIME2
DADA2 denoising with customizable trim/truncation parameters
Phylogeny and core diversity metrics
Alpha rarefaction and optional beta significance testing
Taxonomic classification via pre-trained or custom-trained classifiers
LEfSe-compatible feature table formatting and downstream biomarker analysis

Prerequisites

Ensure the following tools are installed and available in your `PATH`:

zenity (https://help.gnome.org/users/zenity/stable/)
qiime2 (https://qiime2.org/)
bowtie2(http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
biom   (https://github.com/biocore/biom-format)
Python3 with `pandas`
LEfSe Python scripts:
lefse_format_input.py
lefse_run.py
lefse_plot_res.py
lefse_plot_cladogram.py
Pre-trained taxonomic classifier for V3-V4 16S rRNA and associated files at : 10.6084/m9.figshare.29602949


Input Requirements

Place the following files in your **working directory** before running the script:

Paired-end FASTQ files (`*_1.fastq`, `*_2.fastq`)
metadata.tsv (sample metadata with `SampleID` column)
(Optional) Pre-trained classifier `.qza` or reference sequences/taxonomy files for custom training
LEfSe Python scripts listed above

Usage

1. Launch the Script
   Run the script from a terminal:

   bash
   ./amca_qiime2_gui.sh
   

2. Follow Zenity Prompts
   The pipeline will walk you through:

   Selecting the working directory
   Performing host read depletion (optional)
   Importing and denoising sequences
   Running phylogenetic and diversity analysis
   Performing taxonomic classification
   Preparing and running LEfSe analysis

3. Review Output
   Results are saved in the working directory, including `.qza`, `.qzv`, and `.tsv` outputs, LEfSe plots, and classifier files.

Output Summary

QIIME2 artifacts and visualizations
	demux.qza, demux.qzv, rep-seqs.qza, table.qza, etc.
Diversity and phylogeny
	core-metrics-results, alpha-rarefaction.qzv, etc.
Taxonomy
	taxonomy.qzv, taxa-bar-plots.qzv
LEfSe
	lefse_input.tsv, lefse_results.txt, lda_scores.png, cladogram.png

Notes
The script automatically skips previously created output directories to avoid overwriting.
The GUI-based design is ideal for non-technical users or teaching purposes.

Citation
If you use this pipeline in your research, please cite:
Rana, N. et al. (2025). AMCA: A GUI-driven QIIME2 pipeline for microbiome analysis with integrated host depletion and biomarker discovery. (https://github.com/Nitika-Rana/AMCA/).

