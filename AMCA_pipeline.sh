#!/bin/bash
set -euo pipefail

LOGFILE="qiime2_pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

command -v zenity >/dev/null 2>&1 || { echo "Zenity is not installed."; exit 1; }

zenity --info --title="AMCA GUI" --text="Welcome to the QIIME2 based AMCA Microbiome Analysis GUI Workflow.\n\nThis pipeline will guide you through:\n\n\n1. Host depletion (optional)\n\n2. Importing FASTQ and metadata\n\n3. DADA2 denoising\n\n4. Diversity and phylogenetic analysis\n\n5. Taxonomic analysis with Classifier training (optional)\n\n6. LEfSe preprocessing and analysis\n\n\nClick OK to start."

zenity --info --title="Working Directory" --text="\n\nAll necessary files should be in working directory:\n\n\n1.  .fastq raw sequence data\n\n2.  A metadata file with name:metadata.tsv\n\n3.  Custom taxonomic classifier if needed, \n\n      OR reference files to generate a classifier\n\n4.  All LEfSe python scripts\n\n\nNote: Move any previous output directories if run previously performed. The pipeline does not overwrite existing folders and exits if output directories are already present in the working directory.\n\n\nClick OK to start"

# Step 1: Select working directory
workdir=$(zenity --file-selection --directory --title="Select Working Directory")
cd "$workdir" || { zenity --error --text="Invalid directory"; exit 1; }

# Step 2: Host read depletion
host_depleted=false
if zenity --question --title="Host Read Depletion" --text="Would you like to perform host read removal using Bowtie2?\n\nThe program currently performs host-depletion for Human genome with the reference: GRCh38.\n\nThe program would itself download the reference genome pre-built indices."; then
    choice=$(zenity --list --radiolist --title="Bowtie2 Index Option" \
        --text="Choose Bowtie2 index source:" \
        --column="Select" --column="Option" TRUE "Download GRCh38 index" FALSE "Use custom index")

    if [[ "$choice" == "Download GRCh38 index" ]]; then
        zenity --info --text="Downloading GRCh38 Bowtie2 index..."
        wget -O GRCh38_noalt_as.zip https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip 
        unzip -o GRCh38_noalt_as.zip
        INDEX="GRCh38_noalt_as"
    else
        index_path=$(zenity --file-selection --directory --title="Select Bowtie2 Index Directory")
        cd "$index_path"
        INDEX=$(zenity --entry --title="Enter Index Base Name" --text="Enter prefix of .bt2 files:")
        if ! ls "${INDEX}".*.bt2 >/dev/null 2>&1; then
            zenity --error --text="Bowtie2 index files not found. Skipping host read depletion."
            cd "$workdir"
        else
            cd "$workdir"
        fi
    fi

    if [[ -n "${INDEX:-}" && -f "${INDEX}.1.bt2" ]]; then
        mkdir -p host_removed_fastq
        for R1 in *_1.fastq; do
            R2="${R1/_1.fastq/_2.fastq}"
            [[ ! -f "$R2" ]] && continue
            SAMPLE="${R1%_1.fastq}"
            bowtie2 -x "$INDEX" -1 "$R1" -2 "$R2" \
                --sensitive-local --reorder -S /dev/null \
                --un-conc host_removed_fastq/"${SAMPLE}_host_removed.fastq"
            mv host_removed_fastq/"${SAMPLE}_host_removed.fastq.1" "${SAMPLE}_host_removed_R1.fastq"
            mv host_removed_fastq/"${SAMPLE}_host_removed.fastq.2" "${SAMPLE}_host_removed_R2.fastq"
        done
        host_depleted=true
    fi
else
    zenity --info --text="Skipping host read depletion."
fi

# Step 3: Metadata check
[[ ! -f "metadata.tsv" ]] && { zenity --error --text="metadata.tsv not found"; exit 1; }

# Step 4: Manifest creation
zenity --info --title="Step 1: Manifest" --text="Creating manifest.csv from paired-end FASTQ files...\n\nManifest file automates loading of .fastq raw data as a QIIME2 artefact.\n\nPlease make sure the raw .fastq files are in the current working directory"

manifest="manifest.csv"
echo "sample-id,absolute-filepath,direction" > "$manifest"

if [[ "$host_depleted" == true ]]; then
    for file in *_host_removed_R1.fastq; do
        [[ ! -f "$file" ]] && continue
        sample_id=$(basename "$file" "_host_removed_R1.fastq")
        echo "$sample_id,$PWD/$file,forward" >> "$manifest"
    done
    for file in *_host_removed_R2.fastq; do
        [[ ! -f "$file" ]] && continue
        sample_id=$(basename "$file" "_host_removed_R2.fastq")
        echo "$sample_id,$PWD/$file,reverse" >> "$manifest"
    done
else
    for f in *_1.fastq; do
        sid=$(basename "$f" "_1.fastq")
        echo "$sid,$PWD/$f,forward" >> "$manifest"
    done
    for f in *_2.fastq; do
        sid=$(basename "$f" "_2.fastq")
        echo "$sid,$PWD/$f,reverse" >> "$manifest"
    done
fi
zenity --info --title="Manifest file created." --text="\nFile generated: manifest.tsv"

# Step 5: Import and demux
(
echo "10"
echo "# Importing FASTQ files..."
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path "$manifest" --output-path demux.qza --input-format PairedEndFastqManifestPhred33
echo "# File generated: demux.qza"
echo "20"
echo "# Summarizing imported data..."
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
echo "# File generated: demux.qzv\nFile can be visualized at QIIME2View interface on your web browser using https://view.qiime2.org/"
) | zenity --progress --title="QIIME2 Import" --text="Importing data..." --percentage=0 --auto-close
zenity --info --title="QIIME2 input artefact generated" --text="\nFile generated: \ndemux.qza\ndemux.qzv"

# Step 6: DADA2
params=$(zenity --forms --title="DADA2 Parameters" --text="Open QIIME2View interface on your web browser using https://view.qiime2.org/ and visualize demux.qzv.\n\nEnter trimming and truncation values from the 'Interactive Quality Plot' for demux.qzv in browser:\n\nHere \n1.'f' and 'r' stand for forward and reverse reads.\n2.'trim' stands for read length trimming and 'trun' stands for read length truncation.\n\n'0' as an entry is acceptable" --add-entry="trim-left-f" --add-entry="trim-left-r" --add-entry="trunc-len-f" --add-entry="trunc-len-r")
IFS="|" read -r trim_left_f trim_left_r trunc_len_f trunc_len_r <<< "$params"

(
echo "30"
echo "# Running DADA2..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f "$trim_left_f" --p-trim-left-r "$trim_left_r" \
  --p-trunc-len-f "$trunc_len_f" --p-trunc-len-r "$trunc_len_r" \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats-dada2.qza
echo "# Files generated: \n1)rep-seqs.qza\n2)table.qza\n3)stats-dada2.qza"
echo "40"
echo "# Tabulating DADA2 results..."
qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file metadata.tsv
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
echo "# Files generated: \n1)stats-dada2.qzv\n2)table.qzv\n3)rep-seqs.qzv"
) | zenity --progress --title="DADA2 Denoising" --text="Running DADA2..." --percentage=0 --auto-close
zenity --info --title="DADA denoising performed" --text="\nFile generated: \nrep-seqs.qza\ntable.qza\nstats-dada2.qza\nstats-dada2.qzv\ntable.qzv\nrep-seqs.qzv "

# Step 7: Phylogeny & diversity
zenity --info --title="Step 3: Phylogeny & Diversity" --text="Generating phylogenetic tree and running core diversity metrics..."
sampling_depth=$(zenity --entry --title="Sampling Depth for diversity and phylogenetic analysis" --text="Open QIIME2View interface on your web browser using https://view.qiime2.org/ and visualize table.qzv\n\nEnter sampling depth (check 'Interactive Sampling Detail' for table.qzv in browser)\n\nChoose count value so that samples with values less than the selected value are excluded from downstream analysis")
(
echo "60"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
echo "# Files generated: \n1)aligned-rep-seqs.qza\n2)masked-aligned-rep-seqs.qza\n3)unrooted-tree.qza\n4)rooted-tree.qza"
echo "70"
echo "# Diversity analysis ongoing..."
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza --i-table table.qza \
  --p-sampling-depth "$sampling_depth" \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results
echo "# Files generated: \nCore diversity analysis results stored in folder 'core-metrics-results'"
) | zenity --progress --title="Diversity Analysis" --text="Running diversity metrics..." --percentage=0 --auto-close
zenity --info --title="Phylogeny and Diversity performed" --text="\nFile generated: \naligned-rep-seqs.qza\nmasked-aligned-rep-seqs.qza\nunrooted-tree.qza\nrooted-tree.qza\nCore diversity analysis results stored in folder 'core-metrics-results' "

# Step 8: Optional beta significance
while true; do
  column=$(zenity --entry --title="Beta Diversity" --text="Enter metadata column (or 'next' to skip):")
  [[ "$column" == "next" ]] && break
  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column "$column" \
    --o-visualization core-metrics-results/${column}_significance.qzv || zenity --error --text="Failed for column: $column"
done

# Step 9: Alpha rarefaction
zenity --info --title="Computing alpha rarefaction" --text="Add value for the maximum_depth parameter"
max_depth=$(zenity --entry --title="Maximum Depth" --text="Open QIIME2View interface on your web browser using https://view.qiime2.org/ and visualize table.qzv\n\nEnter maximum depth (check Frequency per sample information for table.qzv in browser)\n\nChoose value based on the median frequency")

(
echo "80"
echo "# Running Alpha rarefaction..."
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth "$max_depth" \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
echo "90"
) | zenity --progress --title="Alpha Rarefaction" --text="Running alpha-rarefaction analysis..." --percentage=0 --auto-close

zenity --info --title="Alpha rarefaction performed" --text="\nFile generated:\nalpha-rarefaction.qzv"

# Step 10: Classifier
zenity --info --title="Computing Taxonomic Classification" --text="Analysis will be run on a pre-trained or custom-trained classifier."

if zenity --question --text="Do you want to use a pre-trained classifier?"; then
    while true; do
        classifier_path=$(zenity --file-selection --title="Select Classifier" --file-filter="*.qza")
        if [[ -n "$classifier_path" && -f "$classifier_path" ]]; then
            (
            echo "10"
            echo "# Running classification..."
            qiime feature-classifier classify-sklearn --i-classifier "$classifier_path" --i-reads rep-seqs.qza --o-classification taxonomy.qza
            qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
            qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.tsv --o-visualization taxa-bar-plots.qzv
            echo "100"
            ) | zenity --progress --title="Taxonomic Classification" --text="Classifying with pre-trained classifier..." --percentage=0 --auto-close
            break
        else
            zenity --question --text="Invalid classifier. Retry?" || break
        fi
    done
fi

if [[ -z "${classifier_path:-}" ]]; then
    (
    echo "20"
    echo "# Extracting reads for custom classifier..."
    qiime feature-classifier extract-reads \
      --i-sequences 2022.10.backbone.full-length.fna.qza \
      --p-f-primer CCTACGGGNGGCWGCAG \
      --p-r-primer GACTACHVGGGTATCTAATCC \
      --o-reads ref-seqs.qza

    echo "60"
    echo "# Training custom classifier..."
    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads ref-seqs.qza \
      --i-reference-taxonomy 2022.10.backbone.tax.qza \
      --o-classifier classifier.qza

    echo "# Running classification with custom classifier..."
    qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
    qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
    qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.tsv --o-visualization taxa-bar-plots.qzv
    echo "100"
    ) | zenity --progress --title="Custom Classifier Training & Taxonomy" --text="Training and classifying with custom classifier..." --percentage=0 --auto-close

    classifier_path="classifier.qza"
fi

zenity --info --title="Taxonomic Classification Done" --text="Files generated:\n- taxonomy.qza\n- taxonomy.qzv\n- taxa-bar-plots.qzv"

# Step 11: LEfSe Preparation
taxa_level=$(zenity --entry --title="LEfSe Pre-analysis" --text="Enter taxonomic level (1-7):\n\n1) Kingdom\n2) Phylum\n3) Class\n4) Order\n5) Family\n6) Genus\n7) Species")

(
echo "20"
echo "# Collapsing features to taxonomic level $taxa_level..."
qiime taxa collapse --i-table table.qza --i-taxonomy taxonomy.qza --p-level "$taxa_level" --output-dir exported-table

echo "50"
echo "# Exporting collapsed table..."
qiime tools export --input-path exported-table/collapsed_table.qza --output-path exported-table-with-taxonomy

echo "70"
echo "# Converting BIOM to TSV..."
biom convert -i exported-table-with-taxonomy/feature-table.biom -o exported-table-with-taxonomy/otu-table.tsv --to-tsv

echo "100"
) | zenity --progress --title="LEfSe Table Prep" --text="Preparing LEfSe table..." --percentage=0 --auto-close

zenity --info --title="LEfSe Table Ready" --text="Files generated:\n- exported-table/\n- exported-table-with-taxonomy/otu-table.tsv"

# Step 12a: Format for LEfSe
zenity --info --title="Reformatting OTU Table" \
  --text="Preparing LEfSe input table from exported OTU table and metadata..."

# Step 12b: Select OTU cleaned TSV
tail -n +2 exported-table-with-taxonomy/otu-table.tsv > cleaned_otu_table.tsv
otu_file=$(zenity --file-selection --title="Select OTU cleaned TSV file" --file-filter="TSV files | *.tsv")
[ -z "$otu_file" ] && zenity --error --text="No OTU file selected." && exit 1

# Step 12c: Select metadata TSV
meta_file=$(zenity --file-selection --title="Select metadata TSV file" --file-filter="TSV files | *.tsv")
[ -z "$meta_file" ] && zenity --error --text="No metadata file selected." && exit 1

# Step 12d: Extract column names from metadata (except SampleID)
columns=$(head -n 1 "$meta_file" | tr '\t' '\n' | grep -v '^#\?SampleID$')

# Prompt user to select Class column
lefse_col=$(zenity --list --title="Select Class Column" --column="Metadata Columns" $columns)
[ -z "$lefse_col" ] && zenity --error --text="No class column selected." && exit 1

# Step 12e: Run Python preprocessing
(
echo "10"
echo "# Cleaning and transposing OTU table for LEfSe..."

python3 - <<EOF
import pandas as pd
import sys

otu_file = "$otu_file"
meta_file = "$meta_file"
class_column = "$lefse_col"

# Load OTU table
otu = pd.read_csv(otu_file, sep='\t')
if otu.columns[0].strip() == '#OTU ID':
    otu.rename(columns={otu.columns[0]: 'SampleID'}, inplace=True)
elif otu.columns[0].strip().lower() == 'feature':
    otu.rename(columns={otu.columns[0]: 'SampleID'}, inplace=True)

# Transpose
otu.set_index('SampleID', inplace=True)
otu = otu.T
otu.index.name = 'SampleID'
otu.reset_index(inplace=True)

# Load metadata
meta = pd.read_csv(meta_file, sep='\t', dtype=str)
meta.columns = meta.columns.str.replace('^#', '', regex=True)

# Validate columns
if 'SampleID' not in meta.columns:
    sys.stderr.write("ERROR: 'SampleID' column not found in metadata.tsv\n")
    sys.exit(1)
if class_column not in meta.columns:
    sys.stderr.write(f"ERROR: Column '{class_column}' not found in metadata.tsv\n")
    sys.exit(1)

# Extract and rename class column
meta = meta[['SampleID', class_column]].rename(columns={class_column: 'Class'})

# Merge and reorder
merged = pd.merge(otu, meta, on='SampleID', how='left')
cols = ['SampleID', 'Class'] + [col for col in merged.columns if col not in ['SampleID', 'Class']]
merged = merged[cols]

# Save result
output_file = "lefse_input.tsv"
merged.to_csv(output_file, sep='\t', index=False)
EOF


echo "100"
) | zenity --progress --title="LEfSe Input Preparation" --text="Processing OTU + metadata tables..." --percentage=0 --auto-close

# Confirm output
if [ -f lefse_input.tsv ]; then
    zenity --info --title="LEfSe Input Ready" --text="LEfSe input file saved as:\nlefse_input.tsv"
else
    zenity --error --text="LEfSe input generation failed."
    exit 1
fi


# Step 13: LEfSe analysis
(
echo "20"
python lefse_format_input.py lefse_input.tsv formatted_input.txt -f c -c 2 -u 1
echo "50"
python lefse_run.py formatted_input.txt lefse_results.txt
echo "80"
python lefse_plot_res.py lefse_results.txt lda_scores.png --dpi 600
python lefse_plot_cladogram.py lefse_results.txt cladogram.png --dpi 600
echo "100"
) | zenity --progress --title="Running LEfSe" --text="Performing LEfSe analysis..." --percentage=0 --auto-close

zenity --info --title="LEfSe Analysis Done" --text="Files generated:\n- lefse_results.txt\n- lda_scores.png\n- cladogram.png"

# Step 14: Completion
zenity --info --title="Pipeline Complete" --text="Workflow complete.\nResults saved in Working directory"
