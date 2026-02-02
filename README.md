# TagmentationAnalysis_Rotation
This repsitory contains a modular pipeline and Jupyter notebook for analyzing tagmentation data. It supports quality control, trimming, mapping, insertion-site calling, batch processing, sequence-logo generation, and exploratory visualization.

These tools are designed to be adaptable to different experimental setups and reference genomes. 

---

## Repository Contents

### Scripts

**`01302026_filtered_map_junction_insertions.py`**  
Single-sample end-to-end pipeline:
- FASTQ read-level QC using Phred scores
- Donor/transposon trimming with `cutadapt`
- Mapping with `bwa mem`
- BAM sorting and indexing with `samtools`
- Per-read insertion coordinate calling
- Optional site-level filtering by minimum read count

Outputs per-read TSVs and (optionally) filtered TSVs.

---

**`bulk_from_excel_map_junction_insertions.py`**  
Batch automation of the single-sample pipeline using an Excel sheet.  
Additional behavior:
- Automatically removes samples with zero QC-passing reads
- Automatically removes samples with no insertion sites after filtering
- Produces:
  - `filtered_tsv_manifest.xlsx` (kept samples)
  - `removed_samples.xlsx` (discarded samples + reasons)

---

**`logo_from_manifest.py`**  
Generates insertion-centered sequence logos and enriched k-mer tables from TSV outputs:
- Extracts ±N bp windows from the reference genome
- Reverse-complements minus-strand windows for consistent orientation
- Produces frequency or information-content (“bits”) logos
- Outputs enriched k-mers (default: 6-mers)

---

### Notebook

**`01292026_TagmentationAnalysis.ipynb`**  
Exploratory visualization and comparison of replicate TSVs:
- Genome-wide insertion density plots
- Replicate overlays
- Strand bias summaries
- Hotspot tables
- Optional sgRNA binding site highlighting

---

## Features

- Read-level QC using Phred score thresholds
- Perfect-match genome filtering (`NM == 0`)
- Strand-aware insertion coordinate definition
- Site-level noise removal by read support
- Batch processing via Excel manifests
- Automatic removal of low-quality samples
- Sequence-logo generation with orientation correction
- Enriched k-mer discovery
- Genome-wide density plots and hotspot identification

---
## Getting Started

### Installation

Clone this repository and create the conda environment:

```bash
conda create -n integration_env -c conda-forge -c bioconda python=3.10 pandas numpy matplotlib jupyter openpyxl pysam logomaker biopython cutadapt bwa samtools
conda activate integration_env

## Launch the Notebook
jupyter notebook

If you encounter an architecture-related or user-site conflict error:

PYTHONNOUSERSITE=1 conda run -n integration_env jupyter notebook

Then navigate to this repository and open:

01292026_TagmentationAnalysis.ipynb

Usage
Single Sample Pipeline
python 01302026_filtered_map_junction_insertions.py \
  --fastq sample.fastq \
  --ref reference.fasta \
  --outdir output_dir \
  --donor-seq ACTG... \
  --donor-side 5p \
  --min-read-count 100 \
  --min-mapq 30

Batch Processing from Excel

Excel sheet should include columns such as:

fastq

ref

outdir

donor-seq

donor-side

min-read-count

min-mapq

qc-lowq-thresh

qc-max-lowq-frac

Run:

python bulk_from_excel_map_junction_insertions.py \
  --xlsx inputs.xlsx \
  --sheet Sheet1 \
  --report-dir bulk_output

Generate Sequence Logos
python logo_from_manifest.py \
  --manifest-xlsx bulk_output/filtered_tsv_manifest.xlsx \
  --outdir logos \
  --logo-type bits \
  --logo-name insertion_logo


Outputs include:

Base frequency tables

PNG/PDF logo figures

Enriched k-mer tables

Dependencies

Installed automatically via the conda environment:

python 3.10

pandas

numpy

matplotlib

jupyter

openpyxl

pysam

logomaker

biopython

cutadapt

bwa

samtools

Notes

Reference FASTA record IDs must match the ref column in TSV outputs.

Large FASTA files are loaded into memory during logo generation.

Ensure consistent coordinate systems (0-based vs 1-based) when integrating with other tools.

External tools (bwa, samtools, cutadapt) must be available in PATH via the conda environment.
