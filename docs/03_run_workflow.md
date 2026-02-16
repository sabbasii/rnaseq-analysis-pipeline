## Step 3 — Run the RNA-seq Workflow

This section explains how to run the RNA-seq processing workflow, from quality control through alignment and gene counting.

The workflow is organized under the `workflow/` directory and uses a shared configuration file to ensure reproducibility and consistent paths across all steps.

## Workflow Overview

The workflow consists of three main steps:

| Step | Script                         | Purpose                                   |
|------|---------------------------------|-------------------------------------------|
| 01   | `01_qc_fastqc_multiqc.sh`     | Quality control of FASTQ files           |
| 02   | `02_align_star.sh`            | Align reads to the genome using STAR     |
| 03   | `03_count_featurecounts.sh`   | Generate gene-level count matrix         |

## Shared Configuration

All workflow scripts use a shared configuration file:

```text
workflow/00_config.sh

This file defines all project paths and settings, ensuring consistency and reproducibility across all workflow steps.

## Workflow Directory Structure

```text
workflow/
├── 00_config.sh
├── 01_qc_fastqc_multiqc.sh
├── 02_align_star.sh
└── 03_count_featurecounts.sh

## Step 0 — Configuration (`00_config.sh`)

This file defines all core project paths and settings required by the workflow.


- Project root directory (`REPO_ROOT`)
- Input data locations
- Reference genome paths
- Output directories
- Thread count for parallel processing
- Helper validation functions

This configuration file is automatically sourced by all workflow scripts to ensure consistent and reproducible execution.

---


| Variable          | Description                          |
|-------------------|--------------------------------------|
| `INPUT_DIR`       | Base directory for all input data    |
| `FASTQ_DIR`       | Location of processed FASTQ files   |
| `META_DIR`        | Metadata directory                  |
| `REF_DIR`         | Reference genome directory          |
| `STAR_INDEX_DIR`  | STAR genome index directory         |
| `RESULTS_DIR`     | Base directory for workflow outputs |
| `BAM_DIR`         | Directory for alignment BAM files   |
| `COUNT_DIR`       | Directory for gene count outputs    |
| `THREADS`         | Number of CPU threads to use        |

---

**Important:**  
You do **NOT** run this file directly. It is automatically loaded by each workflow script.

## Step 1 — Create Output Directories and Make Scripts Executable

Run the following commands from the repository root:

```bash
# Create the `workflow/` and `results/` directories
mkdir -p workflow results
# Makes all workflow scripts executable so they can be run directly
chmod +x workflow/*.sh

---

## Step 2 — Run Quality Control

**Script:**

```text
workflow/01_qc_fastqc_multiqc.sh

This step performs:
+ FastQC on all FASTQ files
+ MultiQC summary report

**Run:**

```sh
bash workflow/01_qc_fastqc_multiqc.sh
```

**Output:**

```text
results/qc/fastqc/
results/qc/multiqc/
results/logs/01_qc_fastqc_multiqc.log

**Key report:**
```text
results/qc/multiqc/multiqc_report.html

Open this file in a web browser to inspect sequencing data quality.

---

## Step 3 — Align Reads Using STAR

**Script:**

```text
workflow/02_align_star.sh

This step:

+ Reads sample information from:

```text
metadata/samples.tsv

+ Aligns paired FASTQ files using STAR
+ Produces coordinate-sorted BAM files

**Run:**

```sh
bash workflow/02_align_star.sh
```

**Output:**

```text
results/align/bam/*.bam
results/align/star_logs/
results/logs/02_align_star.log

Each BAM file is named using the following format:

```text
MouseID_Group_TimePoint_Aligned.out.bam

This naming convention preserves key experimental information:
+ MouseID — unique identifier for each animal
+ Group — experimental condition (e.g., Sham, Stroke, NeuDep)
+ TimePoint — sampling time point
+ Aligned.out.bam — indicates STAR-aligned, coordinate-sorted BAM file

Example:

```text
527297-2_Sham_3hPost-Reperfusion_Aligned.out.bam

---

## Step 4 — Generate Gene Count Matrix

**Script:**

```text
workflow/03_count_featurecounts.sh

This step:
+ Reads all BAM files
+ Uses the genome annotation (.gtf)
+ Produces a gene-by-sample count matrix

**Run:**

```sh
bash workflow/03_count_featurecounts.sh
```

**Output:**

```text
results/counts/gene_counts.tsv
results/logs/03_count_featurecounts.log

**This file is the primary input for downstream differential expression analysis.**

---

## Full Run Order

Run the workflow scripts in the following order:

```sh
bash workflow/01_qc_fastqc_multiqc.sh
bash workflow/02_align_star.sh
bash workflow/03_count_featurecounts.sh
```

**This sequence performs quality control, read alignment, and gene count generation.**

## Output Directory Structure

After workflow completion, the `results/` directory will have the following structure:

```text
results/
├── qc/
│   ├── fastqc/
│   └── multiqc/
│
├── align/
│   ├── bam/
│   └── star_logs/
│
├── counts/
│   └── gene_counts.tsv
│
└── logs/

---

## Thread Control (Optional)

Default number of threads:

```bash
THREADS=16

Override the thread count at runtime:

```sh
THREADS=32 bash workflow/02_align_star.sh
```

This allows you to adjust CPU usage depending on available system resources.

---

## Requirements and Installation

### Required Software

The following tools must be installed before running the RNA-seq workflow:

- FastQC  
- MultiQC  
- STAR  
- featureCounts (part of the Subread package)  

---

## Installation Instructions

### Option 1 — Ubuntu / Debian

Install using the system package manager:

```bash
sudo apt update

# Install FastQC
sudo apt install -y fastqc

# Install MultiQC
sudo apt install -y multiqc

# Install STAR
sudo apt install -y star

# Install featureCounts (Subread package)
sudo apt install -y subread

### Option 2 — Conda (Recommended for reproducibility)
Conda provides isolated environments and ensures reproducible software versions.

If Conda is not installed, first follow the Conda installation guide:
```text
docs/install_conda.md

Then create the RNA-seq environment:

```sh
conda create -n rnaseq -y
conda activate rnaseq

conda install -c bioconda -c conda-forge fastqc multiqc star subread -y
```

## Verify Installation
Run the following commands:

```sh
fastqc --version
multiqc --version
STAR --version
featureCounts -v

```

Each command should print a version number without errors.