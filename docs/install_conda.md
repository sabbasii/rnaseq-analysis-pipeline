# Install Conda (Anaconda) on Linux

This document explains how to install Conda using the full Anaconda distribution.

Conda is used to manage software environments and ensure reproducible RNA-seq analysis.

---

## What is Conda?

Conda is an environment and package manager that allows you to:

- Install bioinformatics tools safely
- Avoid software version conflicts
- Create isolated environments
- Ensure reproducibility across systems

---

## Step 1 — Download Anaconda Installer

Navigate to your home directory:

```bash
cd ~

Download the installer:

```sh
wget https://repo.anaconda.com/archive/Anaconda3-2025.06-1-Linux-x86_64.sh
```

---
## Step 2 — Run the Installer

```sh
bash Anaconda3-2025.06-1-Linux-x86_64.sh
```

Follow the prompts:
+ Press Enter to scroll through the license
+ Type yes to accept the license
+ Press Enter to accept the default install location
+ Type yes when asked to initialize Conda

---
## Step 3 — Activate Conda

Reload your shell configuration:

```sh
source ~/.bashrc
```

---
## Step 4 — Verify Installation
Check Conda version:

```sh
conda --version
```

## Step 5 — Create the `rnaseq` Conda Environment

Create the environment:

```bash
conda create -n rnaseq -y

**Optional: Disable Automatic Base Environment Activation**

Run the following command:

```bash
conda config --set auto_activate_base false

Verify the setting

Run:

```bash
conda config --show auto_activate

Expected output:
```text
auto_activate: false

+ This ensures Conda environments are activated only when explicitly requested.

You can activate the RNA-seq environment when needed using:

```sh
conda activate rnaseq
```

---
### Step 5.1 — Install and Verify Core RNA-seq Tools

```sh
conda install -c bioconda -c conda-forge fastqc multiqc star subread -y
```

Verify installation:

```sh
fastqc --version
multiqc --version
STAR --version
featureCounts -v
```

These tools are used for:
+ FastQC → QC
+ MultiQC → QC summary
+ STAR → alignment
+ featureCounts → gene counting

### Step 5.2 — Install R and Core RNA-seq Analysis Packages

```sh
conda install -c bioconda -c conda-forge \
  r-base r-essentials \
  bioconductor-deseq2 bioconductor-edger bioconductor-limma \
  r-tidyverse r-data.table r-ggplot2 r-pheatmap -y
```

Verify:

```sh
R --version
```

---
### Step 5.3 — Install Python and Common Analysis Packages

```sh
conda install -c conda-forge \
  python pip \
  numpy pandas scipy matplotlib seaborn statsmodels \
  jupyterlab -y
```

Verify:

```sh
python --version
pip --version
```

---
## Step 6 — Remove Installer (Optional)

Clean up installer file:

```sh
rm Anaconda3-2025.06-1-Linux-x86_64.sh
```

---
## Installation Location

By default, Conda installs to:
```text
~/anaconda3/

---
### Troubleshooting

If Conda is not recognized:

```sh
source ~/.bashrc
```

**After completing this guide, Conda will be fully installed and ready for RNA-seq analysis.**