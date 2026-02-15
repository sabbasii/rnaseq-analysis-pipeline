## Reference Genome Setup

This pipeline requires a reference genome FASTA, gene annotation GTF, and a STAR index. These files are not stored in the repository because they are very large (10–40 GB) and must be generated locally to ensure reproducibility.

Instead, follow the steps below to download and prepare the reference files.

## Overview

We will download:

| Component      | Description                     | Size   |
|----------------|---------------------------------|--------|
| Genome FASTA  | DNA sequence of the genome     | ~3 GB  |
| Annotation GTF| Gene annotation                | ~200 MB|
| STAR index    | Alignment index (generated locally) | ~30 GB |

## Reference Genome

Reference used in this pipeline:

| Field      | Value              |
|------------|--------------------|
| Species    | _Mus musculus_    |
| Assembly   | GRCm39            |
| Source     | NCBI RefSeq       |
| Accession  | GCF_000001635.27 |

**Official source:**

https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/

---

## Step 1 — Create Reference Directories

Run from the repository root:

```sh
mkdir -p reference/genome
mkdir -p reference/star_index
```

Resulting structure:

```text
reference/
├── genome/
└── star_index/

---

## Step 2 — Download Genome FASTA

Navigate to the genome directory:
```bash
cd reference/genome

Download the FASTA file:

```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
```

Uncompress the file:

You should now have:

```text
reference/genome/GCF_000001635.27_GRCm39_genomic.fna
```

---

## Step 3 — Download Annotation GTF

Still inside:

reference/genome/

Download the annotation file:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz

```

Uncompress the file:
```text
gunzip GCF_000001635.27_GRCm39_genomic.gtf.gz

You should now have:
```text
reference/genome/GCF_000001635.27_GRCm39_genomic.gtf

---

## Step 4 — Verify Downloads

Run:

```bash
ls -lh reference/genome
```

Expected output:
```text
GCF_000001635.27_GRCm39_genomic.fna
GCF_000001635.27_GRCm39_genomic.gtf

---

Step 5 — Build STAR index

The STAR index enables fast RNA-seq alignment.

Run from repository root:

```bash
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir reference/star_index \
  --genomeFastaFiles reference/genome/GCF_000001635.27_GRCm39_genomic.fna \
  --sjdbGTFfile reference/genome/GCF_000001635.27_GRCm39_genomic.gtf \
  --sjdbOverhang 100
```

### Explanation

| Parameter              | Meaning                                      |
|------------------------|----------------------------------------------|
| `--runThreadN`        | Number of CPU threads to use                |
| `--genomeDir`         | Output directory where STAR index is saved  |
| `--genomeFastaFiles`  | Path to the genome FASTA file               |
| `--sjdbGTFfile`       | Path to the gene annotation GTF file        |
| `--sjdbOverhang`      | Read length minus 1 (use 100 for 101 bp reads) |

---

## Step 6 — Verify STAR index build

After STAR finishes, verify that the index was created successfully.

Run:

```bash
ls -lh reference/star_index
```

You should see output similar to:
```text
Genome
SA
SAindex
chrLength.txt
chrName.txt
chrStart.txt
genomeParameters.txt
sjdbInfo.txt
sjdbList.out.tab

## What These Files Are

These files together form the STAR genome index required for RNA-seq alignment. STAR uses them to rapidly locate where each read belongs in the genome.

### Brief Description of Key Files

| File                  | Description                                                |
|-----------------------|------------------------------------------------------------|
| `Genome`             | Binary representation of the genome sequence              |
| `SA`                 | Suffix array used for fast alignment search              |
| `SAindex`            | Index of the suffix array for rapid access               |
| `chrName.txt`        | Names of chromosomes in the genome                       |
| `chrLength.txt`      | Length of each chromosome                                |
| `chrStart.txt`       | Starting position of each chromosome in the index        |
| `genomeParameters.txt` | STAR parameters used during index generation             |
| `sjdbInfo.txt`       | Splice junction database information                     |
| `sjdbList.out.tab`   | List of splice junctions extracted from annotation       |

## Expected Directory Structure

Your `reference` directory should now look like:

```text
reference/
├── genome/
│   ├── GCF_000001635.27_GRCm39_genomic.fna
│   └── GCF_000001635.27_GRCm39_genomic.gtf
│
└── star_index/
    ├── Genome
    ├── SA
    ├── SAindex
    ├── chrLength.txt
    ├── chrName.txt
    ├── chrStart.txt
    ├── genomeParameters.txt
    ├── sjdbInfo.txt
    └── sjdbList.out.tab

Quickly check:

```bash
du -sh reference/star_index
```

Expected size:
```text
~25–35 GB

🟢 The STAR index is complete and can be used by downstream alignment steps.