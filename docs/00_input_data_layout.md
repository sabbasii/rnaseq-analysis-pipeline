### Input Data Layout

#### Overview

The pipeline expects all sequencing input files to be organized under the input_data/ directory.
This directory is not tracked by Git because it contains very large files (FASTQ, MD5 checksums).

```sh
input_data/
├── raw_fastq/
├── processed_fastq/
└── md5/
```

These paths correspond to variables defined in:

```text
workflow/00_config.sh

- `RAW_FASTQ_DIR` = input_data/raw_fastq
- `FASTQ_DIR` = input_data/processed_fastq
- `MD5_DIR` = input_data/md5

---

## raw_fastq/

### Purpose

Stores original sequencing files exactly as delivered by the sequencing facility.

Characteristics:

- Original filenames preserved exactly
- Contains instrument run metadata, barcode indices, and lane information
- Large compressed sequencing files (`.fastq.gz`)
- May include checksum files alongside FASTQs

Example:

_NS.LH00487_0037.005.NEBNext_dual_i7_100---NEBNext_dual_i5_100.Box1_A4_R1.fastq.gz_
_NS.LH00487_0037.005.NEBNext_dual_i7_100---NEBNext_dual_i5_100.Box1_A4_R2.fastq.gz_


### Notes

- These files are considered immutable reference copies.
- Never rename, modify, or delete files in this directory.
- Pipeline scripts do not directly use this directory.
- This directory serves as the permanent archive of sequencing delivery.

---

## Converting raw_fastq → processed_fastq

### Purpose

Create standardized FASTQ filenames suitable for pipeline workflows while preserving original files unchanged.

Sequencing facilities generate filenames containing run identifiers and barcode metadata. The pipeline requires simplified filenames based on biological sample position.

Example transformation:

_NS.LH00487_0037.005.NEBNext_dual_i7_100---NEBNext_dual_i5_100.Box1_A4_R1.fastq.gz_
→
_Box1_A4_R1.fastq.gz_

The standardized files are stored in:
```text
input_data/processed_fastq/


---

### Script used

Renaming is performed using:
```text
workflow/00_rename_fastq.sh

This script:

- Reads files from `$RAW_FASTQ_DIR`
- Extracts the BoxPosition identifier
- Copies files into `$FASTQ_DIR`
- Preserves raw files unchanged
- Ensures compatibility with downstream pipeline steps

---

`processed_fastq/`

Purpose:
Stores standardized FASTQ files used directly by pipeline workflows.

Naming Convention:

```sh
<BoxPosition>_R1.fastq.gz
<BoxPosition>_R2.fastq.gz
```

**Source:**

Files in this directory are generated from:
```text
input_data/raw_fastq/

using:
```text
workflow/00_rename_fastq.sh

**Important**  
All pipeline QC, alignment, and counting steps operate exclusively on:
```text
input_data/processed_fastq/

This directory is the primary input location for the pipeline.

---

`md5/`

Purpose:
Stores checksum files for verifying integrity of sequencing delivery.

Example:

*.fastq.gz.md5

---

### Relationship to Metadata

FASTQ filenames must match entries in:
```text
metadata/samples.tsv

This file defines mapping between sequencing data and biological samples.

---

### Data Flow

Sequencing Facility  
↓  
input_data/raw_fastq/ *(immutable archive)*  
↓  
workflow/00_rename_fastq.sh  
↓  
input_data/processed_fastq/ *(pipeline input)*  
↓  
QC → Alignment → Counting → Analysis

---
### Data Flow

<div align="center">

<span style="font-size: 20px; line-height: 1.8;">

Sequencing Facility  
↓  
input_data/raw_fastq/ *(immutable archive)*  
↓  
workflow/00_rename_fastq.sh  
↓  
input_data/processed_fastq/ *(pipeline input)*  
↓  
QC → Alignment → Counting → Analysis  

</span>

</div>