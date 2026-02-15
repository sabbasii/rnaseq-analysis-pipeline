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

`raw_fastq/`

Purpose:
Stores original sequencing files exactly as delivered by the sequencing facility.

- Original filenames preserved
- Contains run metadata, barcode indices, lane information
- Large compressed FASTQ files (.fastq.gz)
- MD5 checksum files may exist alongside FASTQs

Example:

```sh
NS.LH00487_0037.005.NEBNext_dual_i7_100---NEBNext_dual_i5_100.Box1_A4_R1.fastq.gz
NS.LH00487_0037.005.NEBNext_dual_i7_100---NEBNext_dual_i5_100.Box1_A4_R2.fastq.gz
```

Notes
+ These files are considered immutable reference copies.
+ Pipeline scripts do not directly use this folder for alignment.

---

`processed_fastq/`

Purpose:
Stores standardized FASTQ files used directly by pipeline workflows.

Naming Convention:  

```sh
<BoxPosition>_R1.fastq.gz
<BoxPosition>_R2.fastq.gz
```

Source:

- Generated from raw_fastq/ using renaming scripts
- May include manually recovered or re-downloaded samples

**Important**   
All pipeline alignment and QC steps operate on this directory.

---

`md5/`

Purpose:
Stores checksum files for verifying integrity of sequencing delivery.

Example:

*.fastq.gz.md5

---

### Relationship to Metadata

File naming must match entries in:

```sh
metadata/samples.tsv
```

Where:

- `BoxPosition` = FASTQ filename prefix
- Used to map sequencing reads → biological sample → experimental group → timepoint