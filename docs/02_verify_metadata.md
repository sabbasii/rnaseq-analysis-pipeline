## Verify Metadata Integrity and Consistency

Run the following commands to inspect and validate the metadata file and ensure consistency with FASTQ files:

```bash
# Check metadata file exists and view first lines
head metadata/samples.tsv

# Count total number of samples in metadata
wc -l metadata/samples.tsv

# Check unique experimental groups
cut -f3 metadata/samples.tsv | sort | uniq

# Check unique time points
cut -f4 metadata/samples.tsv | sort | uniq

# Verify specific sample exists
grep Box3_B7 metadata/samples.tsv

# Count FASTQ samples (R1 files)
ls input_data/processed_fastq/*_R1.fastq.gz | wc -l

# Compare with metadata count
wc -l metadata/samples.tsv

### Expected result:
The number of FASTQ samples should match the number of entries in metadata/samples.tsv, confirming that each sequencing file has a corresponding metadata record.