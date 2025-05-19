# 16S DADA2 to Phyloseq standard pipelines for MMs Lab

A collection of R scripts used to analyze various 16S rRNA amplicon data from projects within the MMs lab

- `dada2_raw_pipeline.R` is the standard pipeline used to process raw demultiplexed fastq files from Illumina sequencers. Use this if the sequencing data was all produced in a single run (most cases)

- `dada2_raw_multiple_runs.R` is the pipeline to process raw demultiplexed fastq files from Illumina sequencers. Use this if the sequencing data produced from multiple sequencing runs
