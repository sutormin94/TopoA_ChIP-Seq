# TopoA_ChIP-Seq
Analysis of TopoA binding sites across *E. coli W3110* genome

This repository contains a set of bash and python scripts which have been used for ChIP-Seq data analysis and visualization.


## ChIP-Seq_analysis_pipeline_example.sh

Shell script that makes initial QC of sequencing data, followed by trimming and filtration procedure. 
After post-trimming QC, processed reads are mapped to the reference genome, producing SAM-files which are
converted to BAM, sorted and indexed. Additionally, coverage depth is computed for initial BAm files and 
for ones after removal of PCR-duplicates.

**Requirements:** factqc, trimmomatic, bwa, samtools (1.9 or higher), sra-toolkit, shell

**Input:** Raw reads files (FASTQ), Genome file (FASTA)

**Output:** FastQC reports, SAM files, sorted and indexed BAM files, BED files (coverage depth)


## Bed_to_wig_converter.py

Script takes BED files with coverage depth and converts them to WIG format.

**Requirements:** Python 2 or 3

**Input:** BED files

**Output:** WIG files


## Compute_fold_enrichment.py

Takes two WIG files (for IP and Mock control) and computes by-position fold enrichment. 

**Requirements:** Python 2 or 3

**Input:** WIG files (IP and Mock control)

**Output:** WIG file with Fold Enrichment


## Average_wig_files.py

Takes a set of WIG files (organized as a dictionary) and computes by-position average WIG.

**Requirements:** Python 2 or 3

**Input:** WIG files to be averaged

**Output:** Averaged WIG file


## Return_sequences_under_peaks.py

Takes output of MACS2 for peaks called (NarrowPeak intervals), returns sequences under the peaks as a MFA file,
plots distribution of peaks GC-content in comparison to genome GC-content and the distribution of peaks widths.

**Requirements:** Python 3

**Input:** Peaks coordinates (NarrowPeaks), reference genome (FASTA)

**Output:** MFA with sequences under the peaks, plots


## Normalize_calc_FE.py

Script for advanced data analysis. Not completed yet.