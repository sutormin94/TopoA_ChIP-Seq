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

Takes a set of WIG files (organized as a dictionary) and computes by-position average WIG. Returns correlation matrix 
(file and heatmap) of pair-wize correlations between genomic tracks.

**Requirements:** Python 2 or 3

**Input:** WIG files to be averaged

**Output:** Averaged WIG file, correlation matrix (CSV), correlation matrix heatmap (PNG), clusterized correlation matrix (CSV), clusterized correlation matrix heatmap (PNG)


## Correlate_genome_tracks.py

Takes a set of WIG files (organized as a dictionary) and computes a correlation matrix 
(file and heatmap) of pair-wize correlations between genomic tracks. Creates heatmaps in SFig. 1.

**Requirements:** Python 2 or 3

**Input:** WIG files to be averaged

**Output:** correlation matrix (CSV), correlation matrix heatmap (PNG), clusterized correlation matrix (CSV), clusterized correlation matrix heatmap (PNG)


## Return_sequences_under_peaks_GC_width_FE_distribution.py

Takes output of MACS2 for peaks called (NarrowPeak intervals), returns sequences under the peaks as a MFA file,
plots distribution of peaks GC-content in comparison to genome GC-content, distribution of peaks widths, 
and distribution of peaks fold enrichment. Creates histograms in SFig. 3.

**Requirements:** Python 3

**Input:** Peaks coordinates (NarrowPeaks), reference genome (FASTA)

**Output:** MFA with sequences under the peaks, plots (GC% histogram, peaks width histogram, peaks FE histogram)


## Return_reproducible_peaks.py

Takes a dictionary of narrowPeak files with peaks called by MACS2 for different biological replicas.
Identifies reproducible regions and writes them as a narrowPeak file. Creates heatmaps in SFig. 2.

**Requirements:** Python 2 or 3

**Input:** Peaks coordinates (NarrowPeaks), reference genome (FASTA)

**Output:** Reproducible peaks coordinates (narrowPeak), heatmap (number of shared peaks), heatmap (Jaccardian distance between peak sets), Venn diagram (peak sets overlap)


## FE_for_peaks.py

Takes a dictionary of narrowPeak files with peaks called by MACS2 for different biological replicas or experimental conditions.
Returns average fold enrichment of regions as computed by Compute_fold_enrichment.py and Average_wig_files.py in broadPeak file format.

**Requirements:** Python 2 or 3

**Input:** Peaks coordinates (NarrowPeaks), reference genome (FASTA), WIG files with FE info

**Output:** Peaks with fold enrichment (broadPeak)


## Return_unique_peaks.py

Takes a dictionary of narrowPeak files with peak regions.
Identifies unique regions which do not overlap with other peaks, writes them as a narrowPeak file.
Creates heatmaps in SFig. 2.

**Requirements:** Python 2 or 3

**Input:** Peaks coordinates (NarrowPeaks), reference genome (FASTA)

**Output:** Unique peaks coordinates (narrowPeak), heatmap (number of shared peaks), heatmap (Jaccardian distance between peak sets), Venn diagram (peak sets overlap)


## Peak_overlap_simulation.py

Takes two BroadPeak files of real ChIP-Seq peaks datasets and tests the significance of sets overlay. Uses Monte-Carlo approach to 
place peaks randomly and then samples the overlay. Repetitive sampling gives a distribution for the number of peaks in overlay, which
is then approximated by normal distribution and is used to test the number of peaks in overlay for real datasets.
Creates distibutions in SFig. 6; SFig. 14E,F.

**Requirements:** Python 3

**Input:** Two BroadPeak files with peaks coordinates, genome length, portion of a genome to mask (deletions, etc.)

**Output:** Plot of a distribution of the number of overlaying peaks, p-value, text file with the number of overlaying peaks for every iteration.


## Violin_plots_EcTopoI_vs_RNAP_vs_Gyrase.py

The script tests sets of genomic intervals (Peaks, TUs, BIMEs-1, BIMEs-2, IHF sites, Fis sites, H-NS sites, MatP sites, etc.)
for the enrichment with some continously distributed character CDC (RNApol fold enrichment, score, GC%, etc.) (t-test). 
Makes violin-plots of CDS in intervals vs other sites. Creates violin plots Fig. 2B, F; Fig. 4E; SFig. 5E,F; SFig. 18C,D.

**Requirements:** Python 2 or 3

**Input:** Peaks coordinates (NarrowPeaks), continously distributed character (WIG)

**Output:** Pearson correlation, violin plots


## Regions_features_association.py

The script tests if some continously distributed characters (RNApol fold enrichment, score, GC%, etc.) 
are correlated for a set of genomic intervals (TUs, Peaks, etc.). Creates scatter plots SFig. 5A, B.

**Requirements:** Python 2 or 3

**Input:** Peaks coordinates (NarrowPeaks), continously distributed character (WIG)

**Output:** Pearson correlation, plots CDS1 vs CDS2, violin plots


## FE_over_US_GB_DS.py

Takes wig tracks of different continuous genomic features (GC%, MukB ChiP-Seq, etc.). Computes signal over TUs upstream (US),
downstream (DS) and over TUs bodies (GB). 

**Requirements:** Python 3

**Input:** Files with signal data (WIG), genome annotation (GFF or BroadPeak), regions to be omitted (BroadPeak)

**Output:** WIG files with metagene signal over all TUs, TAB files with average signal for each of TUs, plot of average signal over all TUs, histogram of the signal over TUs


## Plot_signal_over_transcription_units.py

Takes signal over transcriptions units in WIG format generated by FE_over_US_GB_DS.py, plots metagene signal over TUs upstream, downstream and TU bodies.
Creates plots SFig. 5G,H; SFig. 13; SFig. 18E-H; SFig. 19B; SFig. 21II; Fig. C,G,D,H; Fig. 3A,B; Fig. 4D; Fig. 9A,B.

**Requirements:** Python 3

**Input:** WIG files with metagene signal over TUs sets

**Output:** Plots with metagene signal


## FE_over_US_GB_DS_binning_and_statistics.py

Takes wig tracks of different continuous genomic features (GC%, MukB ChiP-Seq, etc.). Computes signal over TUs upstream (US),
downstream (DS) and over TUs bodies (GB). Script is dedicated for noisy data. 
To handle the data it performes data binning for smoothing. Also it keeps signal data for all TUs to contstruct confidential interval.
To plot different features together (e.g. ChIP-Seq data from different expriments or conditions) it normalizes the data.
Creates plots Fig. 2J, K; Fig. 3C; Fig. 2I; Fig. 3D

**Requirements:** Python 3

**Input:** Files with signal data (WIG), genome annotation (GFF or BroadPeak), regions to be omitted (BroadPeak), size of a bin

**Output:** WIG files with metagene signal over all TUs, TAB files with average signal for each of TUs, 
metagen plot of average signal over all TUs, bar plot of average signal in US, DS regions and in TUs bodies


## Compare_signal_US_TSS_GB

Takes wig tracks of different continous genome features (TopoI ChIP-Seq FE, GC%, etc.). Computes relative (X-mean(X)) and normalized ((X-mean(X))/std(X)) signals.
Compares mean signal at upstream (US), transcription start site (TSS), and in transcription unit body (TUB) regions for different ChIP-Seq conditions.
Makes metagene plots without TU body scaling and with signal statistics.

**Requirements:** Python 3

**Input:** WIG files with FE data, TUs annotation (GFF)

**Output:** Barplots comparing US, TSS, TUB; Metagene plot with signal statistics


## FE_over_intergenic_regions.py

Takes wig tracks of different genome features (GC%, MukB ChiP-Seq, etc.). Computes signal in intergenic regions (IGRs, provides as BroadPeak or GFF file), in respect to 
orientation of adjacent genes: ++, +-, -+ or --. 

**Requirements:** Python 3

**Input:** Files with signal data (WIG), IGGs annotation (GFF or BroadPeak), regions to be omitted (e.g. deletions, BroadPeak)

**Output:** WIG files with metagene signal over IGRs classified by adjacent genes orientation, TAB files with average signal for each of IGR, 
plot of average signal over IGR groups, histogram of the signal over IGRs


## Plot_signal_over_intergenic_regions.py

Takes signal over intergenic regions in WIG format generated by FE_over_intergenic_regions.py, plots metagene signal over IGRs.
Creates plot Fig. 4C.

**Requirements:** Python 3

**Input:** WIG files with metagene signal over IGR sets

**Output:** Plots with metagene signal


## Intergenic_regions_promoters_TFs_and_EcTopoI.py

1) Takes set of genes with assigned expression level and returns intergenic regions (IGRs).
Select by length intergenic regins appropriate for further analysis: 50bp<length<1000bp. Creates SFig. 10.
2) Searches for annotated promoters and transcription factor sites in selected IGRs. 
Annotation is taken from RegulonDB.
3) Adds signal info (wig files, e.g. EcTopoI ChIP-Seq FE) to the selected IGRs.
Adds information on whether adjacent genes encode proteins targeting to membrane (from Ecocyc database).
4) Adds information on whether adjacent genes encode proteins targeting to membrane (from PSORT database)
5) Remove anomalous IGRs (near dps gene with extremely high EcTopoI peak or near rRNA genes).
6) Compares EcTopoI FE (RNAP FE, Expression level) between IGRs sets defined by different features.
Makes violin-plots: Fig. 4A, SFig. 10, SFig. 11, SFig. 12, SFig. 19A.

**Requirements:** Python 3

**Input:** TUs annotation (TAB), reference genome (FASTA), promoter sequences (FASTA), TF sites sequences (FASTA), 
files with signal data (WIG), data on genes encoding membrane-targeting protein (TAB)

**Output:** Combine data on IGRs (TAB), violin-plots describing different sets of IGRs


## Run_ChIPMunk.bat

Runs ChipMunk several times to detect motifs in sequences under ChIP-Seq peaks 
(fasta files generated by Return_sequences_under_peaks_GC_width_FE_distribution.py)

**Requirements:** ChipMunk, windows command line

**Input:** MFA with sequences under the peaks

**Output:** ChIPMunk log-files, ChIPMunk output files


## Make_mfa_from_ChIPMunk_output.py

Take ChipMunk output and creates MFA with sequences detected as motif. Creates Logo of the motif using Weblogo
implementation in python. Creates plots in SFig. 15, Fig. 6A.

**Requirements:** Python 3

**Input:** ChIPMunk output

**Output:** MFA with sequences detected as motif, Logo of the motif


## Bar_plot_maker_CFU_Dot_blot_Pull_down_data_quantification.py

Creates barplots for CFU quantification experiments, Dot-blot quantified signals, Pull-down quantified signals.
Creates bar plots in SFig. 7; SFig. 8D; SFig. 17A; SFig. 22.

**Requirements:** Python 3

**Input:** Excell table with raw numbers

**Output:** Bar plots


## Growth_curve_plot.py

Takes culture growth curve data in a form of xlsx table containing series of time-points and corresponding optical density of a culture. 
Plots raw data growth curve and confidential intervals if several replicates provided.
Creates Fig. 7A,D; Fig. 14B.

**Requirements:** Python 3

**Input:** Excell table with time-points and optical density measurments

**Output:** Growth curve plot


## Microscopy_data_analysis.py

Takes raw cell length quantification data and plots distributions, compares datasets means.
Creates Fig. 7B,E.

**Requirements:** Python 3

**Input:** Excell table with cell length measurements

**Output:** Cell length histograms


## MST_curve_plot.py

Plots MST binding curves.
Creates Fig. 6C.

**Requirements:** Python 3

**Input:** Excell table with MST data

**Output:** Plot with binding curves


## qPCR_data_visualization.py

Calculates and plots primers efficiency (for qPCR), plots qPCR data.
Creates SFig. 16A; SFig. 17B; SFig. 23.

**Requirements:** Python 3

**Input:** Excell table with raw qPCR data

**Output:** Primer efficiency plots, qPCR data


## WIG_to_Circos.py

Converts wig file with some genomic feature to Circos-compatible format (bed-like, for Circos).

**Requirements:** Python 3

**Input:** WIG file, bin size, chromosome ID

**Output:** Circos-compatible file




