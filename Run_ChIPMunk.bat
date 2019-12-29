@echo off

echo chipmunk for noRif peaks 2..
java -jar C:\Sutor\Dists\chipmunk.jar s:C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Seq_under_peaks\EcTopoI_noCTD_noRif_rep123_nm_0.001_peaks.fasta 1>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_noRif_rep123_nm_0.001_peaks_motifs_2.txt 2>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_noRif_rep123_nm_0.001_peaks_log_2.txt
echo chipmunk for noRif peaks 3..
java -jar C:\Sutor\Dists\chipmunk.jar s:C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Seq_under_peaks\EcTopoI_noCTD_noRif_rep123_nm_0.001_peaks.fasta 1>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_noRif_rep123_nm_0.001_peaks_motifs_3.txt 2>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_noRif_rep123_nm_0.001_peaks_log_3.txt

echo chipmunk for Rif peaks 2..
java -jar C:\Sutor\Dists\chipmunk.jar s:C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Seq_under_peaks\EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks.fasta 1>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks_motifs_2.txt 2>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks_log_2.txt
echo chipmunk for Rif peaks 3..
java -jar C:\Sutor\Dists\chipmunk.jar s:C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Seq_under_peaks\EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks.fasta 1>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks_motifs_3.txt 2>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks_log_3.txt

echo chipmunk for shared peaks 2..
java -jar C:\Sutor\Dists\chipmunk.jar s:C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Seq_under_peaks\EcTopoI_noCTD_noRif_Rif_rep123_shared_nm_0.001_peaks.fasta 1>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_noRif_Rif_rep123_shared_nm_0.001_peaks_motifs_2.txt 2>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_noRif_Rif_rep123_shared_nm_0.001_peaks_log_2.txt
echo chipmunk for shared peaks 3..
java -jar C:\Sutor\Dists\chipmunk.jar s:C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Seq_under_peaks\EcTopoI_noCTD_noRif_Rif_rep123_shared_nm_0.001_peaks.fasta 1>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_noRif_Rif_rep123_shared_nm_0.001_peaks_motifs_3.txt 2>C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_noRif_Rif_rep123_shared_nm_0.001_peaks_log_3.txt