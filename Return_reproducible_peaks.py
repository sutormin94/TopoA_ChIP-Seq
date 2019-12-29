###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes narrowPeaks files with ChIP-Seq peaks identified in different biological replicas,
#return narrowPeak only with reproducible peaks.
###############################################

#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

#Dictionary of pathes to NarrowPeak file with peaks coordinates (MACS2 output).
#Peaks_data={'noCTD_Rif_1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Old_Data\TopoI_ChIP-Seq\Ec_TopoI_data\Peak_calling\TopoA_rep1_Rif_FC_nm_0.001\TopoA_rep1_Rif_FC_nm_0.001_peaks.narrowPeak",
#            'noCTD_Rif_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Old_Data\TopoI_ChIP-Seq\Ec_TopoI_data\Peak_calling\TopoA_rep2_Rif_FC_nm_0.001\TopoA_rep2_Rif_FC_nm_0.001_peaks.narrowPeak",
#            'noCTD_Rif_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\TopoA_rep3_noCTD_Rif_FC_nm_0.001\TopoA_rep3_noCTD_Rif_FC_nm_0.001_peaks.narrowPeak",
#            }

#Peaks_data={'noCTD_noRif_1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Old_Data\TopoI_ChIP-Seq\Ec_TopoI_data\Peak_calling\TopoA_rep1_FC_nm_0.001\TopoA_rep1_FC_nm_0.001_peaks.narrowPeak",
#            'noCTD_noRif_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Old_Data\TopoI_ChIP-Seq\Ec_TopoI_data\Peak_calling\TopoA_rep2_FC_nm_0.001\TopoA_rep2_FC_nm_0.001_peaks.narrowPeak",
#            'noCTD_noRif_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\TopoA_rep3_noCTD_noRif_FC_nm_0.001\TopoA_rep3_noCTD_noRif_FC_nm_0.001_peaks.narrowPeak",
#            }

Peaks_data={'noCTD_noRif' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks.narrowPeak",
            'noCTD_Rif' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\TopoA_noCTD_Rif_rep123_thr_3_nm_0.001_peaks.narrowPeak",
            }

#Path to the reference genome (e.g. E_coli_w3110_G_Mu.fasta).
Genome="C:\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"
#Threshold for reproducible peaks calling (must not exceed number of replicas).
Threshold=int(2)
#Outpath.
Path_out="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\TopoA_noCTD_Rif_noRif_shared_rep123_thr_2_nm_0.001_peaks.narrowPeak"
#Pics outpath.
Pics_path_out="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\\"
    
    
#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
        genome_id=record.name
    return len(genomefa), genomefa, genome_id

#######
#Opens and reads BED or narrowPeak files.
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar

#######
#Indicate where peaks occures by addition of 1 to these positions to genome-length array.
#######

def Indicate_where_peaks(genome_ar, peaks_ar):
    for peak in peaks_ar:
        for i in range (peak[1]+1-peak[0]):
            genome_ar[peak[0]+i]+=1
    return genome_ar

#######
#Find reproducible regions in genome-length array.
#######  

def Find_rep_peaks(genome_ar, thr):
    peak=0
    rep_peaks_ar=[]
    for i in range(len(genome_ar)):
        if genome_ar[i]<thr and peak==0: #We are not in peak.
            continue
        elif genome_ar[i]>=thr and peak==0: #We are at left peak border.
            peak=1
            current_peak=[i]
            continue
        elif genome_ar[i]>=thr and peak==1: #We are within a peak.
            continue
        elif genome_ar[i]<thr and peak==1: #We are at the right peak border.
            peak=0
            current_peak.append(i)
            rep_peaks_ar.append(current_peak)
            continue
    return rep_peaks_ar
        
#######
#Write reproducible peaks in broadPeak format.
#######   

def write_bed(rep_peaks_ar, chrom_name, outpath):
    fileout=open(outpath, 'w')
    for i in range(len(rep_peaks_ar)):
        fileout.write(chrom_name+'\t'+str(rep_peaks_ar[i][0])+'\t'+str(rep_peaks_ar[i][1])+'\tPeak_'+str(i)+'\t10\t.\t-1\t-1\t-1\n')
    fileout.close()
    return

#######
#Identifies reproducible peaks with threshold (number of samples in which a peak should be present) given.
#######  

def overlap_call(rep_data, thr, genome_length):
    #Create template genome-long array.
    genome_ar=[0]*genome_length
    #Indicates peaks.
    for name, peaks_ar in rep_data.items():
        genome_ar=Indicate_where_peaks(genome_ar, peaks_ar)
    #Identify reproducible peaks.
    Rep_peaks_array=Find_rep_peaks(genome_ar, thr)
    return Rep_peaks_array
    
#######
#Wrapper: takes peaks from different biological replicas,
#Identifies reproducible regions, writes broadPeak file with reproducible peaks.
#######    

def Wrapper(reps_dict, thr, genome_path, outpath, pics_outpath):
    #Reads genome fasta.
    genome_length, genome_seq, chrom_name=read_genome(genome_path)
    #Reads replicas data.
    rep_data={}
    for name, rep_path in reps_dict.items():
        rep_data[name]=deletions_info(rep_path)
    #Create template genome-long array, Indicates peaks, Identify reproducible peaks.
    Rep_peaks_array=overlap_call(rep_data, thr, genome_length)
    #Create Venn diagram represents replicas overlapping.
    plt.figure(figsize=(4,4))   
    keys_list=list(rep_data.keys())
    if len(keys_list)==2:
        venn2(subsets=(len(rep_data[keys_list[0]])-len(Rep_peaks_array), len(rep_data[keys_list[1]])-len(Rep_peaks_array), len(Rep_peaks_array)), set_labels=(keys_list[0], keys_list[1]))
        venn2_circles(subsets=(len(rep_data[keys_list[0]])-len(Rep_peaks_array), len(rep_data[keys_list[1]])-len(Rep_peaks_array), len(Rep_peaks_array)), linestyle='solid')    
        plt.show()
    elif len(keys_list)==3:
        #123-overlap.
        Rep_peaks_array_123=overlap_call(rep_data, 3, genome_length)
        #12-overlap.
        Rep_peaks_array_12=overlap_call({keys_list[0]: rep_data[keys_list[0]], keys_list[1]: rep_data[keys_list[1]]}, 2, genome_length)
        #13-overlap.
        Rep_peaks_array_13=overlap_call({keys_list[0]: rep_data[keys_list[0]], keys_list[2]: rep_data[keys_list[2]]}, 2, genome_length)
        #23-overlap.
        Rep_peaks_array_23=overlap_call({keys_list[1]: rep_data[keys_list[1]], keys_list[2]: rep_data[keys_list[2]]}, 2, genome_length)        
        #Plot.
        Only_1=len(rep_data[keys_list[0]])-len(Rep_peaks_array_12)-len(Rep_peaks_array_13)+len(Rep_peaks_array_123)
        Only_2=len(rep_data[keys_list[1]])-len(Rep_peaks_array_12)-len(Rep_peaks_array_23)+len(Rep_peaks_array_123)
        Only_3=len(rep_data[keys_list[2]])-len(Rep_peaks_array_13)-len(Rep_peaks_array_23)+len(Rep_peaks_array_123)
        Only_12=len(Rep_peaks_array_12)-len(Rep_peaks_array_123)
        Only_13=len(Rep_peaks_array_13)-len(Rep_peaks_array_123)
        Only_23=len(Rep_peaks_array_23)-len(Rep_peaks_array_123)
        Only_123=len(Rep_peaks_array_123)
        print(Only_1, Only_2, Only_12, Only_3, Only_13, Only_23, Only_123)
        venn3(subsets=(Only_1, Only_2, Only_12, Only_3, Only_13, Only_23, Only_123), set_labels=(keys_list[0], keys_list[1], keys_list[2]))
        venn3_circles(subsets=(Only_1, Only_2, Only_12, Only_3, Only_13, Only_23, Only_123), linestyle='solid')    
        plt.show()        
    plt.savefig(pics_outpath+'TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks.png', dpi=400, figsize=(4, 4))
    #Write reproducible peaks.
    write_bed(Rep_peaks_array, chrom_name, outpath)
    return
            
Wrapper(Peaks_data, Threshold, Genome, Path_out, Pics_path_out)       
    