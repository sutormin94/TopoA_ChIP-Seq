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
from matplotlib_venn import venn2, venn2_circles

#Dictionary of pathes to NarrowPeak file with peaks coordinates (MACS2 output).
Peaks_data={'Rif-' : "F:\TopoI_ChIP-Seq\Ec_TopoI_data\Peak_calling\Reproducible_peaks\TopoA_rep12_FC_nm_0.001_peaks.narrowPeak",
            'Rif+' : "F:\TopoI_ChIP-Seq\Ec_TopoI_data\Peak_calling\Reproducible_peaks\TopoA_rep12_Rif_FC_nm_0.001_peaks.narrowPeak",
            }
#Path to the reference genome (e.g. E_coli_w3110_G_Mu.fasta).
Genome="C:\Sutor\science\TopoI_Topo-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"
#Threshold for reproducible peaks calling (must not exceed number of replicas).
Threshold=int(2)
#Outpath.
Path_out="F:\TopoI_ChIP-Seq\Ec_TopoI_data\Peak_calling\Reproducible_peaks\Test_TopoA_min_Rif_plus_Rif_common_0.001_peaks.narrowPeak"
#Pics outpath.
Pics_path_out="F:\TopoI_ChIP-Seq\Ec_TopoI_data\Peak_calling\Reproducible_peaks\\"
    
    
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
    #Create template genome-long array.
    genome_ar=[0]*genome_length
    #Indicates peaks.
    for name, peaks_ar in rep_data.items():
        genome_ar=Indicate_where_peaks(genome_ar, peaks_ar)
    #Identify reproducible peaks.
    Rep_peaks_array=Find_rep_peaks(genome_ar, thr)
    #Create Venn diagram represents replicas overlapping.
    plt.figure(figsize=(4,4))    
    keys_list=list(rep_data.keys())
    venn2(subsets=(len(rep_data[keys_list[0]])-len(Rep_peaks_array), len(rep_data[keys_list[1]])-len(Rep_peaks_array), len(Rep_peaks_array)), set_labels=(keys_list[0], keys_list[1]))
    venn2_circles(subsets=(len(rep_data[keys_list[0]])-len(Rep_peaks_array), len(rep_data[keys_list[1]])-len(Rep_peaks_array), len(Rep_peaks_array)), linestyle='solid')    
    plt.show()
    plt.savefig(pics_outpath+'Test_reproducible_peaks_between_Rif_plus_and_minus.png', dpi=400, figsize=(4, 4))
    #Write reproducible peaks.
    write_bed(Rep_peaks_array, chrom_name, outpath)
    return
            
Wrapper(Peaks_data, Threshold, Genome, Path_out, Pics_path_out)       
    