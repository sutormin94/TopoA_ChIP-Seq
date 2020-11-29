###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes wig file with fold enrichment, finds regions with enrichment > threshold,
#returns BroadPeaks file with peaks.
###############################################

#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles

#Dictionary of pathes to wig file with fold enrichment.
Peaks_data={'RNApol' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Kahramanoglou_RpoB_IP_ME.wig'}
#Threshold for reproducible peaks calling (must not exceed number of replicas).
Threshold=int(750)
#Outpath.
Path_out='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol_Gyrase\RNApol_peaks\RpoB_Kahramanoglou\\'


#######
#Parsing WIG file.
#######

def score_data_parser(inpath, param_name):
    param_file=open(inpath, 'r')
    ar=[]
    for line in param_file:
        line=line.rstrip().split(' ')
        if line[0]=='fixedStep':
            chrom_name=line[1].split('=')[1]      
        if line[0] not in ['track', 'fixedStep']:
            ar.append(float(line[0]))
    param_file.close()
    print('Whole genome average ' + str(param_name) + ' : ' + str(sum(ar)/len(ar)))
    return ar, chrom_name


#######
#Find peaks in fold enrichment array.
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
    print(f'Number of peaks found: {len(rep_peaks_ar)}')
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
#Wrapper: takes wig file with fold enrichment, 
#finds regions with enrichment > threshold now are called peaks,
#returns BroadPeaks file with peaks.
#######    

def Wrapper(reps_dict, thr, outpath):
    rep_data={}
    for name, rep_path in reps_dict.items():
        #Reads FE data.
        FE_data, chrom_name=score_data_parser(rep_path, name)
        #Finds peaks.
        Rep_peaks_array=Find_rep_peaks(FE_data, thr)
        rep_data[name]=Rep_peaks_array
        #Write reproducible peaks.
        write_bed(Rep_peaks_array, chrom_name, f'{outpath}{name}_peaks_threshold_{thr}.BroadPeak')
    return
            
Wrapper(Peaks_data, Threshold, Path_out)       
