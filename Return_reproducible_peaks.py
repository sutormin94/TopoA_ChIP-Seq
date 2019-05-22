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

#Dictionary of pathes to NarrowPeak file with peaks coordinates (MACS2 output).
Peaks_data={'Replic 1' : "Path\Peaks_1.narrowPeak",
            'Replic 2' : "Path\Peaks_2.narrowPeak",
            }
#Path to the reference genome (e.g. E_coli_w3110_G_Mu.fasta).
Genome="Path\Genome.fasta"
#Threshold for reproducible peaks calling (must not exceed number of replicas).
Threshold=int()
#Outpath.
Path_out="Path\Reproducible_peaks.broadPeak"
    
    
#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
        genome_id=record.ID
    return len(genomefa), genome_id

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

def Wrapper(reps_dict, thr, genome_path, outpath):
    #Reads genome fasta.
    genome_length, chrom_name=read_genome(genome_path)
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
    #Write reproducible peaks.
    write_bed(Rep_peaks_array, chrom_name, outpath)
    return
            
Wrapper(Peaks_data, Threshold, Genome, Path_out)       
    