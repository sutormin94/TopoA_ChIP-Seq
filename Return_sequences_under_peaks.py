###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes NarrowPeaks with ChIP-Seq peaks,
#Returns sequences under them,
#Computes GC%, writes multifasta file.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count

#Path to NarrowPeak file with peaks coordinates (MACS2 output).
Peak_data="Path\Peaks.NarrowPeak",
#Path to the reference genome (e.g. E_coli_w3110_G_Mu.fasta).
Genome="Path\Genome.fasta"
#Path to the working directory.
pwd="PWD\\"
#Outpath.
Path_out="PWD\\" + "Seq_under_peaks\\"
if not os.path.exists(Path_out):
    os.makedirs(Path_out)

#Name for output files (plot and mfa).
Name=''
    

#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
    return genomefa

#######
#Opens and reads BED or NarrowPeak files.
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
#Return sub-sequences from reference genome by coordinates provided.
#######

def return_seqs(peaks_ar, genome):
    seq_dict={}
    for peak in peaks_ar:
        seq_dict[peak[0]]=genome[peak[0]:peak[1]]
    return seq_dict

#######
#Writes sequences as mfa file.
#######

def write_seqs(seqs_dict, pathout):
    fileout=open(pathout, 'w')
    for k in seqs_dict:
        fileout.write(f'>{k}\n{seqs_dict[k]}\n')
    fileout.close()
    return

#######
#Plots distribution of peaks GC% in comparision to genome GC%.
#Plots distribution of peaks widths.
#######

def plot_distrib(seq_dict, genome, pathout):
    peaks_gc_values=[]
    peaks_width=[]
    for k in seq_dict:
        peaks_gc_values.append(GC_count(seq_dict[k]))
        peaks_width.append(len(seq_dict[k]))
    mean_len=int(np.mean(peaks_width))-1
    number_of_genome_bin=int(len(genome)/mean_len)-1
    genome_gc_values=[]
    for i in range(number_of_genome_bin):
        genome_gc_values.append(GC_count(genome[mean_len*i:mean_len*(i+1)]))

    GC_genome=GC_count(genome)
    print(f'GC% of the reference genome: {GC_genome}%') 
    print(f'GC% of the reference genome binned: {np.mean(genome_gc_values)}%')  
    print(f'GC% of the peaks regions: {np.mean(peaks_gc_values)}%')
    
    #Plot GC% distributions, peaks width distributions.
    fig=plt.figure(figsize=(15, 5), dpi=100)
    #GC% genome
    plot0=plt.subplot2grid((1,3),(0,0), rowspan=1, colspan=1)     
    plot0.hist(genome_gc_values, color='#ff878b', edgecolor='black', alpha=0.8)
    plot0.annotate(f'Mean genome GC%={round(np.mean(genome_gc_values),2)}%', xy=(0.10, 0.8), xycoords='axes fraction', size=15)
    plot0.annotate(f'Bin width={mean_len}bp', xy=(0.10, 0.7), xycoords='axes fraction', size=15)
    plot0.set_xlabel('Genome bins GC%', size=17)
    plot0.set_ylabel('Number of bins', size=17)
    plot0.set_title('Genome GC%', size=18)
    #GC% peaks
    plot1=plt.subplot2grid((1,3),(0,1), rowspan=1, colspan=1)     
    plot1.hist(peaks_gc_values, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot1.annotate(f'Mean peaks GC%={round(np.mean(peaks_gc_values),2)}%', xy=(0.10, 0.8), xycoords='axes fraction', size=15)
    plot1.annotate(f'Mean peaks width={mean_len}bp', xy=(0.10, 0.7), xycoords='axes fraction', size=15)
    plot1.set_xlabel('Peaks GC%', size=17)
    plot1.set_ylabel('Number of peaks', size=17)
    plot1.set_title('TopoA peaks GC%', size=18)
    #Peaks width
    plot2=plt.subplot2grid((1,3),(0,2), rowspan=1, colspan=1)     
    plot2.hist(peaks_width, color='#BAE85C', edgecolor='black', alpha=0.8)
    plot2.annotate(f'Mean peaks width={mean_len}bp', xy=(0.10, 0.8), xycoords='axes fraction', size=15)
    plot2.annotate(f'Total number of peaks={len(peaks_width)}', xy=(0.10, 0.7), xycoords='axes fraction', size=15)
    plot2.set_xlabel('Peaks width', size=17)
    plot2.set_ylabel('Number of peaks', size=17)
    plot2.set_title('TopoA peaks width', size=18)    
    plt.tight_layout()
    plt.show()
    plt.savefig(pathout, dpi=300, figsize=(10, 4))
    plt.close()    
    return

#######
#Wrapper: reads NarrowPeak input, return sequences by coordinates provided from reference genome,
#writes sequences under peaks to mfa file, plots distributions of GC% of sequences under the peaks and peaks widths.
#######

def wrap_peaks(fasta_inpath, peaks_inpath, name, outpath):
    genome=read_genome(fasta_inpath)
    peaks_data=deletions_info(peaks_inpath)
    seqs=return_seqs(peaks_data, genome)
    plot_distrib(seqs, genome, outpath+name+".png")
    write_seqs(seqs, outpath+name+".fasta")
    return


wrap_peaks(Genome, Peak_data, Name, Path_out)