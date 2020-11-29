###############################################
##Dmitry Sutormin, 2020##
##TopoI ChIP-Seq analysis##

#The script tests sets of genomic intervals (Peaks, TUs, BIMEs-1, BIMEs-2, IHF sites, Fis sites, H-NS sites, MatP sites, etc.)
#for the enrichment with some continously distributed character (RNApol fold enrichment, score, GC%, etc.) (t-test).
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import binom

#######
#Variables to be defined.
#######

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol_Gyrase\\"

#Input: Intervals (e.g. EcTopoI peaks) with additional info (FE, GC%, width, etc.).
path_to_intervals_data_dict={'EcTopoI' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks.narrowPeak',
                             'RpoB' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol_Gyrase\RNApol_peaks\RpoB_Kahramanoglou\RNApol_peaks_threshold_450.BroadPeak',
                             'Gyrase' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Scripts\Gyrase_Topo-seq\Additional_genome_features\Cfx_10mkM_trusted_GCSs_h_s.broadPeak'
                             }

path_to_wig_files_dict={'EcTopoI' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_123.wig',
                        'RpoB' :    'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Kahramanoglou_RpoB_IP_ME.wig',
                        'Gyrase' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\Topo-Seq_data_N3E\EP_Cfx_IP_Mu_10mkM_1_edt_N3E.wig'}

#Output: path to the dir to store output
Outputpath=PWD + "Apr_data"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)
    
    

#######
#Read peaks data.
#######

def read_peaks(path_dict):
    peaks_dict={}
    for name, inpath in path_dict.items():
        print(f'Now is processing {name} peaks')
        filein=open(inpath, 'r')
        peaks_ar=[]        
        if name!='Gyrase':
            for line in filein:
                line=line.rstrip().split('\t')
                start=int(line[1])
                end=int(line[2])
                peaks_ar.append([start, end])
        else:
            for line in filein:
                line=line.rstrip().split('\t')
                start=int(line[1])-50
                end=int(line[2])+50
                peaks_ar.append([start, end])  
        filein.close() 
        
        peaks_dict[name]=peaks_ar
    #print(peaks_dict)
    return peaks_dict


#######
#Read wig data.
#######

def read_wig(path_dict):
    wig_dict={}
    for name, inpath in path_dict.items():
        print(f'Now is processing {name} wig')
        filein=open(inpath, 'r')
        wig_ar=[]    
        for line in filein:
            line=line.rstrip().split(' ')
            if line[0] not in ['track', 'fixedStep']:
                wig_ar.append(float(line[0])) 
        wig_dict[name]=wig_ar
    
    return wig_dict

#######
#Mark peaks.
#######

def mark_peaks_return_subsets(peaks_dict, wig_dict):
    deletions_ar=[[274500, 372148], [793800, 807500], [1199000, 1214000]]
    sets_dict={}
    for wig_name, wig in wig_dict.items():
        for peaks_set_name, peaks_set in peaks_dict.items():
            print(f'{wig_name} masked with {peaks_set_name} peaks')
            mask=[0]*4647454
            for peak in peaks_set:
                mask[peak[0]:(peak[1]+1)]=[1]*((peak[1]+1)-peak[0])
            
            wig_in_peaks=[]
            wig_out_peaks=[]
            for i in range(len(mask)):
                deletion_check=0
                for deletion in deletions_ar:
                    if (i>deletion[0]) & (i<deletion[1]):
                        deletion_check=1
                if deletion_check==0:
                    if wig[i]>0: 
                        if mask[i]==0:
                            wig_out_peaks.append(wig[i])
                        else:
                            wig_in_peaks.append(wig[i])
                    
            sets_dict[f'{wig_name}_wig_in_{peaks_set_name}_peaks']=wig_in_peaks
            sets_dict[f'{wig_name}_wig_out_{peaks_set_name}_peaks']=wig_out_peaks
            
    return sets_dict

#######
#Draw violin plots.
#######

def set_axis_style(ax, labels, pos):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, size=6)
    ax.set_xlim(0.25, max(pos)+0.75)
    return

def draw_violins(sets_dict, pwd):
    
    #Positions of violins.
    pos1=[1, 2]
    
    #Datasets.
    dataset1=[sets_dict['EcTopoI_wig_in_RpoB_peaks'], sets_dict['EcTopoI_wig_out_RpoB_peaks']]

    #Draw violin-plots.
    fig=plt.figure(figsize=(4.5,6), dpi=100)
    plt1=fig.add_subplot(1,1,1) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff7762')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['RpoB\npeaks', 'Other\nsites']
    set_axis_style(plt1, labels, pos1)
    
    yticknames1=np.arange(0, 5, 1)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI FE', size=15)
    plt1.set_ylim(-0.2, 5.1)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    #plt1.set_yscale('log')
    
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(dataset1[0]),2)}', xy=(0.5, 2.2), xycoords='data', size=12, rotation=90)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(dataset1[1]),2)}', xy=(1.5, 2.2), xycoords='data', size=12, rotation=90)        
    
    #Welch t-test.
    Intervals_stat=stats.ttest_ind(dataset1[0], dataset1[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[0]),2)} Mean2={round(np.mean(dataset1[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')

    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Apr_data\EcTopoI_FE_and_RpoB_peaks.png', dpi=400, figsize=(4.5,6))   
    
    return

Peaks_dictionary=read_peaks(path_to_intervals_data_dict)
WIG_dictionary=read_wig(path_to_wig_files_dict)
Sets_dictionary=mark_peaks_return_subsets(Peaks_dictionary, WIG_dictionary)
draw_violins(Sets_dictionary, PWD)