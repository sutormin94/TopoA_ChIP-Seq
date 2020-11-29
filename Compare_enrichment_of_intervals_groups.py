###############################################
##Dmitry Sutormin, 2020##
##TopoI ChIP-Seq analysis##

#The script compares a aignal between different intervals groups (e.g. peaks).
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
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol_Gyrase\Apr_data\\"

#Dictionary of intervals.
Intervals_path_dict={'EcTopoI_Rif_minus_CTD_minus_all' :         'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\\TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks.narrowPeak',
                     'RpoC_Borukhov_all' :                       'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol_Gyrase\RNApol_peaks\RNApol_peaks_threshold_3.BroadPeak',
                     'EcTopoI_Rif_minus_CTD_minus_unique' :      PWD + 'TopoA_Rif_minus_CTD_minus_nm_0.001_peaks_unique_from_RNAP_3_peaks.narrowPeak',
                     'RpoC_Borukhov_unique' :                    PWD + 'RNAP_3_peaks_unique_from_TopoA_Rif_minus_CTD_minus_nm_0.001_peaks.narrowPeak',
                     'EcTopoI_Rif_minus_CTD_minus_shared_RNAP' : PWD + 'TopoA_Rif_minus_CTD_minus_nm_0.001_shared_with_RNAP_3_peaks.narrowPeak',
                     }

#Dictionary of signal data (wig).
Wig_path_dict={'EcTopoI_Rif_minus_CTD_minus' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_123.wig',
               'RpoC_Borukhov' :               'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Borukhov_RpoC_Pol_Sofi_LB_FE.wig'}


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
#Retrive signal by intervals set.
#######

def retrive_signal(wig_data, intervals_set):
    signal_ar=[]
    for peak in intervals_set:
        signal_ar.append(np.mean(wig_data[peak[0]:peak[1]]))
    return signal_ar


#######
#Plot data, compare enrichment of interval sets.
#######

def compare_enrichment_plot(Intervals_dictionary, Wig_dictionary, pathout):
    
    #Get data EcTopoI.
    All_EcTopoI_sites_EcTopoI=retrive_signal(Wig_dictionary['EcTopoI_Rif_minus_CTD_minus'], Intervals_dictionary['EcTopoI_Rif_minus_CTD_minus_all'])
    Unique_EcTopoI_sites_EcTopoI=retrive_signal(Wig_dictionary['EcTopoI_Rif_minus_CTD_minus'], Intervals_dictionary['EcTopoI_Rif_minus_CTD_minus_unique'])
    Shared_EcTopoI_sites_EcTopoI=retrive_signal(Wig_dictionary['EcTopoI_Rif_minus_CTD_minus'], Intervals_dictionary['EcTopoI_Rif_minus_CTD_minus_shared_RNAP'])
    
    EcTopoI_together=[Unique_EcTopoI_sites_EcTopoI, Shared_EcTopoI_sites_EcTopoI]
    
    #Get data RNAP.
    All_RNAP_sites_RNAP=retrive_signal(Wig_dictionary['RpoC_Borukhov'], Intervals_dictionary['RpoC_Borukhov_all'])
    Unique_RNAP_sites_RNAP=retrive_signal(Wig_dictionary['RpoC_Borukhov'], Intervals_dictionary['RpoC_Borukhov_unique'])
    Shared_RNAP_sites_RNAP=retrive_signal(Wig_dictionary['RpoC_Borukhov'], Intervals_dictionary['EcTopoI_Rif_minus_CTD_minus_shared_RNAP'])
    
    RNAP_together=[Unique_RNAP_sites_RNAP, Shared_RNAP_sites_RNAP]
    
    #Plot data.
    fig=plt.figure(figsize=(15, 3), dpi=100)
    
    #EcTopoI peaks FE
    plot2=plt.subplot2grid((1,4),(0,0), rowspan=1, colspan=1) 
    
    bins_FE=np.linspace(min(Unique_EcTopoI_sites_EcTopoI+Shared_EcTopoI_sites_EcTopoI), max(Unique_EcTopoI_sites_EcTopoI+Shared_EcTopoI_sites_EcTopoI), 50)
    plot2.hist(Unique_EcTopoI_sites_EcTopoI, bins_FE, color='#6f89ff', edgecolor='black', alpha=0.5, label='Unique EcTopoI', zorder=3)
    plot2.hist(Shared_EcTopoI_sites_EcTopoI, bins_FE, color='#ffbc47', edgecolor='black', alpha=0.5, label='Shared EcTopoI', zorder=2)
    plot2.set_xlabel('EcTopoI fold enrichment', size=16)
    plot2.set_ylabel('Number of peaks', size=16)
    plot2.set_yscale('log', nonposy='clip')
    #plot2.set_title('Peaks FE\n'+set_name, size=14) 
    plot2.legend(loc='upper right', frameon=False, fontsize=8, handlelength=1)  
    
    plot3=plt.subplot2grid((1,4),(0,1), rowspan=1, colspan=1) 
    
    bins_FE=np.linspace(min(Unique_EcTopoI_sites_EcTopoI+Shared_EcTopoI_sites_EcTopoI), max(Unique_EcTopoI_sites_EcTopoI+Shared_EcTopoI_sites_EcTopoI), 20)
    colors=['#6f89ff', '#ffbc47']
    labels=['Unique EcTopoI', 'Shared EcTopoI']
    plot3.hist(EcTopoI_together, bins_FE, color=colors, edgecolor='black', alpha=0.8, stacked=True, label=labels)
    plot3.set_xlabel('EcTopoI fold enrichment', size=16)
    plot3.set_ylabel('Number of peaks', size=16)
    plot3.set_yscale('log')
    plot3.legend(loc='upper right', frameon=False, fontsize=8, handlelength=1)     
     
    
    #RNAP peaks FE
    plot0=plt.subplot2grid((1,4),(0,2), rowspan=1, colspan=1)  
      
    bins_FE=np.linspace(min(Unique_RNAP_sites_RNAP+Shared_RNAP_sites_RNAP), max(Unique_RNAP_sites_RNAP+Shared_RNAP_sites_RNAP), 50)
    plot0.hist(Unique_RNAP_sites_RNAP, bins_FE, color='#ff67d1', edgecolor='black', alpha=0.5, label='Unique RNAP', zorder=2)
    plot0.hist(Shared_RNAP_sites_RNAP, bins_FE, color='#53f7a5', edgecolor='black', alpha=0.5, label='Shared RNAP', zorder=3)
    plot0.set_xlabel('RNAP fold enrichment', size=16)
    plot0.set_ylabel('Fraction of peaks', size=16)
    plot0.set_yscale('log', nonposy='clip')
    plot0.legend(loc='upper right', frameon=False, fontsize=8, handlelength=1)  
    
    
    plot4=plt.subplot2grid((1,4),(0,3), rowspan=1, colspan=1) 
    
    bins_FE=np.linspace(min(Unique_RNAP_sites_RNAP+Shared_RNAP_sites_RNAP), max(Unique_RNAP_sites_RNAP+Shared_RNAP_sites_RNAP), 20)
    colors=['#ff67d1', '#53f7a5']
    labels=['Unique RNAP', 'Shared RNAP']
    plot4.hist(RNAP_together, bins_FE, color=colors, edgecolor='black', alpha=0.8, stacked=True, label=labels)
    plot4.set_xlabel('RNAP fold enrichment', size=16)
    plot4.set_ylabel('Number of peaks', size=16)
    plot4.set_yscale('log')
    plot4.legend(loc='upper right', frameon=False, fontsize=8, handlelength=1)       
 
    plt.tight_layout()
    plt.show()
    plt.savefig(pathout + 'FE_of_shared_unique_peaks_EcTopoI_RNAP.png', dpi=300, figsize=(15, 3))    
    
    #Welch t-test.
    Intervals_stat=stats.ttest_ind(Unique_EcTopoI_sites_EcTopoI, Shared_EcTopoI_sites_EcTopoI, equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(Unique_EcTopoI_sites_EcTopoI),2)} Mean2={round(np.mean(Shared_EcTopoI_sites_EcTopoI),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    Intervals_stat=stats.ttest_ind(Unique_RNAP_sites_RNAP, Shared_RNAP_sites_RNAP, equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(Unique_RNAP_sites_RNAP),2)} Mean2={round(np.mean(Shared_RNAP_sites_RNAP),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    return

Intervals_data_dictionary=read_peaks(Intervals_path_dict)
Wig_data_dictionary=read_wig(Wig_path_dict)
compare_enrichment_plot(Intervals_data_dictionary, Wig_data_dictionary, PWD)
