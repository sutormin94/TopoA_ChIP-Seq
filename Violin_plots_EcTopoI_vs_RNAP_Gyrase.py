###############################################
##Dmitry Sutormin, 2021##
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
path_to_intervals_data_dict={'EcTopoI'     : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\TopoA_noCTD_noRif_rep346_thr_3_nm_0.001_peaks.narrowPeak',
                             'RpoC'        : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol_Gyrase\RNApol_peaks\RpoC_Borukhov\RNApol_peaks_threshold_3.BroadPeak',
                             'RpoB'        : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol_Gyrase\RNApol_peaks\RpoB_Kahramanoglou\RNApol_peaks_threshold_450.BroadPeak',
                             'Gyrase'      : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Scripts\Gyrase_Topo-seq\Additional_genome_features\Cfx_10mkM_trusted_GCSs_h_s.broadPeak',
                             'EcTopoI Rif' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\TopoA_noCTD_Rif_rep123_thr_3_nm_0.001_peaks.narrowPeak',
                             'RpoC Rif'    : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol_Gyrase\RNApol_peaks\RpoC_Rif_Mooney\Mooney_RpoC_Rif_peaks_threshold_1.5.BroadPeak',
                             'MtbRNAP' :   'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\Mycobacterium_TopoI_Gyrase_RNAP\Peak_calling\\Uplekar_MtbRNAP_peaks_threshold_3.BroadPeak',
                             'MtbGyrase' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\Mycobacterium_TopoI_Gyrase_RNAP\Peak_calling\Ahmed_MtbGyrase_peaks_threshold_2.BroadPeak',                                                          
                             'MsmRNAP' :   'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\Mycobacterium_TopoI_Gyrase_RNAP\Peak_calling\Landick_MsmRNAP_peaks_threshold_4.BroadPeak',
                             'MsmTopoI' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\Mycobacterium_TopoI_Gyrase_RNAP\Peak_calling\Rani_MsmTopoI_FC_peaks_threshold_3.BroadPeak'                             
                             }

path_to_wig_files_dict={'EcTopoI'     : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_346.wig',
                        'RpoC'        : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Borukhov_RpoC_Pol_Sofi_LB_FE.wig',
                        'RpoB'        : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Kahramanoglou_RpoB_IP_ME.wig',
                        'Gyrase'      : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Sutormin_Gyrase_Cfx_10mkM_FE_av.wig',
                        'EcTopoI Rif' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Sutormin_TopA_ChIP_CTD_minus_Rif_plus_FE_av_123.wig',
                        'RpoC Rif'    : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Mooney_RpoC_Rif_FE.wig',
                        'MtbRNAP' :   'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\Mycobacterium_TopoI_Gyrase_RNAP\FE\\Uplekar_MtbRNAP_FE_av.wig',
                        'MtbGyrase' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\Mycobacterium_TopoI_Gyrase_RNAP\FE\Ahmed_Gyrase_MtbRa_nodup_FE_1.wig',                        
                        'MsmRNAP' :   'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\Mycobacterium_TopoI_Gyrase_RNAP\FE\Landick_MsmRNAP_FE_av.wig',
                        'MsmTopoI' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\Mycobacterium_TopoI_Gyrase_RNAP\FE\Rani_MsmTopoI_nodup_FC.wig'                        
                        }

#Output: path to the dir to store output
Outputpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Violin_plots"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)
    
    

#######
#Read peaks data.
#######

def read_peaks(path_dict, interval_name):
    inpath=path_dict[interval_name]
    print(f'Now is processing {interval_name} peaks')
    filein=open(inpath, 'r')
    peaks_ar=[]        
    if interval_name!='Gyrase':
        for line in filein:
            line=line.rstrip().split('\t')
            start=int(line[1])
            end=int(line[2])
            peaks_ar.append([start, end])
    else: #Gyrase data is a Topo-Seq data with the exact positions of cleavage sites. +/-50 bp interval is taken around each cleavage site.
        for line in filein:
            line=line.rstrip().split('\t')
            start=int(line[1])-50
            end=int(line[2])+50
            peaks_ar.append([start, end])  
    filein.close() 

    return peaks_ar


#######
#Read wig data.
#######

def read_wig(path_dict, signal_name):
    inpath=path_dict[signal_name]
    print(f'Now is processing {signal_name} wig')
    filein=open(inpath, 'r')
    wig_ar=[]    
    for line in filein:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            wig_ar.append(float(line[0])) 
    
    return wig_ar

#######
#Mark peaks.
#######

def mark_peaks_return_subsets(peaks_set_name, peaks_set, wig_name, wig, calc_option, deletions_ar):
    mask_len=len(wig)
    sets_dict={}

    print(f'{wig_name} masked with {peaks_set_name} peaks')
    
    if calc_option=='all_positions': #Calculate statistics over all positions.
        mask=[0]*mask_len
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
    
    elif calc_option=='regions': #Calculate statistics over region means.
        wig_in_peaks=[]
        peak_len=[]
        for peak in peaks_set:
            deletion_check=0
            for deletion in deletions_ar:
                if ((peak[0]>deletion[0]) & (peak[0]<deletion[1])) or ((peak[1]>deletion[0]) & (peak[1]<deletion[1])):
                    deletion_check=1
            if deletion_check==0:
                wig_in_peaks.append(np.mean(wig[peak[0]:(peak[1]+1)]))
                peak_len.append(peak[1]-peak[0]+1)
        
        #Exclude peaks and deleted regions.
        mask=[0]*mask_len
        for peak in peaks_set:
            mask[peak[0]:(peak[1]+1)]=[1]*((peak[1]+1)-peak[0])   
        for deletion in deletions_ar:
            mask[deletion[0]:(deletion[1]+1)]=[1]*((deletion[1]+1)-deletion[0])
        
        #Extract and bin the remaining non-masked fraction of a genome by mean width of peaks.
        wig_out_peaks_all_pos=[]
        for i in range(len(mask)):
            if mask[i]==0:
                wig_out_peaks_all_pos.append(wig[i])  
        peak_mean_len=int(np.mean(peak_len))
        wig_out_peaks=[]
        i=0
        while (i+1)*peak_mean_len<len(wig_out_peaks_all_pos):
            wig_out_peaks.append(np.mean(wig_out_peaks_all_pos[i*peak_mean_len:(i+1)*peak_mean_len]))
            i+=1
        wig_out_peaks.append(np.mean(wig_out_peaks_all_pos[i*peak_mean_len:]))
               
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

def draw_violins(dataset1, parameters, pwd):
    
    #Positions of violins.
    pos1=[1, 2]

    #Draw violin-plots.
    fig=plt.figure(figsize=(3.5,6), dpi=100)
    plt1=fig.add_subplot(1,1,1) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.85, showmeans=True, showmedians=True, points=2000)
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
    
    labels=[f'{parameters["interval_name"]}\npeaks', 'Other\nsites']
    set_axis_style(plt1, labels, pos1)
    
    yticknames1_params=parameters['yticknames_pos']
    yticknames1=np.arange(yticknames1_params[0], yticknames1_params[1], yticknames1_params[2])
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel(f'{parameters["signal_name"]} FE', size=15)
    ylim_min, ylim_max=parameters['ylim']
    plt1.set_ylim(ylim_min, ylim_max)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    #plt1.set_yscale('log')
    
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(dataset1[0]),2)}', xy=(0.5, parameters['annotate_pos'][0]), xycoords='data', size=12, rotation=90)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(dataset1[1]),2)}', xy=(1.5, parameters['annotate_pos'][1]), xycoords='data', size=12, rotation=90)        
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}\{parameters["signal_name"]}_FE_and_{parameters["interval_name"]}_peaks.png', dpi=400, figsize=(3.5,6))
    plt.savefig(f'{pwd}\{parameters["signal_name"]}_FE_and_{parameters["interval_name"]}_peaks.svg', dpi=400, figsize=(3.5,6))
    
    return


#######
#Wrapper function.
#######

def violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, parameters, Outpath): 
    
    #Dataset. Specify names of datasets to be compared.
    signal_name=parameters['signal_name']
    interval_name=parameters['interval_name']    
    
    #Get array of intervals (e.g., peaks) and data from wig.
    Peaks_array=read_peaks(path_to_intervals_data_dict, interval_name)
    WIG_array=read_wig(path_to_wig_files_dict, signal_name)
    
    #Mask regions, get data inside and outside of peaks.
    Sets_dictionary=mark_peaks_return_subsets(interval_name, Peaks_array, signal_name, WIG_array, parameters['calc_option'], parameters['del_ar'])
    dataset=[Sets_dictionary[f'{signal_name}_wig_in_{interval_name}_peaks'], Sets_dictionary[f'{signal_name}_wig_out_{interval_name}_peaks']]
    
    #Welch t-test.
    Intervals_stat=stats.ttest_ind(dataset[0], dataset[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset[0]),2)} Mean2={round(np.mean(dataset[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    draw_violins(dataset, parameters, Outpath)
    
    return

Parameters1={'signal_name' : 'RpoC', 'interval_name' : 'EcTopoI', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5], 
             'del_ar' : [[274500, 372148], [793800, 807500], [1199000, 1214000]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters1, Outputpath)

Parameters2={'signal_name' : 'EcTopoI', 'interval_name' : 'RpoC', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
             'del_ar' : [[274500, 372148], [793800, 807500], [1199000, 1214000]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters2, Outputpath)

Parameters3={'signal_name' : 'RpoC Rif', 'interval_name' : 'EcTopoI Rif', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
             'del_ar' : [[274500, 372148], [793800, 807500], [1199000, 1214000]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters3, Outputpath)

Parameters4={'signal_name' : 'EcTopoI Rif', 'interval_name' : 'RpoC Rif', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
             'del_ar' : [[274500, 372148], [793800, 807500], [1199000, 1214000]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters4, Outputpath)

Parameters5={'signal_name' : 'Gyrase', 'interval_name' : 'EcTopoI', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
             'del_ar' : [[274500, 372148], [793800, 807500], [1199000, 1214000]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters5, Outputpath)

Parameters6={'signal_name' : 'EcTopoI', 'interval_name' : 'Gyrase', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
             'del_ar' : [[274500, 372148], [793800, 807500], [1199000, 1214000]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters6, Outputpath)

Parameters7={'signal_name' : 'RpoB', 'interval_name' : 'EcTopoI', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
             'del_ar' : [[274500, 372148], [793800, 807500], [1199000, 1214000]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters7, Outputpath)

Parameters8={'signal_name' : 'EcTopoI', 'interval_name' : 'RpoB', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
             'del_ar' : [[274500, 372148], [793800, 807500], [1199000, 1214000]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters8, Outputpath)

Parameters9={'signal_name' : 'MsmRNAP', 'interval_name' : 'MsmTopoI', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
             'del_ar' : [[0, 0]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters9, Outputpath)

Parameters10={'signal_name' : 'MsmTopoI', 'interval_name' : 'MsmRNAP', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
              'del_ar' : [[0, 0]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters10, Outputpath)

Parameters11={'signal_name' : 'MtbGyrase', 'interval_name' : 'MtbRNAP', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
              'del_ar' : [[0, 0]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters11, Outputpath)

Parameters12={'signal_name' : 'MtbRNAP', 'interval_name' : 'MtbGyrase', 'calc_option' : 'regions', 'yticknames_pos' : [0, 8.5, 2], 'ylim' : [-0.5, 8], 'annotate_pos' : [5, 5],
              'del_ar' : [[0, 0]]}
violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, Parameters12, Outputpath)