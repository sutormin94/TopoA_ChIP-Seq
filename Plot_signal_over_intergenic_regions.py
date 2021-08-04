###############################################
##Dmitry Sutormin, 2021##
##TopoI ChIP-Seq analysis##

#Takes WIG files with information about the cumulative distribution of some 
#protein over intergenic regions (IGRs). Plots this information.

###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
from scipy import stats
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
import matplotlib.patheffects as PathEffects
from matplotlib_venn import venn2, venn3, venn3_circles
from matplotlib import cm as cm
import collections
from collections import OrderedDict
import pandas as pd
from pandas import DataFrame


#Path to the working directory.
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\\'

#Name of the signal to plotted (protein or smth.).
Signal_name='RpoC_Borukhov_Fig5'
#Half-window width will be used to smooth signal.
Sm_window=10
#Dictionary of pathes to input WIG data.
Wig_data_in_dict_IGR={'All no dps gyrase ++' :            PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_Gyrase_Cfx_over_++IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps gyrase +-' :            PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_Gyrase_Cfx_over_+-IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps gyrase -+' :            PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_Gyrase_Cfx_over_-+IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps gyrase --' :            PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_Gyrase_Cfx_over_--IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps gyrase Rif ++' :        PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_Gyrase_Cfx_Rif_over_++IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps gyrase Rif +-' :        PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_Gyrase_Cfx_Rif_over_+-IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps gyrase Rif -+' :        PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_Gyrase_Cfx_Rif_over_-+IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps gyrase Rif --' :        PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_Gyrase_Cfx_Rif_over_--IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps RpoC ++' :              PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_RpoC_Borukhov_over_++IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps RpoC +-' :              PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_RpoC_Borukhov_over_+-IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps RpoC -+' :              PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_RpoC_Borukhov_over_-+IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps RpoC --' :              PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_RpoC_Borukhov_over_--IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps RpoB ++' :              PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_RpoB_Kahramanoglou_over_++IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps RpoB +-' :              PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_RpoB_Kahramanoglou_over_+-IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps RpoB -+' :              PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_RpoB_Kahramanoglou_over_-+IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps RpoB --' :              PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_RpoB_Kahramanoglou_over_--IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps EcTopoI-CTD/-Rif ++' :  PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_TopA_CTD_minus_Rif_minus_av_3_4_6_over_++IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps EcTopoI-CTD/-Rif +-' :  PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_TopA_CTD_minus_Rif_minus_av_3_4_6_over_+-IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps EcTopoI-CTD/-Rif -+' :  PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_TopA_CTD_minus_Rif_minus_av_3_4_6_over_-+IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',
                      'All no dps EcTopoI-CTD/-Rif --' :  PWD + 'Signal_of_TUs_wig\IGR_50_1000_all_no_dps\\Signal_TopA_CTD_minus_Rif_minus_av_3_4_6_over_--IGR_50_1000_all_no_dps_width_200bp_gb_100bp.wig',                          
                      }

#Set type to choose plotting parameters: genes, operons or transcripts.
Set_type="IGR"

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=f'{PWD}\Figures\Plot_combinations\{Signal_name}'
Dir_check_create(Out_path)


#######
#Parses WIG file with FE over TUs.
#######

def wig_FE_over_genes_parsing(name, wigfile):
    print('Now is processing: ' + str(name) + ' ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] in ['track']:
            ww_l=line[2].split('=')[1].rstrip('"').lstrip('"').split('_')
            win_width=int(ww_l[0])
            length=int(ww_l[1])
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    print(f'Window width: {win_width}, length of TU: {length}')
    return NE_values, win_width, length


#######
#Returns smoothed tracks.
#######

def Smoothing(ends, window):
    smoothed=[]
    #Calculating the value for the first position
    sm=0.0
    window_float=float(window)
    sm+=np.mean(ends[:2*window])
    smoothed.append(sm)
    #Calculating values for the part of the array remains
    for i in range(len(ends)-2*window):
        sm+=(ends[i+2*window]-ends[i])/(2*window_float)
        smoothed.append(sm)
    return smoothed


#######
#Plot the signal for different groups of genes together.
#######

def plot_intergenic_FE_TUs_groups(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    TU_sets_v={'All no dps gyrase ++' : 671, 'All no dps gyrase +-' : 320, 'All no dps gyrase -+' : 576, 'All no dps gyrase --' : 727,
               'All no dps gyrase Rif ++' : 671, 'All no dps gyrase Rif +-' : 320, 'All no dps gyrase Rif -+' : 576, 'All no dps gyrase Rif --' : 727,
               'All no dps RpoC ++' : 671, 'All no dps RpoC +-' : 320, 'All no dps RpoC -+' : 576, 'All no dps RpoC --' : 727,
               'All no dps RpoB ++' : 671, 'All no dps RpoB +-' : 320, 'All no dps RpoB -+' : 576, 'All no dps RpoB --' : 727,
               'All no dps EcTopoI-CTD/-Rif ++' : 671, 'All no dps EcTopoI-CTD/-Rif +-' : 320, 'All no dps EcTopoI-CTD/-Rif -+' : 576, 'All no dps EcTopoI-CTD/-Rif --' : 727
               }
    
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=15000
    length=5000
    for name, file in wig_in_dict.items():
        data=wig_FE_over_genes_parsing(name, file)
        win_width=data[1]
        length=data[2]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(positions[0], positions[-1])
       
    
    #Plot FE over genes.
    plt.figure(figsize=(4, 4), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    
    ##Transcription units below.
    if set_type=="IGR":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        #EcTopoI CTD-/Rif-
        plot1.plot(positions, np.array(dict_of_wigs['All no dps EcTopoI-CTD/-Rif ++'])-0.18, linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({TU_sets_v["All no dps EcTopoI-CTD/-Rif ++"]})', zorder=1) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['All no dps EcTopoI-CTD/-Rif +-'])-0.18, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({TU_sets_v["All no dps EcTopoI-CTD/-Rif +-"]})', zorder=2) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42
        plot1.plot(positions, np.array(dict_of_wigs['All no dps EcTopoI-CTD/-Rif -+'])-0.18, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({TU_sets_v["All no dps EcTopoI-CTD/-Rif -+"]})', zorder=3) #Def linewidth=2; #R123 -0.25; R123 -0.25; R123 -0.27
        plot1.plot(positions, np.array(dict_of_wigs['All no dps EcTopoI-CTD/-Rif --'])-0.18, linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({TU_sets_v["All no dps EcTopoI-CTD/-Rif --"]})', zorder=4) #Def linewidth=1; #R23 -0.25; R23 -0.22; R23 -0.24
        ##DNA-gyrase.
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps gyrase ++'])-0.1, linestyle='-', color='#757d8b', linewidth=1, alpha=1, label=f'-> -> ({TU_sets_v["All no dps gyrase ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps gyrase +-'])-0.1, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'-> <- ({TU_sets_v["All no dps gyrase +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps gyrase -+'])-0.1, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'<- -> ({TU_sets_v["All no dps gyrase -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps gyrase --'])-0.1, linestyle='-', color='#e4d1b4', linewidth=1, alpha=1, label=f'<- <- ({TU_sets_v["All no dps gyrase --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        ##DNA-gyrase Rif.
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps gyrase Rif ++'])+0.05, linestyle='-', color='#757d8b', linewidth=1, alpha=1, label=f'-> -> ({TU_sets_v["All no dps gyrase Rif ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps gyrase Rif +-'])+0.05, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'-> <- ({TU_sets_v["All no dps gyrase Rif +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps gyrase Rif -+'])+0.05, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'<- -> ({TU_sets_v["All no dps gyrase Rif -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps gyrase Rif --'])+0.05, linestyle='-', color='#e4d1b4', linewidth=1, alpha=1, label=f'<- <- ({TU_sets_v["All no dps gyrase Rif --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        ##RpoC RNAP.
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps RpoC ++'])-0.22, linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({TU_sets_v["All no dps RpoC ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps RpoC +-'])-0.22, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({TU_sets_v["All no dps RpoC +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps RpoC -+'])-0.22, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({TU_sets_v["All no dps RpoC -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps RpoC --'])-0.22, linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({TU_sets_v["All no dps RpoC --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        #RpoB RNAP.
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps RpoB ++']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({TU_sets_v["All no dps RpoB ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps RpoB +-']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({TU_sets_v["All no dps RpoB +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps RpoB -+']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({TU_sets_v["All no dps RpoB -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All no dps RpoB --']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({TU_sets_v["All no dps RpoB --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        
                    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS/TE'
    ticks_lables[ticks.index(length)]='TE/TS'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE/TS')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1) 
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    #plot1.axhline(1, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=11, frameon=False, markerscale=2, handlelength=0.5, handletextpad=0.2, loc='upper left')    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI fold enrichment', size=20)
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.tight_layout()
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_EcTopoI_no_CTD_no_Rif_346_{win_width}bp_with_IGR_{length}_bp.png', dpi=400, figsize=(4, 4), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(4, 4), dpi=100)
    plot1=plt.subplot(111)
    
    ##Transcription units below.
    if set_type=="IGR":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        #EcTopoI CTD-/Rif-
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps EcTopoI-CTD/-Rif ++'])-0.18, linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({TU_sets_v["All no dps EcTopoI-CTD/-Rif ++"]})', zorder=1) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps EcTopoI-CTD/-Rif +-'])-0.18, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({TU_sets_v["All no dps EcTopoI-CTD/-Rif +-"]})', zorder=2) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps EcTopoI-CTD/-Rif -+'])-0.18, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({TU_sets_v["All no dps EcTopoI-CTD/-Rif -+"]})', zorder=3) #Def linewidth=2; #R123 -0.25; R123 -0.25; R123 -0.27
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps EcTopoI-CTD/-Rif --'])-0.18, linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({TU_sets_v["All no dps EcTopoI-CTD/-Rif --"]})', zorder=4) #Def linewidth=1; #R23 -0.25; R23 -0.22; R23 -0.24
        ##DNA-gyrase. 
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps gyrase ++'])-0.1, linestyle='-', color='#757d8b', linewidth=1, alpha=1, label=f'-> -> ({TU_sets_v["All no dps gyrase ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps gyrase +-'])-0.1, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'-> <- ({TU_sets_v["All no dps gyrase +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps gyrase -+'])-0.1, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'<- -> ({TU_sets_v["All no dps gyrase -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps gyrase --'])-0.1, linestyle='-', color='#e4d1b4', linewidth=1, alpha=1, label=f'<- <- ({TU_sets_v["All no dps gyrase --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        ##DNA-gyrase Rif. 
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps gyrase Rif ++'])+0.05, linestyle='-', color='#757d8b', linewidth=1, alpha=1, label=f'-> -> ({TU_sets_v["All no dps gyrase Rif ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps gyrase Rif +-'])+0.05, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'-> <- ({TU_sets_v["All no dps gyrase Rif +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps gyrase Rif -+'])+0.05, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'<- -> ({TU_sets_v["All no dps gyrase Rif -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps gyrase Rif --'])+0.05, linestyle='-', color='#e4d1b4', linewidth=1, alpha=1, label=f'<- <- ({TU_sets_v["All no dps gyrase Rif --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        ##RpoC RNAP. 
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps RpoC ++'])-0.22, linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({TU_sets_v["All no dps RpoC ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps RpoC +-'])-0.22, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({TU_sets_v["All no dps RpoC +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps RpoC -+'])-0.22, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({TU_sets_v["All no dps RpoC -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps RpoC --'])-0.22, linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({TU_sets_v["All no dps RpoC --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        #RpoB RNAP.       
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps RpoB ++']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({TU_sets_v["All no dps RpoB ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps RpoB +-']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({TU_sets_v["All no dps RpoB +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps RpoB -+']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({TU_sets_v["All no dps RpoB -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All no dps RpoB --']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({TU_sets_v["All no dps RpoB --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS/TE'
    ticks_lables[ticks.index(length)]='TE/TS'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE/TS')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axvline(0, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.axvline(length, color='black', linestyle=':', alpha=0.7, linewidth=1.5)    
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=1)   
    #plot1.axhline(1, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.legend(fontsize=11, frameon=False, markerscale=2, handlelength=0.5, handletextpad=0.2, loc='upper left')    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI fold enrichment', size=20)
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)   
    plt.tight_layout()
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_EcTopoI_no_CTD_no_Rif_346_smoothed_{win_width}bp_with_IGR_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(4, 4), transparent=True)   
    plt.show()
    plt.close()    
    return

plot_intergenic_FE_TUs_groups(Wig_data_in_dict_IGR, Sm_window, Out_path, Signal_name, Set_type)