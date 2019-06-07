###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes WIG files with information about the cumulative distribution of some 
#protein over transcription units (TUs). Plots this information.
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
PWD='F:\Signal_over_TUs'
#Name of the signal to plotted (protein or smth.).
Signal_name='Gyrase +Rif'
#Half-window width will be used to smooth signal.
Sm_window=100
#Dictionary of pathes to input WIG data.
Wig_data_in_dict_operons={'All operons' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\All_operons\Signal_Gyrase Cfx +Rif_over_All_operons_width_15000bp_gb_5000bp.wig',
                  'HEO 144' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\HEO_144\Signal_Gyrase Cfx +Rif_over_HEO_144_width_15000bp_gb_5000bp.wig',
                  'LEO 144' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\LEO_144\Signal_Gyrase Cfx +Rif_over_LEO_144_width_15000bp_gb_5000bp.wig',
                  'LAO 27' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\LAO_27\Signal_Gyrase Cfx +Rif_over_LAO_27_width_15000bp_gb_5000bp.wig',
                  'SAO 27' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\SAO_27\Signal_Gyrase Cfx +Rif_over_SAO_27_width_15000bp_gb_5000bp.wig',
                  }
Wig_data_in_dict_genes={'All genes' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\All_operons\Signal_Gyrase Cfx +Rif_over_All_operons_width_15000bp_gb_5000bp.wig',
                  'HEG 270' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\HEO_144\Signal_Gyrase Cfx +Rif_over_HEO_144_width_15000bp_gb_5000bp.wig',
                  'LEG 270' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\LEO_144\Signal_Gyrase Cfx +Rif_over_LEO_144_width_15000bp_gb_5000bp.wig',
                  }
#Set type to choose plotting parameters: genes_1 for Wig_data_in_dict_genes; operons_1 for Wig_data_in_dict_operons.
Set_type="operons"

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=PWD
Dir_check_create(Out_path)
Dir_check_create(PWD+'\Figures\Plots_together\\'+Signal_name)


#######
#Parses WIG file with FE over TUs.
#######

def wig_FE_over_genes_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
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
    print(win_width, length)
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
#Plot the signal for all groups of genes together.
#######


def plot_FE_all_expression_gg_Rif_no_Rif(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes within sets.
    TU_sets_v={'All genes' : 4119, 'HEG 270' : 269, 'LEG 270' : 270, 'HEG 370' : 369, 'LEG 370' : 370,
               'All operons' : 2327, 'HEO 144' : 143, 'LEO 144' : 144, 'HEO 186' : 185, 'LEO 186' : 186,
               'LAO 27' : 27, 'SAO 27' : 27, 'rRNA_operons' : 7}
    #FE averaged WIG parsing
    dict_of_wigs={}
    #win_width=15000
    #length=5000
    for name, file in wig_in_dict.items():
        print(name, file)
        data=wig_FE_over_genes_parsing(file)
        win_width=data[1]
        length=data[2]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(win_width, length)
       
    
    #Plot FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    ##Genes below
    if set_type=='genes':
        #LEG_270
        plot1.plot(positions, dict_of_wigs['LEG 270'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEG ({TU_sets_v["LEG 270"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['LEG_270_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEG Rif ({TU_sets_v["LEG_270_Rif"]})', zorder=5)   
        #LEG_370
        #plot1.plot(positions, dict_of_wigs['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_genes
        plot1.plot(positions, dict_of_wigs['All genes'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes ({TU_sets_v["All genes"]})', zorder=10)
        #plot1.plot(positions, dict_of_wigs['All_genes_Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All genes Rif ({TU_sets_v["All_genes_Rif"]})', zorder=9)
        #HEG_270
        plot1.plot(positions, dict_of_wigs['HEG 270'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEG no rRNA, tRNA ({TU_sets_v["HEG 270"]})', zorder=8)
        #plot1.plot(positions, dict_of_wigs['HEG_no_tRNA_rRNA_270_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEG no rRNA, tRNA Rif ({TU_sets_v["HEG_no_tRNA_rRNA_270_Rif"]})', zorder=7)
        #HEG_370
        #plot1.plot(positions, dict_of_wigs['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)
    ##Operons below.
    elif set_type=="operons":
        #LEO_144
        plot1.plot(positions, dict_of_wigs['LEO 144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO 144"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
        #LEO_186 or SAO_27
        plot1.plot(positions, dict_of_wigs['SAO 27'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'SAO ({TU_sets_v["SAO 27"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['LEO_186'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEO ({TU_sets_v["LEO_186"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['LEO_186_Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_operons
        plot1.plot(positions, dict_of_wigs['All operons'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons ({TU_sets_v["All operons"]})', zorder=10)
        #plot1.plot(positions, dict_of_wigs['All_operons_Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All operons Rif ({TU_sets_v["All_operons_Rif"]})', zorder=9)
        #HEO_144
        plot1.plot(positions, dict_of_wigs['HEO 144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO 144"]})', zorder=8)
        #plot1.plot(positions, dict_of_wigs['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
        #HEO_186 or rRNA operons or LAO_27
        plot1.plot(positions, dict_of_wigs['LAO 27'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'LAO ({TU_sets_v["LAO 27"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['HEO_186'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEO ({TU_sets_v["HEO_186"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['HEO_186_Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1) 
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{set_name} fold enrichment', size=20)
    plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\Figures\Plots_together\\{set_name}\\{set_name}_FE_over_{set_type}_{win_width}bp_nd_with_body_{length}.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    ##Genes below
    if set_type=='genes':
        #LEG_270
        plot1.plot(positions_sm, dict_of_wigs_sm['LEG 270'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEG ({TU_sets_v["LEG 270"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEG_270_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEG Rif ({TU_sets_v["LEG_270_Rif"]})', zorder=5)   
        #LEG_370
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_genes
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes ({TU_sets_v["All genes"]})', zorder=10)
        #plot1.plot(positions_sm, dict_of_wigs_sm['All_genes_Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All genes Rif ({TU_sets_v["All_genes_Rif"]})', zorder=9)
        #HEG_270
        plot1.plot(positions_sm, dict_of_wigs_sm['HEG 270'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEG no rRNA, tRNA ({TU_sets_v["HEG 270"]})', zorder=8)
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEG_no_tRNA_rRNA_270_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEG no rRNA, tRNA Rif ({TU_sets_v["HEG_no_tRNA_rRNA_270_Rif"]})', zorder=7)
        #HEG_370
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)  
    ##Operons below.
    elif set_type=="operons":
        #LEO_144
        plot1.plot(positions_sm, dict_of_wigs_sm['LEO 144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO 144"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
        #LEO_186 or SAO 27
        plot1.plot(positions_sm, dict_of_wigs_sm['SAO 27'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'SAO ({TU_sets_v["SAO 27"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO_186'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEO ({TU_sets_v["LEO_186"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO_186_Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_operons
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons ({TU_sets_v["All operons"]})', zorder=10)
        #plot1.plot(positions_sm, dict_of_wigs_sm['All_operons_Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All operons Rif ({TU_sets_v["All_operons_Rif"]})', zorder=9)
        #HEO_144
        plot1.plot(positions_sm, dict_of_wigs_sm['HEO 144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO 144"]})', zorder=8)
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
        #HEO_186 or rRNA operons or LAO_27
        plot1.plot(positions_sm, dict_of_wigs_sm['LAO 27'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'LAO ({TU_sets_v["LAO 27"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO_186'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEO ({TU_sets_v["HEO_186"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO_186_Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)      
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)   
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)    
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{set_name} fold enrichment', size=20)
    plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\Figures\Plots_together\\{set_name}\\{set_name}_FE_over_{set_type}_smoothed_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()    
    return


plot_FE_all_expression_gg_Rif_no_Rif(Wig_data_in_dict_operons, Sm_window, Out_path, Signal_name, Set_type)