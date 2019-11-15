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
PWD='F:\Signal_over_TUs\Transcript-based\Figures'
#Name of the signal to plotted (protein or smth.).
Signal_name='DNA-gyrase'
#Half-window width will be used to smooth signal.
Sm_window=100
#Dictionary of pathes to input WIG data.
Wig_data_in_dict_operons={'All operons' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\All_operons_no_dps\\Signal_TopoA -Rif_over_All_operons_no_dps_width_15000bp_gb_5000bp.wig',
                          'All operons Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\All_operons_no_dps\\Signal_TopoA +Rif_over_All_operons_no_dps_width_15000bp_gb_5000bp.wig',
                          'HEO 186' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\HEO_186_no_dps\\Signal_TopoA -Rif_over_HEO_186_no_dps_width_15000bp_gb_5000bp.wig',
                          'HEO 186 Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\HEO_186_no_dps\\Signal_TopoA +Rif_over_HEO_186_no_dps_width_15000bp_gb_5000bp.wig',
                          'LEO 186' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\LEO_186_no_dps\\Signal_TopoA -Rif_over_LEO_186_no_dps_width_15000bp_gb_5000bp.wig',
                          'LEO 186 Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\LEO_186_no_dps\\Signal_TopoA +Rif_over_LEO_186_no_dps_width_15000bp_gb_5000bp.wig',
                          'LAO 27' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\LAO_27\Signal_TopoA +Rif_over_LAO_27_width_15000bp_gb_5000bp.wig',
                          'SAO 27' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\SAO_27\Signal_TopoA +Rif_over_SAO_27_width_15000bp_gb_5000bp.wig',
                          }
Wig_data_in_dict_genes={'All genes' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\All_genes\\Signal_TopoA -Rif_over_All_genes_width_15000bp_gb_5000bp.wig',
                        'All genes Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\All_genes\\Signal_TopoA +Rif_over_All_genes_width_15000bp_gb_5000bp.wig',
                        'HEG 270' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\HEG_270\\Signal_TopoA -Rif_over_HEG_270_width_15000bp_gb_5000bp.wig',
                        'HEG 270 Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\HEG_270\\Signal_TopoA +Rif_over_HEG_270_width_15000bp_gb_5000bp.wig',
                        'LEG 270' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\LEG_270\\Signal_TopoA -Rif_over_LEG_270_width_15000bp_gb_5000bp.wig',
                        'LEG 270 Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_wig\LEG_270\\Signal_TopoA +Rif_over_LEG_270_width_15000bp_gb_5000bp.wig',
                        }
Wig_data_in_dict_transcripts={'All TUs' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_2173\Signal_TopoA -Rif_over_All_TU_2173_width_15000bp_gb_5000bp.wig', 
                              'All TUs Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_2173\Signal_TopoA +Rif_over_All_TU_2173_width_15000bp_gb_5000bp.wig', 
                              'All TUs no tRNA, rRNA' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_249_no_tRNA_rRNA\Signal_TopoA -Rif_over_All_TU_249_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig', 
                              'All TUs no tRNA, rRNA Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_249_no_tRNA_rRNA\Signal_TopoA +Rif_over_All_TU_249_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig', 
                              'HETU no tRNA, rRNA 317' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_317_no_tRNA_rRNA_ompX\Signal_TopoA -Rif_over_HE_TU_317_no_tRNA_rRNA_ompX_width_15000bp_gb_5000bp.wig', 
                              'HETU no tRNA, rRNA 317 Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_317_no_tRNA_rRNA_ompX\Signal_TopoA +Rif_over_HE_TU_317_no_tRNA_rRNA_ompX_width_15000bp_gb_5000bp.wig', 
                              'HETU 323' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_323\Signal_TopoA -Rif_over_HE_TU_323_width_15000bp_gb_5000bp.wig',
                              'HETU 323 Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_323\Signal_TopoA +Rif_over_HE_TU_323_width_15000bp_gb_5000bp.wig',
                              'HETU 321 no ompX' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_TopoA -Rif_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',
                              'HETU 321 no ompX Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_TopoA +Rif_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',                              
                              'LETU 245' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_245\Signal_TopoA -Rif_over_LE_TU_245_width_15000bp_gb_5000bp.wig',
                              'LETU 245 Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_245\Signal_TopoA +Rif_over_LE_TU_245_width_15000bp_gb_5000bp.wig',
                              'LETU 244 no ybiI' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_TopoA -Rif_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',
                              'LETU 244 no ybiI Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_TopoA +Rif_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',
                              'rRNA 7' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\\rRNA_7\Signal_TopoA -Rif_over_rRNA_7_width_15000bp_gb_5000bp.wig',
                              'rRNA 7 Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\\rRNA_7\Signal_TopoA +Rif_over_rRNA_7_width_15000bp_gb_5000bp.wig', 
                              'tRNA 49' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\\tRNA_49\Signal_TopoA -Rif_over_tRNA_49_width_15000bp_gb_5000bp.wig',
                              'tRNA 49 Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\\tRNA_49\Signal_TopoA +Rif_over_tRNA_49_width_15000bp_gb_5000bp.wig'                               
                              }
Wig_data_in_dict_transcripts_RNApol={'All TUs' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_2173\Signal_PolSofi_over_All_TU_2173_width_15000bp_gb_5000bp.wig',         
                                     'HETU 321 no ompX' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_PolSofi_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',               
                                     'LETU 244 no ybiI' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_PolSofi_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig'
                                     }

Wig_data_in_dict_transcripts_gyrase={'All TUs' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_2173\Signal_Gyrase Cfx_over_All_TU_2173_width_15000bp_gb_5000bp.wig', 
                                     'All TUs Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_2173\Signal_Gyrase Cfx +Rif_over_All_TU_2173_width_15000bp_gb_5000bp.wig', 
                                     'All TUs no tRNA, rRNA' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_249_no_tRNA_rRNA\Signal_Gyrase Cfx_over_All_TU_249_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig', 
                                     'All TUs no tRNA, rRNA Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_249_no_tRNA_rRNA\Signal_Gyrase Cfx +Rif_over_All_TU_249_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig',  
                                     'HETU 321 no ompX' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_Gyrase Cfx_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',
                                     'HETU 321 no ompX Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_Gyrase Cfx +Rif_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',                              
                                     'LETU 244 no ybiI' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_Gyrase Cfx_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',
                                     'LETU 244 no ybiI Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_Gyrase Cfx +Rif_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',                           
                                     }
#Set type to choose plotting parameters: genes, operons or transcripts.
Set_type="transcripts"

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=f'{PWD}\Plots_together\{Signal_name}'
Dir_check_create(Out_path)


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
               'All operons' : 2327, 'HEO 144' : 143, 'LEO 144' : 144, 'HEO 186' : 185, 'LEO 186' : 186, 'LAO 27' : 27, 'SAO 27' : 27, 'rRNA_operons' : 7, 
               'All TUs' : 1672, 'All TUs no tRNA, rRNA' : 1634, 'HETU no tRNA, rRNA 317' : 202, 'LETU 249' : 205, 
               'HETU 323' : 201, 'HETU 321 no ompX' : 200, 'LETU 245' : 202, 'LETU 244 no ybiI' : 201,
               'rRNA 7' : 7, 'tRNA 49' : 39}
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=15000
    length=5000
    for name, file in wig_in_dict.items():
        print(name, file)
        data=wig_FE_over_genes_parsing(file)
        win_width=data[1]
        length=data[2]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(win_width, length)
       
    
    #Plot FE over genes.
    plt.figure(figsize=(3, 6), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    ##Genes below
    if set_type=='genes':
        #LEG_270
        plot1.plot(positions, dict_of_wigs['LEG 270'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEG ({TU_sets_v["LEG 270"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['LEG 270 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEG Rif ({TU_sets_v["LEG 270"]})', zorder=5)   
        #LEG_370
        #plot1.plot(positions, dict_of_wigs['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_genes
        plot1.plot(positions, dict_of_wigs['All genes'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes ({TU_sets_v["All genes"]})', zorder=10)
        plot1.plot(positions, dict_of_wigs['All genes Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All genes Rif ({TU_sets_v["All genes"]})', zorder=9)
        #HEG_270
        plot1.plot(positions, dict_of_wigs['HEG 270'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEG ({TU_sets_v["HEG 270"]})', zorder=8)
        plot1.plot(positions, dict_of_wigs['HEG 270 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEG ({TU_sets_v["HEG 270"]})', zorder=7)
        #HEG_370
        #plot1.plot(positions, dict_of_wigs['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)
    ##Operons below.
    elif set_type=="operons":
        #LEO_144
        #plot1.plot(positions, dict_of_wigs['LEO 144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO 144"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
        #LEO_186 or SAO_27
        #plot1.plot(positions, dict_of_wigs['SAO 27'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'SAO ({TU_sets_v["SAO 27"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['LEO 186'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1.5, label=f'LEO ({TU_sets_v["LEO 186"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['LEO 186 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=0.8, label=f'LEO Rif ({TU_sets_v["LEO 186"]})', zorder=5)        
        #All_operons
        plot1.plot(positions, dict_of_wigs['All operons'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons ({TU_sets_v["All operons"]})', zorder=10)
        plot1.plot(positions, dict_of_wigs['All operons Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All operons Rif ({TU_sets_v["All operons"]})', zorder=9)
        #HEO_144
        #plot1.plot(positions, dict_of_wigs['HEO 144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO 144"]})', zorder=8)
        #plot1.plot(positions, dict_of_wigs['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
        #HEO_186 or rRNA operons or LAO_27
        #plot1.plot(positions, dict_of_wigs['LAO 27'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'LAO ({TU_sets_v["LAO 27"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
        plot1.plot(positions, dict_of_wigs['HEO 186'], linestyle='-', color='#b08642', linewidth=1.5, alpha=1.5, label=f'HEO ({TU_sets_v["HEO 186"]})', zorder=4)
        plot1.plot(positions, dict_of_wigs['HEO 186 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=0.8, label=f'HEO Rif ({TU_sets_v["HEO 186"]})', zorder=3)  
        ##Transcription units below.
    elif set_type=="transcripts":
        #LETU
        #plot1.plot(positions, dict_of_wigs['LETU 245'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU ({TU_sets_v["LETU 245"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['LETU 245 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LETU Rif ({TU_sets_v["LETU 245"]})', zorder=5)                 
        #LETU, no ybiI
        plot1.plot(positions, dict_of_wigs['LETU 244 no ybiI'], linestyle='-', color='#757d8b', linewidth=4, alpha=1, label=f'LETU ({TU_sets_v["LETU 244 no ybiI"]})', zorder=6) #Def linewidth=1.5
        plot1.plot(positions, dict_of_wigs['LETU 244 no ybiI Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LETU Rif ({TU_sets_v["LETU 244 no ybiI"]})', zorder=5) #Def linewidth=1     
        #All_TUs
        plot1.plot(positions, dict_of_wigs['All TUs'], linestyle='-', color='#333738', linewidth=4, alpha=0.8, label=f'All TUs ({TU_sets_v["All TUs"]})', zorder=10) #Def linewidth=2
        plot1.plot(positions, dict_of_wigs['All TUs Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs Rif ({TU_sets_v["All TUs"]})', zorder=9) #Def linewidth=1
        #All_TUs no tRNA, rRNA.
        #plot1.plot(positions, dict_of_wigs['All TUs no tRNA, rRNA'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All TUs no tRNA, rRNA ({TU_sets_v["All TUs no tRNA, rRNA"]})', zorder=10)
        #plot1.plot(positions, dict_of_wigs['All TUs no tRNA, rRNA Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs no tRNA, rRNA Rif ({TU_sets_v["All TUs no tRNA, rRNA"]})', zorder=9)
        #HETU
        #plot1.plot(positions, dict_of_wigs['HETU 323'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HETU ({TU_sets_v["HETU 323"]})', zorder=8)
        #plot1.plot(positions, dict_of_wigs['HETU 323 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU 323"]})', zorder=7)        
        #HETU, no ompX
        plot1.plot(positions, dict_of_wigs['HETU 321 no ompX'], linestyle='-', color='#b08642', linewidth=4, alpha=0.8, label=f'HETU ({TU_sets_v["HETU 321 no ompX"]})', zorder=8) #Def linewidth=1.5
        plot1.plot(positions, dict_of_wigs['HETU 321 no ompX Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU 321 no ompX"]})', zorder=7) #Def linewidth=1
        #rRNA
        #plot1.plot(positions, dict_of_wigs['rRNA 7'], linestyle='-', color='#bec1cb', linewidth=1.5, alpha=0.8, label=f'rRNA ({TU_sets_v["rRNA 7"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions, dict_of_wigs['rRNA 7 Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=1, label=f'rRNA Rif ({TU_sets_v["rRNA 7"]})', zorder=7) #Def linewidth=1
        #tRNA
        #plot1.plot(positions, dict_of_wigs['tRNA 49'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=0.8, label=f'tRNA ({TU_sets_v["tRNA 49"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions, dict_of_wigs['tRNA 49 Rif'], linestyle='--', color='#FFC000', linewidth=1, alpha=1, label=f'tRNA Rif ({TU_sets_v["tRNA 49"]})', zorder=7) #Def linewidth=1   
                     
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks+list(range(-100, 100, 10)))
    #plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1) 
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.set_yticks([1], minor='True')
    plot1.set_xlim(-100, 100)
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)
    #plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{set_name} fold enrichment', size=20)
    plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_{win_width}bp_nd_with_body_{length}_bp_no_dps_closer_look.png', dpi=400, figsize=(10, 6))  #Def size - 10, 6; Closer look - 3, 6
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
        plot1.plot(positions_sm, dict_of_wigs_sm['LEG 270 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEG Rif ({TU_sets_v["LEG 270"]})', zorder=5)   
        #LEG_370
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_genes
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes ({TU_sets_v["All genes"]})', zorder=10)
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All genes Rif ({TU_sets_v["All genes"]})', zorder=9)
        #HEG_270
        plot1.plot(positions_sm, dict_of_wigs_sm['HEG 270'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEG ({TU_sets_v["HEG 270"]})', zorder=8)
        plot1.plot(positions_sm, dict_of_wigs_sm['HEG 270 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEG Rif ({TU_sets_v["HEG 270"]})', zorder=7)
        #HEG_370
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)  
    ##Operons below.
    elif set_type=="operons":
        #LEO_144
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO 144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO 144"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
        #LEO_186 or SAO 27
        #plot1.plot(positions_sm, dict_of_wigs_sm['SAO 27'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'SAO ({TU_sets_v["SAO 27"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['LEO 186'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1.5, label=f'LEO ({TU_sets_v["LEO 186"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['LEO 186 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=0.8, label=f'LEO Rif ({TU_sets_v["LEO 186"]})', zorder=5)        
        #All_operons
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons ({TU_sets_v["All operons"]})', zorder=10)
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All operons Rif ({TU_sets_v["All operons"]})', zorder=9)
        #HEO_144
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO 144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO 144"]})', zorder=8)
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
        #HEO_186 or rRNA operons or LAO_27
        #plot1.plot(positions_sm, dict_of_wigs_sm['LAO 27'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'LAO ({TU_sets_v["LAO 27"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
        plot1.plot(positions_sm, dict_of_wigs_sm['HEO 186'], linestyle='-', color='#b08642', linewidth=1.5, alpha=1.5, label=f'HEO ({TU_sets_v["HEO 186"]})', zorder=4)
        plot1.plot(positions_sm, dict_of_wigs_sm['HEO 186 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=0.8, label=f'HEO Rif ({TU_sets_v["HEO 186"]})', zorder=3)  
    ##Transcription units below.
    elif set_type=="transcripts":
        #LETU
        #plot1.plot(positions_sm, dict_of_wigs_sm['LETU 245'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU ({TU_sets_v["LETU 245"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['LETU 245 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LETU Rif ({TU_sets_v["LETU 245"]})', zorder=5)                 
        #LETU, no ybiI
        plot1.plot(positions_sm, dict_of_wigs_sm['LETU 244 no ybiI'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU ({TU_sets_v["LETU 244 no ybiI"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['LETU 244 no ybiI Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LETU Rif ({TU_sets_v["LETU 244 no ybiI"]})', zorder=5)         
        #All_TUs
        plot1.plot(positions_sm, dict_of_wigs_sm['All TUs'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All TUs ({TU_sets_v["All TUs"]})', zorder=10)
        plot1.plot(positions_sm, dict_of_wigs_sm['All TUs Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs Rif ({TU_sets_v["All TUs"]})', zorder=9)
        #All_TUs no tRNA, rRNA.
        #plot1.plot(positions_sm, dict_of_wigs_sm['All TUs no tRNA, rRNA'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All TUs no tRNA, rRNA ({TU_sets_v["All TUs no tRNA, rRNA"]})', zorder=10)
        #plot1.plot(positions_sm, dict_of_wigs_sm['All TUs no tRNA, rRNA Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs no tRNA, rRNA Rif ({TU_sets_v["All TUs no tRNA, rRNA"]})', zorder=9)
        #HETU
        #plot1.plot(positions_sm, dict_of_wigs_sm['HETU 323'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HETU ({TU_sets_v["HETU 323"]})', zorder=8)
        #plot1.plot(positions_sm, dict_of_wigs_sm['HETU 323 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU 323"]})', zorder=7)        
        #HETU, no ompX
        plot1.plot(positions_sm, dict_of_wigs_sm['HETU 321 no ompX'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HETU ({TU_sets_v["HETU 321 no ompX"]})', zorder=8)
        plot1.plot(positions_sm, dict_of_wigs_sm['HETU 321 no ompX Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU 321 no ompX"]})', zorder=7)
        #rRNA
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7'], linestyle='-', color='#bec1cb', linewidth=1.5, alpha=0.8, label=f'rRNA ({TU_sets_v["rRNA 7"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=1, label=f'rRNA Rif ({TU_sets_v["rRNA 7"]})', zorder=7) #Def linewidth=1
        #tRNA
        #plot1.plot(positions_sm, dict_of_wigs_sm['tRNA 49'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=0.8, label=f'tRNA ({TU_sets_v["tRNA 49"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, dict_of_wigs_sm['tRNA 49 Rif'], linestyle='--', color='#FFC000', linewidth=1, alpha=1, label=f'tRNA Rif ({TU_sets_v["tRNA 49"]})', zorder=7) #Def linewidth=1   
        
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
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{set_name} fold enrichment', size=20)
    plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_smoothed_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp_no_dps.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()    
    return


plot_FE_all_expression_gg_Rif_no_Rif(Wig_data_in_dict_transcripts_gyrase, Sm_window, Out_path, Signal_name, Set_type)