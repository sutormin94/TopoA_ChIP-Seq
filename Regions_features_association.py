###############################################
##Dmitry Sutormin, 2021##
##TopoI ChIP-Seq analysis##

#The script tests if some continously distributed characters (RNApol fold enrichment, score, GC%, etc.) 
#are correlated for a set of genomic intervals (TUs, Peaks, etc.).
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
from scipy.stats import spearmanr

#######
#Variables to be defined.
#######

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Signal_over_TUs\Representative_transcripts\Signal_of_TUs_tab\All_TUs_1672\\"

#Input: Intervals (e.g. EcTopoI peaks) with additional info (FE, GC%, width, etc.).
path_to_intervals_data=PWD + "Signal_over_TUs_All_representative_TUs.txt"

#Output: path to the dir to store output
Outputpath=PWD + "Figures\\"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)

#######
#Correlation of RpoC and RpoB fold enrichements for TUs.
#######

def factors_association_RpoC_RpoB(Intervals_info, out_plot_path):
    
    ##Filtering.
    factor_1=Intervals_info[(Intervals_info['RpoC_Borukhov_FE_GB']>0) & (Intervals_info['RpoB_Kahramanoglou_FE_GB']>0)]['RpoC_Borukhov_FE_GB']
    factor_2=Intervals_info[(Intervals_info['RpoC_Borukhov_FE_GB']>0) & (Intervals_info['RpoB_Kahramanoglou_FE_GB']>0)]['RpoB_Kahramanoglou_FE_GB']   
    
    ##Factor changes.
    factor_1=np.log2(factor_1)
    factor_2=np.log2(factor_2)
    
    ##Linnear fitting of log data.
    fit=np.polyfit(factor_1, factor_2, 1)
    print(fit)
    fit_fn=np.poly1d(fit) 
    xdata=np.linspace(min(factor_1), max(factor_1), 50)
    
    ##Correlation.
    pearson_cor=scipy.stats.pearsonr(factor_1, factor_2)
    print(f'Paerson correlation (RpoC log-fold enrichment, RpoC log-fold enrichment) for TUs: {pearson_cor}') 
    spearman_cor=scipy.stats.spearmanr(factor_1, factor_2)
    print(f'Spearman correlation (RpoC log-fold enrichment, RpoC log-fold enrichment) for TUs: {spearman_cor}') 
    
    ##Plot data.
    fig, plot=plt.subplots(1,1,figsize=(3,3), dpi=100)
    plot.scatter(factor_1, factor_2, s=2)
    plot.plot(xdata, fit_fn(xdata), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot.annotate(f'Spearman cor coef: {np.round(spearman_cor[0],3)}\np-value: {"{:.1e}".format(spearman_cor[1])}', xy=(-3, 10), xycoords='data', size=7)
    plot.set_xlabel('RpoC fold enrichment, log2', size=12)
    plot.set_ylabel('RpoB coverage depth, log2', size=12)
    plot.spines["top"].set_visible(False)
    plot.spines["right"].set_visible(False)
    plt.legend(fontsize=8, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc='upper left')
    plt.tight_layout()
    plt.show()
    
    plt.savefig(out_plot_path+'TEST_RpoC_log_enrichment_vs_RpoB_log_enrichment_for_TUs.png', dpi=300)    
    return


#######
#Correlation of RpoC fold enrichment and transcription level of TUs determined by RNA-Seq.
#######

def factors_association_RpoC_RNA_Seq(Intervals_info, out_plot_path):
    
    ##Filtering.
    factor_1=Intervals_info[(Intervals_info['RpoC_Borukhov_FE_GB']>0) & (Intervals_info['RNA_Seq_FE_GB']>0)]['RpoC_Borukhov_FE_GB']
    factor_2=Intervals_info[(Intervals_info['RpoC_Borukhov_FE_GB']>0) & (Intervals_info['RNA_Seq_FE_GB']>0)]['RNA_Seq_FE_GB']   
    
    ##Factor changes.
    factor_1=np.log2(factor_1)
    factor_2=np.log2(factor_2)
    
    ##Linnear fitting of log data.
    fit=np.polyfit(factor_1, factor_2, 1)
    print(fit)
    fit_fn=np.poly1d(fit) 
    xdata=np.linspace(min(factor_1), max(factor_1), 50)
    
    ##Correlation.
    pearson_cor=scipy.stats.pearsonr(factor_1, factor_2)
    print(f'Paerson correlation (RpoC log-fold enrichment, log transcription level) for TUs: {pearson_cor}') 
    spearman_cor=scipy.stats.spearmanr(factor_1, factor_2)
    print(f'Spearman correlation (RpoC log-fold enrichment, log transcription level) for TUs: {spearman_cor}') 
    
    ##Plot data.
    fig, plot=plt.subplots(1,1,figsize=(3,3), dpi=100)
    plot.scatter(factor_1, factor_2, s=2)
    plot.plot(xdata, fit_fn(xdata), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot.annotate(f'Spearman cor coef: {np.round(spearman_cor[0],3)}\np-value: {"{:.1e}".format(spearman_cor[1])}', xy=(-3, 10), xycoords='data', size=7)
    plot.set_xlabel('RpoC fold enrichment, log2', size=12)
    plot.set_ylabel('Transcription level, log2 FPKM', size=12)
    plot.spines["top"].set_visible(False)
    plot.spines["right"].set_visible(False)
    plt.legend(fontsize=8, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc='upper left')
    plt.tight_layout()
    plt.show()
    
    plt.savefig(out_plot_path+'TEST_RpoC_log_enrichment_vs_Transcription_log_level_for_TUs.png', dpi=300)    
    return


#######
#Correlation of EcTopoI and EcTopoI fold enrichments for TUs.
#######

def factors_association_EcTopoI_RpoC(Intervals_info, out_plot_path):
    
    ##Filtering.
    factor_1=Intervals_info[(Intervals_info['RpoC_Borukhov_FE_GB']>0) & (Intervals_info['TopA_CTD_minus_Rif_minus_av_1_2_3_FE_GB']>0)]['TopA_CTD_minus_Rif_minus_av_1_2_3_FE_GB']
    factor_2=Intervals_info[(Intervals_info['RpoC_Borukhov_FE_GB']>0) & (Intervals_info['TopA_CTD_minus_Rif_minus_av_1_2_3_FE_GB']>0)]['RpoC_Borukhov_FE_GB']
       
    
    ##Factor changes.
    factor_1=np.log2(factor_1)
    factor_2=np.log2(factor_2)
    
    ##Linnear fitting of log data.
    fit=np.polyfit(factor_1, factor_2, 1)
    print(fit)
    fit_fn=np.poly1d(fit) 
    xdata=np.linspace(min(factor_1), max(factor_1), 50)
    
    ##Correlation.
    pearson_cor=scipy.stats.pearsonr(factor_1, factor_2)
    print(f'Paerson correlation (EcTopoI log-fold enrichment, RpoC log-fold enrichment) for TUs: {pearson_cor}') 
    spearman_cor=scipy.stats.spearmanr(factor_1, factor_2)
    print(f'Spearman correlation (EcTopoI log-fold enrichment, RpoC log-fold enrichment) for TUs: {spearman_cor}') 
    
    ##Plot data.
    fig, plot=plt.subplots(1,1,figsize=(3,3), dpi=100)
    plot.scatter(factor_1, factor_2, s=2)
    plot.plot(xdata, fit_fn(xdata), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot.annotate(f'Spearman cor coef: {np.round(spearman_cor[0],3)}\np-value: {"{:.1e}".format(spearman_cor[1])}', xy=(-0.5, 4), xycoords='data', size=7)
    plot.set_xlabel('EcTopoI fold enrichment, log2', size=12)
    plot.set_ylabel('RpoC fold enrichment, log2', size=12)
    plot.spines["top"].set_visible(False)
    plot.spines["right"].set_visible(False)
    plt.legend(fontsize=8, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc='upper left')
    plt.tight_layout()
    plt.show()
    
    plt.savefig(out_plot_path+'TEST_EcTopoI_log_enrichment_vs_RpoC_log_enrichment_for_TUs.png', dpi=300)    
    return


#######
#Functions wrapper.
#######

def func_wrapper(intervals_sets_path, outpath):
    
    ##Read data
    Intervals_info=pd.read_csv(intervals_sets_path, sep='\t', header=0, index_col=False, dtype={'Start' : np.int64, 'End' : np.int64})
    
    ##EcTopoI FE vs RpoC FE.
    factors_association_EcTopoI_RpoC(Intervals_info, outpath)
    
    ##RpoC FE vs RpoB FE.
    factors_association_RpoC_RpoB(Intervals_info, outpath)
    
    ##RpoC FE vs RNA-Seq data.
    factors_association_RpoC_RNA_Seq(Intervals_info, outpath)
    
    return

func_wrapper(path_to_intervals_data, Outputpath)

print('Script ended its work succesfully!') 