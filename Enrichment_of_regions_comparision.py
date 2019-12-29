###############################################
##Dmitry Sutormin, 2019##
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
import scipy
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import binom

#######
#Variables to be defined.
#######

#Input: Intervals (e.g. EcTopoI peaks) (BroadPeak or NarrowPeaks).
path_to_intervals_sets={'RNApol': "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol\RNApol_peaks\RNApol_peaks_threshold_3.BroadPeak",
                        'EcTopoI': "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\TopoA_noCTD_noRif_rep123_thr_2_nm_0.001_peaks.narrowPeak"
                        }

#Input: Continously distributed value (e.g. RNApol distribution or EcTopoI distribution) (WIG).
cont_characters_dict={'RNApol' : "C:\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\Pol_Sofi_LB_w3110_for_Mu.wig",
                      'EcTopoI' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_minus_average_FE_1_2_3.wig",
                      'GC%': "C:\Sutor\Science\Gyrase_Topo-Seq\Scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig",
                      'EcTopoI ms' : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_scanning\ChIP-Munk\Rep12_thr_0.001\EcTopoI_motif_w3110_scanned_both.wig"                      
                      }

#Output: path to the dir to store output
Outputpath="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_vs_RNApol\\Nov_data\\"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)
    

#######
#BroadPeak-formalized intervals parsing (NAPs sites or BIMEs, etc) and filtering intervals that are not deleted.
#######

def broadpeak_pars(intervals_sets_path):
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    intervals_sets_dict={}
    for k, v in intervals_sets_path.items():
        filein=open(v, 'r')
        ar=[]        
        for line in filein:
            line=line.rstrip().split('\t')
            int_start=int(line[1])
            int_end=int(line[2])
            del_check=0
            for j in range(len(deletions)):
                if deletions[j][1]>int_start>deletions[j][0] or deletions[j][1]>int_end>deletions[j][0]:
                    del_check=1
            if del_check==0:
                ar.append([int_start, int_end]) #Interval start, interval end        
        intervals_sets_dict[k]=ar
        print("Number of " + str(k) + " regions: " + str(len(ar)))
        filein.close()
    return intervals_sets_dict


#######
#Parsing WIG file.
#######

def score_data_parser(inpath_dict):
    cont_char_dict={}
    for param_name, inpath in inpath_dict.items():
        param_file=open(inpath, 'r')
        ar=[]
        for line in param_file:
            line=line.rstrip().split(' ')
            if line[0] not in ['track', 'fixedStep']:
                ar.append(float(line[0]))
        param_file.close()
        print('Whole genome average ' + str(param_name) + ' : ' + str(sum(ar)/len(ar)))
        cont_char_dict[param_name]=ar 
    return cont_char_dict

#######
#Mask deletions.
#######

def mask_array(ar, regions_to_mask):
    #Mask deletions or smth in FE array.
    maska=[0]*len(ar)
    for deletion in regions_to_mask:
        del_len=deletion[1]-deletion[0]
        for i in range(del_len):
            maska[deletion[0]+i]=1
    ar_masked=np.ma.masked_array(ar, mask=maska)  
    ar_non_del_only=list(ar_masked[~ar_masked.mask])
    print(len(ar), len(ar_non_del_only))
    return ar_non_del_only

#######
#Return lengths of intervals.
#######

def width_of_intervals(intervals_sets_dict, name_of_intervals):
    intervals=intervals_sets_dict[name_of_intervals]
    intervals_len=[]
    for interval in intervals:
        intervals_len.append(interval[1]-interval[0])
    return intervals_len

#######
#Return continous characteristics of intervals, mask deletions and intervals regions.
#######

def cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals, name_of_cont_char):
    intervals=intervals_sets_dict[name_of_intervals]
    cont_char=cont_char_dict[name_of_cont_char]
    #Return enrichment of intervals.
    intervals_param_ar_tg=[]
    intervals_param_ar_sp=[]
    for interval in intervals:
        intervals_param_ar_tg+=cont_char[interval[0]:interval[1]]   
        intervals_param_ar_sp.append(np.mean(cont_char[interval[0]:interval[1]]))
    
    #Prepare continuos character without deleted regions and regions defined by intervals. 
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    regions_to_mask=deletions+intervals
    cont_char_masked=mask_array(cont_char, regions_to_mask)
    return intervals_param_ar_tg, cont_char_masked, intervals_param_ar_sp

#######
#Compare RNApol occupation of BIMEs and overall RNApol occupation.
#######

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels)+1))
    ax.set_xticklabels(labels, size=35)
    ax.set_xlim(0.25, len(labels)+0.75)
    return

def Intervals_occupation(intervals_param_ar, cont_char_masked, name_of_intervals, name_of_cont_char, params, outpath):
    #Plot it.
    pos=[1,2]
    dataset=[intervals_param_ar, cont_char_masked]
    #print(len(pos))
    #print(len(dataset))
    #Violin plots for whole-genome signal vs signal of intervals (e.g. peaks).
    fig=plt.figure(figsize=(5.5,10), dpi=100)
    plt1=fig.add_subplot(1,1,1) 
    violins=plt1.violinplot(dataset, positions=pos, widths=0.9, showmeans=True, showmedians=True, points=200)
    #print(violins)
    for vio in violins['bodies']:
        vio.set_facecolor('#ff7762')
        vio.set_edgecolor('black')
        vio.set_alpha(1)
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
    vbars.set_linewidth(2)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    labels=[f'{name_of_intervals}\nsites', 'Other\nsites']
    set_axis_style(plt1, labels)
    yticknames1=np.arange(params[0], params[1], params[2]) #EcTopoI 0, 35, 1, GC% 0, 100, 10, RNApol 0, 21, 2
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylim(params[3], params[4]) #EcTopoI -0.2, 4, GC% 10, 77 RNApol -1, 21
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=35)   
    plt1.annotate(f'Mean {name_of_cont_char}\nsignal={round(np.mean(intervals_param_ar),2)}', xy=(params[5], params[6]), xycoords='data', size=35, rotation=90) #EcTopoI 0.4, 3.8, GC% 0.3, 35, RNApol 0.4, 11
    plt1.annotate(f'Mean {name_of_cont_char}\nsignal={round(np.mean(cont_char_masked),2)}', xy=(params[7], params[6]), xycoords='data', size=35, rotation=90) #EcTopoI 1.4, 3.8, GC% 1.3, 35, RNApol 1.4, 11
    
    Intervals_stat=stats.ttest_ind(cont_char_masked, intervals_param_ar)
    print(f'\nT-test for all genome sites vs {name_of_intervals} sites {name_of_cont_char} FE means\n' + 'p-value=' + str(Intervals_stat[1]) +'\n' + 't-statistic=' + str(Intervals_stat[0]) + '\n')    
    #plt.show()
    plt.savefig(f'{outpath}{name_of_cont_char}_signal_of_{name_of_intervals}_sites_and_genome.png', dpi=400, figsize=(5.5, 10)) 
    plt.close()    
    return

#######
#2D plot.
#######

def plot_2D(intervals_char_1, intervals_char_2, char_name_1, char_name_2, intervals_name, method, outpath):
    ##Fitting.
    if method=='lin':
        #Linear fitting of linear data.
        fit=np.polyfit(intervals_char_1, intervals_char_2, 1)
        print(fit)
        fit_fn=np.poly1d(fit)  
    elif method=='log':
        #Linnear fitting of log data.
        fit=np.polyfit(intervals_char_1, np.log10(intervals_char_2), 1)
        print(fit)
        fit_fn=np.poly1d(fit)     
    
    ##Pearson correlation.
    if method=='lin':
        #Pearson linear data.
        pearson_cor=scipy.stats.pearsonr(intervals_char_1, intervals_char_2)
        print(f'Paerson correlation ({char_name_1}, {char_name_2}) for {intervals_name} {pearson_cor}')
    elif method=='log':
        #Pearson log data.
        pearson_cor=scipy.stats.pearsonr(intervals_char_1, np.log10(intervals_char_2))
        print(f'Paerson correlation ({char_name_1}, {char_name_2}) log for {intervals_name} {pearson_cor}')    
      
    #Plot data.
    fig, ax = plt.subplots()
    if method=='lin':
        ax.plot(intervals_char_1, intervals_char_2, 'ro')
    elif method=='log':
        ax.plot(intervals_char_1, np.log10(intervals_char_2), 'bo')
    ax.plot(intervals_char_1, fit_fn(intervals_char_1), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3))) 
    ax.annotate(f'Pearson correlation=\n{round(pearson_cor[0], 3)}', xy=(0.7, 0.2), xycoords='axes fraction', size=9)
    ax.set_xlabel(char_name_1)
    ax.set_ylabel(char_name_2)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_title(f'{char_name_1} vs {char_name_2} for {intervals_name} peaks')
    #plt.show()
    if method=='lin':
        plt.savefig(f'{outpath}\{char_name_1}_vs_{char_name_2}_for_{intervals_name}_peaks.png')
    elif method=='log':
        plt.savefig(f'{outpath}\{char_name_1}_vs_{char_name_2}_log_for_{intervals_name}_peaks.png')
    return

#######
#Functions wrapper.
#######

def func_wrapper(intervals_sets_path_dict, cont_char_path_dict, outpath):
    intervals_sets_dict=broadpeak_pars(intervals_sets_path_dict)
    cont_char_dict=score_data_parser(cont_char_path_dict)
    #Process intervals.
    name_of_intervals_1='EcTopoI'
    intervals_width_1=width_of_intervals(intervals_sets_dict, name_of_intervals_1)
    name_of_intervals_2='RNApol'
    intervals_width_2=width_of_intervals(intervals_sets_dict, name_of_intervals_2)    
    #Process continuos data.
    name_of_cont_char_1='EcTopoI ms'
    name_of_cont_char_2='RNApol'
    name_of_cont_char_3='EcTopoI'
    name_of_cont_char_4='GC%'
    #EcTopoI intervals.
    intervals_param_ar_tg_1, cont_char_masked_1, intervals_param_ar_sp_1=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals_1, name_of_cont_char_1) 
    intervals_param_ar_tg_2, cont_char_masked_2, intervals_param_ar_sp_2=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals_1, name_of_cont_char_2)
    intervals_param_ar_tg_3, cont_char_masked_3, intervals_param_ar_sp_3=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals_1, name_of_cont_char_3)
    intervals_param_ar_tg_4, cont_char_masked_4, intervals_param_ar_sp_4=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals_1, name_of_cont_char_4)
    #RNApol intervals.
    intervals_param_ar_tg_2_1, cont_char_masked_2_1, intervals_param_ar_sp_2_1=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals_2, name_of_cont_char_1) 
    intervals_param_ar_tg_2_2, cont_char_masked_2_2, intervals_param_ar_sp_2_2=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals_2, name_of_cont_char_2)
    intervals_param_ar_tg_2_3, cont_char_masked_2_3, intervals_param_ar_sp_2_3=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals_2, name_of_cont_char_3)
    intervals_param_ar_tg_2_4, cont_char_masked_2_4, intervals_param_ar_sp_2_4=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals_2, name_of_cont_char_4)    
    ##2D plots.
    #EcTopoI intervals.
    #EcTopoI ms vs GC%
    plot_2D(intervals_param_ar_sp_1, intervals_param_ar_sp_4, name_of_cont_char_1, name_of_cont_char_4, name_of_intervals_1, 'log', outpath+'GC_vs_EcTopoI_vs_RNApol\\')
    #EcTopoI ms vs RNApol
    plot_2D(intervals_param_ar_sp_1, intervals_param_ar_sp_2, name_of_cont_char_1, name_of_cont_char_2, name_of_intervals_1, 'log', outpath+'GC_vs_EcTopoI_vs_RNApol\\') 
    #EcTopoI ms vs EcTopoI
    plot_2D(intervals_param_ar_sp_1, intervals_param_ar_sp_3, name_of_cont_char_1, name_of_cont_char_3, name_of_intervals_1, 'log', outpath+'GC_vs_EcTopoI_vs_RNApol\\')     
    #EcTopoI vs RNApol
    plot_2D(intervals_param_ar_sp_3, intervals_param_ar_sp_2, name_of_cont_char_3, name_of_cont_char_2, name_of_intervals_1, 'log', outpath+'GC_vs_EcTopoI_vs_RNApol\\') 
    #EcTopoI vs EcTopoI peaks width.
    plot_2D(intervals_param_ar_sp_3, intervals_width_1, name_of_cont_char_3, 'Peaks width', name_of_intervals_1, 'log', outpath+'GC_vs_EcTopoI_vs_RNApol\\')     
    #RNApol intervals.
    #EcTopoI ms vs GC%
    plot_2D(intervals_param_ar_sp_2_1, intervals_param_ar_sp_2_4, name_of_cont_char_1, name_of_cont_char_4, name_of_intervals_2, 'lin', outpath+'GC_vs_EcTopoI_vs_RNApol\\')
    #EcTopoI ms vs RNApol
    plot_2D(intervals_param_ar_sp_2_1, intervals_param_ar_sp_2_2, name_of_cont_char_1, name_of_cont_char_2, name_of_intervals_2, 'lin', outpath+'GC_vs_EcTopoI_vs_RNApol\\') 
    #EcTopoI ms vs EcTopoI
    plot_2D(intervals_param_ar_sp_2_1, intervals_param_ar_sp_2_3, name_of_cont_char_1, name_of_cont_char_3, name_of_intervals_2, 'lin', outpath+'GC_vs_EcTopoI_vs_RNApol\\')     
    #EcTopoI vs RNApol
    plot_2D(intervals_param_ar_sp_2_3, intervals_param_ar_sp_2_2, name_of_cont_char_3, name_of_cont_char_2, name_of_intervals_2, 'lin', outpath+'GC_vs_EcTopoI_vs_RNApol\\') 
    #EcTopoI vs EcTopoI peaks width.
    plot_2D(intervals_param_ar_sp_2_3, intervals_width_2, name_of_cont_char_3, 'Peaks width', name_of_intervals_2, 'lin', outpath+'GC_vs_EcTopoI_vs_RNApol\\')         
    ##Violin plots.
    #EcTopoI intervals.
    Intervals_occupation(intervals_param_ar_tg_2, cont_char_masked_2, name_of_intervals_1, name_of_cont_char_2, [0, 21, 2, -1, 21, 0.35, 11, 1.35], outpath)   
    #RNApol intervals.
    Intervals_occupation(intervals_param_ar_tg_2_3, cont_char_masked_2_3, name_of_intervals_2, name_of_cont_char_3, [0, 5, 1, -0.2, 4.5, 0.35, 2.1, 1.35], outpath)
    return

func_wrapper(path_to_intervals_sets, cont_characters_dict, Outputpath)

print('Script ended its work succesfully!') 