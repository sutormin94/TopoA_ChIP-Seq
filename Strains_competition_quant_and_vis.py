###############################################
##Dmitry Sutormin, 2021##
##Strains competition experiments quantification and visualization##

#Takes raw qPCR data.
#Tracks strains ratio dynamics.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import pearsonr, binom, ttest_ind

#######
#Import data for strains competition.
#######

#Path to the raw data.
BW25113_strains_compet="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\BW25113_topA_mutants\Strains_competiton\Competition_all_results.xlsx"
#Name of a worksheet.
WS_name_BW25113_delta11="wt_vs_delta11_tr"
#Mutant name.
Mutant_name='delta11'
#Path to the output plots.
Outpath_BW25113_delta11="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\BW25113_topA_mutants\Strains_competiton\BW25113_topA_delta11_vs_wt_tr.svg"


def strains_comp_analysis(pathin, worksheet, mut_name, outpath):
    
    #Read growth curves data.
    competition_data=pd.read_excel(pathin, sheet_name=worksheet, header=0, index_col=0)  
    print(competition_data)
    
    #Average technical replicates.
    Time_points_ar=[]
    for i in range(6):
        tech_replicates_data=competition_data.loc[:, [f'{i}_1', f'{i}_2', f'{i}_3']]
        tech_replicates_mean=tech_replicates_data.mean(axis=1)
        competition_data[i]=tech_replicates_mean
        Time_points_ar.append(i)
    
    print(competition_data)
    
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    
    #Calculate ratio topA_mut/wt.
    Mut_wt_ratio_dict={}
    
    for i in range(2):
        for j in range(3):
            #Get amplification parameters.
            len_wt=competition_data.loc[f'wt_{i+1}_{j+1}', 'Amplicon length']
            lam_wt=competition_data.loc[f'wt_{i+1}_{j+1}', 'Primers efficiency']
            len_mut=competition_data.loc[f'{mut_name}_{i+1}_{j+1}', 'Amplicon length']
            lam_mut=competition_data.loc[f'{mut_name}_{i+1}_{j+1}', 'Primers efficiency']   
            #Get qPCR Ct data.
            data_wt=competition_data.loc[f'wt_{i+1}_{j+1}', Time_points_ar]
            data_mut=competition_data.loc[f'{mut_name}_{i+1}_{j+1}', Time_points_ar]
            #Calculate mut/wt ratio.
            mut_wt_ratio=(len_wt*(lam_wt**data_wt))/(len_mut*(lam_mut**data_mut))
            #Scale by ratio in a starter mixture.
            initial_ratio=mut_wt_ratio[0]
            mut_wt_ratio_norm=mut_wt_ratio/initial_ratio
            #Visualize competition data.
            plot_1.plot(Time_points_ar, mut_wt_ratio_norm, '.-', label=f'Replicate {i+1}_{j+1}')
            #Keep data.
            Mut_wt_ratio_dict[f'{i+1}_{j+1}']=mut_wt_ratio_norm
            
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1)
    plot_1.spines["left"].set_linewidth(1)     
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('time, days', size=16)
    plt.ylabel(r'$topA\Delta11$/wt ratio', size=16)
    plt.legend(fontsize=10, frameon=False, markerscale=1, handlelength=1, handletextpad=0.3, loc='upper right')
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, figsize=(3, 2))  
    
    #Keep ratio data in a dataframe.
    Mut_wt_ratio_dataframe=pd.DataFrame.from_dict(Mut_wt_ratio_dict, orient='index', columns=Time_points_ar)   
    print(Mut_wt_ratio_dataframe)
    
    #Welch t-test.
    Intervals_stat=stats.ttest_ind(Mut_wt_ratio_dataframe.loc[:, 0], Mut_wt_ratio_dataframe.loc[:, 5], equal_var=False)
    print(f'\nT-test mut/wt ratio means\nMean1={round(np.mean(Mut_wt_ratio_dataframe.loc[:, 0]),2)} Mean2={round(np.mean(Mut_wt_ratio_dataframe.loc[:, 5]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    
    return

strains_comp_analysis(BW25113_strains_compet, WS_name_BW25113_delta11, Mutant_name, Outpath_BW25113_delta11)