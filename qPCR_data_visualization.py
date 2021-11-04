###############################################
##Dmitry Sutormin, 2021##
##qPCR data visualization##

#Takes table with Ct data and primers efficiency and makes barpolts.
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
from scipy.stats import pearsonr
from scipy.stats import binom


#################
### Primers calibration data analysis.
#################

#Plot data.
def qPCR_primers_calibration(dataframe, Primers_list, suptitle_text, outpath):
    
    #Read data.
    Conc_data=dataframe.loc['Template DNA concentration, log10', :].tolist()
    print(f'Primer pairs names: {Primers_list}')    
    Num_of_datasets=len(Primers_list)
    print(f'Number of primer pairs: {Num_of_datasets}')
    
    
    #Determine the dimensions of the picture.
    Num_columns=3
    Num_rows=int(Num_of_datasets/Num_columns)
    
    if (Num_of_datasets%Num_columns)!=0:
        Num_rows+=1
    if Num_rows==0:
        Num_columns=Num_of_datasets
    
    #Plot data.
    fig, plot_av=plt.subplots(Num_rows,Num_columns,figsize=(4*Num_columns,2.2*Num_rows), dpi=100)
    fig.suptitle(suptitle_text, size=15)
    
    for i in range(Num_of_datasets):
        
        #Plot points.
        Primers_pair=Primers_list[i]
        print(f'Plotting calibration of the {Primers_pair} pair')
        qC_data=dataframe.loc[Primers_pair, :]
        print(f'Primer pair position: {i}, Primer pair row number: {int(i/Num_columns)}, Primer pair column number: {i%Num_columns}')
        plot_av[int(i/Num_columns), i%Num_columns].scatter(Conc_data, qC_data.tolist(), s=2, color='k', edgecolors='black', linewidth=0.2, alpha=1, zorder=1000) 
        print(Conc_data, qC_data.tolist())
        
        #Filter data, remove nan.
        qC_data_ar=qC_data.tolist()
        Conc_data_filtered=[]
        qC_data_ar_filtered=[]
        for j in range(len(Conc_data)):
            if np.isnan(qC_data_ar[j])==False:
                Conc_data_filtered.append(Conc_data[j])
                qC_data_ar_filtered.append(qC_data_ar[j])
        
        #Linear fitting of calibration data.
        fit=np.polyfit(Conc_data_filtered, qC_data_ar_filtered, 1)
        print(fit)
        fit_fn=np.poly1d(fit) 
        plot_av[int(i/Num_columns), i%Num_columns].plot(Conc_data, fit_fn(Conc_data), '--b', linewidth=0.5, label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3))) 
        primer_effectiveness=10**(-1/fit[0])
        plot_av[int(i/Num_columns), i%Num_columns].annotate(f"$\lambda$={round(primer_effectiveness, 2)}", xy=(0.7, 0.5), xycoords='axes fraction', size=16)
        plot_av[int(i/Num_columns), i%Num_columns].set_ylabel('Cq', size=15)
        plot_av[int(i/Num_columns), i%Num_columns].set_xlabel('Relative template concentration', size=10)
        plot_av[int(i/Num_columns), i%Num_columns].set_xticks(np.unique(Conc_data), minor=False)
        plot_av[int(i/Num_columns), i%Num_columns].set_xticklabels([1, 10, 100, 1000], rotation=0, size=12)
        plot_av[int(i/Num_columns), i%Num_columns].set_title(f"Primers pair {Primers_pair}", size=12)
        plot_av[int(i/Num_columns), i%Num_columns].legend()
    
    #Remove empty plots.  
    if Num_of_datasets%Num_columns!=0:
        for j in range(Num_columns-(Num_of_datasets%Num_columns)):
            fig.delaxes(plot_av[Num_rows-1,Num_columns-1-j])
        
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    plt.savefig(outpath, dpi=300, size=(12,13))
    plt.close()

    return

Input_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
PC_data=pd.read_excel(Input_table, sheet_name='Sup_Table_14', header=None, index_col=0)
print(PC_data)

#topA_primers_set=['topA_1', 'topA_2', 'topA_3', 'topA_4', 'topA_5', 'topA_6', 'topA_7']
#qPCR_primers_calibration(PC_data, topA_primers_set, 'topA primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\TopA_primers_calibration.png")
#
#RpmH_primers_set=['rpmH_1', 'rpmH_2', 'rpmH_3', 'rpmH_4', 'rpmH_5', 'rpmH_6', 'rpmH_7', 'rpmH_8']
#qPCR_primers_calibration(PC_data, RpmH_primers_set, 'rpmH-rnpA-yidD-yidC cluster primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\RpmH_primers_calibration.png")
#
#BW_comp_set=['BW25113_wt', 'delta11', 'delta11_SPA', 'delta14', 'delta14_SPA']
#qPCR_primers_calibration(PC_data, BW_comp_set, 'BW25113 topA mutants competition primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\BW25113_topA_mut_primers_calibration.png")




#################
### ChIP-qPCR data analysis.
#################

#Compute fold enrichment.
def FE_calc_pair(ar_IP, ar_Mock, primers_eff, length):

    #Calculate FE for data points.
    FE_ar=[]
    for i in range(len(ar_IP)):
        FE_Rep_ar=[]
        for j in range(len(ar_IP[i])):
            FE=primers_eff[i]**(ar_Mock[i][j]-ar_IP[i][j])
            FE_Rep_ar.append(FE) 
        FE_ar.append(FE_Rep_ar) 
        
    return FE_ar

#Count mean of replicates.
def count_means(list_of_lists):
    list_of_means=[]
    
    for some_list in list_of_lists:
        List_mean=np.mean(some_list)
        list_of_means.append(List_mean)
    
    return list_of_means

#Compute fold enrichment mean and standard deviation.
def FE_mean_and_std(ar):
    Points_mean=[]
    Points_std=[]
    for points_set in ar:
        points_set_mean=np.mean(points_set)
        points_set_std=stats.sem(points_set)*1.96 #0.95 confidence interval for mean
        Points_mean.append(points_set_mean)
        Points_std.append(points_set_std)
    
    return Points_mean, Points_std

#Plot data.
def ChIP_qPCR_data_analysis(dataframe, outpath):
    
    ###
    ##Plot all Cts.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(4,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['dps', 'potF', 'dif', 'H_2394']
    print(len(Conditions))
    
    X_coords=[1,2,3,4]
    print(len(X_coords))
    
    Efficiency=dataframe.loc[:, 'Primers efficiency'].tolist()
    Length=dataframe.loc[:, 'Amplicon length'].tolist()    
    
    IP_data_means=[]
    Mock_data_means=[]
    IP_Rif_data_means=[]
    Mock_Rif_data_means=[]    
    for pair in Conditions:
        IP_data_means.append(dataframe.loc[pair, ['IP_-CTD_-Rif_1', 'IP_-CTD_-Rif_2', 'IP_-CTD_-Rif_3']].dropna().tolist())
        Mock_data_means.append(dataframe.loc[pair, ['Mock_-CTD_-Rif_1', 'Mock_-CTD_-Rif_2', 'Mock_-CTD_-Rif_3']].dropna().tolist())
        IP_Rif_data_means.append(dataframe.loc[pair, ['IP_-CTD_+Rif_1', 'IP_-CTD_+Rif_2', 'IP_-CTD_+Rif_3']].dropna().tolist())
        Mock_Rif_data_means.append(dataframe.loc[pair, ['Mock_-CTD_+Rif_1', 'Mock_-CTD_+Rif_2', 'Mock_-CTD_+Rif_3']].dropna().tolist())
        
    #Calculate fold enrichment of IP over mock samples.
    FE_data=FE_calc_pair(IP_data_means, Mock_data_means, Efficiency, Length)
    FE_Rif_data=FE_calc_pair(IP_Rif_data_means, Mock_Rif_data_means, Efficiency, Length)
    print(FE_data)
    
    #Count mean of replicates.
    FE_data_means=count_means(FE_data)
    FE_Rif_data_means=count_means(FE_Rif_data)
    print(FE_data_means)
    
    #Normalize by minimal FE value.
    FE_data=[in_ar/min(FE_data_means) for in_ar in FE_data]
    FE_Rif_data=[in_ar/min(FE_Rif_data_means) for in_ar in FE_Rif_data]
    print(FE_data)
    
    #Compute mean and standard deviation for fold enrichment.
    FE_mean, FE_points_std=FE_mean_and_std(FE_data)
    FE_Rif_mean, FE_Rif_points_std=FE_mean_and_std(FE_Rif_data)
    print(FE_mean)
    print(FE_points_std)
    
    #Plot data.
    print(np.asarray(X_coords)-0.2)
    plot_av.bar(np.asarray(X_coords)-0.2, FE_mean, yerr=FE_points_std, error_kw=dict(lw=1, capsize=3, capthick=1), width=0.3, color='#89d8fa', edgecolor='k', linewidth=0.6, label='-CTD/-Rif')
    plot_av.bar(np.asarray(X_coords)+0.2, FE_Rif_mean, yerr=FE_Rif_points_std, error_kw=dict(lw=1, capsize=3, capthick=1), width=0.3, color='#e8d844', edgecolor='k', linewidth=0.6, label='-CTD/+Rif')
    
    for i in range(len(FE_data)):
        for j in range(len(FE_data[i])):
            plot_av.plot(X_coords[i]-0.2, FE_data[i][j], 'ko', markersize=1)
    
    for i in range(len(FE_Rif_data)):
        for j in range(len(FE_Rif_data[i])):
            plot_av.plot(X_coords[i]+0.2, FE_Rif_data[i][j], 'ko', markersize=1)
    
    plot_av.set_ylabel('Fold enrichment\nIP/Mock', size=15)
    plot_av.set_yticks([1,10,100])
    plot_av.set_yticklabels([1,10,100], rotation=0, size=13)     
    plot_av.set_xticks(X_coords)
    plot_av.set_xticklabels(Conditions, rotation=0, size=13)   
    plot_av.tick_params(axis='x', which='major', pad=0.5)
    plot_av.set_yscale('log')
    plot_av.set_xlim([0,4.7])
    
    plt.legend(loc='upper right', frameon=False)
    plt.tight_layout()
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(12,3))
    
    return
    
#Path to qPCR data table.
ChIP_qPCR_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
ChIP_qPCR_data=pd.read_excel(ChIP_qPCR_table, sheet_name='Sup_Table_17', header=0, index_col=0)
print(ChIP_qPCR_data)    

ChIP_qPCR_data_analysis(ChIP_qPCR_data, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\EcTopoI_ChIP_qPCR_plot.png")




#################
### CTD overexpression effects. 3' RNA decay of long transcripts.
#################

#Compute fold enrichment.
def FE_calc(ar, primers_eff, length):

    #Calculate FE for data points.
    FE_ar=[]
    for i in range(len(ar)):
        replicas_data=[]
        for j in range(len(ar[i])):
            FE=(length[-1]*(primers_eff[-1]**ar[-1][j]))/(length[i]*(primers_eff[i]**ar[i][j]))
            replicas_data.append(FE)
        FE_ar.append(replicas_data) 
        
    return FE_ar


#Plot data.
def CTD_ind_RNA_decay_qPCR(dataframe, rows_ar, outpath):
    
    ###
    ##Plot all Cts.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(12,3), dpi=100)
    
    #Prepare x axis, extract data.
    Conditions=rows_ar
    
    X_coords_c=np.asarray(dataframe.loc[rows_ar, 'Distance'].tolist())
    X_coords_l=X_coords_c-25
    X_coords_r=X_coords_c+25
    
    Efficiency=dataframe.loc[rows_ar, 'Effectiveness'].tolist()
    Length=dataframe.loc[rows_ar, 'Length'].tolist()
    
    wt_data_points=[]
    IPTG_data_points=[]
    for pair in Conditions:
        wt_data_points.append(dataframe.loc[pair, ['wt1', 'wt2', 'wt3']].tolist())
        IPTG_data_points.append(dataframe.loc[pair, ['IPTG1', 'IPTG2', 'IPTG3']].tolist())
    
    #Compute fold enrichment for mean data points.
    FE_wt=FE_calc(wt_data_points, Efficiency, Length)
    FE_IPTG=FE_calc(IPTG_data_points, Efficiency, Length)
    
    #Compute mean and standard deviation for fold enrichment.
    Points_mean_wt, Points_std_wt=FE_mean_and_std(FE_wt)
    Points_mean_IPTG, Points_std_IPTG=FE_mean_and_std(FE_IPTG)
        
    #Plot data.
    plot_av.bar(X_coords_l, Points_mean_wt, yerr=Points_std_wt, error_kw=dict(lw=1, capsize=3, capthick=1), width=50, color='#b2e69a', edgecolor='k', linewidth=0.6, label='-IPTG')
    plot_av.bar(X_coords_r, Points_mean_IPTG, yerr=Points_std_IPTG, error_kw=dict(lw=1, capsize=3, capthick=1), width=50, color='#89d8fa', edgecolor='k', linewidth=0.6, label='+IPTG 1mM')
    
    plot_av.plot(X_coords_l, FE_wt, 'ko', markersize=1)
    plot_av.plot(X_coords_r, FE_IPTG, 'ko', markersize=1)
    
    plot_av.set_ylabel('Fold enrichment', size=20)
    plot_av.set_xticks(X_coords_c, minor=True)
    plot_av.set_xticklabels(Conditions, minor=True, rotation=0, size=12)  
    plot_av.set_xticks(range(0, max(X_coords_c)+200, 500), minor=False)
    plot_av.set_xticklabels(range(0, max(X_coords_c)+200, 500), minor=False, rotation=90, size=7) 
    
    plot_av.tick_params(axis='x', which='minor', pad=1.5)
    plot_av.tick_params(axis='x', which='major', pad=15)
    
    plot_av.set_yticks(np.arange(0, max(Points_mean_wt)+0.2, 0.25), minor=False)
    plot_av.set_yticklabels(np.arange(0, max(Points_mean_wt)+0.2, 0.25), minor=False, rotation=0, size=12)  
    plot_av.spines["top"].set_visible(False)
    plot_av.spines["right"].set_visible(False)  
    plot_av.spines["bottom"].set_linewidth(1.5)
    plot_av.spines["left"].set_linewidth(1.5)    
    
    #Place legend outside of a graph. Taken from: https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    box=plot_av.get_position()
    plot_av.set_position([box.x0, box.y0, box.width * 0.95, box.height])    
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.75), frameon=False)
    
    plt.tight_layout(rect=[0,0,0.95,1])
    plt.show()
    plt.savefig(outpath, dpi=300, size=(12,3))
    
    return

#Path to qPCR data table.
Input_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
Input_dec=pd.read_excel(Input_table, sheet_name='Sup_Table_16', header=0, index_col=0)
print(Input_dec)

#topA_points=['topA7', 'topA6', 'topA5', 'topA4', 'topA3', 'topA2', 'topA1']
#CTD_ind_RNA_decay_qPCR(Input_dec, topA_points, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\TopA_CTD_expression_effect_qPCR.svg")
#rpmH_points=['rpmH8', 'rpmH7', 'rpmH6', 'rpmH5', 'rpmH4', 'rpmH3', 'rpmH2', 'rpmH1']
#CTD_ind_RNA_decay_qPCR(Input_dec, rpmH_points, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\RpmH_CTD_expression_effect_qPCR.svg")




#################
### BW25113 wt and topA mutants co-cultivation competition.
#################

def strains_comp_analysis(pathin, worksheet, wt_set_name, mut_name, outpath):
    
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
        
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    
    #Calculate ratio topA_mut/wt.
    Mut_wt_ratio_dict={}
    
    for i in range(int(len(competition_data.index)/4)):
        #Get amplification parameters.
        len_wt=competition_data.loc[f'{wt_set_name}_{i+1}', 'Amplicon length']
        lam_wt=competition_data.loc[f'{wt_set_name}_{i+1}', 'Primers efficiency']
        len_mut=competition_data.loc[f'{mut_name}_{i+1}', 'Amplicon length']
        lam_mut=competition_data.loc[f'{mut_name}_{i+1}', 'Primers efficiency']   
        #Get qPCR Ct data.
        data_wt=competition_data.loc[f'{wt_set_name}_{i+1}', Time_points_ar]
        data_mut=competition_data.loc[f'{mut_name}_{i+1}', Time_points_ar]
        #Calculate mut/wt ratio.
        mut_wt_ratio=(len_wt*(lam_wt**data_wt))/(len_mut*(lam_mut**data_mut))
        #Scale by ratio in a starter mixture.
        initial_ratio=mut_wt_ratio[0]
        mut_wt_ratio_norm=mut_wt_ratio/initial_ratio
        #Visualize competition data.
        plot_1.plot(Time_points_ar, mut_wt_ratio_norm, '.-', label=f'Replicate {i+1}')
        #Keep data.
        Mut_wt_ratio_dict[f'{i+1}']=mut_wt_ratio_norm
        
    #Keep ratio data in a dataframe.
    Mut_wt_ratio_dataframe=pd.DataFrame.from_dict(Mut_wt_ratio_dict, orient='index', columns=Time_points_ar)   
    
    #Welch t-test.
    Intervals_stat=stats.ttest_ind(Mut_wt_ratio_dataframe.loc[:, 0], Mut_wt_ratio_dataframe.loc[:, 5], equal_var=False)
    print(f'\nT-test mut/wt ratio means\nMean {wt_set_name}={round(np.mean(Mut_wt_ratio_dataframe.loc[:, 0]),2)} Mean {mut_name}={round(np.mean(Mut_wt_ratio_dataframe.loc[:, 5]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    plot_1.bar(0, np.mean(Mut_wt_ratio_dataframe.loc[:, 0]), yerr=stats.sem(Mut_wt_ratio_dataframe.loc[:, 0]), error_kw=dict(lw=0.7, capsize=2, capthick=0.7), align='center', width=0.3, color='blue', alpha=0.4, edgecolor='k', linewidth=0.6, zorder=0)
    plot_1.bar(5, np.mean(Mut_wt_ratio_dataframe.loc[:, 5]), yerr=stats.sem(Mut_wt_ratio_dataframe.loc[:, 5]), error_kw=dict(lw=0.7, capsize=2, capthick=0.7), align='center', width=0.3, color='red', alpha=0.4, edgecolor='k', linewidth=0.6, zorder=1)
     
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1)
    plot_1.spines["left"].set_linewidth(1)     
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('time, days', size=16)
    plt.ylabel(r'$topA\Delta$' + mut_name.lstrip("delta") + '/wt ratio', size=16)
    plt.legend(fontsize=10, frameon=False, markerscale=1, handlelength=1, handletextpad=0.3, loc='upper right')
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, figsize=(3, 2)) 
    plt.close()

    return

#Path to the raw data.
BW25113_strains_compet="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
#Name of a worksheet.
WS_name_BW25113="Sup_Table_15"

#Path to the output plots.
Outpath_BW25113_delta11="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\TEST_BW25113_topA_delta11_vs_wt.svg"
#strains_comp_analysis(BW25113_strains_compet, WS_name_BW25113, 'wt_1', 'delta11', Outpath_BW25113_delta11)

#Path to the output plots.
Outpath_BW25113_delta14="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\TEST_BW25113_topA_delta14_vs_wt.svg"
#strains_comp_analysis(BW25113_strains_compet, WS_name_BW25113, 'wt_2', 'delta14', Outpath_BW25113_delta14)