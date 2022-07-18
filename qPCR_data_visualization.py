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
import math
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

Input_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables\Sutormin_et_al.,2022,Supplementary_Tables.xlsx"
PC_data=pd.read_excel(Input_table, sheet_name='Table S7', header=None, index_col=0)
print(PC_data)

topA_primers_set=['topA_1', 'topA_2', 'topA_3', 'topA_4', 'topA_5', 'topA_6', 'topA_7']
qPCR_primers_calibration(PC_data, topA_primers_set, 'topA primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\TopA_primers_calibration.png")

RpmH_primers_set=['rpmH_1', 'rpmH_2', 'rpmH_3', 'rpmH_4', 'rpmH_5', 'rpmH_6', 'rpmH_7', 'rpmH_8']
qPCR_primers_calibration(PC_data, RpmH_primers_set, 'rpmH-rnpA-yidD-yidC cluster primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\RpmH_primers_calibration.png")

BW_comp_set=['BW25113_wt', 'delta11', 'delta11_SPA', 'delta14', 'delta14_SPA']
qPCR_primers_calibration(PC_data, BW_comp_set, 'BW25113 topA mutants competition primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\BW25113_topA_mut_primers_calibration.png")

sup_rel_primers_set=['rho', 'rnhB', 'rnhA', 'topB', 'parE', 'parC', 'gyrB', 'gyrA']
qPCR_primers_calibration(PC_data, sup_rel_primers_set, 'gyrAB parCE rnhAB topB rho primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\gyrAB_parCE_rnhAB_topB_rho_primers_calibration.png")




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
        points_set_mean=np.nanmean(points_set)
        points_set_std=np.nanstd(points_set)
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
ChIP_qPCR_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables\Sutormin_et_al.,2022,Supplementary_Tables.xlsx"
ChIP_qPCR_data=pd.read_excel(ChIP_qPCR_table, sheet_name='Table S12', header=0, index_col=0)
print(ChIP_qPCR_data)    

ChIP_qPCR_data_analysis(ChIP_qPCR_data, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\EcTopoI_ChIP_qPCR_plot.png")



##############
### EcTopoI 14kDa CTD overexpression effects. RpoC ChIP-qPCR data visualization, RNAP enrichement over long transcripts.
##############

#Compute fold enrichment.
def FE_calc_norm(ar_mock, ar_IP, primers_eff, length):

    #Calculate FE for data points IP/Mock.
    FE_ar=[]
    for i in range(len(ar_mock)):
        replicas_data=[]
        for j in range(len(ar_mock[i])):
            if (not math.isnan(ar_mock[i][j])) and (not math.isnan(ar_mock[i][j])):
                FE=(length[i]*(primers_eff[i]**ar_mock[i][j]))/(length[i]*(primers_eff[i]**ar_IP[i][j]))
                replicas_data.append(FE)
            else:
                replicas_data.append(np.nan)
        FE_ar.append(replicas_data) 
    
    #Normalize FE by FE in control site.
    FE_ar_norm=[]
    for i in range(len(FE_ar)):
        replicas_data_norm=[]
        for j in range(len(FE_ar[i])):
            if (not math.isnan(FE_ar[i][j])) and (not math.isnan(FE_ar[-1][j])):
                FE_norm=FE_ar[i][j]/FE_ar[-1][j]
                replicas_data_norm.append(FE_norm)
            else:
                replicas_data_norm.append(np.nan)
        FE_ar_norm.append(replicas_data_norm)
        
    return FE_ar_norm


#Plot data.
def topA_rpmH_ChIP_qPCR(dataframe, rows_ar, bar_width, outpath):
    
    ###
    ##Plot all Cts.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(12,3), dpi=100)
    
    #Prepare x axis, extract data.
    Conditions=rows_ar
    
    X_coords_c_base=np.asarray(dataframe.loc[rows_ar, 'Distance'].tolist())
    X_coords_c=X_coords_c_base+int(bar_width/2)+1
    X_coords_l=X_coords_c_base-int(bar_width/2)
    
    Efficiency=dataframe.loc[rows_ar, 'Effectiveness'].tolist()
    Length=dataframe.loc[rows_ar, 'Length'].tolist()
    
    wt_mock_data_points=[]
    wt_IP_data_points=[]
    CTD_mock_data_points=[]
    CTD_IP_data_points=[]
    
    for pair in Conditions:
        #-CTD/-Rif
        wt_mock_data_points.append(dataframe.loc[pair, ['wt_mock_1', 'wt_mock_2', 'wt_mock_3',  'wt_mock_4']].tolist())
        wt_IP_data_points.append(dataframe.loc[pair, ['wt_IP_1', 'wt_IP_2', 'wt_IP_3', 'wt_IP_4']].tolist())        
        #+CTD/-Rif
        CTD_mock_data_points.append(dataframe.loc[pair, ['CTD_mock_1',  'CTD_mock_2', 'CTD_mock_3', 'CTD_mock_4']].tolist())
        CTD_IP_data_points.append(dataframe.loc[pair, ['CTD_IP_1',  'CTD_IP_2', 'CTD_IP_3', 'CTD_IP_4']].tolist())     
    
    #Compute fold enrichment for mean data points.
    FE_norm_wt  = FE_calc_norm(wt_mock_data_points,  wt_IP_data_points,  Efficiency, Length)
    FE_norm_CTD = FE_calc_norm(CTD_mock_data_points, CTD_IP_data_points, Efficiency, Length)
    
    #Compute mean and standard deviation for fold enrichment.
    Points_mean_wt,  Points_std_wt  =FE_mean_and_std(FE_norm_wt)
    Points_mean_CTD, Points_std_CTD =FE_mean_and_std(FE_norm_CTD)
    
    #Compare means, stat.
    cleanedFE_norm_wt=[]
    for primer_pair in FE_norm_wt:
        parimer_pair_cleaned=[]
        for value in primer_pair:
            if math.isnan(value) != True:
                parimer_pair_cleaned.append(value)
        cleanedFE_norm_wt.append(parimer_pair_cleaned)
        
    cleanedFE_norm_CTD=[]
    for primer_pair in FE_norm_CTD:
        parimer_pair_cleaned=[]
        for value in primer_pair:
            if math.isnan(value) != True:
                parimer_pair_cleaned.append(value)
        cleanedFE_norm_CTD.append(parimer_pair_cleaned)
                
    Data_points=[cleanedFE_norm_wt, cleanedFE_norm_CTD]
    print(cleanedFE_norm_wt)
    print(cleanedFE_norm_CTD)
    for i in range(len(Data_points[0])):
        print(f'T-test for {Conditions[i]} : {stats.ttest_ind(Data_points[0][i], Data_points[1][i])}')    
        
    #Plot data.
    print('\n\n\n\n')
    print(X_coords_l)
    print(Points_mean_wt)
    print(FE_norm_wt)
    print('\n\n\n\n')
    plot_av.bar(X_coords_l, Points_mean_wt,  yerr=Points_std_wt,  error_kw=dict(lw=1, capsize=3, capthick=1), width=bar_width, color='#b2e69a', edgecolor='k', linewidth=0.6, label='-CTD')
    plot_av.bar(X_coords_c, Points_mean_CTD, yerr=Points_std_CTD, error_kw=dict(lw=1, capsize=3, capthick=1), width=bar_width, color='#89d8fa', edgecolor='k', linewidth=0.6, label='+CTD')    
    
    plot_av.plot(X_coords_l, FE_norm_wt, 'ko', markersize=1)
    plot_av.plot(X_coords_c, FE_norm_CTD, 'ko', markersize=1)
    
    plot_av.set_ylabel('Fold enrichment', size=20)
    plot_av.set_xticks(X_coords_c_base, minor=True)
    plot_av.set_xticklabels(Conditions, minor=True, rotation=0, size=12)  
    plot_av.set_xticks(range(0, max(X_coords_c_base)+200, 500), minor=False)
    plot_av.set_xticklabels(range(0, max(X_coords_c_base)+200, 500), minor=False, rotation=90, size=7) 
    
    plot_av.tick_params(axis='x', which='minor', pad=1.5)
    plot_av.tick_params(axis='x', which='major', pad=15) 
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
    plt.savefig(f'{outpath}.png', dpi=300, size=(12,3))
    plt.savefig(f'{outpath}.svg', dpi=300, size=(12,3))
    
    return

#Path to qPCR data table.
Input_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables\Sutormin_et_al.,2022,Supplementary_Tables.xlsx"
Input_dec=pd.read_excel(Input_table, sheet_name='Table S9', header=0, index_col=0)
print(Input_dec)

topA_points=['topA7', 'topA6', 'topA5', 'topA4', 'topA3', 'topA2', 'topA1']
topA_rpmH_ChIP_qPCR(Input_dec, topA_points, 50, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\TopA_CTD_Rif_RNAP_effect_ChIP_qPCR_topA1_norm_no_Rif")
rpmH_points=['rpmH8', 'rpmH7', 'rpmH6', 'rpmH5', 'rpmH4', 'rpmH3', 'rpmH2', 'rpmH1']
topA_rpmH_ChIP_qPCR(Input_dec, rpmH_points, 50, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\RpmH_CTD_Rif_RNAP_effect_ChIP_qPCR_all_reps_rpmH1_norm_no_Rif")



#################
### E. coli BW25113 wt and topA mutants co-cultivation competition.
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
BW25113_strains_compet="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables\Sutormin_et_al.,2022,Supplementary_Tables.xlsx"
#Name of a worksheet.
WS_name_BW25113="Table S16"

#Path to the output plots.
Outpath_BW25113_delta11="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\TEST_BW25113_topA_delta11_vs_wt.svg"
#strains_comp_analysis(BW25113_strains_compet, WS_name_BW25113, 'wt_1', 'delta11', Outpath_BW25113_delta11)

#Path to the output plots.
Outpath_BW25113_delta14="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\TEST_BW25113_topA_delta14_vs_wt.svg"
#strains_comp_analysis(BW25113_strains_compet, WS_name_BW25113, 'wt_2', 'delta14', Outpath_BW25113_delta14)



#################
### rnhAB gyrAB parCE topB expression response on EcTopoI Y319F mutant overexpression.
#################

#Compute fold enrichment.
def FE_calc_qPCR(ar, primers_eff, length):

    #Calculate FE for data points.
    FE_ar=[]
    for i in range(len(ar)):
        replicas_data=[]
        for j in range(len(ar[i])):
            FE=(length[-1]*(primers_eff[-1]**ar[-1][j]))/(length[i]*(primers_eff[i]**ar[i][j]))
            replicas_data.append(FE)
        FE_ar.append(replicas_data) 
        
    return FE_ar


#Path to qPCR data table.
data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables\Sutormin_et_al.,2022,Supplementary_Tables.xlsx"
data_tab=pd.read_excel(data_table, sheet_name='Table S8', header=0, index_col=0)
print(data_tab)


#Plot data.
def qPCR_Y319F_expression_response(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(6.5,3), dpi=100)
    
    #Prepare x axis.
    Genes_studied=['$\it{gyrA}$', '$\it{gyrB}$', '$\it{topB}$', '$\it{parC}$', '$\it{parE}$', '$\it{rnhA}$', '$\it{rnhB}$', '$\it{rho}$']
    print(len(Genes_studied))
    
    X_coords_main=[1.9,4.9,7.9,10.9,13.9,16.9,19.9, 22.9]
    
    Efficiency=dataframe.loc[:, 'Effectiveness'].tolist()
    Length=dataframe.loc[:, 'Length'].tolist()
    
    Conditions_rho_gyrAB=['gyrA2', 'gyrB1', 'rho2_R1']
    dataframe_rho_gyrAB=dataframe.loc[Conditions_rho_gyrAB, :]
    print(dataframe_rho_gyrAB)
    Conditions_rho_topB_parE_parC=['topB1', 'parC2', 'parE1', 'rho2_R2']
    dataframe_rho_topB_parE_parC=dataframe.loc[Conditions_rho_topB_parE_parC, :]
    print(dataframe_rho_topB_parE_parC)
    Conditions_rho_rnhAB=['rnhA1', 'rnhB1', 'rho2_R3']
    dataframe_rho_rnhAB=dataframe.loc[Conditions_rho_rnhAB, :]
    print(dataframe_rho_rnhAB)
    
    #rho and gyrAB
    Efficiency_rho_gyrAB=dataframe_rho_gyrAB.loc[:, 'Effectiveness'].tolist()
    Length_rho_gyrAB=dataframe_rho_gyrAB.loc[:, 'Length'].tolist()
    GFP_m_rho_gyrAB=[]
    GFP_p_rho_gyrAB=[]
    Y319F_m_rho_gyrAB=[]
    Y319F_p_rho_gyrAB=[]    
    for pair in Conditions_rho_gyrAB:
        GFP_m_rho_gyrAB.append(dataframe_rho_gyrAB.loc[pair, ['GFP-1', 'GFP-2', 'GFP-3', 'GFP-4', 'GFP-5', 'GFP-6']].tolist())
        GFP_p_rho_gyrAB.append(dataframe_rho_gyrAB.loc[pair, ['GFP+1', 'GFP+2', 'GFP+3', 'GFP+4', 'GFP+5', 'GFP+6']].tolist())
        Y319F_m_rho_gyrAB.append(dataframe_rho_gyrAB.loc[pair, ['Y319F-1', 'Y319F-2', 'Y319F-3', 'Y319F-4', 'Y319F-5', 'Y319F-6']].tolist())
        Y319F_p_rho_gyrAB.append(dataframe_rho_gyrAB.loc[pair, ['Y319F+1', 'Y319F+2', 'Y319F+3', 'Y319F+4', 'Y319F+5', 'Y319F+6']].tolist())        
    
    #Compute normalized fold enrichment for mean data points.
    Norm_GFP_m_rho_gyrAB=FE_calc_qPCR(GFP_m_rho_gyrAB, Efficiency_rho_gyrAB, Length_rho_gyrAB)
    Norm_GFP_p_rho_gyrAB=FE_calc_qPCR(GFP_p_rho_gyrAB, Efficiency_rho_gyrAB, Length_rho_gyrAB)
    Norm_Y319F_m_rho_gyrAB=FE_calc_qPCR(Y319F_m_rho_gyrAB, Efficiency_rho_gyrAB, Length_rho_gyrAB)
    Norm_Y319F_p_rho_gyrAB=FE_calc_qPCR(Y319F_p_rho_gyrAB, Efficiency_rho_gyrAB, Length_rho_gyrAB)
    print(Norm_GFP_m_rho_gyrAB, Norm_GFP_p_rho_gyrAB)
    
    #rho, topB, parE, and parC
    Efficiency_rho_topB_parE_parC=dataframe_rho_topB_parE_parC.loc[:, 'Effectiveness'].tolist()
    Length_rho_topB_parE_parC=dataframe_rho_topB_parE_parC.loc[:, 'Length'].tolist()
    GFP_m_rho_topB_parE_parC=[]
    GFP_p_rho_topB_parE_parC=[]
    Y319F_m_rho_topB_parE_parC=[]
    Y319F_p_rho_topB_parE_parC=[]      
    for pair in Conditions_rho_topB_parE_parC:
        GFP_m_rho_topB_parE_parC.append(dataframe_rho_topB_parE_parC.loc[pair, ['GFP-1', 'GFP-2', 'GFP-3', 'GFP-4', 'GFP-5', 'GFP-6']].tolist())
        GFP_p_rho_topB_parE_parC.append(dataframe_rho_topB_parE_parC.loc[pair, ['GFP+1', 'GFP+2', 'GFP+3', 'GFP+4', 'GFP+5', 'GFP+6']].tolist())
        Y319F_m_rho_topB_parE_parC.append(dataframe_rho_topB_parE_parC.loc[pair, ['Y319F-1', 'Y319F-2', 'Y319F-3', 'Y319F-4', 'Y319F-5', 'Y319F-6']].tolist())
        Y319F_p_rho_topB_parE_parC.append(dataframe_rho_topB_parE_parC.loc[pair, ['Y319F+1', 'Y319F+2', 'Y319F+3', 'Y319F+4', 'Y319F+5', 'Y319F+6']].tolist())          
    
    #Compute normalized fold enrichment for mean data points.
    Norm_GFP_m_rho_topB_parE_parC=FE_calc_qPCR(GFP_m_rho_topB_parE_parC, Efficiency_rho_topB_parE_parC, Length_rho_topB_parE_parC)
    Norm_GFP_p_rho_topB_parE_parC=FE_calc_qPCR(GFP_p_rho_topB_parE_parC, Efficiency_rho_topB_parE_parC, Length_rho_topB_parE_parC)
    Norm_Y319F_m_rho_topB_parE_parC=FE_calc_qPCR(Y319F_m_rho_topB_parE_parC, Efficiency_rho_topB_parE_parC, Length_rho_topB_parE_parC)
    Norm_Y319F_p_rho_topB_parE_parC=FE_calc_qPCR(Y319F_p_rho_topB_parE_parC, Efficiency_rho_topB_parE_parC, Length_rho_topB_parE_parC)    
    print(Norm_GFP_m_rho_topB_parE_parC, Norm_GFP_p_rho_topB_parE_parC)
    
    #rho and rnhAB
    Efficiency_rho_rnhAB=dataframe_rho_rnhAB.loc[:, 'Effectiveness'].tolist()
    Length_rho_rnhAB=dataframe_rho_rnhAB.loc[:, 'Length'].tolist()
    GFP_m_rho_rnhAB=[]
    GFP_p_rho_rnhAB=[]
    Y319F_m_rho_rnhAB=[]
    Y319F_p_rho_rnhAB=[]      
    for pair in Conditions_rho_rnhAB:
        GFP_m_rho_rnhAB.append(dataframe_rho_rnhAB.loc[pair, ['GFP-1', 'GFP-2', 'GFP-3', 'GFP-4', 'GFP-5', 'GFP-6']].tolist())
        GFP_p_rho_rnhAB.append(dataframe_rho_rnhAB.loc[pair, ['GFP+1', 'GFP+2', 'GFP+3', 'GFP+4', 'GFP+5', 'GFP+6']].tolist())
        Y319F_m_rho_rnhAB.append(dataframe_rho_rnhAB.loc[pair, ['Y319F-1', 'Y319F-2', 'Y319F-3', 'Y319F-4', 'Y319F-5', 'Y319F-6']].tolist())
        Y319F_p_rho_rnhAB.append(dataframe_rho_rnhAB.loc[pair, ['Y319F+1', 'Y319F+2', 'Y319F+3', 'Y319F+4', 'Y319F+5', 'Y319F+6']].tolist())        
    
    #Compute normalized fold enrichment for mean data points.
    Norm_GFP_m_rho_rnhAB=FE_calc_qPCR(GFP_m_rho_rnhAB, Efficiency_rho_rnhAB, Length_rho_rnhAB)
    Norm_GFP_p_rho_rnhAB=FE_calc_qPCR(GFP_p_rho_rnhAB, Efficiency_rho_rnhAB, Length_rho_rnhAB)
    Norm_Y319F_m_rho_rnhAB=FE_calc_qPCR(Y319F_m_rho_rnhAB, Efficiency_rho_rnhAB, Length_rho_rnhAB)
    Norm_Y319F_p_rho_rnhAB=FE_calc_qPCR(Y319F_p_rho_rnhAB, Efficiency_rho_rnhAB, Length_rho_rnhAB)    
    print(Norm_GFP_m_rho_rnhAB, Norm_GFP_p_rho_rnhAB)    
    
    #Compute mean and standard deviation for fold enrichment.
    Norm_GFP_m=Norm_GFP_m_rho_gyrAB[:-1]+Norm_GFP_m_rho_topB_parE_parC[:-1]+Norm_GFP_m_rho_rnhAB
    Norm_GFP_p=Norm_GFP_p_rho_gyrAB[:-1]+Norm_GFP_p_rho_topB_parE_parC[:-1]+Norm_GFP_p_rho_rnhAB
    Norm_Y319F_m=Norm_Y319F_m_rho_gyrAB[:-1]+Norm_Y319F_m_rho_topB_parE_parC[:-1]+Norm_Y319F_m_rho_rnhAB
    Norm_Y319F_p=Norm_Y319F_p_rho_gyrAB[:-1]+Norm_Y319F_p_rho_topB_parE_parC[:-1]+Norm_Y319F_p_rho_rnhAB
    
    Points_mean_GFP_m, Points_std_GFP_m=FE_mean_and_std(Norm_GFP_m)
    Points_mean_GFP_p, Points_std_GFP_p=FE_mean_and_std(Norm_GFP_p)
    Points_mean_Y319F_m, Points_std_Y319F_m=FE_mean_and_std(Norm_Y319F_m)
    Points_mean_Y319F_p, Points_std_Y319F_p=FE_mean_and_std(Norm_Y319F_p)  
    
    print(len(Norm_GFP_m), len(Norm_GFP_p), len(Norm_Y319F_m), len(Norm_Y319F_p))
    print(len(Points_mean_GFP_m), len(Points_mean_GFP_p), len(Points_mean_Y319F_m), len(Points_mean_Y319F_p))
    
    #Compare conditions.
    All_data_ar=[Norm_GFP_m, Norm_GFP_p, Norm_Y319F_m, Norm_Y319F_p]
    Experimental_conditions=['GFP-', 'GFP+', 'topA Y319F-', 'topA Y319F+']
    for k in range(len(All_data_ar[0])):
        for i in range(len(All_data_ar)):
            for j in range(len(All_data_ar)):
                if j>i:
                    qPCR_stat=stats.ttest_ind(All_data_ar[i][k], All_data_ar[j][k], equal_var=False, nan_policy='omit')
                    print(f'Test difference between: {Experimental_conditions[i]} and {Experimental_conditions[j]} for gene {Genes_studied[k]}')
                    print(f'Sample size: {len(All_data_ar[i][k])}, Sample size: {len(All_data_ar[j][k])}')   
                    print(f'\nT-test FE Mean1={round(np.mean(All_data_ar[i][k]),3)}; Mean2={round(np.mean(All_data_ar[j][k]),3)}\np-value={qPCR_stat[1]}\nt-statistic={qPCR_stat[0]}\n')                  
                    
        
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa']*len(Norm_GFP_m)
    print(len(Colors))
    
    X_coords_GFP_m=np.array(X_coords_main)-0.9
    X_coords_GFP_p=np.array(X_coords_main)-0.3
    X_coords_Y319F_m=np.array(X_coords_main)+0.3
    X_coords_Y319F_p=np.array(X_coords_main)+0.9
    
    #Plot data.
    Bars_GFP_m=plot_av.bar(X_coords_GFP_m, Points_mean_GFP_m, yerr=Points_std_GFP_m, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color='#b2e69a', edgecolor='k', linewidth=0.6)
    Bars_GFP_p=plot_av.bar(X_coords_GFP_p, Points_mean_GFP_p, yerr=Points_std_GFP_p, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color='#f598b8', edgecolor='k', linewidth=0.6)
    Bars_Y319F_m=plot_av.bar(X_coords_Y319F_m, Points_mean_Y319F_m, yerr=Points_std_Y319F_m, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color='#f5ab87', edgecolor='k', linewidth=0.6)
    Bars_Y319F_p=plot_av.bar(X_coords_Y319F_p, Points_mean_Y319F_p, yerr=Points_std_Y319F_p, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color='#89d8fa', edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_GFP_m, Norm_GFP_m, 'ko', markersize=1) 
    plot_av.plot(X_coords_GFP_p, Norm_GFP_p, 'ko', markersize=1) 
    plot_av.plot(X_coords_Y319F_m, Norm_Y319F_m, 'ko', markersize=1) 
    plot_av.plot(X_coords_Y319F_p, Norm_Y319F_p, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Genes_studied, rotation=0, size=14)  
    plot_av.set_yticklabels([0,1,2,3,4,5], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 5])
    
    plt.legend((Bars_GFP_m[0],Bars_GFP_p[0],Bars_Y319F_m[0],Bars_Y319F_p[0]), ('$\it{gfp}$-', '$\it{gfp}$+', '$\it{topA}$ Y319F-', '$\it{topA}$ Y319F+'), fontsize=14, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(6.5,3))   

    return

qPCR_Y319F_expression_response(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\qPCR\\rnhAB_gyrAB_parCE_topB_response_on_Y319F_expression_qPCR.png")
