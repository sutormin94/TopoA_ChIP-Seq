###############################################
##Dmitry Sutormin, 2021##
##Growth curves plot##

#Plots growth curve data with standard error intervals.
###############################################

#######
#Packages to be imported.
#######

import random as rd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import scipy
from scipy.stats import norm
import pandas as pd


#########
## Import data for BW25113 topA mutants: topA delta11, topA delta14, topA delta30.
#########

#Path to the raw data.
Growth_curves_data_BW25113_mutants="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
#Name of a worksheet.
WS_name_BW25113_mutants="Sup_Table_12"
#Path to the output plots.
Outpath_BW25113_mutants="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Growth_curves\BW25113_topA_mutants_growth_curve.png"
Outpath1_BW25113_mutants="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Growth_curves\BW25113_topA_mutants_growth_rate_doubling_time.png"


def prepare_data(data_frame, columns_list):
    #Take data.
    replicates_data=data_frame.loc[:, columns_list]
    #Add mean.
    replicates_data['Mean']=replicates_data.loc[:, columns_list].mean(axis=1, skipna=True) 
    #Add standard error of mean.
    replicates_data['CI95']=replicates_data.loc[:, columns_list].std(axis=1, skipna=True)*1.96/np.sqrt(len(columns_list))
    #Add confidential intervals.
    replicates_data['Mean_up']=replicates_data['Mean']+replicates_data['CI95']
    replicates_data['Mean_dn']=replicates_data['Mean']-replicates_data['CI95']    
    return replicates_data


def Plot_growth_curve_BW25113_mutants(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    WT=['BW25113_wt1', 'BW25113_wt2', 'BW25113_wt3']
    topA66=['BW25113_delta11_1', 'BW25113_delta11_2', 'BW25113_delta11_3']
    delta14=['BW25113_delta14_1', 'BW25113_delta14_2', 'BW25113_delta14_3']
    delta30=['BW25113_delta30_1', 'BW25113_delta30_2', 'BW25113_delta30_3']
    
    Sets=[WT, topA66, delta14, delta30]
    
    All_dataframes=[]
    for deta_set in Sets:
        All_dataframes.append(prepare_data(gc_data, deta_set))

    Colors_list=['#e07a5f', '#3d405b', '#81b29a', '#f2cc8f']
    Fill_color_list=['#e07a5f', '#3d405b', '#81b29a', '#f2cc8f']
    Label_list=[r'$\it{wt}$', r'$topA\Delta11$', r'$topA\Delta14$', r'$topA\Delta30$']
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    for i in range(len(All_dataframes)):
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean, Colors_list[i], linewidth=2, label=Label_list[i])
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_up, Colors_list[i], linewidth=1, alpha=0.3)
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_dn, Colors_list[i], linewidth=1, alpha=0.3)
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Colors_list[i], alpha=0.3, linewidth=0.2)
        
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1.5)
    plot_1.spines["left"].set_linewidth(1.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('time, min', size=16)
    plt.ylabel('OD$_{600}$', size=16)
    plt.legend(fontsize=10, frameon=False, markerscale=2, handlelength=1, handletextpad=0.3)
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, figsize=(3, 2))
    
    return

Plot_growth_curve_BW25113_mutants(Growth_curves_data_BW25113_mutants, WS_name_BW25113_mutants, Outpath_BW25113_mutants)


def growth_rate_BW25113_mutants(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    WT=['BW25113_wt1', 'BW25113_wt2', 'BW25113_wt3']
    topA66=['BW25113_delta11_1', 'BW25113_delta11_2', 'BW25113_delta11_3']
    delta14=['BW25113_delta14_1', 'BW25113_delta14_2', 'BW25113_delta14_3']
    delta30=['BW25113_delta30_1', 'BW25113_delta30_2', 'BW25113_delta30_3']
    
    #Set time points.
    Time_points=[20, 40, 60, 75, 90, 110, 130, 150]
    
    Sets=[WT, topA66, delta14, delta30]
    
    All_dataframes=[]
    All_dataframes_fits=[]
    for data_set in Sets:
        
        #Select and log transform data.
        Get_data=prepare_data(gc_data, data_set)
        Get_data=Get_data.loc[Time_points, :]
        Get_data=Get_data.apply(np.log, axis=0)*(1/np.log(2))
        All_dataframes.append(Get_data)
        
        #Fit data.
        k=[]
        b=[]
        for replicate in data_set:
            print(replicate)
            fit=np.polyfit(Get_data.index.tolist(), list(Get_data[replicate]), 1)
            print(fit)
            k.append(fit[0])
            b.append(fit[1])
        
        k_mean=np.mean(k)
        b_mean=np.mean(b)
        
        All_dataframes_fits.append([k, k_mean, b, b_mean]) 

    Colors_list=['#e07a5f', '#3d405b', '#81b29a', '#f2cc8f']
    Fill_color_list=['#e07a5f', '#3d405b', '#81b29a', '#f2cc8f']
    Label_list=[r'$\it{wt}$', r'$topA\Delta11$', r'$topA\Delta14$', r'$topA\Delta30$']   
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    for i in range(len(All_dataframes)):
        fit=[All_dataframes_fits[i][1], All_dataframes_fits[i][3]]
        fit_fn=np.poly1d(fit)
        plot_1.plot(All_dataframes[i].index, fit_fn(All_dataframes[i].index), color=Colors_list[i], linestyle='dashed', linewidth=1) 
        division_time=1/fit[0]  
        print(f'Doubling time: {division_time} min')
        
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean, Colors_list[i], linewidth=2, label=Label_list[i] + r' $\tau$=' + f'{int(division_time)} min')
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_up, Colors_list[i], linewidth=1, alpha=0.3)
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_dn, Colors_list[i], linewidth=1, alpha=0.3)
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Colors_list[i], alpha=0.3, linewidth=0.2)
           
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1.5)
    plot_1.spines["left"].set_linewidth(1.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('time, min', size=16)
    plt.ylabel('log$_{2}$(OD$_{600}$)', size=16)
    plt.legend(fontsize=10, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc="lower right", bbox_to_anchor=(0.58,0.61))
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, figsize=(3, 2))    
    
    return

growth_rate_BW25113_mutants(Growth_curves_data_BW25113_mutants, WS_name_BW25113_mutants, Outpath1_BW25113_mutants)



#########
## Import data for overexpression of EcTopoI G116S M320V from the pBAD33 plasmid.
#########

#Path to the raw data.
Growth_curves_data_pBAD33_topA_mut="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
#Name of a worksheet.
WS_name_pBAD33_topA_mut="Sup_Table_11"
#Path to the output plots.
Outpath_pBAD33_topA_mut="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Growth_curves\DY330_pBAD33_EcTopoI_G116S_M320V_induction_growth_curve.png"


def Plot_growth_curve_pBAD33_topA_mut(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    no_Ara=['EcTopoI_mut_Replicate_1', 'EcTopoI_mut_Replicate_2', 'EcTopoI_mut_Replicate_3']
    plus_Ara=['EcTopoI_mut_Replicate_1_Ara', 'EcTopoI_mut_Replicate_2_Ara', 'EcTopoI_mut_Replicate_3_Ara']
    
    Sets=[no_Ara, plus_Ara]
    
    All_dataframes=[]
    for deta_set in Sets:
        All_dataframes.append(prepare_data(gc_data, deta_set))

    Colors_list=['#000000', '#ACCFFF']
    Fill_color_list=['#000000', '#0000FF']
    Label_list=['pBAD33 EcTopoI\nG116S M320V\n-Ara', 'pBAD33 EcTopoI\nG116S M320V\n+Ara 10mM']
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    for i in range(len(All_dataframes)):
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean, Colors_list[i], linewidth=2, label=Label_list[i])
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_up, Fill_color_list[i], linewidth=1, alpha=0.3)
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_dn, Fill_color_list[i], linewidth=1, alpha=0.3)
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Fill_color_list[i], alpha=0.3, linewidth=0.2)
        
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1.5)
    plot_1.spines["left"].set_linewidth(1.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('time, min', size=16)
    plt.ylabel('OD$_{600}$', size=16)
    plt.legend(fontsize=10, frameon=False, markerscale=2, handlelength=1, handletextpad=0.3)
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, figsize=(3, 2))
    
    return

Plot_growth_curve_pBAD33_topA_mut(Growth_curves_data_pBAD33_topA_mut, WS_name_pBAD33_topA_mut, Outpath_pBAD33_topA_mut)



#########
## Import data for overexpression of EcTopoI from the pCA25 plasmid.
#########

#Path to the raw data.
Growth_curves_data_pCA25_topA="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
#Name of a worksheet.
WS_name_pCA25_topA="Sup_Table_11"
#Path to the output plots.
Outpath_pCA25_topA="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Growth_curves\DY330_EcTopoI_induction_growth_curve.png"


def Plot_growth_curve_pCA25_topA(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    no_Ara=['EcTopoI_Replicate_1', 'EcTopoI_Replicate_2', 'EcTopoI_Replicate_3']
    plus_Ara=['EcTopoI_Replicate_1_IPTG', 'EcTopoI_Replicate_2_IPTG', 'EcTopoI_Replicate_3_IPTG']
    
    Sets=[no_Ara, plus_Ara]
    
    All_dataframes=[]
    for deta_set in Sets:
        All_dataframes.append(prepare_data(gc_data, deta_set))

    Colors_list=['#000000', '#FF4500']
    Fill_color_list=['#000000', '#FF4500']
    Label_list=['pCA25 topA -IPTG', 'pCA25 topA +IPTG 1mM']
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    for i in range(len(All_dataframes)):
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean, Colors_list[i], linewidth=2, label=Label_list[i])
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_up, Fill_color_list[i], linewidth=1, alpha=0.3)
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_dn, Fill_color_list[i], linewidth=1, alpha=0.3)
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Fill_color_list[i], alpha=0.3, linewidth=0.2)
        
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1.5)
    plot_1.spines["left"].set_linewidth(1.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('time, min', size=16)
    plt.ylabel('OD$_{600}$', size=16)
    plt.legend(fontsize=10, frameon=False, markerscale=2, handlelength=1, handletextpad=0.3)
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, figsize=(3, 2))
    
    return

Plot_growth_curve_pCA25_topA(Growth_curves_data_pCA25_topA, WS_name_pCA25_topA, Outpath_pCA25_topA)



#########
## Import data for overexpression of EcTopoI CTD from the pCA25 plasmid.
#########

#Path to the raw data.
Growth_curves_data_pCA25_EcTopoI_CTD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
#Name of a worksheet.
WS_name_pCA25_EcTopoI_CTD="Sup_Table_11"
#Path to the output plots.
Outpath_pCA25_EcTopoI_CTD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Growth_curves\DY330_CTD_induction_growth_curve.png"


def Plot_growth_curve_pCA25_CTD(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    no_Ara=['CTD_Replicate_1', 'CTD_Replicate_2', 'CTD_Replicate_3']
    plus_Ara=['CTD_Replicate_1_IPTG', 'CTD_Replicate_2_IPTG', 'CTD_Replicate_3_IPTG']
    
    Sets=[no_Ara, plus_Ara]
    
    All_dataframes=[]
    for deta_set in Sets:
        All_dataframes.append(prepare_data(gc_data, deta_set))

    Colors_list=['#000000', '#0000FF']
    Fill_color_list=['#000000', '#0000FF']
    Label_list=['pCA25 14kDa CTD -IPTG', 'pCA25 14kDa CTD +IPTG 1mM']
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    for i in range(len(All_dataframes)):
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean, Colors_list[i], linewidth=2, label=Label_list[i])
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_up, Fill_color_list[i], linewidth=1, alpha=0.3)
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_dn, Fill_color_list[i], linewidth=1, alpha=0.3)
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Fill_color_list[i], alpha=0.3, linewidth=0.2)
        
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1.5)
    plot_1.spines["left"].set_linewidth(1.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('time, min', size=16)
    plt.ylabel('OD$_{600}$', size=16)
    plt.legend(fontsize=10, frameon=False, markerscale=2, handlelength=1, handletextpad=0.3)
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, figsize=(3, 2))
    
    return

Plot_growth_curve_pCA25_CTD(Growth_curves_data_pCA25_EcTopoI_CTD, WS_name_pCA25_EcTopoI_CTD, Outpath_pCA25_EcTopoI_CTD)



#########
## Import data for overexpression of GFP from the pCA25 plasmid.
#########

#Path to the raw data.
Growth_curves_data_pCA25_gfp="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables.xlsx"
#Name of a worksheet.
WS_name_pCA25_gfp="Sup_Table_11"
#Path to the output plots.
Outpath_pCA25_gfp="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Growth_curves\DY330_GFP_induction_growth_curve.png"


def Plot_growth_curve_pCA25_gfp(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    no_Ara=['GFP_Replicate_1', 'GFP_Replicate_2', 'GFP_Replicate_3']
    plus_Ara=['GFP_Replicate_1_IPTG', 'GFP_Replicate_2_IPTG', 'GFP_Replicate_3_IPTG']
    
    Sets=[no_Ara, plus_Ara]
    
    All_dataframes=[]
    for deta_set in Sets:
        All_dataframes.append(prepare_data(gc_data, deta_set))

    Colors_list=['#000000', '#008000']
    Fill_color_list=['#000000', '#008000']
    Label_list=['pCA25 GFP -IPTG', 'pCA25 GFP +IPTG 1mM']
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    for i in range(len(All_dataframes)):
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean, Colors_list[i], linewidth=2, label=Label_list[i])
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_up, Fill_color_list[i], linewidth=1, alpha=0.3)
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_dn, Fill_color_list[i], linewidth=1, alpha=0.3)
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Fill_color_list[i], alpha=0.3, linewidth=0.2)
        
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1.5)
    plot_1.spines["left"].set_linewidth(1.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('time, min', size=16)
    plt.ylabel('OD$_{600}$', size=16)
    plt.legend(fontsize=10, frameon=False, markerscale=2, handlelength=1, handletextpad=0.3)
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, figsize=(3, 2))
    
    return

Plot_growth_curve_pCA25_gfp(Growth_curves_data_pCA25_gfp, WS_name_pCA25_gfp, Outpath_pCA25_gfp)