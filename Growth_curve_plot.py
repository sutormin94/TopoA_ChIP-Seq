###############################################
##Dmitry Sutormin, 2020##
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


#######
#Import data for BW25113 topA mutants: topA delta11, topA delta14, topA delta30.
#######

#Path to the raw data.
Growth_curves_data_BW25113_mutants="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\BW25113_topA_mutants\Growth_curves_BW25113_topA_mutants.xlsx"
#Name of a worksheet.
WS_name_BW25113_mutants="Raw_data"
#Path to the output plots.
Outpath_BW25113_mutants="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\BW25113_topA_mutants\BW25113_topA_mutants_growth_curve.png"
Outpath1_BW25113_mutants="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\BW25113_topA_mutants\BW25113_topA_mutants_growth_rate_doubling_time.png"


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
    WT=['wt1', 'wt2', 'wt3']
    topA66=['delta11_1', 'delta11_2', 'delta11_3']
    delta14=['delta14_1', 'delta14_2', 'delta14_3']
    delta30=['delta30_1', 'delta30_2', 'delta30_3']
    
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
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Colors_list[i], alpha=0.3, linewidth=0)
        
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


#Plot_growth_curve_BW25113_mutants(Growth_curves_data_BW25113_mutants, WS_name_BW25113_mutants, Outpath_BW25113_mutants)


def growth_rate_BW25113_mutants(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    WT=['wt1', 'wt2', 'wt3']
    topA66=['delta11_1', 'delta11_2', 'delta11_3']
    delta14=['delta14_1', 'delta14_2', 'delta14_3']
    delta30=['delta30_1', 'delta30_2', 'delta30_3']
    
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
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Colors_list[i], alpha=0.3, linewidth=0)
        
        
        
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

#growth_rate_BW25113_mutants(Growth_curves_data_BW25113_mutants, WS_name_BW25113_mutants, Outpath1_BW25113_mutants)


#######
#Import data for BW25113 wt and topA delta14 responses on CTD o.e.
#######

#Path to the raw data.
Growth_curves_data_BW25113_o_e_CTD_wt_vs_delta14="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\Growth_curves\Growth_curves_and_topology.xlsx"
#Name of a worksheet.
WS_name_BW25113_o_e_CTD="BW25113_oeCTD_wt_vs_topAdelta14"
#Path to the output plots.
Outpath_BW25113_o_e_CTD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\Growth_curves\Growth_curves_plots\BW25113_wt_vs_topA_delta14_CTD_o_e_response_growth_curve.png"

def Plot_growth_curve_BW25113_mutants(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    WT=['wt_R1', 'wt_R2', 'wt_R3']
    WT_IPTG=['wt_IPTG_R1', 'wt_IPTG_R2', 'wt_IPTG_R3']
    delta14=['delta14_R1', 'delta14_R2', 'delta14_R3']
    delta14_IPTG=['delta14_IPTG_R1', 'delta14_IPTG_R2', 'delta14_IPTG_R3']
    
    Sets=[WT, WT_IPTG, delta14, delta14_IPTG]
    
    All_dataframes=[]
    for deta_set in Sets:
        All_dataframes.append(prepare_data(gc_data, deta_set))

    Colors_list=['#e07a5f', '#3d405b', '#81b29a', '#f2cc8f']
    Fill_color_list=['#e07a5f', '#3d405b', '#81b29a', '#f2cc8f']
    Label_list=[r'$\it{wt}$', r'$\it{wt}$'+' 1 mM IPTG', r'$topA\Delta14$', r'$topA\Delta14$'+'\n1 mM IPTG']
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    for i in range(len(All_dataframes)):
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean, Colors_list[i], linewidth=2, label=Label_list[i])
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_up, Colors_list[i], linewidth=1, alpha=0.3)
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_dn, Colors_list[i], linewidth=1, alpha=0.3)
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Colors_list[i], alpha=0.3, linewidth=0)
        
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


#Plot_growth_curve_BW25113_mutants(Growth_curves_data_BW25113_o_e_CTD_wt_vs_delta14, WS_name_BW25113_o_e_CTD, Outpath_BW25113_o_e_CTD)


#######
#Import data for overexpression of EcTopoI G116S M320V from pBAD33 plasmid.
#######

#Path to the raw data.
Growth_curves_data_pBAD33_topA_mut="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Poisoned_TopoI_mutant\Growth_curves\Growth_curves.xlsx"
#Name of a worksheet.
WS_name_pBAD33_topA_mut="DY330_pBAD33_TopoI_G116S_M320Vt"
#Path to the output plots.
Outpath_pBAD33_topA_mut="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Poisoned_TopoI_mutant\Growth_curves\DY330_pBAD33_EcTopoI_G116S_M320V_induction_growth_curve_2tr.png"


def Plot_growth_curve_pBAD33_topA_mut(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    no_Ara=['Replicate_1', 'Replicate_2', 'Replicate_3']
    plus_Ara=['Replicate_1_Ara', 'Replicate_2_Ara', 'Replicate_3_Ara']
    
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
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Fill_color_list[i], alpha=0.3, linewidth=0)
        
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

#Plot_growth_curve_pBAD33_topA_mut(Growth_curves_data_pBAD33_topA_mut, WS_name_pBAD33_topA_mut, Outpath_pBAD33_topA_mut)



#######
#Import data for overexpression of EcTopoI from pCA25 plasmid.
#######

#Path to the raw data.
Growth_curves_data_pCA25_topA="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\Growth_curves_and_topology.xlsx"
#Name of a worksheet.
WS_name_pCA25_topA="27_02_21_pCA25_EcTopoI"
#Path to the output plots.
Outpath_pCA25_topA="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\Growth_curves_plots\DY330_EcTopoI_induction_growth_curve.png"


def Plot_growth_curve_pCA25_topA(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    no_Ara=['Replicate_1', 'Replicate_2', 'Replicate_3']
    plus_Ara=['Replicate_1_IPTG', 'Replicate_2_IPTG', 'Replicate_3_IPTG']
    
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
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Fill_color_list[i], alpha=0.3, linewidth=0)
        
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

#Plot_growth_curve_pCA25_topA(Growth_curves_data_pCA25_topA, WS_name_pCA25_topA, Outpath_pCA25_topA)


#######
#Import data for E. coli DY330 topA-SPA and DY330 rpoC-TAP.
#######

#Path to the raw data.
Growth_curves_data="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\Growth_curves\Growth_curves_and_topology.xlsx"
#Name of a worksheet.
WS_name="DY330_topA_SPA_vs_rpoC_TAP"
#Path to the output plots.
Outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\Growth_curves\Growth_curves_plots\DY330_topA_SPA_vs_rpoC_TAP_growth_curve.png"


def Plot_growth_curve_DY330_topA_vs_rpoC(data_inpath, sheetname, outpath):
    
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    rpoC=['rpoC_TAP_1', 'rpoC_TAP_2', 'rpoC_TAP_3']
    topA=['topA_SPA_1', 'topA_SPA_2', 'topA_SPA_3']
    
    Sets=[rpoC, topA]
    
    All_dataframes=[]
    for deta_set in Sets:
        All_dataframes.append(prepare_data(gc_data, deta_set))

    Colors_list=['#e07a5f', '#3d405b']
    Fill_color_list=['#e07a5f', '#3d405b']
    Label_list=[r'$\it{rpoC-TAP}$', r'$\it{topA-SPA}$']
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    for i in range(len(All_dataframes)):
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean, Colors_list[i], linewidth=2, label=Label_list[i])
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_up, Colors_list[i], linewidth=1, alpha=0.3)
        plot_1.plot(All_dataframes[i].index, All_dataframes[i].Mean_dn, Colors_list[i], linewidth=1, alpha=0.3)
        plot_1.fill_between(All_dataframes[i].index, All_dataframes[i].Mean_dn, All_dataframes[i].Mean_up, color=Colors_list[i], alpha=0.3, linewidth=0)
        
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


Plot_growth_curve_DY330_topA_vs_rpoC(Growth_curves_data, WS_name, Outpath)