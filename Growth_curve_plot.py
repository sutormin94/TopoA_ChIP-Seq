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
#Import data.
#######

#Path to the raw data.
Growth_curves_data="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Poisoned_TopoI_mutant\Growth_curves\Growth_curves.xlsx"

#Name of a worksheet.
WS_name="DY330pBAD33_TopoI_G116S_M320Vtr"

#Path to the output plots.
Outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Poisoned_TopoI_mutant\Growth_curves\DY330_pBAD33_EcTopoI_G116S_M320V_induction_growth_curve_trunk_time.png"


def Plot_growth_curve(data_inpath, sheetname, outpath):
    #Read growth curves data.
    gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Retrive non induced data and induced data.
    Non_induced_cols=['Replicate_1', 'Replicate_2', 'Replicate_3']
    Induced_cols=['Replicate_1_Ara', 'Replicate_2_Ara', 'Replicate_3_Ara']
    
    NI_data=gc_data.loc[:, Non_induced_cols]
    I_data=gc_data.loc[:, Induced_cols]
    
    #Add mean column.
    NI_data['Mean']=NI_data.loc[:, Non_induced_cols].mean(axis=1, skipna=True) 
    I_data['Mean']=I_data.loc[:, Induced_cols].mean(axis=1, skipna=True)  
    
    #Add 95% two-tail confident interval range.
    NI_data['CI95']=NI_data.loc[:, Non_induced_cols].std(axis=1, skipna=True)*1.96/np.sqrt(len(Non_induced_cols))
    I_data['CI95']=I_data.loc[:, Induced_cols].std(axis=1, skipna=True)*1.96/np.sqrt(len(Induced_cols))
    NI_data['Mean_up']=NI_data['Mean']+NI_data['CI95']
    NI_data['Mean_dn']=NI_data['Mean']-NI_data['CI95']
    I_data['Mean_up']=I_data['Mean']+I_data['CI95']
    I_data['Mean_dn']=I_data['Mean']-I_data['CI95']    
    
    #Plot data.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    plot_1.plot(NI_data.index, NI_data.Mean, 'k', linewidth=2, label="pBAD33 topA_mut -Ara")
    plot_1.plot(NI_data.index, NI_data.Mean_up, 'k', linewidth=1, alpha=0.3)
    plot_1.plot(NI_data.index, NI_data.Mean_dn, 'k', linewidth=1, alpha=0.3)
    plot_1.fill_between(NI_data.index, NI_data.Mean_dn, NI_data.Mean_up, color='k', alpha=0.3, linewidth=0)
    plot_1.plot(I_data.index, I_data.Mean, '#accfff', linewidth=2, label="pBAD33 topA_mut +Ara 10mM")  
    plot_1.plot(I_data.index, I_data.Mean_up, '#accfff', linewidth=1, alpha=0.3)
    plot_1.plot(I_data.index, I_data.Mean_dn, '#accfff', linewidth=1, alpha=0.3)
    plot_1.fill_between(I_data.index, I_data.Mean_dn, I_data.Mean_up, color='b', alpha=0.3, linewidth=0)
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


Plot_growth_curve(Growth_curves_data, WS_name, Outpath)