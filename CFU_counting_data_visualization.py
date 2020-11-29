###############################################
##Dmitry Sutormin, 2020##
##CFU counting data visualization##

#Takes adtaframe with CFU counts and makes barplots.
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
### Input data.
#################


#Path to CFU counts.
Data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\TopoA_14kDa_CTD_expression\CFU_counting_data.xlsx"
Exp1=pd.read_excel(Data_table, sheet_name='Experiment_1', header=0, index_col=0)
Exp2=pd.read_excel(Data_table, sheet_name='Experiment_2', header=0, index_col=0)
print(Exp1)
print(Exp2)



#################
### Make barplots.
#################


#Plot data.
def CFU_counting_exp_1(dataframe, outpath):
    
    fig, plot_av=plt.subplots(1,1,figsize=(6,3), dpi=100)
    
    #Prepare x axis, extract data.
    Conditions=dataframe.columns.tolist()
    print(Conditions)
    
    X_coords_no=[1,5,9]
    X_coords_Glc=[2,6,10]
    X_coords_IPTG=[3,7,11]
    X_coords=X_coords_no+X_coords_Glc+X_coords_IPTG
    X_coords.sort()
    print(X_coords)
    
    Mean_CFU_number=dataframe.mean(axis=0).tolist()
    Mean_CFU_number_no=[Mean_CFU_number[0], Mean_CFU_number[3], Mean_CFU_number[6]]
    Mean_CFU_number_Glc=[Mean_CFU_number[1], Mean_CFU_number[4], Mean_CFU_number[7]]
    Mean_CFU_number_IPTG=[Mean_CFU_number[2], Mean_CFU_number[5], Mean_CFU_number[8]]
    
    CFU_number=[]
    for exp in Conditions:
        CFU_number.append(dataframe.loc[:, exp].tolist())
    CFU_number_no=[CFU_number[0], CFU_number[3], CFU_number[6]]
    CFU_number_Glc=[CFU_number[1], CFU_number[4], CFU_number[7]]
    CFU_number_IPTG=[CFU_number[2], CFU_number[5], CFU_number[8]]
        
    #Plot data.
    plot_av.bar(X_coords_no, Mean_CFU_number_no, width=1, color='#b2e69a', edgecolor='k', linewidth=0.6, label='-Glc/-IPTG')
    plot_av.bar(X_coords_Glc, Mean_CFU_number_Glc, width=1, color='#f598b8', edgecolor='k', linewidth=0.6, label='Glc 0.5%/-IPTG')
    plot_av.bar(X_coords_IPTG, Mean_CFU_number_IPTG, width=1, color='#89d8fa', edgecolor='k', linewidth=0.6, label='-Glc/IPTG 1mM')
    
    plot_av.plot(X_coords_no, CFU_number_no, 'ko', markersize=1)
    plot_av.plot(X_coords_Glc, CFU_number_Glc, 'ko', markersize=1)
    plot_av.plot(X_coords_IPTG, CFU_number_IPTG, 'ko', markersize=1)
    
    plot_av.set_ylabel('CFU number, log', size=20)
    plot_av.set_xticks(X_coords_Glc, minor=False)
    plot_av.set_xticklabels(['-\nplasmid', 'pCA25\nGFP', 'pCA25\n14kDa CTD'], minor=False, rotation=0, size=12)  
      
    plot_av.set_yscale('log')
    
    #Place legend outside of a graph. Stolen from: https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    box=plot_av.get_position()
    plot_av.set_position([box.x0, box.y0, box.width * 0.95, box.height])    
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.75))
    
    plt.tight_layout(rect=[0,0,0.95,1])
    plt.show()
    plt.savefig(outpath, dpi=300, size=(6,3))
        

    return

#CFU_counting_exp_1(Exp1, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\TopoA_14kDa_CTD_expression\Experiment_1_IPTG_nc_IPTG_in_plate\Experiment_1_CFU_counting.png")



#Plot data.
def CFU_counting_exp_2(dataframe, outpath):
    
    fig, plot_av=plt.subplots(1,1,figsize=(9,3), dpi=100)
    
    #Prepare x axis, extract data.
    Conditions=dataframe.columns.tolist()
    print(Conditions)
    
    X_coords_Glc=[1,4,7,10, 13, 16]
    X_coords_IPTG=[2,5,8,11, 14, 17]
    X_coords=X_coords_Glc+X_coords_IPTG
    X_coords.sort()
    print(X_coords)
    
    Mean_CFU_number=dataframe.mean(axis=0).tolist()
    Mean_CFU_number_Glc=[Mean_CFU_number[0], Mean_CFU_number[2], Mean_CFU_number[4], Mean_CFU_number[6], Mean_CFU_number[8], Mean_CFU_number[10]]
    Mean_CFU_number_IPTG=[Mean_CFU_number[1], Mean_CFU_number[3], Mean_CFU_number[5], Mean_CFU_number[7], Mean_CFU_number[9], Mean_CFU_number[11]]
    
    CFU_number=[]
    for exp in Conditions:
        CFU_number.append(dataframe.loc[:, exp].tolist())
    CFU_number_Glc=[CFU_number[0], CFU_number[2], CFU_number[4], CFU_number[6], CFU_number[8], CFU_number[10]]
    CFU_number_IPTG=[CFU_number[1], CFU_number[3], CFU_number[5], CFU_number[7], CFU_number[9], CFU_number[11]]
        
    #Plot data.
    plot_av.bar(X_coords_Glc, Mean_CFU_number_Glc, width=1, color='#f598b8', edgecolor='k', linewidth=0.6, label='Glc 0.5%/-IPTG')
    plot_av.bar(X_coords_IPTG, Mean_CFU_number_IPTG, width=1, color='#89d8fa', edgecolor='k', linewidth=0.6, label='-Glc/IPTG 1mM')
    
    plot_av.plot(X_coords_Glc, CFU_number_Glc, 'ko', markersize=1)
    plot_av.plot(X_coords_IPTG, CFU_number_IPTG, 'ko', markersize=1)
    
    plot_av.set_ylabel('CFU number, log', size=20)
    plot_av.set_xticks(np.array(X_coords_Glc)+0.5, minor=False)
    plot_av.set_xticklabels(['- plasmid', 'pCA25 GFP', 'pSRK TopA', 'pCA25 GFP\npSRK TopA', 'pCA25 14kDa CTD', 'pCA25 14kDa CTD\npSRK TopA'], minor=False, rotation=0, size=7)  
      
    plot_av.set_yscale('log')
    
    #Place legend outside of a graph. Stolen from: https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    box=plot_av.get_position()
    plot_av.set_position([box.x0, box.y0, box.width * 0.95, box.height])    
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.75))
    
    plt.tight_layout(rect=[0,0,0.95,1])
    plt.show()
    plt.savefig(outpath, dpi=600, size=(9,3))
        

    return

CFU_counting_exp_2(Exp2, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\TopoA_14kDa_CTD_expression\Experiment_2_IPTG_in_plate\Experiment_2_CFU_counting.png")