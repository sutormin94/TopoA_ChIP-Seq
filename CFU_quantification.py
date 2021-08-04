###############################################
##Dmitry Sutormin, 2021##
##CFU & Dot-blot quantification & Pull-down quantification##

#Plots CFU number with barplots.
#Plots Dot-blot signal barplots.
#Plots pull-down signal barplots.
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


#################
### BW25113 topA mutants CFU counting.
#################


#Path to CFU data table.
data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\BW25113_topA_mutants\CFU_counting\BW25113_CFU_counting.xlsx"
data_tab=pd.read_excel(data_table, sheet_name='Sheet1', header=0, index_col=0)
print(data_tab)

def compute_standard_error(column):
    std_err=(np.std(column)/np.sqrt(len(column)))*1.96
    return std_err

#Plot data.
def count_CFU_BW25113(dataframe, outpath):
    
    ###
    ##Plot all CFU.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(2.2,3.5), dpi=100)
    
    #Prepare x axis.
    Conditions=[r'$wt$', r'$topA\Delta11$', r'$topA\Delta14$', r'$topA\Delta30$']
    #print(len(Conditions))
    
    X_coords=[1,2,3,4]
    #print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar_1=dataframe.loc[:, 'BW25113_wt']
    Data_points_ar_2=dataframe.loc[:, 'BW25113_delta11']
    Data_points_ar_3=dataframe.loc[:, 'BW25113_delta14']
    Data_points_ar_4=dataframe.loc[:, 'BW25113_delta30']
    #print(Data_points_ar_1, Data_points_ar_2, Data_points_ar_3, Data_points_ar_4)
    
    #Compare means, stat.
    Data_points=[Data_points_ar_1, Data_points_ar_2, Data_points_ar_3, Data_points_ar_4]
    #print(Data_points)
    for i in range(len(Data_points)-1):
        print(stats.ttest_ind(Data_points[0].dropna(), Data_points[i+1].dropna()))
    
    #Prepare mean CFU.
    Mean_CFU=dataframe.mean(axis=0)
    print(Mean_CFU)
    
    #Prepare data for error bars.
    Std_err_CFU=dataframe.apply(compute_standard_error, axis=0)
    #print(Std_err_CFU)
    
    #Set colors for bars.
    Colors=['#ff4564', '#ff99af', '#4589ff', '#99c9ff',]
    #print(Colors)
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i].tolist())
        Data_points_ar=Data_points_ar+Data_points[i].tolist()
    
    #print(len(X_coords_for_points))
    #print(len(Data_points_ar))
    
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_CFU, yerr=Std_err_CFU, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.7, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('CFU number', size=17)
    plot_av.set_xticks(X_coords)
    plot_av.set_xticklabels(Conditions, rotation=90, size=10)  
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_yscale('log')
    plot_av.set_ylim([1e8, 7e9])
    
    #plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), (r'wt', r'topA$\Delta$11', r'topA$\Delta$14', r'topA$\Delta$30'), fontsize=14, ncol=3, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, loc='upper left')
    plt.tight_layout()
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(6.5,3))

    return

#count_CFU_BW25113(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\BW25113_topA_mutants\CFU_counting\BW25113_topA_mutants_CFU_counting.png")



#Path to CFU data table.
#data_table_CTD_st="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\TopoA_14kDa_CTD_expression_CFU_count\CFU_counting_data.xlsx"
#data_tab_CTD_st=pd.read_excel(data_table_CTD_st, sheet_name='Experiment_3', header=0, index_col=0)
#print(data_tab_CTD_st)

#Plot data.
def count_CFU_CTD_short_term(dataframe, outpath):
    
    ###
    ##Plot all CFU.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(2,3.5), dpi=100)
    
    #Prepare x axis.
    Conditions=[ 'pCA25\n14kDa CTD']
    #print(len(Conditions))
    
    X_coords=[1,2]
    print(f'Len of X_coords: {len(X_coords)}')
    
    X_coords_major=[1.5]
    print(f'Len of X_coords_major: {len(X_coords_major)}')
    
    #Prepare data for bars.
    Data_points_ar_1=dataframe.loc[:, 'CTD']
    Data_points_ar_2=dataframe.loc[:, 'CTD IPTG']
    #print(Data_points_ar_1, Data_points_ar_2, Data_points_ar_3, Data_points_ar_4)
    
    #Compare means, stat.
    Data_points=[Data_points_ar_1, Data_points_ar_2]
    #print(Data_points)
    for i in range(len(Data_points)-1):
        print(stats.ttest_ind(Data_points[0].dropna(), Data_points[i+1].dropna()))
    
    #Prepare mean CFU.
    Mean_CFU=dataframe.mean(axis=0)
    print(f'Len of Mean_CFU: {len(Mean_CFU)}')
    
    #Prepare data for error bars.
    Std_err_CFU=dataframe.apply(compute_standard_error, axis=0)
    print(f'Len of Std_err_CFU: {len(Std_err_CFU)}')
    
    #Set colors for bars.
    Colors=['#ff4564', '#ff99af']
    print(f'Len of Colors: {len(Colors)}')
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i].tolist())
        Data_points_ar=Data_points_ar+Data_points[i].tolist()
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_CFU, yerr=Std_err_CFU, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.7, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('CFU number', size=17)
    plot_av.set_xticks(X_coords_major)
    plot_av.set_xticklabels(Conditions, rotation=0, size=10)  
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_yscale('log')
    plot_av.set_ylim([1e7, 5e9])
    
    plt.legend((Bars[0],Bars[1]), (r'-IPTG', r'+IPTG 1mM'), fontsize=9.5, ncol=1, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, loc='upper left')
    plt.tight_layout()
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(2,3.5))

    return

#count_CFU_CTD_short_term(data_tab_CTD_st, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Properties_of_CTD\TopoA_14kDa_CTD_expression_CFU_count\Experiment_3_short_term_oe_1h_plate_no_IPTG\CTD_effect_short_term_1h.png")



#################
### S9.6 dot-blot quantification; R-loops quantification. Experiment with EcTopoI Y319F overexpression.
#################

#Path to Dot-blot data table.
#data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\R-loop_detection\Results\Dot_blot_S9_6\Dot_blot_quantification.xlsx"
#data_tab=pd.read_excel(data_table, sheet_name='Blot_07_09_20_plotting', header=0, index_col=0)
#print(data_tab)

#Plot data.
def plot_dot_blot(dataframe, outpath):
    
    ###
    ##Plot dot-blot.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(2.8,3.5), dpi=100)
    
    #Prepare x axis.
    Conditions=['EcTopoI\nY319F\n-HI', 'EcTopoI\nY319F\n+HI', 'GFP\n-HI', 'GFP\n+HI']
    #print(len(Conditions))
    
    X_coords=[1,1.5, 2.5,3, 4,4.5, 5.5,6]
    X_coords_major=[1.25, 2.75, 4.25, 5.75]
    print(X_coords)
    
    #Prepare data for bars.
    Data_points=[dataframe.loc['Y319F-IPTG', '+III/-HI'], dataframe.loc['Y319F+IPTG', '+III/-HI'],
                 dataframe.loc['GFP-IPTG', '+III/-HI'], dataframe.loc['GFP+IPTG', '+III/-HI'],
                 dataframe.loc['Y319F-IPTG', '+III/+HI'], dataframe.loc['Y319F+IPTG', '+III/+HI'],
                 dataframe.loc['GFP-IPTG', '+III/+HI'], dataframe.loc['GFP+IPTG', '+III/+HI']]
    print(Data_points)
    
    #Set colors for bars.
    Colors=['#99c9ff', '#ff99af', '#99c9ff', '#ff99af','#99c9ff', '#ff99af', '#99c9ff', '#ff99af']
    print(Colors)
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Data_points, align='center', width=0.5, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.set_ylabel('Pixel intensity', size=17)
    plot_av.set_xticks(X_coords_major)
    plot_av.set_xticklabels(Conditions, rotation=0, size=8)  
    plot_av.tick_params(axis='x', which='major', pad=5)
    
    plt.legend((Bars[0],Bars[1]), ('-IPTG', '+1 mM IPTG'), fontsize=12, ncol=1, frameon=False, markerscale=1, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, loc='upper right')
    plt.tight_layout()
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(6.5,3))

    return

#plot_dot_blot(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\R-loop_detection\Results\Dot_blot_S9_6\Dot_blot_07_09_20_quantification.png")



#################
### Pull-down experiments quantification. EcTopoI pull-down experiments. 
#################

#Path to pull-down data table.
data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Pull_down_experiments\Travin_EcTopoI_coIP_RNAP_detection\Pull_down_experiments_qantification.xlsx"
data_tab=pd.read_excel(data_table, sheet_name='Final_data_background_corrected', header=0, index_col=0)
print(data_tab)

#Plot data.
def plot_pull_down(dataframe, outpath):
    
    ###
    ##Plot pull-down.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(2,2.5), dpi=100)
    
    #Prepare x axis.
    Conditions=['pCA25\nGFP', 'pCA25\n14kDa CTD']
    #print(len(Conditions))
    
    X_coords=[1,1.5, 2.5,3]
    X_coords_major=[1.25, 2.75]
    print(X_coords)
    
    #Prepare data for bars.
    Data_points=[dataframe.loc['GFP -IPTG', 'RNAP/EcTopoI'],   dataframe.loc['GFP +IPTG', 'RNAP/EcTopoI'],
                 dataframe.loc['CTD -IPTG', 'RNAP/EcTopoI'],   dataframe.loc['CTD +IPTG', 'RNAP/EcTopoI']]
    print(Data_points)
    
    #Set colors for bars.
    Colors=['#99c9ff', '#ff99af', '#99c9ff', '#ff99af']
    print(Colors)
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Data_points, align='center', width=0.45, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.set_ylabel('RNAP/EcTopoI\npixel intensity', size=17)
    plot_av.set_xticks(X_coords_major)
    plot_av.set_xticklabels(Conditions, rotation=0, size=8)  
    plot_av.tick_params(axis='x', which='major', pad=5)
    
    plt.legend((Bars[0],Bars[1]), ('-IPTG', '+1 mM IPTG'), fontsize=12, ncol=1, frameon=False, markerscale=1, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, loc='upper right')
    plt.tight_layout()
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(2,2.5))

    return

plot_pull_down(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Pull_down_experiments\Travin_EcTopoI_coIP_RNAP_detection\EcTopoI_pull_down_quantification_background_corr.svg")

