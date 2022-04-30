###############################################
##Dmitry Sutormin, 2021##
##Microscopy data analysis##

#Analysis of cell length distribution.
###############################################

#######
#Packages to be imported.
#######

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns 
from scipy.stats import norm, shapiro, normaltest, mannwhitneyu


#######
#Import data.
#######

#Path to raw data table.
Cell_length_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables\Supplementary_Tables.xlsx"
#Name of a worksheet.
WS_name_cell_length="Table S13"


#######
#Fit data.
#######

def fitted_data(data):
    avg=np.mean(data)
    var=np.var(data)
    
    pdf_x=np.linspace(0,10,100)
    pdf_y=1.0/np.sqrt(2*np.pi*var)*np.exp(-0.5*(pdf_x-avg)**2/var)
    return pdf_x, pdf_y


##########
## Read cell length data and plot for E. coli BW25113 topA mutants.
#########

def plot_cell_length_dist_BW25113_topA_mut(data_inpath, sheetname, out_plot_path):
    
    cl_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Import data.
    BW25113_wt=cl_data.loc[:, 'BW25113 wt length, mkm'].dropna()
    BW25113_topA_delta11=cl_data.loc[:, 'BW25113 topA delta11 length, mkm'].dropna()
    BW25113_topA_delta14=cl_data.loc[:, 'BW25113 topA delta14 length, mkm'].dropna()
    BW25113_topA_delta30=cl_data.loc[:, 'BW25113 topA delta30 length, mkm'].dropna()
    
    #Prepare data.
    Data_dict={'topA_wt' : [BW25113_wt, "#f5d81a", r'$wt$'], 'topA_delta11' : [BW25113_topA_delta11, "#5bff7c", r'$topA\Delta11$'],
               'topA_delta14' : [BW25113_topA_delta14, "#878787", r'$topA\Delta14$'], 'topA_delta30' : [BW25113_topA_delta30, "#212DA7", r'$topA\Delta30$']}
    
    #Compare the data statistically.
    #Test for normality.
    print('\nTest if data is normally distributed using normaltest')
    p_value_thr=1e-3
    for dataset_name, dataset_data in Data_dict.items():
        statistics, p_value=normaltest(dataset_data[0])
        if p_value<p_value_thr:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean length={round(np.mean(dataset_data[0]),2)}um: {statistics}, {p_value} - not normal')
        else:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean length={round(np.mean(dataset_data[0]),2)}um: {statistics}, {p_value} - normal')
            
    #Compare means. Use non-parametric Mann-Whitney U test. Null hypothesis: sub-populations came from the same general population - means are the same.
    print('\nCompare datasets means by mannwhitneyu')
    for dataset_name_1, dataset_data_1 in Data_dict.items():
        for dataset_name_2, dataset_data_2 in Data_dict.items():
            if dataset_name_1!=dataset_name_2:
                statistics, p_value=mannwhitneyu(dataset_data_1[0], dataset_data_2[0], alternative='two-sided')
                if p_value<p_value_thr:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)}um): {statistics}, {p_value} - means are different')
                else:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)}um): {statistics}, {p_value} - means are equal')
    
    #Calculate fraction of long cells.    
    print('\nDataset name\tNumber of long cells\tTotal number of cells\tFraction of long cells')
    wt_mean=np.mean(Data_dict['topA_wt'][0])
    for dataset_name, dataset_data in Data_dict.items():
        Number_of_long_cells=dataset_data[0][dataset_data[0]>2*wt_mean].count()
        Total_number_of_cells=len(dataset_data[0])
        Fraction_of_long_cells=float(Number_of_long_cells)/Total_number_of_cells
        print(dataset_name, Number_of_long_cells, Total_number_of_cells, Fraction_of_long_cells)
    
    #Plot data.
    xticks_ar=[1,2,4,6,8,10.1]
    xtickslabels_ar=[1,2,4,6,8,'>10']
    bins_ar=list(np.linspace(0,10.2,50))
    fig, plot_1=plt.subplots(1,1,figsize=(4,2.5), dpi=100)   
    
    for dataset_name, dataset_data in Data_dict.items():
        clipped_array=np.clip(dataset_data[0], bins_ar[0], bins_ar[-1])
        weights=np.ones_like(clipped_array)/(len(clipped_array)) #Taken from https://stackoverflow.com/questions/42481698/probability-density-histogram-with-matplotlib-doesnt-make-sense     
        plot_1.hist(clipped_array, bins=bins_ar, weights=weights, color=dataset_data[1], alpha=0.3, label=f'{dataset_data[2]} ({len(dataset_data[0])})')
        plot_1.hist(clipped_array, bins=bins_ar, weights=weights, histtype=u'step', color=dataset_data[1])
       
    plot_1.axvline(2*wt_mean, color='black', linestyle='--', linewidth=1)
    plot_1.set_xticks(xticks_ar)
    plot_1.set_xticklabels(xtickslabels_ar)
    plot_1.set_xlabel(r'Cell length, $\mu$m', size=12)
    plot_1.set_ylabel('Fraction', size=12)
    plot_1.set_xlim([0.5,10.5])
    
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)
    
    plt.legend(fontsize=8, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3)
    plt.tight_layout()
    plt.show()
    
    plt.savefig(out_plot_path, dpi=300)
    
    return

#Path to output plot.
Plot_outpath_BW25113_topA_mut="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Microscopy_cell_length\Cells_length_measurement_distribution_BW25113_topA_mutants.png"

plot_cell_length_dist_BW25113_topA_mut(Cell_length_path, WS_name_cell_length, Plot_outpath_BW25113_topA_mut)



##########
## Read cell length data and plot for EcTopoI 14kDa CTD anf GFP overexpression from pCA24 plasmid in E. coli DY330.
#########

def plot_cell_length_dist_GFP_CTD_oe(data_inpath, sheetname, out_plot_path):
    
    cl_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Import data.
    GFP_minus=cl_data.loc[:, 'GFP- length, mkm'].dropna()
    GFP_plus=cl_data.loc[:, 'GFP+ length, mkm'].dropna()
    CTD_minus=cl_data.loc[:, 'CTD- length, mkm'].dropna()
    CTD_plus=cl_data.loc[:, 'CTD+ length, mkm'].dropna() 
    
    #Prepare data.
    Data_dict={'GFP_minus' : [GFP_minus, "#f5d81a", 'pCA24 GFP -IPTG'], 'GFP_plus' : [GFP_plus, "#5bff7c", 'pCA24 GFP +IPTG 1mM'],
               'CTD_minus' : [CTD_minus, "#878787", 'pCA24 CTD\n-IPTG'], 'CTD_plus' : [CTD_plus, "#212DA7", 'pCA24 CTD\n+IPTG 1mM']}
    
    #Compare the data statistically.
    #Test for normality.
    print('\nTest if data is normally distributed using normaltest')
    p_value_thr=1e-3
    for dataset_name, dataset_data in Data_dict.items():
        statistics, p_value=normaltest(dataset_data[0])
        if p_value<p_value_thr:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean length={round(np.mean(dataset_data[0]),2)}um: {statistics}, {p_value} - not normal')
        else:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean length={round(np.mean(dataset_data[0]),2)}um: {statistics}, {p_value} - normal')
            
    #Compare means. Use non-parametric Mann-Whitney U test. Null hypothesis: sub-populations came from the same general population - means are the same.
    print('\nCompare datasets means by mannwhitneyu')
    for dataset_name_1, dataset_data_1 in Data_dict.items():
        for dataset_name_2, dataset_data_2 in Data_dict.items():
            if dataset_name_1!=dataset_name_2:
                statistics, p_value=mannwhitneyu(dataset_data_1[0], dataset_data_2[0], alternative='two-sided')
                if p_value<p_value_thr:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)}um): {statistics}, {p_value} - means are different')
                else:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)}um): {statistics}, {p_value} - means are equal')
    
    #Calculate fraction of long cells. 
    print('\nDataset name\tNumber of long cells\tTotal number of cells\tFraction of long cells')
    wt_mean=np.mean(Data_dict['GFP_minus'][0])
    for dataset_name, dataset_data in Data_dict.items():
        Number_of_long_cells=dataset_data[0][dataset_data[0]>2*wt_mean].count()
        Total_number_of_cells=len(dataset_data[0])
        Fraction_of_long_cells=float(Number_of_long_cells)/Total_number_of_cells
        print(dataset_name, Number_of_long_cells, Total_number_of_cells, Fraction_of_long_cells)
    
    #Plot data.
    xticks_ar=[1,2,4,6,8,10.1]
    xtickslabels_ar=[1,2,4,6,8,'>10']
    bins_ar=list(np.linspace(0,10.2,50))
    fig, plot_1=plt.subplots(1,1,figsize=(4,2.5), dpi=100)   
    
    for dataset_name, dataset_data in Data_dict.items():
        clipped_array=np.clip(dataset_data[0], bins_ar[0], bins_ar[-1])
        weights=np.ones_like(clipped_array)/(len(clipped_array)) #Taken from https://stackoverflow.com/questions/42481698/probability-density-histogram-with-matplotlib-doesnt-make-sense     
        plot_1.hist(clipped_array, bins=bins_ar, weights=weights, color=dataset_data[1], alpha=0.3, label=f'{dataset_data[2]} ({len(dataset_data[0])})')
        plot_1.hist(clipped_array, bins=bins_ar, weights=weights, histtype=u'step', color=dataset_data[1])
       
    plot_1.set_xticks(xticks_ar)
    plot_1.set_xticklabels(xtickslabels_ar)
    plot_1.set_xlabel(r'Cell length, $\mu$m', size=12)
    plot_1.set_ylabel('Fraction', size=12)
    plot_1.set_xlim([0.5,10.5])
    
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)
    
    plt.legend(fontsize=7.5, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc='upper right')
    plt.tight_layout()
    plt.show()
    
    plt.savefig(out_plot_path, dpi=300)
    
    return

#Path to output plot.
Plot_outpath_GFP_CTD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\Microscopy_cell_length\Cells_length_measurement_distribution_DY330_GFP_CTD_oe.png"

plot_cell_length_dist_GFP_CTD_oe(Cell_length_path, WS_name_cell_length, Plot_outpath_GFP_CTD)