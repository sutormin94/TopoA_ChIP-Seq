###############################################
##Dmitry Sutormin, 2020##
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
Cell_length_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\BW25113_topA_mutants\Microscopy\Cell_length_distribution.xlsx"

#Path to output plot.
Plot_outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\BW25113_topA_mutants\Microscopy\Cells_length_measurement_distribution_topA_mutants.png"


#######
#Fit data.
#######


def fitted_data(data):
    avg=np.mean(data)
    var=np.var(data)
    
    pdf_x=np.linspace(0,10,100)
    pdf_y=1.0/np.sqrt(2*np.pi*var)*np.exp(-0.5*(pdf_x-avg)**2/var)
    return pdf_x, pdf_y


#######
#Read data and plot.
#######

def plot_cell_length_dist(data_table_path, out_plot_path):
    
    #Import data.
    GFP_minus=pd.read_excel(data_table_path, sheet_name='topA_wt')['Length']
    GFP_plus=pd.read_excel(data_table_path, sheet_name='topA_delta11')['Length']
    CTD_minus=pd.read_excel(data_table_path, sheet_name='topA_delta14')['Length']
    CTD_plus=pd.read_excel(data_table_path, sheet_name='topA_delta30')['Length']
    
    #Prepare data.
    Data_dict={'topA_wt' : [GFP_minus, "#f5d81a", r'$wt$'], 'topA_delta11' : [GFP_plus, "#5bff7c", r'$topA\Delta11$'],
               'topA_delta14' : [CTD_minus, "#878787", r'$topA\Delta14$'], 'topA_delta30' : [CTD_plus, "#212DA7", r'$topA\Delta30$']}
    
    #Compare the data statistically.
    #Test for normality.
    p_value_thr=1e-3
    for dataset_name, dataset_data in Data_dict.items():
        statistics, p_value=normaltest(dataset_data[0])
        if p_value<p_value_thr:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean length={round(np.mean(dataset_data[0]),2)}um: {statistics}, {p_value} - not normal')
        else:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean length={round(np.mean(dataset_data[0]),2)}um: {statistics}, {p_value} - normal')
            
    #Compare means. Use non-parametric Mann-Whitney U test. Null hypothesis: sub-populations came from the same general population - means are the same.
    for dataset_name_1, dataset_data_1 in Data_dict.items():
        for dataset_name_2, dataset_data_2 in Data_dict.items():
            if dataset_name_1!=dataset_name_2:
                statistics, p_value=mannwhitneyu(dataset_data_1[0], dataset_data_2[0], alternative='two-sided')
                if p_value<p_value_thr:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)}um): {statistics}, {p_value} - means are different')
                else:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)}um): {statistics}, {p_value} - means are equal')
    
    #Calculate fraction of long cells.    
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

plot_cell_length_dist(Cell_length_path, Plot_outpath)