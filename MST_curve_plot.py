###############################################
##Dmitry Sutormin, 2020##
##Growth curves plot##

#Plots MST binding curves visualization.
###############################################

#######
#Packages to be imported.
#######

import random as rd
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit
import pandas as pd



#######
#Import data.
#######

#Path to the raw data.
MST_curves_data="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Binding_motif_validation\MST\MST_raw_data.xlsx"

#Name of a worksheet.
WS_names_dict={'Consensus 2' : "Consensus_short_raw_data",
               'Consensus 1' : "Consensus_long_raw_data",
               'Poly-T' : "PolyT_raw_data",
               'Random' : "Random_raw_data",}

#Path to the output plots.
Outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Binding_motif_validation\MST\\"

def binding_curve(L0, A0, Kd):
    return 0.5*((A0+L0+Kd)-np.sqrt(((A0+L0+Kd)**2)-4*A0*L0))/A0


def Plot_MST_curve(data_inpath, sheetnames_dict, outpath):
    
    #Create plot.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    #Prepare array of colors.
    Colors={'Consensus 1' : "#1D9659", 'Consensus 2' : "#E5001F", 'Poly-T' : "#040073", 'Random' : "#A25AE6"}
    #Binding constants.
    #Kd_dict={'Consensus 1' : 7.1246E-08, 'Consensus 2' : 2.107E-08, 'Poly-T' : 6.1697E-08, 'Random' : 3.2083E-07}
    Kd_dict={'Consensus 1' : 6.6E-08, 'Consensus 2' : 2.1E-08, 'Poly-T' : 5.3E-08, 'Random' : 3.06E-07}
    Kd_dict1={'Consensus 1' : 6.6E-08, 'Consensus 2' : 2.1E-08, 'Poly-T' : 5.3E-08, 'Random' : 2.06E-07}
    Kd_std={'Consensus 1' : 2.4E-08, 'Consensus 2' : 0.7E-08, 'Poly-T' : 2.6E-08, 'Random' : 1.1E-08}
    A0_dict={'Consensus 1' : 1.25E-09, 'Consensus 2' : 0.5E-09, 'Poly-T' : 0.5E-09, 'Random' : 0.5E-09}
    
    Unbound={'Consensus 1' : 36.021, 'Consensus 2' : 39.455, 'Poly-T' : 40.99, 'Random' : 36.146}
    Bound={'Consensus 1' : 21.939, 'Consensus 2' : 21.045, 'Poly-T' : 23.291, 'Random' : 20.425}
    
    #Iterate datasheets.
    for name, sheetname in sheetnames_dict.items():
    
        #Read growth curves data.
        gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=None)
        
        Sat=Unbound[name]-Bound[name]
        
        #Add 95% two-tail confident interval range for mean.
        gc_data['Rel_Mean']=(gc_data['Mean']-Unbound[name])*(-1)/Sat
        gc_data['CI95']=((gc_data['STD']*1.96)/np.sqrt(gc_data['Replics']))/Sat
        gc_data['Rel_Mean_up']=gc_data['Rel_Mean']+gc_data['CI95']
        gc_data['Rel_Mean_dn']=gc_data['Rel_Mean']-gc_data['CI95']
        
        conc_range=np.linspace(2E-11,10E-6, 10000)
    
        #Plot data.
        plot_1.plot(gc_data.Concentration, gc_data.Rel_Mean, marker='o', color=Colors[name], linewidth=0, label=name, markersize=2)
        plot_1.plot(conc_range, binding_curve(conc_range, A0_dict[name], Kd_dict[name]), '--', color=Colors[name], linewidth=1, label=f'fit: Kd={int(Kd_dict1[name]*1E9)}$\pm${int(Kd_std[name]*1E9)} nM')
        #plot_1.plot(gc_data.Concentration, gc_data.Rel_Mean_up, '-', color=Colors[name], linewidth=1, alpha=0.3)
        #plot_1.plot(gc_data.Concentration, gc_data.Rel_Mean_dn, '-', color=Colors[name], linewidth=1, alpha=0.3)
        plot_1.fill_between(gc_data.Concentration, gc_data.Rel_Mean_dn, gc_data.Rel_Mean_up, color=Colors[name], alpha=0.3, linewidth=0)
    
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1.5)
    plot_1.spines["left"].set_linewidth(1.5)
    plt.xlim([1E-11,10E-6])
    plt.xscale('log')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('EcTopoI concentration, mol/L', size=16)
    plt.ylabel('Fraction bound', size=16)
    plt.legend(fontsize=8.5, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, labelspacing=0.35, loc="lower right", bbox_to_anchor=(0.53,0.35))
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath+'MST_binding_curves_cor_Kd.png', dpi=300, figsize=(3, 2))
    plt.close()
    return

def Plot_MST_curve_1(data_inpath, sheetnames_dict, outpath):
    
    #Create plot.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    #Prepare array of colors.
    Colors={'Consensus 1' : "#1D9659", 'Consensus 2' : "#E5001F", 'Poly-T' : "#040073", 'Random' : "#A25AE6"}
    #Binding constants.
    #Kd_dict={'Consensus 1' : 7.1246E-08, 'Consensus 2' : 2.107E-08, 'Poly-T' : 6.1697E-08, 'Random' : 3.2083E-07}
    Kd_dict={'Consensus 1' : 6.6E-08, 'Consensus 2' : 2.1E-08, 'Poly-T' : 5.3E-08, 'Random' : 3.06E-07}
    Kd_dict1={'Consensus 1' : 6.6E-08, 'Consensus 2' : 2.1E-08, 'Poly-T' : 5.3E-08, 'Random' : 2.06E-07}
    Kd_std={'Consensus 1' : 2.4E-08, 'Consensus 2' : 0.7E-08, 'Poly-T' : 2.6E-08, 'Random' : 1.1E-08}
    A0_dict={'Consensus 1' : 1.25E-09, 'Consensus 2' : 0.5E-09, 'Poly-T' : 0.5E-09, 'Random' : 0.5E-09}
    
    Unbound={'Consensus 1' : 36.021, 'Consensus 2' : 39.455, 'Poly-T' : 40.99, 'Random' : 36.146}
    Bound={'Consensus 1' : 21.939, 'Consensus 2' : 21.045, 'Poly-T' : 23.291, 'Random' : 20.425}
    
    #Iterate datasheets.
    for name, sheetname in sheetnames_dict.items():
    
        #Read growth curves data.
        gc_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=None)
        
        Sat=Unbound[name]-Bound[name]
        
        #Add 95% two-tail confident interval range for mean.
        gc_data['Rel_Mean']=(gc_data['Mean']-Unbound[name])*(-1)/Sat
        gc_data['CI95']=((gc_data['STD']*1.96)/np.sqrt(gc_data['Replics']))/Sat
        gc_data['Rel_Mean_up']=gc_data['Rel_Mean']+gc_data['CI95']
        gc_data['Rel_Mean_dn']=gc_data['Rel_Mean']-gc_data['CI95']
        
        conc_range=np.linspace(2E-11,10E-6, 10000)
    
        #Plot data.
        plot_1.plot(gc_data.Concentration, gc_data.Rel_Mean, marker='o', color=Colors[name], linewidth=0, label=f'{name} Kd={int(Kd_dict1[name]*1E9)}$\pm${int(Kd_std[name]*1E9)} nM', markersize=2)
        plot_1.plot(conc_range, binding_curve(conc_range, A0_dict[name], Kd_dict[name]), '--', color=Colors[name], linewidth=1)
        plot_1.errorbar(gc_data.Concentration, gc_data.Rel_Mean, yerr=gc_data['CI95'], marker='o', color=Colors[name], ecolor=Colors[name], elinewidth=0.5, capsize=2, markersize=3, linewidth=0)
    
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)  
    plot_1.spines["bottom"].set_linewidth(1.5)
    plot_1.spines["left"].set_linewidth(1.5)
    plt.xlim([1E-11,10E-6])
    plt.xscale('log')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('EcTopoI concentration, mol/L', size=16)
    plt.ylabel('Fraction bound', size=16)
    plt.legend(fontsize=8.5, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, labelspacing=0.35, loc="lower right", bbox_to_anchor=(0.70,0.63))
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath+'MST_binding_curves_cor_Kd_1.png', dpi=300, figsize=(3, 2))
    plt.close()
    return


Plot_MST_curve(MST_curves_data, WS_names_dict, Outpath)
Plot_MST_curve_1(MST_curves_data, WS_names_dict, Outpath)