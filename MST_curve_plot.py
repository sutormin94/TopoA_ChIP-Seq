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
MST_curves_data="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuscript\Supplementary_Tables\Supplementary_Tables.xlsx"
#Name of a worksheet.
MST_WS_name="Table S3"

#Name of a worksheet.
WS_names_ar=['Consensus F', 'Consensus R', 'Random F', 'Random R', 'Poly-T', 'Poly-A']

#Path to the output plots.
Outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts_test\MST\\"

def binding_curve(L0, A0, Kd):
    return 0.5*((A0+L0+Kd)-np.sqrt(((A0+L0+Kd)**2)-4*A0*L0))/A0


def Plot_MST_curve(data_inpath, sheetname, ws_names_ar, outpath):
    
    #Read growth curves data.
    MST_data=pd.read_excel(data_inpath, sheet_name=sheetname, skiprows=[lambda x: x in range(11)], header=11, index_col=None)   
    print(MST_data)
    
    #Create plot.
    fig, plot_1=plt.subplots(1,1,figsize=(4,3), dpi=100)
    #Prepare array of colors.
    Mfc_colors={'Consensus F' : "#E5001F", 'Poly-T' : "#040073", 'Random F' : "#A25AE6", 'Consensus R' : "none", 'Poly-A' : "none", 'Random R' : "none"}
    Mec_colors={'Consensus F' : "#E5001F", 'Poly-T' : "#040073", 'Random F' : "#A25AE6", 'Consensus R' : "#E5001F", 'Poly-A' : "#040073", 'Random R' : "#A25AE6"}
    #Binding curves parameters.
    Kd_dict={'Consensus F' : 2.1E-08, 'Poly-T' : 6.2E-08, 'Random F' : 3.21E-07, 'Consensus R' : 2.7E-08, 'Poly-A' : 4.0E-08, 'Random R' : 3.28E-07}
    Kd_std={'Consensus F' : 0.7E-08, 'Poly-T' : 2.2E-08, 'Random F' : 1.1E-08, 'Consensus R' : 1.5E-09, 'Poly-A' : 2.6E-09, 'Random R' : 3.8E-08}
    A0_dict={'Consensus F' : 0.5E-09, 'Poly-T' : 0.5E-09, 'Random F' : 0.5E-09, 'Consensus R' : 0.5E-09, 'Poly-A' : 0.5E-09, 'Random R' : 0.5E-09}
    Unbound={'Consensus F' : 39.455, 'Poly-T' : 40.99, 'Random F' : 36.146, 'Consensus R' : 836.52, 'Poly-A' : 836.78, 'Random R' : 882.86}
    Bound={'Consensus F' : 21.045, 'Poly-T' : 23.291, 'Random F' : 20.425, 'Consensus R' : 931.37, 'Poly-A' : 946.91, 'Random R' : 941.45}
    
    #Iterate datasheets.
    for dataset_name in ws_names_ar:
    
        Sat=Unbound[dataset_name]-Bound[dataset_name]
        
        #Add 95% two-tail confident interval range for mean.
        MST_data[dataset_name + ' rel_Mean']=(MST_data[dataset_name + ' MST signal, Mean']-Unbound[dataset_name])*(-1)/Sat
        MST_data[dataset_name + ' CI95']=((MST_data[dataset_name + ' MST signal, STD']*1.96)/np.sqrt(MST_data[dataset_name + ' Number of replicates']))/Sat
        MST_data[dataset_name + ' rel_Mean_up']=MST_data[dataset_name + ' rel_Mean']+MST_data[dataset_name + ' CI95']
        MST_data[dataset_name + ' rel_Mean_dn']=MST_data[dataset_name + ' rel_Mean']-MST_data[dataset_name + ' CI95']
        
        conc_range=np.linspace(2E-11,10E-6, 10000)
    
        #Plot data.
        plot_1.plot(MST_data[dataset_name + ' EcTopoI concentration, M'], MST_data[dataset_name + ' rel_Mean'], marker='o', mfc=Mfc_colors[dataset_name], mec=Mec_colors[dataset_name], linewidth=0, label=f'{dataset_name} Kd={int(Kd_dict[dataset_name]*1E9)}$\pm${int(Kd_std[dataset_name]*1E9)} nM', markersize=2)
        plot_1.plot(conc_range, binding_curve(conc_range, A0_dict[dataset_name], Kd_dict[dataset_name]), '--', color=Mec_colors[dataset_name], linewidth=0.7)
        plot_1.errorbar(MST_data[dataset_name + ' EcTopoI concentration, M'], MST_data[dataset_name + ' rel_Mean'], yerr=MST_data[dataset_name + ' CI95'], marker='o', color=Mfc_colors[dataset_name], ecolor=Mec_colors[dataset_name], elinewidth=0.5, capsize=2, markersize=2, linewidth=0)
    
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
    plt.savefig(outpath+'MST_binding_curves.png', dpi=300, figsize=(4,3))
    #plt.close()
    
    return

Plot_MST_curve(MST_curves_data, MST_WS_name, WS_names_ar, Outpath)