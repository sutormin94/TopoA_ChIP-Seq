###############################################
##Dmitry Sutormin, 2020##
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
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import binom



#################
### Primers calibration data analysis.
#################


#Path to qPCR data table.
#topA_pc_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\EcTopoI_RNAP_uncoupling_effects_qPCR\\rnhAB_expression_analysis_qPCR\gyrAB_parCE_rnhAB_topB_calibration_16-22_09_20\gyrAB_parCE_rnhAB_topB_calibration_16-22_09_2020.xlsx"
#rpmH_pc_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\qPCR_results\CTD_effect_qPCR_primers_calibration\RpmH_primers_calibration_18_03_20_admin_2020-03-18 21-12-35_CT018440 -  Quantification Plate View Results.xlsx"
#topA_pc=pd.read_excel(topA_pc_table, sheet_name='Processed_data', header=0, index_col=0)
#rpmH_pc=pd.read_excel(rpmH_pc_table, sheet_name='Processed_data_rpmH', header=0, index_col=0)
#print(topA_pc)
#print(rpmH_pc)


#Plot data.
def qPCR_primers_calibration(dataframe, suptitle_text, outpath):
    fig, plot_av=plt.subplots(3,3,figsize=(12,6), dpi=100)
    fig.suptitle(suptitle_text, size=15)
    
    #Prepare x axis.
    Conc_data=dataframe.loc['Concentration', :].tolist()
    
    #Plot data.
    Num_of_datasets=len(dataframe.index.tolist())-1
    print(Num_of_datasets)
    Primers_list=dataframe.index.tolist()[1:]
    print(Primers_list)
    for i in range(Num_of_datasets):
        #Points.
        Primers_pair=Primers_list[i]
        print(Primers_pair)
        qC_data=dataframe.loc[Primers_pair, :]
        print(Primers_pair, qC_data)
        print(i, int(i/3), i%3)
        plot_av[int(i/3), i%3].scatter(Conc_data, qC_data.tolist(), s=2, color='k', edgecolors='black', linewidth=0.2, alpha=1, zorder=1000) 
        print(Conc_data, qC_data.tolist())
        
        #Filter data, remove nan.
        qC_data_ar=qC_data.tolist()
        Conc_data_filtered=[]
        qC_data_ar_filtered=[]
        for j in range(len(Conc_data)):
            if np.isnan(qC_data_ar[j])==False:
                Conc_data_filtered.append(Conc_data[j])
                qC_data_ar_filtered.append(qC_data_ar[j])
                
        print(Conc_data_filtered, qC_data_ar_filtered)
        #Linear fitting of calibration data.
        fit=np.polyfit(Conc_data_filtered, qC_data_ar_filtered, 1)
        print(fit)
        fit_fn=np.poly1d(fit) 
        plot_av[int(i/3), i%3].plot(Conc_data, fit_fn(Conc_data), '--b', linewidth=0.5, label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3))) 
        primer_effectiveness=10**(-1/fit[0])
        plot_av[int(i/3), i%3].annotate(f"$\lambda$={round(primer_effectiveness, 2)}", xy=(0.7, 0.5), xycoords='axes fraction', size=16)

        plot_av[int(i/3), i%3].set_ylabel('qC', size=20)
        plot_av[int(i/3), i%3].set_xticks(np.unique(Conc_data), minor=False)
        plot_av[int(i/3), i%3].set_xticklabels([1, 10, 100, 1000], rotation=0, size=12)
        #plot_av[int(i/3), i%3].set_xscale('log')
        plot_av[int(i/3), i%3].set_title(f"Primers pair {Primers_pair}", size=12)
        plot_av[int(i/3), i%3].legend()
        
    fig.delaxes(plot_av[2,2])
    #fig.delaxes(plot_av[2,1])
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    plt.savefig(outpath, dpi=300, size=(12,13))

    return

#qPCR_primers_calibration(topA_pc, 'rho, gyrAB, parCE, topB, rnhAB primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\EcTopoI_RNAP_uncoupling_effects_qPCR\\rnhAB_expression_analysis_qPCR\gyrAB_parCE_rnhAB_topB_calibration_16-22_09_20\\gyrAB_parCE_rnhAB_topB_primers_calibration_09_20.png")
#qPCR_primers_calibration(rpmH_pc, 'rpmH primers calibration', "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\qPCR_results\CTD_effect_qPCR_primers_calibration\RpmH_primers_calibration_18_03_20.png")


#################
### CTD expression analysis.
### !!!!!RETURN TO ITS ORIGINAL CODE FOR CTD!!!!
#################


#Path to qPCR data table.
#CTD_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\rnhAB_17_09_20_1\\rnhAB_17_09_20_1_admin_2020-09-17 17-44-05_Alpha -  Quantification Plate View Results.xlsx"
#CTD=pd.read_excel(CTD_table, sheet_name='Processed_data_rnhAB1', header=0, index_col=0)
#print(CTD)


#Plot data.
def qPCR_CTD_expression(dataframe, outpath):
    
    ###
    ##Plot all Cts.
    ###
    
    fig, plot_av=plt.subplots(1,3,figsize=(20,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['rho2', 'rho1', 'cysG', 'rhhB2', 'rnhB1', 'rnhA1',
                'rho2', 'rho1', 'cysG', 'rhhB2', 'rnhB1', 'rnhA1',
                'rho2', 'rho1', 'cysG', 'rhhB2', 'rnhB1', 'rnhA1',
                'rho2', 'rho1', 'cysG', 'rhhB2', 'rnhB1', 'rnhA1']
    print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,3.4,4,
              6.4,7,7.6,8.2,8.8,9.4,
              11.8,12.4,13,13.6,14.2,14.8,
              17.2,17.8,18.4,19,19.6,20.2]
    print(len(X_coords))
    
    Mean_Ct=dataframe.loc[:, 'GFP+_av'].tolist() + dataframe.loc[:, 'GFP-_av'].tolist() + dataframe.loc[:, 'Y319F+_av'].tolist() + dataframe.loc[:, 'Y319F-_av'].tolist()
    print(len(Mean_Ct))
    
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa', '#e8d844', '#94ebc7',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa', '#e8d844', '#94ebc7',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa', '#e8d844', '#94ebc7',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa', '#e8d844', '#94ebc7']
    print(len(Colors))
    
    Data_points=[dataframe.loc['rho2', ['GFP+1', 'GFP+2']].tolist(), dataframe.loc['rho1', ['GFP+1', 'GFP+2']].tolist(), dataframe.loc['cysG', ['GFP+1', 'GFP+2']].tolist(), dataframe.loc['rnhB2', ['GFP+1', 'GFP+2']].tolist(), dataframe.loc['rnhB1', ['GFP+1', 'GFP+2']].tolist(), dataframe.loc['rnhA1', ['GFP+1', 'GFP+2']].tolist(), 
                 dataframe.loc['rho2', ['GFP-1', 'GFP-2', 'GFP-3']].tolist(), dataframe.loc['rho1', ['GFP-1', 'GFP-2', 'GFP-3']].tolist(), dataframe.loc['cysG', ['GFP-1', 'GFP-3']].tolist(), dataframe.loc['rnhB2', ['GFP-1', 'GFP-2', 'GFP-3']].tolist(), dataframe.loc['rnhB1', ['GFP-1', 'GFP-2', 'GFP-3']].tolist(), dataframe.loc['rnhA1', ['GFP-1', 'GFP-2', 'GFP-3']].tolist(), 
                 dataframe.loc['rho2', ['Y319F+1', 'Y319F+2', 'Y319F+3']].tolist(), dataframe.loc['rho1', ['Y319F+1', 'Y319F+2', 'Y319F+3']].tolist(), dataframe.loc['cysG', ['Y319F+1', 'Y319F+2', 'Y319F+3']].tolist(), dataframe.loc['rnhB2', ['Y319F+1', 'Y319F+3']].tolist(), dataframe.loc['rnhB1', ['Y319F+1', 'Y319F+3']].tolist(), dataframe.loc['rnhA1', ['Y319F+1', 'Y319F+3']].tolist(), 
                 dataframe.loc['rho1', ['Y319F-1', 'Y319F-2']].tolist(), dataframe.loc['cysG', ['Y319F-1', 'Y319F-2']].tolist(), dataframe.loc['rnhB2', ['Y319F-1', 'Y319F-2']].tolist(), dataframe.loc['rnhB1', ['Y319F-1', 'Y319F-2']].tolist(), dataframe.loc['rnhA1', ['Y319F-1', 'Y319F-2']].tolist(), ]
    
    #Plot data.
    plot_av[0].bar(X_coords, Mean_Ct, width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    #plot_av[0].plot(X_coords, Data_points, 'ko', markersize=1)
    plot_av[0].set_ylabel('Ct', size=20)
    plot_av[0].set_xticks(X_coords)
    plot_av[0].set_xticklabels(Conditions, rotation=90, size=8)   
    plot_av[0].tick_params(axis='x', which='major', pad=0.5)
    plot_av[0].set_ylim([0, 32])
    plot_av[0].annotate('pCA25 GFP\nIPTG 1mM', xy=(1, 25), xycoords='data', size=9)
    plot_av[0].annotate('pCA25 GFP\nGlc 0.5%', xy=(6.4, 25), xycoords='data', size=9)
    plot_av[0].annotate('pCA25 Y319F\nIPTG 1mM', xy=(11.8, 25), xycoords='data', size=9)
    plot_av[0].annotate('pCA25 Y319F\nGlc 0.5%', xy=(17.2, 25), xycoords='data', size=9)
    
    plt.tight_layout()
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(12,3))
    
    print('Here')
    

    ###
    ##Plot fold change. CTD over rpmH1.
    ###
    
    #Extract CTD and reference data - rpmH1.
    rpmH1_Data_points=[dataframe.loc['rpmH1', ['wt1', 'wt2', 'wt3']].tolist(), dataframe.loc['rpmH1', ['Glc1', 'Glc2', 'Glc3']].tolist(), dataframe.loc['rpmH1', ['IPTG1', 'IPTG2', 'IPTG3']].tolist()]
    CTD_Data_points=[dataframe.loc['CTD', ['wt1', 'wt2', 'wt3']].tolist(), dataframe.loc['CTD', ['Glc1', 'Glc2', 'Glc3']].tolist(), dataframe.loc['CTD', ['IPTG1', 'IPTG2', 'IPTG3']].tolist()]
    
    #Calculate fold-enrichment.
    FE_av=[]
    FE_data_points=[]
    for i in range(len(rpmH1_Data_points)):
        data_points=[]
        for j in range(len(rpmH1_Data_points[i])):
            delta_Ct=rpmH1_Data_points[i][j]-CTD_Data_points[i][j]
            FE=2**delta_Ct
            data_points.append(FE)
        FE_mean=np.mean(data_points)
        FE_av.append(FE_mean)
        FE_data_points.append(data_points)
    
    Conditions=['-\nplasmid', 'pCA25 CTD\nGlc 0.5%', 'pCA25 CTD\nIPTG 1mM']
    
    X_coords=[1,2,3]
    
        
    plot_av[1].bar(X_coords, FE_av, width=0.6, color='#e62552', edgecolor='k', linewidth=0.6)
    plot_av[1].plot(X_coords, FE_data_points, 'ko', markersize=1)
    plot_av[1].set_ylabel('Fold enrichment\nCTD/rpmH1', size=16)
    plot_av[1].set_xticks(X_coords)
    plot_av[1].set_xticklabels(Conditions, rotation=0, size=12)     
    plot_av[1].set_yscale('log')   
    
    
    ###
    ##Plot fold change. CTD over topA_7.
    ###
    
    #Extract CTD and reference data - rpmH1.
    topA_7_Data_points=[dataframe.loc['7top', ['wt1', 'wt2', 'wt3']].tolist(), dataframe.loc['7top', ['Glc1', 'Glc2', 'Glc3']].tolist(), dataframe.loc['7top', ['IPTG1', 'IPTG2', 'IPTG3']].tolist()]
    CTD_Data_points=[dataframe.loc['CTD', ['wt1', 'wt2', 'wt3']].tolist(), dataframe.loc['CTD', ['Glc1', 'Glc2', 'Glc3']].tolist(), dataframe.loc['CTD', ['IPTG1', 'IPTG2', 'IPTG3']].tolist()]
    
    #Calculate fold-enrichment.
    top7_FE_av=[]
    top7_FE_data_points=[]
    for i in range(len(topA_7_Data_points)):
        data_points=[]
        for j in range(len(topA_7_Data_points[i])):
            delta_Ct=topA_7_Data_points[i][j]-CTD_Data_points[i][j]
            FE=2**delta_Ct
            data_points.append(FE)
        FE_mean=np.mean(data_points)
        top7_FE_av.append(FE_mean)
        top7_FE_data_points.append(data_points)
    
    Conditions=['-\nplasmid', 'pCA25 CTD\nGlc 0.5%', 'pCA25 CTD\nIPTG 1mM']
    
    X_coords=[1,2,3]
    
        
    plot_av[2].bar(X_coords, top7_FE_av, width=0.6, color='#59a9de', edgecolor='k', linewidth=0.6)
    plot_av[2].plot(X_coords, top7_FE_data_points, 'ko', markersize=1)
    plot_av[2].set_ylabel('Fold enrichment\nCTD/topA7', size=16)
    plot_av[2].set_xticks(X_coords)
    plot_av[2].set_xticklabels(Conditions, rotation=0, size=12)     
    plot_av[2].set_yscale('log')       
    
    
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, size=(12,3))
   

    #return

#qPCR_CTD_expression(CTD, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\rnhAB_17_09_20_1\rnhAB_expression_qPCR_1.png")


#################
### 3' RNA decay.
#################


#Path to qPCR data table.
#rpmH_dec_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\EcTopoI_RNAP_uncoupling_effects_qPCR\CTD_overexpression_qPCR_results\CTD_effect_21_03_20\RpmH_CTD_effect_21_03_20_admin_2020-03-18 21-12-35_CT018440 -  Quantification Plate View Results.xlsx"
#topA_dec_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\EcTopoI_RNAP_uncoupling_effects_qPCR\CTD_overexpression_qPCR_results\CTD_effect_21_03_20\TopoA_CTD_effect_admin_2020-03-18 21-06-35_CT013337 -  Quantification Plate View Results.xlsx"
#topA_dec=pd.read_excel(topA_dec_table, sheet_name='Data_processed_topA', header=0, index_col=0)
#rpmH_dec=pd.read_excel(rpmH_dec_table, sheet_name='Data_processed_rpmH', header=0, index_col=0)
#print(topA_dec)
#print(rpmH_dec)


#Compute fold enrichment.
def FE_calc(ar, primers_eff, length):

    #Calculate FE for data points.
    FE_ar=[]
    for i in range(len(ar)):
        replicas_data=[]
        for j in range(len(ar[i])):
            FE=(length[-1]*(primers_eff[-1]**ar[-1][j]))/(length[i]*(primers_eff[i]**ar[i][j]))
            replicas_data.append(FE)
        FE_ar.append(replicas_data) 
        
    return FE_ar

#Compute fold enrichment mean and standard deviation.
def FE_mean_and_std(ar):
    Points_mean=[]
    Points_std=[]
    for points_set in ar:
        points_set_mean=np.mean(points_set)
        points_set_std=np.std(points_set)
        Points_mean.append(points_set_mean)
        Points_std.append(points_set_std)
    
    return Points_mean, Points_std


#Plot data.
def qPCR_CTD_effect(dataframe, outpath):
    
    ###
    ##Plot all Cts.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(12,3), dpi=100)
    
    #Prepare x axis, extract data.
    Conditions=dataframe.index.tolist()
    
    X_coords_c=np.asarray(dataframe.loc[:, 'Distance'].tolist())
    X_coords_l=X_coords_c-25
    X_coords_r=X_coords_c+25
    
    Efficiency=dataframe.loc[:, 'Effectiveness'].tolist()
    
    Length=dataframe.loc[:, 'Length'].tolist()
    
    wt_data_points=[]
    Glc_data_points=[]
    IPTG_data_points=[]
    for pair in Conditions:
        wt_data_points.append(dataframe.loc[pair, ['wt1', 'wt2', 'wt3']].tolist())
        Glc_data_points.append(dataframe.loc[pair, ['Glc1', 'Glc2', 'Glc3']].tolist())
        IPTG_data_points.append(dataframe.loc[pair, ['IPTG1', 'IPTG2', 'IPTG3']].tolist())
    
    #Compute fold enrichment for mean data points.
    FE_wt=FE_calc(wt_data_points, Efficiency, Length)
    FE_Glc=FE_calc(Glc_data_points, Efficiency, Length)
    FE_IPTG=FE_calc(IPTG_data_points, Efficiency, Length)
    
    #Compute mean and standard deviation for fold enrichment.
    Points_mean_wt, Points_std_wt=FE_mean_and_std(FE_wt)
    Points_mean_Glc, Points_std_Glc=FE_mean_and_std(FE_Glc)
    Points_mean_IPTG, Points_std_IPTG=FE_mean_and_std(FE_IPTG)
        
    #Plot data.
    #plot_av.bar(X_coords_l, Points_mean_wt, yerr=Points_std_wt, error_kw=dict(lw=1, capsize=3, capthick=1), width=50, color='#b2e69a', edgecolor='k', linewidth=0.6, label='- plasmid')
    plot_av.bar(X_coords_l, Points_mean_Glc, yerr=Points_std_Glc, error_kw=dict(lw=1, capsize=3, capthick=1), width=50, color='#f598b8', edgecolor='k', linewidth=0.6, label='pCA25 CTD\nGlc 0.5%')
    plot_av.bar(X_coords_r, Points_mean_IPTG, yerr=Points_std_IPTG, error_kw=dict(lw=1, capsize=3, capthick=1), width=50, color='#89d8fa', edgecolor='k', linewidth=0.6, label='pCA25 CTD\nIPTG 1mM')
    
    #plot_av.plot(X_coords_l, FE_wt, 'ko', markersize=1)
    plot_av.plot(X_coords_l, FE_Glc, 'ko', markersize=1)
    plot_av.plot(X_coords_r, FE_IPTG, 'ko', markersize=1)
    
    plot_av.set_ylabel('Fold enrichment', size=20)
    plot_av.set_xticks(X_coords_c, minor=True)
    plot_av.set_xticklabels(Conditions, minor=True, rotation=0, size=12)  
    plot_av.set_xticks(range(0, max(X_coords_c)+200, 500), minor=False)
    plot_av.set_xticklabels(range(0, max(X_coords_c)+200, 500), minor=False, rotation=90, size=7) 
    
    plot_av.tick_params(axis='x', which='minor', pad=1.5)
    plot_av.tick_params(axis='x', which='major', pad=15)
    
    plot_av.set_yticks(np.arange(0, max(Points_mean_wt)+0.2, 0.25), minor=False)
    plot_av.set_yticklabels(np.arange(0, max(Points_mean_wt)+0.2, 0.25), minor=False, rotation=0, size=12)  
    plot_av.spines["top"].set_visible(False)
    plot_av.spines["right"].set_visible(False)  
    plot_av.spines["bottom"].set_linewidth(1.5)
    plot_av.spines["left"].set_linewidth(1.5)    
    
    
    #Place legend outside of a graph. Stolen from: https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    box=plot_av.get_position()
    plot_av.set_position([box.x0, box.y0, box.width * 0.95, box.height])    
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.75), frameon=False)
    
    plt.tight_layout(rect=[0,0,0.95,1])
    plt.show()
    plt.savefig(outpath, dpi=300, size=(12,3))
        

    return

#qPCR_CTD_effect(topA_dec, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\EcTopoI_RNAP_uncoupling_effects_qPCR\CTD_overexpression_qPCR_results\TopA_CTD_expression_effect_qPCR.png")
#qPCR_CTD_effect(rpmH_dec, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\EcTopoI_RNAP_uncoupling_effects_qPCR\CTD_overexpression_qPCR_results\RpmH_CTD_expression_effect_qPCR.png")



#################
### TopA exression in response to CTD over-expression and GFP over-expression.
#################


#Path to qPCR data table.
#TopA_vs_rpmH_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\EcTopoI_RNAP_uncoupling_effects_qPCR\CTD_overexpression_qPCR_results\CTD_induction_detection_23_03_20\CTD_induction_detection_admin_2020-03-23 21-23-44_CT013337 -  Quantification Plate View Results.xlsx"
#TopA_vs_rpmH=pd.read_excel(TopA_vs_rpmH_table, sheet_name='Processed_data_topA_expression', header=0, index_col=0)
#print(TopA_vs_rpmH)


#Plot data.
def qPCR_CTD_expression(dataframe, outpath):
    
    ###
    ##Plot fold change. topA1 or topA7 over rpmH1.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(3,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['plasmid -', '$\it{gfp}$ -', '$\it{gfp}$ +', '$\it{CTD}$ -', '$\it{CTD}$ +']
    print(len(Conditions))
    
    X_coords=[1,2.2,2.8,4,4.6]
    X_coords_main=[1,2.2,2.8,4,4.6]
    print(len(X_coords))    
    
    Colors=['#ff9fe5', '#9febcc', '#37b024', '#f7b783', '#fa6e52']
    
    #Prepare data for bars.
    Ct_points_ar=['wt_av', 'GFP-av', 'GFP+av', 'CTD-av', 'CTD+av']
    Mean_Ct_topA1=dataframe.loc['topA1', Ct_points_ar].tolist()
    print(len(Mean_Ct_topA1))
    
    Mean_Ct_topA7=dataframe.loc['topA7', Ct_points_ar].tolist()
    print(len(Mean_Ct_topA7))    
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['wt_std', 'GFP-std', 'GFP+std', 'CTD-std', 'CTD+std']
    Errors_Ct_topA1=dataframe.loc['topA1', Data_errors_ar].tolist()
    print(len(Errors_Ct_topA1))  
    
    Errors_Ct_topA7=dataframe.loc['topA7', Data_errors_ar].tolist()
    print(len(Errors_Ct_topA7))        

    #Extract topA1 and topA7 data.
    topA1_Data_points=[dataframe.loc['topA1', ['wt1', 'wt2', 'wt3']].tolist(), dataframe.loc['topA1', ['GFP-1', 'GFP-2', 'GFP-3']].tolist(), dataframe.loc['topA1', ['GFP+1', 'GFP+2', 'GFP+3']].tolist(), dataframe.loc['topA1', ['CTD-1', 'CTD-2', 'CTD-3']].tolist(), dataframe.loc['topA1', ['CTD+1', 'CTD+2', 'CTD+3']].tolist()]
    topA7_Data_points=[dataframe.loc['topA7', ['wt1', 'wt2', 'wt3']].tolist(), dataframe.loc['topA7', ['GFP-1', 'GFP-2', 'GFP-3']].tolist(), dataframe.loc['topA7', ['GFP+1', 'GFP+2', 'GFP+3']].tolist(), dataframe.loc['topA7', ['CTD-1', 'CTD-2', 'CTD-3']].tolist(), dataframe.loc['topA7', ['CTD+1', 'CTD+2', 'CTD+3']].tolist()]
    
    #Plot data points.
    plot_av.bar(X_coords, Mean_Ct_topA7, yerr=Errors_Ct_topA7, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords, topA7_Data_points, 'ko', markersize=1)
    plot_av.set_ylabel('Normalized value\n(topA/rpmH)', size=16)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=40, size=12)    
    plot_av.tick_params(axis='x', which='major', pad=-3)
    plot_av.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], rotation=0, size=12)  
    plot_av.set_ylim([0, 1])
    
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=300, size=(3,3))

    return

#qPCR_CTD_expression(TopA_vs_rpmH, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\EcTopoI_RNAP_uncoupling_effects_qPCR\CTD_overexpression_qPCR_results\CTD_induction_detection_23_03_20\TopA7_response_on_CTD_GFP_induction_qPCR.png")


#################
### rnhAB expression analysis.
#################


#Path to qPCR data table.
#rnhAB_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\rnhAB_17_09_20_1\\rnhAB_17_09_20_1_admin_2020-09-17 17-44-05_Alpha -  Quantification Plate View Results.xlsx"
#rnhAB=pd.read_excel(rnhAB_table, sheet_name='Final_data_norm_by_rho2', header=0, index_col=0)
#print(rnhAB)


#Plot data.
def qPCR_rnhAB_expression(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(5,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{rho}$', '$\it{rnhA}$', '$\it{rnhB}$ 1', '$\it{rnhB}$ 2']
    print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8]
    X_coords_main=[1.9,4.9,7.9,10.9]
    print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['GFP-', 'GFP+', 'Y319F-', 'Y319F+']
    Mean_Ct=dataframe.loc['rho', Data_points_ar].tolist() + dataframe.loc['rnhA1', Data_points_ar].tolist() + dataframe.loc['rnhB1', Data_points_ar].tolist() + dataframe.loc['rnhB2', Data_points_ar].tolist()
    print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['GFP-_std', 'GFP+_std', 'Y319F-_std', 'Y319F+_std']
    Errors_Ct=dataframe.loc['rho', Data_errors_ar].tolist() + dataframe.loc['rnhA1', Data_errors_ar].tolist() + dataframe.loc['rnhB1', Data_errors_ar].tolist() + dataframe.loc['rnhB2', Data_errors_ar].tolist()
    print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    print(len(Colors))
    
    #Prepare data for real data points.
    GFP_minus_labs=['GFP-1', 'GFP-2', 'GFP-3', 'GFP-4', 'GFP-5', 'GFP-6']
    GFP_plus_labs=['GFP+1', 'GFP+2', 'GFP+3', 'GFP+4']  
    Y319F_minus_labs=['Y319F-1', 'Y319F-2', 'Y319F-3', 'Y319F-4', 'Y319F-5']
    Y319F_plus_labs=['Y319F+1', 'Y319F+2', 'Y319F+3', 'Y319F+4', 'Y319F+5']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['rho', GFP_minus_labs].tolist(), dataframe.loc['rho', GFP_plus_labs].tolist(), dataframe.loc['rho', Y319F_minus_labs].tolist(), dataframe.loc['rho', Y319F_plus_labs].tolist(),
                 dataframe.loc['rnhA1', GFP_minus_labs].tolist(), dataframe.loc['rnhA1', GFP_plus_labs].tolist(), dataframe.loc['rnhA1', Y319F_minus_labs].tolist(), dataframe.loc['rnhA1', Y319F_plus_labs].tolist(),
                 dataframe.loc['rnhB1', GFP_minus_labs].tolist(), dataframe.loc['rnhB1', GFP_plus_labs].tolist(), dataframe.loc['rnhB1', Y319F_minus_labs].tolist(), dataframe.loc['rnhB1', Y319F_plus_labs].tolist(),
                 dataframe.loc['rnhB2', GFP_minus_labs].tolist(), dataframe.loc['rnhB2', GFP_plus_labs].tolist(), dataframe.loc['rnhB2', Y319F_minus_labs].tolist(), dataframe.loc['rnhB2', Y319F_plus_labs].tolist()]
    print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    print(len(X_coords_for_points))
    print(len(Data_points_ar))
    
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    plot_av.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 1.25])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('$\it{gfp}$-', '$\it{gfp}$+', '$\it{topA}$ Y319F-', '$\it{topA}$ Y319F+'), fontsize=11.3, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(5,3))
    

    return

#qPCR_rnhAB_expression(rnhAB, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\rnhAB_17_09_20_1\\rnhAB_expression_qPCR_1.png")



#################
### gyrAB parC expression analysis.
#################


#Path to qPCR data table.
#gyrAB_parCE_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\gyrAB_parC_22_09_20_1\\admin_2020-09-22 23-15-07_Alpha -  Quantification Plate View Results.xlsx"
#gyrAB_parCE=pd.read_excel(gyrAB_parCE_table, sheet_name='Final_data_gyrAB_parCE_rho2_nm', header=0, index_col=0)
#print(gyrAB_parCE)


#Plot data.
def qPCR_gyrAB_parCE_expression(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(6,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{rho}$', '$\it{gyrA}$ 1', '$\it{gyrA}$ 2', '$\it{gyrB}$ 1', '$\it{gyrB}$ 2', '$\it{parC}$ 1']
    print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8,
              13,13.6,14.2,14.8,
              16,16.6,17.2,17.8]
    X_coords_main=[1.9,4.9,7.9,10.9,13.9,16.9]
    print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['GFP-', 'GFP+', 'Y319F-', 'Y319F+']
    Mean_Ct=dataframe.loc['rho', Data_points_ar].tolist() + dataframe.loc['gyrA1', Data_points_ar].tolist() + dataframe.loc['gyrA2', Data_points_ar].tolist() + dataframe.loc['gyrB1', Data_points_ar].tolist() + dataframe.loc['gyrB2', Data_points_ar].tolist() + dataframe.loc['parC1', Data_points_ar].tolist()
    print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['GFP-_std', 'GFP+_std', 'Y319F-_std', 'Y319F+_std']
    Errors_Ct=dataframe.loc['rho', Data_errors_ar].tolist() + dataframe.loc['gyrA1', Data_errors_ar].tolist() + dataframe.loc['gyrA2', Data_errors_ar].tolist() + dataframe.loc['gyrB1', Data_errors_ar].tolist() + dataframe.loc['gyrB2', Data_errors_ar].tolist() + dataframe.loc['parC1', Data_errors_ar].tolist()
    print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    print(len(Colors))
    
    #Prepare data for real data points.
    GFP_minus_labs=['GFP-1', 'GFP-2', 'GFP-3', 'GFP-4', 'GFP-5', 'GFP-6']
    GFP_plus_labs=['GFP+1', 'GFP+2', 'GFP+3', 'GFP+4', 'GFP+5', 'GFP+6']  
    Y319F_minus_labs=['Y319F-1', 'Y319F-2', 'Y319F-3', 'Y319F-4', 'Y319F-5', 'Y319F-6']
    Y319F_plus_labs=['Y319F+1', 'Y319F+2', 'Y319F+3', 'Y319F+4', 'Y319F+5', 'Y319F+6']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['rho', GFP_minus_labs].tolist(), dataframe.loc['rho', GFP_plus_labs].tolist(), dataframe.loc['rho', Y319F_minus_labs].tolist(), dataframe.loc['rho', Y319F_plus_labs].tolist(),
                 dataframe.loc['gyrA1', GFP_minus_labs].tolist(), dataframe.loc['gyrA1', GFP_plus_labs].tolist(), dataframe.loc['gyrA1', Y319F_minus_labs].tolist(), dataframe.loc['gyrA1', Y319F_plus_labs].tolist(),
                 dataframe.loc['gyrA2', GFP_minus_labs].tolist(), dataframe.loc['gyrA2', GFP_plus_labs].tolist(), dataframe.loc['gyrA2', Y319F_minus_labs].tolist(), dataframe.loc['gyrA2', Y319F_plus_labs].tolist(),
                 dataframe.loc['gyrB1', GFP_minus_labs].tolist(), dataframe.loc['gyrB1', GFP_plus_labs].tolist(), dataframe.loc['gyrB1', Y319F_minus_labs].tolist(), dataframe.loc['gyrB1', Y319F_plus_labs].tolist(),
                 dataframe.loc['gyrB2', GFP_minus_labs].tolist(), dataframe.loc['gyrB2', GFP_plus_labs].tolist(), dataframe.loc['gyrB2', Y319F_minus_labs].tolist(), dataframe.loc['gyrB2', Y319F_plus_labs].tolist(),
                 dataframe.loc['parC1', GFP_minus_labs].tolist(), dataframe.loc['parC1', GFP_plus_labs].tolist(), dataframe.loc['parC1', Y319F_minus_labs].tolist(), dataframe.loc['parC1', Y319F_plus_labs].tolist()]
    print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    print(len(X_coords_for_points))
    print(len(Data_points_ar))
    
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    plot_av.set_yticklabels([0, 1, 2, 3, 4, 5, 6, 7], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 7])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('$\it{gfp}$-', '$\it{gfp}$+', '$\it{topA}$ Y319F-', '$\it{topA}$ Y319F+'), fontsize=11.3, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(6,3))
    

    return

#qPCR_gyrAB_parCE_expression(gyrAB_parCE, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\gyrAB_parC_22_09_20_1\\gyrAB_parC_expression_qPCR.png")


#################
### parCE topB expression analysis.
#################


#Path to qPCR data table.
#parCE_topB_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\parCE_topB_23_09_20_1\\parCE_topB_23_09_20_1_admin_2020-09-23 22-23-18_Alpha -  Quantification Plate View Results.xlsx"
#parCE_topB=pd.read_excel(parCE_topB_table, sheet_name='Final_data_parCE_topB_nm_rho2', header=0, index_col=0)
#print(parCE_topB)


#Plot data.
def qPCR_parCE_topB_expression(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(6,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{rho}$', '$\it{parC}$ 2', '$\it{parE}$ 1', '$\it{parE}$ 2', '$\it{topB}$ 1', '$\it{topB}$ 2']
    print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8,
              13,13.6,14.2,14.8,
              16,16.6,17.2,17.8]
    X_coords_main=[1.9,4.9,7.9,10.9,13.9,16.9]
    print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['GFP-', 'GFP+', 'Y319F-', 'Y319F+']
    Mean_Ct=dataframe.loc['rho', Data_points_ar].tolist() + dataframe.loc['parC2', Data_points_ar].tolist() + dataframe.loc['parE1', Data_points_ar].tolist() + dataframe.loc['parE2', Data_points_ar].tolist() + dataframe.loc['topB1', Data_points_ar].tolist() + dataframe.loc['topB2', Data_points_ar].tolist()
    print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['GFP-_std', 'GFP+_std', 'Y319F-_std', 'Y319F+_std']
    Errors_Ct=dataframe.loc['rho', Data_errors_ar].tolist() + dataframe.loc['parC2', Data_errors_ar].tolist() + dataframe.loc['parE1', Data_errors_ar].tolist() + dataframe.loc['parE2', Data_errors_ar].tolist() + dataframe.loc['topB1', Data_errors_ar].tolist() + dataframe.loc['topB2', Data_errors_ar].tolist()
    print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    print(len(Colors))
    
    #Prepare data for real data points.
    GFP_minus_labs=['GFP-1', 'GFP-2', 'GFP-3', 'GFP-4', 'GFP-5', 'GFP-6']
    GFP_plus_labs=['GFP+1', 'GFP+2', 'GFP+3', 'GFP+4', 'GFP+5', 'GFP+6']  
    Y319F_minus_labs=['Y319F-1', 'Y319F-2', 'Y319F-3', 'Y319F-4', 'Y319F-5', 'Y319F-6']
    Y319F_plus_labs=['Y319F+1', 'Y319F+2', 'Y319F+3', 'Y319F+4', 'Y319F+5', 'Y319F+6']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['rho', GFP_minus_labs].tolist(), dataframe.loc['rho', GFP_plus_labs].tolist(), dataframe.loc['rho', Y319F_minus_labs].tolist(), dataframe.loc['rho', Y319F_plus_labs].tolist(),
                 dataframe.loc['parC2', GFP_minus_labs].tolist(), dataframe.loc['parC2', GFP_plus_labs].tolist(), dataframe.loc['parC2', Y319F_minus_labs].tolist(), dataframe.loc['parC2', Y319F_plus_labs].tolist(),
                 dataframe.loc['parE1', GFP_minus_labs].tolist(), dataframe.loc['parE1', GFP_plus_labs].tolist(), dataframe.loc['parE1', Y319F_minus_labs].tolist(), dataframe.loc['parE1', Y319F_plus_labs].tolist(),
                 dataframe.loc['parE2', GFP_minus_labs].tolist(), dataframe.loc['parE2', GFP_plus_labs].tolist(), dataframe.loc['parE2', Y319F_minus_labs].tolist(), dataframe.loc['parE2', Y319F_plus_labs].tolist(),
                 dataframe.loc['topB1', GFP_minus_labs].tolist(), dataframe.loc['topB1', GFP_plus_labs].tolist(), dataframe.loc['topB1', Y319F_minus_labs].tolist(), dataframe.loc['topB1', Y319F_plus_labs].tolist(),
                 dataframe.loc['topB2', GFP_minus_labs].tolist(), dataframe.loc['topB2', GFP_plus_labs].tolist(), dataframe.loc['topB2', Y319F_minus_labs].tolist(), dataframe.loc['topB2', Y319F_plus_labs].tolist()]
    print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    print(len(X_coords_for_points))
    print(len(Data_points_ar))
    
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    plot_av.set_yticklabels([0, 0.2, 0.4, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 2])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('$\it{gfp}$-', '$\it{gfp}$+', '$\it{topA}$ Y319F-', '$\it{topA}$ Y319F+'), fontsize=11.3, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(6,3))
    

    return

#qPCR_parCE_topB_expression(parCE_topB, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\parCE_topB_23_09_20_1\\parCE_topB_expression_qPCR.png")


#################
### rnhAB gyrAB parCE topB expression analysis. All data together.
#################


#Path to qPCR data table.
#data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\rnhAB_expression_analysis.xlsx"
#data_tab=pd.read_excel(data_table, sheet_name='qPCR_data_rho2_norm', header=0, index_col=0)
#print(data_tab)


#Plot data.
def qPCR_expression(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(12,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{rho}$', '$\it{rnhA}$', '$\it{rnhB}$ 1', '$\it{rnhB}$ 2', '$\it{gyrA}$ 1', '$\it{gyrA}$ 2', '$\it{gyrB}$ 1', '$\it{gyrB}$ 2', '$\it{parC}$ 1', '$\it{parC}$ 2', '$\it{parE}$ 1', '$\it{parE}$ 2', '$\it{topB}$ 1', '$\it{topB}$ 2']
    print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8,
              13,13.6,14.2,14.8,
              16,16.6,17.2,17.8,
              19, 19.6, 20.2, 20.8,
              22, 22.6, 23.2, 23.8,
              25, 25.6, 26.2, 26.8,
              28, 28.6, 29.2, 29.8,
              31, 31.6, 32.2, 32.8,
              34, 34.6, 35.2, 35.8,
              37, 37.6, 38.2, 38.8,
              40, 40.6, 41.2, 41.8]
    X_coords_main=[1.9,4.9,7.9,10.9,13.9,16.9, 19.9, 22.9, 25.9, 28.9, 31.9, 34.9, 37.9, 40.9]
    print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['GFP-', 'GFP+', 'Y319F-', 'Y319F+']
    Mean_Ct=dataframe.loc['rho', Data_points_ar].tolist() + dataframe.loc['rnhA1', Data_points_ar].tolist() + dataframe.loc['rnhB1', Data_points_ar].tolist() + dataframe.loc['rnhB2', Data_points_ar].tolist() + dataframe.loc['gyrA1', Data_points_ar].tolist() + dataframe.loc['gyrA2', Data_points_ar].tolist() + dataframe.loc['gyrB1', Data_points_ar].tolist() + dataframe.loc['gyrB2', Data_points_ar].tolist() + dataframe.loc['parC1', Data_points_ar].tolist() + dataframe.loc['parC2', Data_points_ar].tolist() + dataframe.loc['parE1', Data_points_ar].tolist() + dataframe.loc['parE2', Data_points_ar].tolist() + dataframe.loc['topB1', Data_points_ar].tolist() + dataframe.loc['topB2', Data_points_ar].tolist()
    print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['GFP-_std', 'GFP+_std', 'Y319F-_std', 'Y319F+_std']
    Errors_Ct=dataframe.loc['rho', Data_errors_ar].tolist() + dataframe.loc['rnhA1', Data_errors_ar].tolist() + dataframe.loc['rnhB1', Data_errors_ar].tolist() + dataframe.loc['rnhB2', Data_errors_ar].tolist() + dataframe.loc['gyrA1', Data_errors_ar].tolist() + dataframe.loc['gyrA2', Data_errors_ar].tolist() + dataframe.loc['gyrB1', Data_errors_ar].tolist() + dataframe.loc['gyrB2', Data_errors_ar].tolist() + dataframe.loc['parC1', Data_errors_ar].tolist() + dataframe.loc['parC2', Data_errors_ar].tolist() + dataframe.loc['parE1', Data_errors_ar].tolist() + dataframe.loc['parE2', Data_errors_ar].tolist() + dataframe.loc['topB1', Data_errors_ar].tolist() + dataframe.loc['topB2', Data_errors_ar].tolist()
    print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',        
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    print(len(Colors))
    
    #Prepare data for real data points.
    GFP_minus_labs=['GFP-1', 'GFP-2', 'GFP-3', 'GFP-4', 'GFP-5', 'GFP-6']
    GFP_plus_labs=['GFP+1', 'GFP+2', 'GFP+3', 'GFP+4', 'GFP+5', 'GFP+6']  
    Y319F_minus_labs=['Y319F-1', 'Y319F-2', 'Y319F-3', 'Y319F-4', 'Y319F-5', 'Y319F-6']
    Y319F_plus_labs=['Y319F+1', 'Y319F+2', 'Y319F+3', 'Y319F+4', 'Y319F+5', 'Y319F+6']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['rho', GFP_minus_labs].tolist(), dataframe.loc['rho', GFP_plus_labs].tolist(), dataframe.loc['rho', Y319F_minus_labs].tolist(), dataframe.loc['rho', Y319F_plus_labs].tolist(),
                 dataframe.loc['rnhA1', GFP_minus_labs].tolist(), dataframe.loc['rnhA1', GFP_plus_labs].tolist(), dataframe.loc['rnhA1', Y319F_minus_labs].tolist(), dataframe.loc['rnhA1', Y319F_plus_labs].tolist(),
                 dataframe.loc['rnhB1', GFP_minus_labs].tolist(), dataframe.loc['rnhB1', GFP_plus_labs].tolist(), dataframe.loc['rnhB1', Y319F_minus_labs].tolist(), dataframe.loc['rnhB1', Y319F_plus_labs].tolist(),
                 dataframe.loc['rnhB2', GFP_minus_labs].tolist(), dataframe.loc['rnhB2', GFP_plus_labs].tolist(), dataframe.loc['rnhB2', Y319F_minus_labs].tolist(), dataframe.loc['rnhB2', Y319F_plus_labs].tolist(),
                 dataframe.loc['gyrA1', GFP_minus_labs].tolist(), dataframe.loc['gyrA1', GFP_plus_labs].tolist(), dataframe.loc['gyrA1', Y319F_minus_labs].tolist(), dataframe.loc['gyrA1', Y319F_plus_labs].tolist(),
                 dataframe.loc['gyrA2', GFP_minus_labs].tolist(), dataframe.loc['gyrA2', GFP_plus_labs].tolist(), dataframe.loc['gyrA2', Y319F_minus_labs].tolist(), dataframe.loc['gyrA2', Y319F_plus_labs].tolist(),
                 dataframe.loc['gyrB1', GFP_minus_labs].tolist(), dataframe.loc['gyrB1', GFP_plus_labs].tolist(), dataframe.loc['gyrB1', Y319F_minus_labs].tolist(), dataframe.loc['gyrB1', Y319F_plus_labs].tolist(),  
                 dataframe.loc['gyrB2', GFP_minus_labs].tolist(), dataframe.loc['gyrB2', GFP_plus_labs].tolist(), dataframe.loc['gyrB2', Y319F_minus_labs].tolist(), dataframe.loc['gyrB2', Y319F_plus_labs].tolist(),
                 dataframe.loc['parC1', GFP_minus_labs].tolist(), dataframe.loc['parC1', GFP_plus_labs].tolist(), dataframe.loc['parC1', Y319F_minus_labs].tolist(), dataframe.loc['parC1', Y319F_plus_labs].tolist(),                 
                 dataframe.loc['parC2', GFP_minus_labs].tolist(), dataframe.loc['parC2', GFP_plus_labs].tolist(), dataframe.loc['parC2', Y319F_minus_labs].tolist(), dataframe.loc['parC2', Y319F_plus_labs].tolist(),
                 dataframe.loc['parE1', GFP_minus_labs].tolist(), dataframe.loc['parE1', GFP_plus_labs].tolist(), dataframe.loc['parE1', Y319F_minus_labs].tolist(), dataframe.loc['parE1', Y319F_plus_labs].tolist(),
                 dataframe.loc['parE2', GFP_minus_labs].tolist(), dataframe.loc['parE2', GFP_plus_labs].tolist(), dataframe.loc['parE2', Y319F_minus_labs].tolist(), dataframe.loc['parE2', Y319F_plus_labs].tolist(),
                 dataframe.loc['topB1', GFP_minus_labs].tolist(), dataframe.loc['topB1', GFP_plus_labs].tolist(), dataframe.loc['topB1', Y319F_minus_labs].tolist(), dataframe.loc['topB1', Y319F_plus_labs].tolist(),
                 dataframe.loc['topB2', GFP_minus_labs].tolist(), dataframe.loc['topB2', GFP_plus_labs].tolist(), dataframe.loc['topB2', Y319F_minus_labs].tolist(), dataframe.loc['topB2', Y319F_plus_labs].tolist()]
    print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    print(len(X_coords_for_points))
    print(len(Data_points_ar))
    
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    plot_av.set_yticklabels([0,1,2,3,4,5], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 6])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('$\it{gfp}$-', '$\it{gfp}$+', '$\it{topA}$ Y319F-', '$\it{topA}$ Y319F+'), fontsize=14, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(12,3))
    

    return

#qPCR_expression(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\rnhAB_gyrAB_parCE_topB_expression_qPCR.png")


#################
### rnhAB gyrAB parCE topB expression analysis. All data together, polished.
#################


#Path to qPCR data table.
#data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\rnhAB_expression_analysis.xlsx"
#data_tab=pd.read_excel(data_table, sheet_name='qPCR_data_rho2_norm', header=0, index_col=0)
#print(data_tab)


#Plot data.
def qPCR_expression_2(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(6.5,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{gyrA}$', '$\it{gyrB}$', '$\it{parC}$', '$\it{parE}$', '$\it{topB}$', '$\it{rnhA}$', '$\it{rnhB}$']
    print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8,
              13,13.6,14.2,14.8,
              16,16.6,17.2,17.8,
              19, 19.6, 20.2, 20.8,]
    X_coords_main=[1.9,4.9,7.9,10.9,13.9,16.9,19.9]
    print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['GFP-', 'GFP+', 'Y319F-', 'Y319F+']
    Mean_Ct=dataframe.loc['gyrA2', Data_points_ar].tolist() + dataframe.loc['gyrB1', Data_points_ar].tolist() + dataframe.loc['parC2', Data_points_ar].tolist() + dataframe.loc['parE1', Data_points_ar].tolist() + dataframe.loc['topB1', Data_points_ar].tolist() + dataframe.loc['rnhA1', Data_points_ar].tolist() + dataframe.loc['rnhB1', Data_points_ar].tolist()
    print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['GFP-_std', 'GFP+_std', 'Y319F-_std', 'Y319F+_std']
    Errors_Ct=dataframe.loc['gyrA2', Data_errors_ar].tolist() + dataframe.loc['gyrB1', Data_errors_ar].tolist() + dataframe.loc['parC2', Data_errors_ar].tolist() + dataframe.loc['parE1', Data_errors_ar].tolist() + dataframe.loc['topB1', Data_errors_ar].tolist() + dataframe.loc['rnhA1', Data_errors_ar].tolist() + dataframe.loc['rnhB1', Data_errors_ar].tolist()
    print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    print(len(Colors))
    
    #Prepare data for real data points.
    GFP_minus_labs=['GFP-1', 'GFP-2', 'GFP-3', 'GFP-4', 'GFP-5', 'GFP-6']
    GFP_plus_labs=['GFP+1', 'GFP+2', 'GFP+3', 'GFP+4', 'GFP+5', 'GFP+6']  
    Y319F_minus_labs=['Y319F-1', 'Y319F-2', 'Y319F-3', 'Y319F-4', 'Y319F-5', 'Y319F-6']
    Y319F_plus_labs=['Y319F+1', 'Y319F+2', 'Y319F+3', 'Y319F+4', 'Y319F+5', 'Y319F+6']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['gyrA2', GFP_minus_labs].tolist(), dataframe.loc['gyrA2', GFP_plus_labs].tolist(), dataframe.loc['gyrA2', Y319F_minus_labs].tolist(), dataframe.loc['gyrA2', Y319F_plus_labs].tolist(),
                 dataframe.loc['gyrB1', GFP_minus_labs].tolist(), dataframe.loc['gyrB1', GFP_plus_labs].tolist(), dataframe.loc['gyrB1', Y319F_minus_labs].tolist(), dataframe.loc['gyrB1', Y319F_plus_labs].tolist(),                  
                 dataframe.loc['parC2', GFP_minus_labs].tolist(), dataframe.loc['parC2', GFP_plus_labs].tolist(), dataframe.loc['parC2', Y319F_minus_labs].tolist(), dataframe.loc['parC2', Y319F_plus_labs].tolist(),
                 dataframe.loc['parE1', GFP_minus_labs].tolist(), dataframe.loc['parE1', GFP_plus_labs].tolist(), dataframe.loc['parE1', Y319F_minus_labs].tolist(), dataframe.loc['parE1', Y319F_plus_labs].tolist(),
                 dataframe.loc['topB1', GFP_minus_labs].tolist(), dataframe.loc['topB1', GFP_plus_labs].tolist(), dataframe.loc['topB1', Y319F_minus_labs].tolist(), dataframe.loc['topB1', Y319F_plus_labs].tolist(),
                 dataframe.loc['rnhA1', GFP_minus_labs].tolist(), dataframe.loc['rnhA1', GFP_plus_labs].tolist(), dataframe.loc['rnhA1', Y319F_minus_labs].tolist(), dataframe.loc['rnhA1', Y319F_plus_labs].tolist(),
                 dataframe.loc['rnhB1', GFP_minus_labs].tolist(), dataframe.loc['rnhB1', GFP_plus_labs].tolist(), dataframe.loc['rnhB1', Y319F_minus_labs].tolist(), dataframe.loc['rnhB1', Y319F_plus_labs].tolist(),]
    print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    print(len(X_coords_for_points))
    print(len(Data_points_ar))
    
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    plot_av.set_yticklabels([0,1,2,3,4,5], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 5])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('$\it{gfp}$-', '$\it{gfp}$+', '$\it{topA}$ Y319F-', '$\it{topA}$ Y319F+'), fontsize=14, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(6.5,3))
    

    return

#qPCR_expression_2(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\CTD_and_transcription_qPCR\\rnhAB_expression_analysis_qPCR\\rnhAB_gyrAB_parCE_topB_expression_qPCR_nr.png")


#################
### EcTopoI Topo-Seq qPCR.
#################


#Path to qPCR data table.
data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Poisoned_TopoI_mutant\qPCR\EcTopoI_Topo-Seq_25_11_20\\EcTopoI_Topo_Seq_pBAD_admin_2020-11-25 22-50-23_Alpha -  Quantification Plate View Results.xlsx"
data_tab=pd.read_excel(data_table, sheet_name='Data_for_plot', header=0, index_col=0)
print(data_tab)


#Plot data.
def qPCR_Topo_Seq(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(6.5,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['Replicate 1', 'Replicate 2']
    print(len(Conditions))
    
    X_coords=[1,1.6, 2.3,2.9, 3.6,4.2,
              5.4,6, 6.7,7.3, 8,8.6,]
    X_coords_main=[2.6,7,]
    print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar_1=['1_IP-/1_Mock-_av', '1_IP+/1_Mock+_av']
    Data_points_ar_2=['2_IP-/2_Mock-_av', '2_IP+/2_Mock+_av']
    Mean_Ct=dataframe.loc['L2394', Data_points_ar_1].tolist() + dataframe.loc['potF', Data_points_ar_1].tolist() + dataframe.loc['dps', Data_points_ar_1].tolist() + dataframe.loc['L2394', Data_points_ar_2].tolist() + dataframe.loc['potF', Data_points_ar_2].tolist() + dataframe.loc['dps', Data_points_ar_2].tolist()
    print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar_1=['1_IP-/1_Mock-_ster', '1_IP+/1_Mock+_ster']
    Data_errors_ar_2=['2_IP-/2_Mock-_ster', '2_IP+/2_Mock+_ster']
    Errors_Ct=dataframe.loc['L2394', Data_errors_ar_1].tolist() + dataframe.loc['potF', Data_errors_ar_1].tolist() + dataframe.loc['dps', Data_errors_ar_1].tolist() + dataframe.loc['L2394', Data_errors_ar_2].tolist() + dataframe.loc['potF', Data_errors_ar_2].tolist() + dataframe.loc['dps', Data_errors_ar_2].tolist()
    print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#ff4564', '#ff99af', '#4589ff', '#99c9ff', '#fff045', '#fff599', 
            '#ff4564', '#ff99af', '#4589ff', '#99c9ff', '#fff045', '#fff599',]
    print(len(Colors))
    
    #Set texture for bars.
    Hatch_ar=['', 'x', '', 'x', '', 'x',
              '', 'x', '', 'x', '', 'x',]
    
    #Prepare data for real data points.
    IP1_min=['1_IP_-_1/1_Mock_-_1', '1_IP_-_2/1_Mock_-_2', '1_IP_-_3/1_Mock_-_3']
    IP1_pl=['1_IP_+_1/1_Mock_+_1', '1_IP_+_2/1_Mock_+_2', '1_IP_+_3/1_Mock_+_3']  
    IP2_min=['2_IP_-_1/2_Mock_-_1', '2_IP_-_2/2_Mock_-_2', '2_IP_-_3/2_Mock_-_3']
    IP2_pl=['2_IP_+_1/2_Mock_+_1', '2_IP_+_2/2_Mock_+_2', '2_IP_+_3/2_Mock_+_3'] 
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['L2394', IP1_min].tolist(), dataframe.loc['L2394', IP1_pl].tolist(), dataframe.loc['potF', IP1_min].tolist(), dataframe.loc['potF', IP1_pl].tolist(), dataframe.loc['dps', IP1_min].tolist(), dataframe.loc['dps', IP1_pl].tolist(),
                 dataframe.loc['L2394', IP2_min].tolist(), dataframe.loc['L2394', IP2_pl].tolist(), dataframe.loc['potF', IP2_min].tolist(), dataframe.loc['potF', IP2_pl].tolist(), dataframe.loc['dps', IP2_min].tolist(), dataframe.loc['dps', IP2_pl].tolist(),
                 ]
    print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    print(len(X_coords_for_points))
    print(len(Data_points_ar))
    
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6, hatch=Hatch_ar)
    for bar, hatch in zip(Bars, Hatch_ar):  # loop over bars and hatches to set hatches in correct order
        bar.set_hatch(hatch)   #Taken from: https://stackoverflow.com/questions/55826167/matplotlib-assigning-different-hatch-to-bars
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    #plot_av.set_yticklabels([0,1,2,3], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 1.9])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3],Bars[4],Bars[5]), ('L2394 -IP', 'L2394 +IP', 'potF -IP', 'potF +IP', 'dps -IP', 'dps +IP'), fontsize=14, ncol=3, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, loc='upper left')
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(outpath, dpi=300, size=(6.5,3))
    

    return

qPCR_Topo_Seq(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Poisoned_TopoI_mutant\qPCR\EcTopoI_Topo-Seq_25_11_20\EcTopoI_Topo-Seq_pBAD33_topA_mut_25_11_20.png")
