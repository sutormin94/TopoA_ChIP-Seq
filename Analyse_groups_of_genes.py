###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes narrowPeaks files with ChIP-Seq peaks identified in different biological replicas,
#return narrowPeak only with reproducible peaks.
###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Path to the all-containing table.
Data_path="F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Signal_over_TUs_and_regulonDB_info_eq_len_All_genes_membrane_syns_TF.xlsx"
#Path to the output histogram.
Outpath="F:\TopoI_ChIP-Seq\Ec_TopoI_data\Figures\Membrane_Promoter_complexity_and_FE\Membrane_loc_promoter_complexity_and_FE.png"



#######
#Plot distributions.
#######

def plot_signal_dist(S1, S2, S3, S4, S1a, S2a, S3a, S4a, outpath):
    all_values=pd.concat([S1, S2, S3, S4, S1a, S2a, S3a, S4a], ignore_index=True)
    min_FE=min(all_values)
    max_FE=max(all_values)
    #Plot distribution of S1.
    fig=plt.figure(figsize=(15, 5), dpi=100)
    bins0=np.arange(min_FE, max_FE, 0.2)
    plot0=plt.subplot2grid((2,4),(0,0), rowspan=1, colspan=1)
    plot0.hist(S1, bins0, color='#ff878b', edgecolor='black', alpha=0.8)
    plot0.annotate(f'Mean FE={round(np.mean(S1),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot0.set_yscale('log')
    plot0.set_xlabel('Fold enrichment', size=13)
    plot0.set_ylabel('Number of genes', size=13)
    plot0.set_title('All genes -Rif', size=15)  
    #Plot distribution of S1a.
    plot0a=plt.subplot2grid((2,4),(1,0), rowspan=1, colspan=1)
    plot0a.hist(S1a, bins0, color='#ff878b', edgecolor='black', alpha=0.8)
    plot0a.annotate(f'Mean FE={round(np.mean(S1a),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot0a.set_yscale('log')
    plot0a.set_xlabel('Fold enrichment', size=13)
    plot0a.set_ylabel('Number of genes', size=13)
    plot0a.set_title('All genes +Rif', size=15)   
    
    #Plot distribution of S2.  
    plot1=plt.subplot2grid((2,4),(0,1), rowspan=1, colspan=1)     
    plot1.hist(S2, bins0, color='#ffce91', edgecolor='black', alpha=0.5)
    plot1.annotate(f'Mean FE={round(np.mean(S2),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot1.set_yscale('log')
    plot1.set_xlabel('Fold enrichment', size=13)
    plot1.set_ylabel('Number of genes', size=13)
    plot1.set_title('Membrane proteins genes \n-Rif', size=15) 
    #Plot distribution of S2a.
    plot1a=plt.subplot2grid((2,4),(1,1), rowspan=1, colspan=1)     
    plot1a.hist(S2a, bins0, color='#ffce91', edgecolor='black', alpha=0.5)
    plot1a.annotate(f'Mean FE={round(np.mean(S2a),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot1a.set_yscale('log')
    plot1a.set_xlabel('Fold enrichment', size=13)
    plot1a.set_ylabel('Number of genes', size=13)
    plot1a.set_title('Membrane proteins genes \n+Rif', size=15) 
    
    #Plot distribution of Ded FE values over IG.
    plot2=plt.subplot2grid((2,4),(0,2), rowspan=1, colspan=1) 
    plot2.hist(S3, bins0, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot2.annotate(f'Mean FE={round(np.mean(S3),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot2.set_yscale('log')
    plot2.set_xlabel('Fold enrichment', size=13)
    plot2.set_ylabel('Number of genes', size=13)
    plot2.set_title('Complex promoter genes \n-Rif', size=15)  
    #Plot distribution of Ded FE values over GB.
    plot2a=plt.subplot2grid((2,4),(1,2), rowspan=1, colspan=1) 
    plot2a.hist(S3a, bins0, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot2a.annotate(f'Mean FE={round(np.mean(S3a),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot2a.set_yscale('log')
    plot2a.set_xlabel('Fold enrichment', size=13)
    plot2a.set_ylabel('Number of genes', size=13)
    plot2a.set_title('Complex promoter genes \n+Rif', size=15)      
    
    #Plot distribution of Ded FE values over IG.
    plot3=plt.subplot2grid((2,4),(0,3), rowspan=1, colspan=1) 
    plot3.hist(S4, bins0, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot3.annotate(f'Mean FE={round(np.mean(S4),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot3.set_yscale('log')
    plot3.set_xlabel('Fold enrichment', size=13)
    plot3.set_ylabel('Number of genes', size=13)
    plot3.set_title('Complex promoter genes \n-Rif', size=15)  
    #Plot distribution of Ded FE values over GB.
    plot3a=plt.subplot2grid((2,4),(1,3), rowspan=1, colspan=1) 
    plot3a.hist(S4a, bins0, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot3a.annotate(f'Mean FE={round(np.mean(S4a),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot3a.set_yscale('log')
    plot3a.set_xlabel('Fold enrichment', size=13)
    plot3a.set_ylabel('Number of genes', size=13)
    plot3a.set_title('Complex promoter genes \n+Rif', size=15)      
        
    
    plt.tight_layout()
    #plt.show()
    plt.savefig(outpath, dpi=300, figsize=(15, 10))
    plt.close()     
    return


#######
#Wrap data.
#######

def wrapper_gene_groups(pathin, pathout):
    #Read data.
    Data_frame=pd.read_excel(pathin, header=0)
    #All genes.
    All_noRif=Data_frame[Data_frame['TopoA -Rif_FE_GB'].notnull()]['TopoA -Rif_FE_GB']
    All_Rif=Data_frame[Data_frame['TopoA +Rif_FE_GB'].notnull()]['TopoA +Rif_FE_GB']
    #Return data based on membrane localization.
    In_membrane_dataframe=Data_frame.loc[Data_frame['Membrane_localization']!='-']
    Membrane_TopoA_noRif=In_membrane_dataframe[In_membrane_dataframe['TopoA -Rif_FE_GB'].notnull()]['TopoA -Rif_FE_GB']
    Membrane_TopoA_Rif=In_membrane_dataframe[In_membrane_dataframe['TopoA +Rif_FE_GB'].notnull()]['TopoA +Rif_FE_GB']
    #Return data based on promoter complexity.
    TF_known_dataframe=Data_frame[Data_frame['Number_of_TF_sites'].notnull()]
    Complex_promoter_TopoA_noRif=TF_known_dataframe[TF_known_dataframe['TopoA -Rif_FE_GB'].notnull()]['TopoA -Rif_FE_GB']
    Complex_promoter_TopoA_Rif=TF_known_dataframe[TF_known_dataframe['TopoA +Rif_FE_GB'].notnull()]['TopoA +Rif_FE_GB']
    #Select data based on membrane localization & promoter complexity.
    In_mem_complex_prom_dataframe=In_membrane_dataframe[In_membrane_dataframe['Number_of_TF_sites'].notnull()]
    IMCP_TopoA_noRif=In_mem_complex_prom_dataframe[In_mem_complex_prom_dataframe['TopoA -Rif_FE_GB'].notnull()]['TopoA -Rif_FE_GB']
    IMCP_TopoA_Rif=In_mem_complex_prom_dataframe[In_mem_complex_prom_dataframe['TopoA +Rif_FE_GB'].notnull()]['TopoA +Rif_FE_GB']    
    
    print(min(All_noRif), np.mean(All_noRif), max(All_noRif))
    plot_signal_dist(All_noRif, Membrane_TopoA_noRif, Complex_promoter_TopoA_noRif, IMCP_TopoA_noRif, 
                     All_Rif, Membrane_TopoA_Rif, Complex_promoter_TopoA_Rif, IMCP_TopoA_Rif, 
                     pathout)
    return

wrapper_gene_groups(Data_path, Outpath)