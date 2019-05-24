###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes wig files for IP and mock control,
#Performes normalization on the number of reads,
#Calculates fold enrichment of IP over the mock control.
#Calculates correlation between tracks.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import scipy
import pandas as pd
from pandas import DataFrame

#Path to the working directory.
PWD='F:\Signal_over_TUs'
#Dictionary of input TAB data files (information about the signal over TUs).
TUs_data_dict={'TopoA -Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\TopoA -Rif_over_All genes_15000bp.txt',
               'TopoA +Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\TopoA +Rif_over_All genes_15000bp.txt',
               'TopoIV' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\TopoIV Cfx_over_All genes_15000bp.txt',
               'Gyrase' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\Gyrase Cfx_over_All genes_15000bp.txt',
               'MukB' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\MukB_over_All genes_15000bp.txt',
               'RpoB' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\RpoB_over_All genes_15000bp.txt',
               'H-NS' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\HNS_over_All genes_15000bp.txt',
               'Fis' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\Fis_over_All genes_15000bp.txt',
               'GC' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\GC_over_All genes_15000bp.txt',
               'Expression' : 'C:\Sutor\science\TopoI_Topo-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\DOOR_Mu_del_cor_genes_expression.txt',
               'MatP' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\MatP_over_All genes_15000bp.txt',
               'PolSofi' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\PolSofi_over_All genes_15000bp.txt',
               'RpoS' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\RpoS_over_All genes_15000bp.txt',
               }
#Main table (to be added to).
Main_table='F:\TopoI_ChIP-Seq\Ec_TopoI_data\FE_of_genes\noRif_Rif_Ded_FE_over_Genes_5000bp.xlsx'
#Set of genes.
Set_of_genes='All_genes'

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=PWD
Dir_check_create(Out_path)
Dir_check_create(PWD+'\Signal_of_TUs_tab_all\\'+Set_of_genes)


#######
#Read .tab-file with TopoA FE over TUs.
#######

def read_FE_tables_combine_together(data_dict, main_table, path_out):
    #Read input dataframes.
    Dict_of_dataset={}
    for dataset_name, dataset in data_dict.items():
        data_dataframe=pd.read_csv(dataset, sep='\t')
        print(data_dataframe.iloc[: , 5:6])
        #Dict_of_dataset[dataset_name]=
        
    #Read main table.
    #main_data=pd.read_excel(main_table, header=1)
    #Iteratively merge new data with main table.
    #for dataset_name, dataset in Dict_of_dataset.items():
    #    noRif_Rif_Ded_df_exp_trRNA=pd.merge(noRif_Rif_Ded_df_exp, tRNA_rRNA_all_df)
        
    #noRif_df=pd.read_csv(Path_in+noRif_table_name, sep='\t')
    #Rif_df=pd.read_csv(Path_in+Rif_table_name, sep='\t')
    #Ded_df=pd.read_csv(Path_in+Ded_table_name, sep='\t')
    #Combune into new dataframe.
    #Data={'Gene_name' : noRif_df['Gene_name'], 'Strand' : noRif_df['Strand'], 'Start' : noRif_df['Start'], 'End' : noRif_df['End'],
    #      'TopoA_noRif_FE_US' : noRif_df['TopoA_FE_US'], 'TopoA_noRif_FE_GB' : noRif_df['TopoA_FE_GB'], 'TopoA_noRif_FE_DS' : noRif_df['TopoA_FE_DS'], 
    #      'TopoA_Rif_FE_US' : Rif_df['TopoA_FE_US'], 'TopoA_Rif_FE_GB' : Rif_df['TopoA_FE_GB'], 'TopoA_Rif_FE_DS' : Rif_df['TopoA_FE_DS'],
    #      'TopoA_Ded_FE_US' : Ded_df['TopoA_FE_US'], 'TopoA_Ded_FE_GB' : Ded_df['TopoA_FE_GB'], 'TopoA_Ded_FE_DS' : Ded_df['TopoA_FE_DS']
    #      }
    #noRif_Rif_Ded_df=pd.DataFrame(Data, columns=['Gene_name', 'Strand', 'Start', 'End', 'TopoA_noRif_FE_US', 'TopoA_noRif_FE_GB', 'TopoA_noRif_FE_DS',
    #                                            'TopoA_Rif_FE_US', 'TopoA_Rif_FE_GB', 'TopoA_Rif_FE_DS', 'TopoA_Ded_FE_US', 'TopoA_Ded_FE_GB', 'TopoA_Ded_FE_DS'])
    
    #Read expression data.
    #Expression=pd.read_csv(Expression_data, sep='\t')
    #print(Expression)
    #Make some changes over expression data.
    #for index, row in Expression.iterrows():
    #    #print(index, row)
    #    Genes_ID_changed=row['Gene_name'].rstrip(';')
    #    Expression_dot=row['Expression'].replace(',', '.')        
    #    Expression.at[index, 'Gene_name']=Genes_ID_changed
    #    Expression.at[index, 'Expression']=Expression_dot
    #Combine TopoA FE and expression data.
    #Expression=Expression.astype({'Expression': np.float64})
    #noRif_Rif_Ded_df_exp=pd.merge(noRif_Rif_Ded_df, Expression)
    #Write new dataframe.
    #noRif_Rif_Ded_df_exp_trRNA.to_csv(f'{path_out}FE_of_genes\\noRif_Rif_Ded_FE_Exp_tRNA_rRNA_over_{TU_set_name}_15000bp.txt', sep='\t', index=False)    
    return

read_FE_tables_combine_together(TUs_data_dict, Main_table, Out_path)