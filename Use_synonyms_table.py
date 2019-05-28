###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Reads table of E. coli W3110 genes synonyms, identifies synonyms o request.
#Additionally, returns information about membrane localisation of query gene.
###############################################

#######
#Packages to be imported.
#######

import manage_synonyms
from manage_synonyms import Synonyms as syn
import pandas as pd

#######
#Variables to be defined.
#######

#Path to the synonyms table.
Path_to_syns_table="F:\E_coli_membrane_proteins\Tables_of_synonyms\E_coli_W3110_based_table_of_synonyms_membrane_info.txt"

#List (table) of genes to be ajusted with synonyms table.
Path_to_genes_list="F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Signal_over_TUs_and_regulonDB_info_eq_len_All_genes.txt"


Dataframe_of_syns=pd.read_csv(Path_to_syns_table, sep='\t', index_col=0, header=None).T #Genes become columns, not rows.
Some_genes_list=pd.read_csv(Path_to_genes_list, sep='\t', header=0)

syn.return_info('nalA', Dataframe_of_syns)

#######
#Add membrane assignment to some table.
#######

def add_membrane_info(some_genes_dataframe, dataframe_of_syns, path_out):
    Addit_info_dict={'Gene_name': [], 'Membrane_localization': [], 'Synonyms' : []}
    for index, row in some_genes_dataframe.iterrows():
        gene_name=row['Gene_name']
        mem_and_syns=syn.return_info(gene_name, dataframe_of_syns)
        if mem_and_syns!=-1:
            Addit_info_dict['Gene_name'].append(gene_name)
            Addit_info_dict['Membrane_localization'].append(mem_and_syns[1])
            Addit_info_dict['Synonyms'].append(mem_and_syns[2])
        else:
            Addit_info_dict['Gene_name'].append('-')
            Addit_info_dict['Membrane_localization'].append('-')
            Addit_info_dict['Synonyms'].append('-')            
    Addit_info_df=pd.DataFrame(Addit_info_dict, columns=['Gene_name', 'Membrane_localization', 'Synonyms'])
    
    List_and_add_info=pd.merge(some_genes_dataframe, Addit_info_df, how='left', on='Gene_name')
    List_and_add_info.to_csv(path_out, sep='\t', index=False)   
    return


#######
#.
#######