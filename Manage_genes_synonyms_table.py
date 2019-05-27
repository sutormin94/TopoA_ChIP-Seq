###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Reads table of E. coli W3110 genes synonyms, identifies synonyms o request.
#Additionally, returns information about membrane localisation of query gene.
###############################################

#######
#Packages to be imported.
#######

import pandas as pd

#######
#Variables to be defined.
#######

#Path to the synonyms table.
Path_to_syns_table="F:\E_coli_membrane_proteins\Tables_of_synonyms\E_coli_W3110_based_table_of_synonyms_membrane_info.txt"

#List (table) of genes to be ajusted with synonyms table.
Path_to_genes_list="F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Signal_over_TUs_and_regulonDB_info_eq_len_All_genes.txt"

def read_synonyms_table(inpath):
    syns_dataframe=pd.read_csv(inpath, sep='\t', index_col=0, header=None).T #Genes become columns, not rows.
    return syns_dataframe

Dataframe_of_syns=read_synonyms_table(Path_to_syns_table)

Some_genes_list=pd.read_csv(Path_to_genes_list, sep='\t', header=0)

def return_info(query, syns_dataframe):
    Membrane_assignment=syns_dataframe[query].dropna()[1]
    Synonyms=syns_dataframe[query].dropna()[1:].tolist()
    return Membrane_assignment, Synonyms

return_info('ftsQ', Dataframe_of_syns)

Addit_info_dict={'Gene_name': [], 'Membrane_localization': [], 'Synonyms' : []}
for index, row in Some_genes_list.iterrows():
    gene_name=row['Gene_name']
    mem_and_syns=return_info(gene_name, Dataframe_of_syns)
    Addit_info_dict['Gene_name'].append(gene_name)
    Addit_info_dict['Membrane_localization'].append(mem_and_syns[0])
    Addit_info_dict['Synonyms'].append(mem_and_syns[1])
Addit_info_df=pd.DataFrame(Addit_info_dict, columns=['Gene_name', 'Membrane_localization', 'Synonyms'])

List_and_add_info=pd.merge(Some_genes_list, Addit_info_df, on='Gene_name')
List_and_add_info.to_csv('F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Signal_over_TUs_and_regulonDB_info_eq_len_All_genes_membrane_syns.txt', sep='\t', index=False)    