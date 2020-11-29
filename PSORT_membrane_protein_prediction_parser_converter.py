###############################################
##Dmitry Sutormin, 2020##
##TopoA ChIP-Seq analysi##

####
#Reads results PSORT prediction of membrane proteins.
#Takes DNA sequence of classified genes from MG1655 genome.
#Converts to W3110 MuSGS genome by BLAST.
#Assign features to W3110 genes.
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import re
import matplotlib.pyplot as plt
import os
import pandas as pd

#######
#Data to be used.
#######

#Source genome.
MG1655_genome_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Genomes\E_coli_K-12_MG1655_NC_000913.3.fasta"

#PSORT data.
PSORT_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_membrane_proteins\Databases\PSORT_3_database_E_coli.xlsx"

#W3110 genes.
W3110_genes_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\Intergenic_regions_NO_DPS_NO_rRNA_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_RNAP_signal.xlsx"

#PWD.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_membrane_proteins\PSORT_to_w3110\\"


#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
    return genomefa

#######
#Read PSORT data, return genes multifasta, blast.
#######

def read_PSORT_blast(PSORT_path, source_genome_path, pwd):
    #Read PSORT data.
    PSORT_data=pd.read_excel(PSORT_path, header=0, sheet_name='psortdb-results')
    #Read source genome.
    source_genome=read_genome(source_genome_path)
    
    #Write genes multifasta.
    genes=open(pwd+'MG1655_PSORT_genes.fasta', 'w')
    j=0
    for i in PSORT_data.index:
        #Retrive gene coordinates.
        gene_coordinates=PSORT_data.loc[i, 'Genome_Location']
        gene_coordinates=gene_coordinates.lstrip('c').split('..')
        gene_start=int(gene_coordinates[0])
        gene_end=int(gene_coordinates[1])
        
        #Retrive gene score.
        gene_score=PSORT_data.loc[i, 'Final_Score']
        
        #Retrive gene localization prediction.
        gene_localization=PSORT_data.loc[i, 'Final_Localization']
        gene_localization=re.sub(' ', '', gene_localization)
        
        gene_sequence=source_genome[gene_start:gene_end]
        
        #Write gene sequence.
        genes.write(f'>{gene_start}_{gene_end}_{gene_score}_{gene_localization}\n{gene_sequence}\n')
        j+=1
    genes.close()
    
    print(f'Number of genes from source genome: {j}')
    
    #Create BLAST DB.
    Make_w3110_db=NcbimakeblastdbCommandline(dbtype="nucl", input_file=pwd+'Genome\E_coli_w3110_G_Mu.fasta')    
    print('Making blast database: ' + str(Make_w3110_db))
    Make_w3110_db()
    
    #Blast source genome gene sequences.
    Genes_blast=NcbiblastnCommandline(query=pwd+'MG1655_PSORT_genes.fasta', db=pwd+'Genome\E_coli_w3110_G_Mu.fasta', out=pwd+'MG1655_genes_in_W3110_blast.txt', outfmt=6)  
    print('Run blast of gene sequences: ' + str(Genes_blast))
    Genes_blast()    
    
    return


#read_PSORT_blast(PSORT_data_path, MG1655_genome_path, PWD)


def read_blast_convert(pwd, target_genome_path):
    
    #Read gene sequences blast results.
    Gblast_output=open(pwd+'MG1655_genes_in_W3110_blast.txt', 'r')
    Genes_dict={}
    j=0
    i=0
    for line in Gblast_output:
        j+=1
        line=line.rstrip().split('\t')
        if line[0] in Genes_dict:
            continue
        else:
            Genes_dict[line[0]]=[min([int(line[8]), int(line[9])]), max([int(line[8]), int(line[9])])]
            i+=1
        
    Gblast_output.close()   
    
    print(f'Number of BLAST hits: {j}; Number of unique hits: {i}')
    
    #Read W3110 genes data.
    Target_genome_data=pd.read_excel(target_genome_path, header=0, index_col=0, sheet_name='Intergenic_regions_info')
    
    #Find corresponding genes.
    G1_score_ar=[]
    G1_Loc_ar=[]
    G2_score_ar=[]
    G2_Loc_ar=[]
    for i in Target_genome_data.index:
        G1_start=Target_genome_data.loc[i, 'G1_start']
        G1_end=Target_genome_data.loc[i, 'G1_end']
        G2_start=Target_genome_data.loc[i, 'G2_start']
        G2_end=Target_genome_data.loc[i, 'G2_end']   
        
        check_G1=0
        check_G2=0
        find=[0,0]
        for gene_info, coordinates in Genes_dict.items():
            #For G1.
            if (abs(G1_start-coordinates[0])<10) & (abs(G1_end-coordinates[1])<10):
                check_G1+=1
                
                gene_info=gene_info.split('_')
                score=gene_info[2]
                localization=gene_info[3]
                
                if check_G1<2:
                
                    G1_score_ar.append(score)
                    G1_Loc_ar.append(localization)
                
                else:
                    print(G1_start, G1_end, gene_info, coordinates)                   
                
                find[0]+=1   
            #For G2.
            if (abs(G2_start-coordinates[0])<10) & (abs(G2_end-coordinates[1])<10):
                check_G2+=1
                
                gene_info=gene_info.split('_')
                score=gene_info[2]
                localization=gene_info[3]
                
                if check_G2<2:
                
                    G2_score_ar.append(score)
                    G2_Loc_ar.append(localization)   
                    
                else:
                    print(G2_start, G2_end, gene_info, coordinates)                
                
                find[1]+=1
        
        if check_G1==0:
            G1_score_ar.append(0)
            G1_Loc_ar.append('-')
        if check_G2==0:
            G2_score_ar.append(0)
            G2_Loc_ar.append('-')
        if find!=[1,1]:
            print(find)
        if len(G1_score_ar)!=len(G2_score_ar):
            print('Smth wrong!')
    
    print(len(Target_genome_data.index), len(G1_score_ar), len(G1_Loc_ar), len(G2_score_ar), len(G2_Loc_ar))    
    
    #Add new data to dataframe.
    Target_genome_data['G1_PSORT_score']=G1_score_ar
    Target_genome_data['G1_PSORT_localization']=G1_Loc_ar
    Target_genome_data['G2_PSORT_score']=G2_score_ar
    Target_genome_data['G2_PSORT_localization']=G2_Loc_ar
    
    #Write updated dataframe.
    Target_genome_data.to_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_RNAP_signal_PSORT.xlsx', header=True, index=True, sheet_name='Intergenic_regions_info')
    return

read_blast_convert(PWD, W3110_genes_path)