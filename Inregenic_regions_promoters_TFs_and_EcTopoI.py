###############################################
##Dmitry Sutormin, 2021##
##EcTopoI ChIP-Seq analysis##

# 1) Takes set of genes with assigned expression level and returns intergenic regions (IGRs).
# Select by length intergenic regins appropriate for further analysis: 50bp<length<1000bp. Creates SFig. 10.
# 2) Searches for annotated promoters and transcription factor sites in selected IGRs. 
# Annotation is taken from RegulonDB.
# 3) Adds signal info (wig files, e.g. EcTopoI ChIP-Seq FE) to the selected IGRs.
# Adds information on whether adjacent genes encode proteins targeting to membrane (Ecocyc).
# 4) Adds information on whether adjacent genes encode proteins targeting to membrane (PSORT)
# 5) Remove anomalous IGRs (near dps gene with extremely high EcTopoI peak or near rRNA genes).
# 6) Compares EcTopoI FE (RNAP FE, Expression level) between IGRs sets defined by different features.
# Makes violin-plots: Fig. 4A-C, SFig. 11, SFig. 12, SFig. 14A, SFig. 15.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import Bio
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import pandas as pd
from pandas import DataFrame
import math
import re


##############
##Part 1.
##Read expression data, identify intergenic regions, return sequences of intergenic regions.
##############

#Path to expression level of genes.
Genes_EL="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\DY330_RNA-Seq_genes_EP_del_cor.txt"
#Path to genome sequence of a strain you are working with (.fasta).
Genome_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"
#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\\"


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
#Read expression data, identify intergenic regions, return sequences of intergenic regions.
#######

def identify_IG_regions(filein_path, genome_path, fileout_path):
    
    #Read data.
    EL_dataframe=pd.read_csv(filein_path, sep='\t', header=0, index_col=0)
    EL_dataframe.sort_values('Start', axis=0, ascending=True, inplace=True)
    print(EL_dataframe.head(10))
    
    #Length of IG regions.
    Start=0
    End=0
    Name='yjtD'
    Strand='+'
    Expression=2.5320008333937762
    IG_data=[]
    IG_len_ar=[]
    for index, row in EL_dataframe.iterrows():
        print(row['Gene_name'], row['Start'], row['End'], row['Strand'], row['Expression_E']) 
        IG_len=row['Start']-End
        Intergenic_region_info=[Name, Start, End, Strand, Expression, IG_len, row['Gene_name'], row['Start'], row['End'], row['Strand'], row['Expression_E']]
        IG_len_ar.append(IG_len)
        IG_data.append(Intergenic_region_info)
        Start=row['Start']
        End=row['End']
        Name=row['Gene_name']
        Strand=row['Strand']
        Expression=row['Expression_E']
        
    print(min(IG_len_ar), max(IG_len_ar))
    
    #Plot distribution of intergenic regions length.
    plt.hist(IG_len_ar, bins=10000)
    plt.axvline(50, linestyle='--', linewidth=1, color='black')
    plt.axvline(1000, linestyle='--', linewidth=1, color='black')
    plt.yscale('log')
    plt.xlim([-10, 2000])
    plt.xlabel('Length of intergenic region', size=15)
    plt.ylabel('Number of intergenic regions, log', size=15)
    plt.show()
    plt.savefig(fileout_path+'Distribution_of_intergenic_regions_length.png', dpi=300)
    
    #Read genome sequence.
    Genome_fasta=read_genome(genome_path)
    
    #Sort data by the length of IG regions and write. Also write sequences of IG regions.
    IG_data_sorted=sorted(IG_data, key = lambda x: int(x[5]))
    fileout_IG_data=open(fileout_path+'DY330_intergenic_regions_all_TEST.txt', 'w')
    fileout_IG_seqs=open(fileout_path+'DY330_intergenic_regions_sequences_all_TEST.fasta', 'w')
    for ele in IG_data_sorted:
        fileout_IG_data.write(f'{ele[0]}\t{ele[1]}\t{ele[2]}\t{ele[3]}\t{ele[4]}\t{ele[5]}\t{ele[6]}\t{ele[7]}\t{ele[8]}\t{ele[9]}\t{ele[10]}\n')
        fileout_IG_seqs.write(f'>{ele[0]}_{ele[1]}_{ele[2]}_{ele[3]}_{ele[4]}_{ele[5]}_{ele[6]}_{ele[7]}_{ele[8]}_{ele[9]}_{ele[10]}\n{Genome_fasta[ele[2]:ele[7]]}\n')
    fileout_IG_data.close()
    fileout_IG_seqs.close()
    return

#identify_IG_regions(Genes_EL, Genome_path, PWD)



##############
##Part 2.
##Read promoters data from RegulonDB, blast in a database constructed from a collection of intergenic sequences (see Part 1).
##Read TFs data from RegulonDB, blast in a database constructed from a collection of intergenic sequences (see Part 1).
##############

#Data on promoters from RegulonDB (for E. coli MG1655).
Promoters_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\RegulonDB_PromoterSet.txt"
#Data on transcription factor binding sites from RegulonDB (for E. coli MG1655).
Transcription_Factors_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\RegulonDB_BindingSiteSet.txt"

#######
#Read promoters data. 
#######

def read_promoters(PD_path, P_seq_path):
    PD_file=open(PD_path, 'r')
    PS_file=open(P_seq_path+'RegulonDB_promoters_sequences.fasta', 'w')
    for line in PD_file:
        if line[0]!='#':
            line=line.rstrip().split('\t')
            Pname=line[1]
            Pstrand=line[2]
            Psigma=line[4]
            Psequence=line[5]
            Pevidence=line[7]
            #Write promoter sequence and data.
            PS_file.write(f'>{Pname}_{Pstrand}_{Psigma}_{Pevidence}\n{Psequence}\n')
    PD_file.close()
    PS_file.close()
    return


#######
#Read TFs data. 
#######

def read_TFs(TF_path, TF_seq_path):
    TF_file=open(TF_path, 'r')
    TFS_file=open(TF_seq_path+'RegulonDB_TF_sites_sequences.fasta', 'w')
    TF_start=-1
    TF_end=-1
    for line in TF_file:
        if line[0]!='#':
            line=line.rstrip().split('\t')
            if (int(line[3])!=0) & (int(line[4])!=0) & (len(line)==14):
                print(len(line))
                if (int(line[3])!=TF_start) & (int(line[4])!=TF_end):
                    print(line)
                    TFname=line[1]
                    TF_start=int(line[3])
                    TF_end=int(line[4])
                    TFstrand=line[5]
                    TFeffect=line[8]
                    TFsequence=line[11]
                    TFevidence=line[13]
                    #Write TF sequence and data.
                    TFS_file.write(f'>{TFname}_{TFstrand}_{TFeffect}_{TFevidence}\n{TFsequence}\n')
    TF_file.close()
    TFS_file.close()
    return


#######
#Create database from a collection of intergenic regions sequences. 
#Blast promoters sequences and TF sites sequences in the database.
#######

def db_create_and_blast(PWD):
    
    #Create blast database.
    Make_IG_sequences_db=NcbimakeblastdbCommandline(dbtype="nucl", input_file=PWD+'DY330_intergenic_regions_sequences_filtrated_50_1000bp.fasta')    
    print('Making blast database: ' + str(Make_IG_sequences_db))
    Make_IG_sequences_db()
    
    #Blast promoters.
    Promoters_blast=NcbiblastnCommandline(query=PWD+'RegulonDB_promoters_sequences.fasta', db=PWD+'DY330_intergenic_regions_sequences_filtrated_50_1000bp.fasta', out=PWD+'RegulonDB_promoters_sequences_blast_results_IG_50_1000bp.txt', outfmt=6)  
    print('Run blast of promoters sequences: ' + str(Promoters_blast))
    Promoters_blast()
    
    #Blast TF sites.
    TFS_blast=NcbiblastnCommandline(query=PWD+'RegulonDB_TF_sites_sequences.fasta', db=PWD+'DY330_intergenic_regions_sequences_filtrated_50_1000bp.fasta', out=PWD+'RegulonDB_TF_sites_sequences_blast_results_IG_50_1000bp.txt', outfmt=6)  
    print('Run blast of TF sites sequences: ' + str(TFS_blast))
    TFS_blast()      
    
    
    return

#read_promoters(Promoters_data_path, PWD)
#read_TFs(Transcription_Factors_data_path, PWD)
#db_create_and_blast(PWD)



##############
##Part 3.
##Read results of TF sites and promoters blast.
##Combine data about intergenic regions and assign them with a fold enrichment for EcTopoI.
##Add information on whether adjacent genes encode proteins targeting to membrane.
##############

#Dictionary with .wig files containing continuous signal to be assigned with the set of IGRs.
EcTopoI_data_dict={'-Rif-CTD' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_minus_average_FE_3_4_6.wig',
                   '+Rif-CTD' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_plus_average_FE_1_2_3.wig',
                   '-Rif+CTD' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_minus_average_FE.wig',
                   '+Rif+CTD' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_plus_average_FE.wig',
                   'RNAP' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\Pol_Sofi_LB_w3110_for_Mu.wig'}

#######
#Read promoters and TF sites blast output and combine.
#######

def read_blast_output_combine(pwd):
    #Read promoter sequences blast results.
    Pblast_output=open(pwd+'RegulonDB_promoters_sequences_blast_results_IG_50_1000bp_polished.txt', 'r')
    IG_info_dict={}
    for line in Pblast_output:
        line=line.rstrip().split('\t')
        if line[1] in IG_info_dict:
            IG_info_dict[line[1]]['Promoter_info'].append(line[0])
        else:
            IG_info_dict[line[1]]={'Promoter_info' : [line[0]], 'TF_info' : []}
    Pblast_output.close()
    
    #Read TF sites sequences blast results.
    TFblast_output=open(pwd+'RegulonDB_TF_sites_sequences_blast_results_IG_50_1000bp.txt', 'r')
    for line in TFblast_output:
        line=line.rstrip().split('\t')
        if line[1] in IG_info_dict:
            IG_info_dict[line[1]]['TF_info'].append(line[0])
        else:
            IG_info_dict[line[1]]={'Promoter_info' : [], 'TF_info' : [line[0]]}    
    TFblast_output.close()
    return IG_info_dict


#######
#Parsing WIG file.
#######

def score_data_parser(inpath, param_name):
    param_file=open(inpath, 'r')
    ar=[]
    for line in param_file:
        line=line.rstrip().split(' ')
        if line[0]=='fixedStep':
            chrom_name=line[1].split('=')[1]      
        if line[0] not in ['track', 'fixedStep']:
            ar.append(float(line[0]))
    param_file.close()
    print('Whole genome average ' + str(param_name) + ' : ' + str(sum(ar)/len(ar)))
    return ar, chrom_name


#######
#Read synonims and genes, encoding membrane proteins.
#######

def mem_prot_genes(pwd):
    mem_prot_data=open(pwd+'E_coli_W3110_based_table_of_synonyms_membrane_info.txt', 'r')
    
    mem_prot_genes_ar=[]
    for line in mem_prot_data:
        line=line.rstrip().split('\t')
        if line[1]!='-':
            mem_prot_genes_ar.append(line)
    
    mem_prot_data.close()
    return mem_prot_genes_ar


#######
#Split and write TF and promoters data, assign IGs with EcTopoI signal.
#######

def write_IG_data(IG_info_dict, pwd, EcTopoI_dict, mem_prot_genes_data):
    #Read EcTopoI signal.
    EcTopoI_data={}
    for name, path in EcTopoI_dict.items():
        wig_data, chr_id=score_data_parser(path, name)
        EcTopoI_data[name]=wig_data
    
    IG_output=open(pwd+'Intergenic_regions_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal.txt', 'w')
    IG_output.write('ID\tG1_name\tG1_start\tG1_end\tG1_strand\tG1_expression\tIG_length\tG2_name\tG2_start\tG2_end\tG2_strand\tG2_expression\t'+
                    'G1_membrane\tG2_membrane\t'+
                    'Promoters_names\tPromoters_sigmas\tPromoters_confidence\tPromoters_number\t'+
                    'TF_names\tTF_effects\tTF_confidence\tTF_number\t'+
                    'EcTopoI-Rif-CTD_FE\tEcTopoI+Rif-CTD_FE\tEcTopoI-Rif+CTD_FE\tEcTopoI+Rif+CTD_FE\tRNAP\n')
    #Split IG data and write.
    ID=0
    for IG_info, data in IG_info_dict.items():
        IG_info=IG_info.split('_')
        G1_name=IG_info[0]
        G2_name=IG_info[6]
        IG_start=int(IG_info[2])
        IG_end=int(IG_info[7])
        
        IG_output.write(f'{ID}\t')
        
        #Info about ajucent genes.
        for ele in IG_info:
            IG_output.write(f'{ele}\t')
        
        #Membrane genes.
        G1_mem, G2_mem=0, 0
        for cand_gene in mem_prot_genes_data:
            for synonim in cand_gene:
                if G1_name==synonim:
                    G1_mem=cand_gene[1]
                if G2_name==synonim:
                    G2_mem=cand_gene[1]
        if G1_mem==0:
            IG_output.write(f'-\t')
        else:
            IG_output.write(f'{G1_mem}\t')
        if G2_mem==0:
            IG_output.write(f'-\t')
        else:
            IG_output.write(f'{G2_mem}\t')
                
        #List promoter names.
        names_string=''
        for TSS in data['Promoter_info']:
            TSS=TSS.split('_')
            names_string=names_string+TSS[0]+';'
        names_string=names_string.rstrip(';')
        IG_output.write(f'{names_string}\t')
        
        #List promoter's sigma factor.
        sigma_string=''
        for TSS in data['Promoter_info']:
            TSS=TSS.split('_')
            sigma_string=sigma_string+TSS[2]+';'
        sigma_string=sigma_string.rstrip(';')
        IG_output.write(f'{sigma_string}\t')
        
        #List promoter's confidence.
        Confidence_string=''
        for TSS in data['Promoter_info']:
            TSS=TSS.split('_')
            Confidence_string=Confidence_string+TSS[3]+';'
        Confidence_string=Confidence_string.rstrip(';')    
        IG_output.write(f'{Confidence_string}\t')
        
        IG_output.write(f'{len(data["Promoter_info"])}\t')
        
        #List TF names.
        TF_names_string=''
        for TF in data['TF_info']:
            TF=TF.split('_')
            TF_names_string=TF_names_string+TF[0]+';'
        TF_names_string=TF_names_string.rstrip(';')
        IG_output.write(f'{TF_names_string}\t')        
            
        #List TF effects.
        TF_effect_string=''
        for TF in data['TF_info']:
            TF=TF.split('_')
            TF_effect_string=TF_effect_string+TF[2]+';'
        TF_effect_string=TF_effect_string.rstrip(';')
        IG_output.write(f'{TF_effect_string}\t')   
        
        #List TF confidence.
        TF_conf_string=''
        for TF in data['TF_info']:
            TF=TF.split('_')
            TF_conf_string=TF_conf_string+TF[3]+';'
        TF_conf_string=TF_conf_string.rstrip(';')
        IG_output.write(f'{TF_conf_string}\t')          
           
        IG_output.write(f'{len(data["TF_info"])}\t')
        
        #Assign EcTopoI and RNAP signal.
        for name, data in EcTopoI_data.items():
            Ec_signal=np.mean(data[IG_start:IG_end])
            IG_output.write(f'{Ec_signal}\t')
        IG_output.write(f'\n')
        
        ID+=1
   
    IG_output.close()
    return

#IG_info_dictionary=read_blast_output_combine(PWD)
#Membrane_proteins_data=mem_prot_genes(PWD)
#write_IG_data(IG_info_dictionary, PWD, EcTopoI_data_dict, Membrane_proteins_data)



##############
##Part 4.
##Add information from PSORT database (membrane-targeting proteins).
##############

#Source genome to extract sequences of genes from.
MG1655_genome_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Genomes\E_coli_K-12_MG1655_NC_000913.3.fasta"
#PSORT data.
PSORT_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_membrane_proteins\Databases\PSORT_3_database_E_coli.xlsx"
#W3110 genes.
W3110_genes_path_txt="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\Intergenic_regions_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal.txt"
W3110_genes_path_xlsx="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\Intergenic_regions_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal.xlsx"


#######
#Converts .txt file containing IGSs info (TAB) to .xlsx
#######

def convert_csv_to_xlsx(pathin, pathout):
    
    input_dataframe=pd.read_csv(pathin, sep='\t', header=0, index_col=False)
    input_dataframe.to_excel(pathout, header=True, index=False, sheet_name='Intergenic_regions_info')
    
    return


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


#######
#Read blast data, add to IGR info, return new extended .xlsx.
#######

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
        G1_start=int(Target_genome_data.loc[i, 'G1_start'])
        G1_end=int(Target_genome_data.loc[i, 'G1_end'])
        G2_start=int(Target_genome_data.loc[i, 'G2_start'])
        G2_end=int(Target_genome_data.loc[i, 'G2_end'])
        
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
    Target_genome_data.to_excel(pwd+'Intergenic_regions_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=True, index=True, sheet_name='Intergenic_regions_info')
    return

#convert_csv_to_xlsx(W3110_genes_path_txt, W3110_genes_path_xlsx)
#read_PSORT_blast(PSORT_data_path, MG1655_genome_path, PWD)
#read_blast_convert(PWD, W3110_genes_path_xlsx)



##############
##Part 5.
##Remove IGRs with anomalous signals (extremely high EcTopoI peaks near dps gene; or rRNA/tRNA genes with extremely high expression levels).
##############


#Genes to be removed.
Genes_to_remove={'G1_name' : ['rrsG', 'rrsA', 'rrsB', 'rrsC', 'rrsD', 'dps', 'ybhB', 'ybhC'],
                 'G2_name' : ['rrsE', 'rrsH']}

def remove_anomalous_regions(pwd, genes_to_remove_dict):
    
    #Read input dataframe.
    IGR_data=pd.read_excel(pwd+'Intergenic_regions_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index=0, sheet_name='Intergenic_regions_info')
    #Remove rows containing undesired gene names.
    IGR_data_filtered=IGR_data[(~IGR_data['G1_name'].isin(genes_to_remove_dict['G1_name'])) & (~IGR_data['G2_name'].isin(genes_to_remove_dict['G2_name']))]
    #Write a resultant dataframe.
    IGR_data_filtered.to_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=True, index=False, sheet_name='Intergenic_regions_info')
    return

#remove_anomalous_regions(PWD, Genes_to_remove)



##############
##Part 6.
##IGRs data analysis and visualization.
##############

PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\\"

#######
#Read final table, test gene orientation factor.
#######

def set_axis_style(ax, labels, pos):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, size=6)
    ax.set_xlim(0.25, max(pos)+0.75)
    return

def read_test_gene_orientation(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']
    
    #Classify gene pairs by orientation.
    Genes_plus_plus=input_data[(input_data['G1_strand']=='+') & (input_data['G2_strand']=='+')]
    Genes_plus_minus=input_data[(input_data['G1_strand']=='+') & (input_data['G2_strand']=='-')]
    Genes_minus_plus=input_data[(input_data['G1_strand']=='-') & (input_data['G2_strand']=='+')]
    Genes_minus_minus=input_data[(input_data['G1_strand']=='-') & (input_data['G2_strand']=='-')]
    
    print(Genes_plus_plus.shape, Genes_minus_minus.shape, Genes_plus_minus.shape, Genes_minus_plus.shape)

    #Plot EcTopoI enrichment.
    pos1=[1, 2, 3, 4, 7, 8, 9, 10]
    dataset1=[Genes_plus_plus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_minus_minus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_plus_minus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),  Genes_minus_plus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),
              Genes_plus_plus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_minus_minus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_plus_minus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(),  Genes_minus_plus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    pos2=[1, 2, 3, 4]
    dataset2=[Genes_plus_plus.loc[:, 'RNAP'].tolist(), Genes_minus_minus.loc[:, 'RNAP'].tolist(), Genes_plus_minus.loc[:, 'RNAP'].tolist(),  Genes_minus_plus.loc[:, 'RNAP'].tolist()]
    dataset3=[Genes_plus_plus.loc[:, 'Cumulative_expression'].tolist(), Genes_minus_minus.loc[:, 'Cumulative_expression'].tolist(), Genes_plus_minus.loc[:, 'Cumulative_expression'].tolist(),  Genes_minus_plus.loc[:, 'Cumulative_expression'].tolist()]
    
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(9,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i<4:
            violins['bodies'][i].set_facecolor('#ff7762')
        else:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['->->', '<-<-', '-><-', '<-->', '->->', '<-<-', '-><-', '<-->']
    set_axis_style(plt1, labels, pos1)
    
    yticknames1=np.arange(0.1, 40, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.1, 50)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f'   {Genes_plus_plus.shape[0]}',   xy=(0.5, 28), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'   {Genes_minus_minus.shape[0]}', xy=(1.5, 28), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_plus_minus.shape[0]}', xy=(2.5, 3), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'   {Genes_minus_plus.shape[0]}',  xy=(3.5, 29), xycoords='data', size=15, rotation=0) 
    
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.35), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.32), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(2.5, 0.45), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_minus_plus.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.31), xycoords='data', size=12, rotation=0)      
    
    plt1.annotate(f'   {Genes_plus_plus.shape[0]}',   xy=(6.5, 23), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'   {Genes_minus_minus.shape[0]}', xy=(7.5, 12.5), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_plus_minus.shape[0]}', xy=(8.5, 3.2), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'   {Genes_minus_plus.shape[0]}',  xy=(9.5, 11), xycoords='data', size=15, rotation=0) 
    
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(6.5, 0.45), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(7.5, 0.42), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(8.5, 0.55), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_minus_plus.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(9.5, 0.45), xycoords='data', size=12, rotation=0)     
    
    #Test EcTopoI FE difference between groups of IGRs.
    for i in range(4):
        for j in range(4):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j], dataset1[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i]),2)} Mean2={round(np.mean(dataset1[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    for i in range(4):
        for j in range(4):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j+4], dataset1[i+4], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i+4]),2)} Mean2={round(np.mean(dataset1[j+4]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')   
                
                
    #RNAP fold enrichment.
    #Draw violin-plots.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['->->', '<-<-', '-><-', '<-->']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.01, 60, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.01, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f'     {Genes_plus_plus.shape[0]}', xy=(0.5, 30), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'   {Genes_minus_minus.shape[0]}', xy=(1.5, 30), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'    {Genes_plus_minus.shape[0]}', xy=(2.5, 12), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'    {Genes_minus_plus.shape[0]}', xy=(3.5, 30), xycoords='data', size=15, rotation=0)
    
    plt2.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "RNAP"].tolist()),2)}', xy=(0.5, 0.055), xycoords='data', size=12, rotation=0)
    plt2.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"RNAP"].tolist()),2)}', xy=(1.5, 0.05), xycoords='data', size=12, rotation=0)
    plt2.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "RNAP"].tolist()),2)}', xy=(2.5, 0.25), xycoords='data', size=12, rotation=0)  
    plt2.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_plus.loc[:, "RNAP"].tolist()),2)}', xy=(3.5, 0.05), xycoords='data', size=12, rotation=0)      
    
    #Test RNAP FE difference between groups of IGRs.
    for i in range(4):
        for j in range(4):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset2[j], dataset2[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[i]),2)} Mean2={round(np.mean(dataset2[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    #Expression level.
    #Draw violin-plots.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['->->', '<-<-', '-><-', '<-->']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.0002, 8000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.0002, 8000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log') 
    
    plt3.annotate(f'     {Genes_plus_plus.shape[0]}', xy=(0.5, 1500), xycoords='data', size=15, rotation=0)
    plt3.annotate(f'   {Genes_minus_minus.shape[0]}', xy=(1.5, 2100), xycoords='data', size=15, rotation=0)
    plt3.annotate(f'    {Genes_plus_minus.shape[0]}', xy=(2.5, 30), xycoords='data', size=15, rotation=0)
    plt3.annotate(f'    {Genes_minus_plus.shape[0]}', xy=(3.5, 800), xycoords='data', size=15, rotation=0)

    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "Cumulative_expression"].tolist()),1)}', xy=(0.5, 0.0035), xycoords='data', size=12, rotation=0)
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.0025), xycoords='data', size=12, rotation=0) 
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(2.5, 0.12), xycoords='data', size=12, rotation=0)  
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_plus.loc[:, "Cumulative_expression"].tolist()),2)}', xy=(3.5, 0.0015), xycoords='data', size=12, rotation=0)   
        
    
    #Test Expression level difference between groups of IGRs.
    for i in range(4):
        for j in range(4):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset3[j], dataset3[i], equal_var=False)
                print(f'\nT-test Expression level means\nMean1={round(np.mean(dataset3[i]),2)} Mean2={round(np.mean(dataset3[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')                 
                
           
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_Genes_orientation_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.png', dpi=400, figsize=(4.2, 6.5))         
    
    return

#read_test_gene_orientation(PWD)


#######
#Read final table, test genes-encode-membrane-targeting-protein factor effect on EcTopo signal.
#Data from GO terms (Ecocyc).
#######

def read_test_membrane_orientation_GO_terms(pwd, DB_type):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']
    
    #Classify gene pairs by membranness of the proteins they encode.
    if DB_type=='Ecocyc':
        Genes_plus_plus=input_data[(input_data['G1_membrane']!='-') & (input_data['G2_membrane']!='-')]
        Genes_plus_minus=input_data[((input_data['G1_membrane']=='-') & (input_data['G2_membrane']!='-')) | ((input_data['G1_membrane']!='-') & (input_data['G2_membrane']=='-'))]
        Genes_minus_minus=input_data[(input_data['G1_membrane']=='-') & (input_data['G2_membrane']=='-')]
    elif DB_type=='PSORT':
        input_data['G1_PSORT_mem']=(input_data['G1_PSORT_localization']=='CytoplasmicMembrane') | (input_data['G1_PSORT_localization']=='OuterMembrane')
        input_data['G2_PSORT_mem']=(input_data['G2_PSORT_localization']=='CytoplasmicMembrane') | (input_data['G2_PSORT_localization']=='OuterMembrane')
    
        Genes_plus_plus=input_data[(input_data['G1_PSORT_mem']==True) & (input_data['G2_PSORT_mem']==True)]
        Genes_plus_minus=input_data[((input_data['G1_PSORT_mem']==True) & (input_data['G2_PSORT_mem']==False)) | ((input_data['G1_PSORT_mem']==False) & (input_data['G2_PSORT_mem']==True))]
        Genes_minus_minus=input_data[(input_data['G1_PSORT_mem']==False) & (input_data['G2_PSORT_mem']==False) ]       
    
    print(Genes_plus_plus.shape, Genes_minus_minus.shape, Genes_plus_minus.shape)

    #Plot EcTopoI enrichment.
    pos1=[1, 2, 3, 5, 6, 7]
    dataset1=[Genes_plus_plus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_plus_minus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_minus_minus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), 
              Genes_plus_plus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_plus_minus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_minus_minus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    pos2=[1, 2, 3]
    dataset2=[Genes_plus_plus.loc[:, 'RNAP'].tolist(), Genes_plus_minus.loc[:, 'RNAP'].tolist(), Genes_minus_minus.loc[:, 'RNAP'].tolist()]    
    dataset3=[Genes_plus_plus.loc[:, 'Cumulative_expression'].tolist(), Genes_plus_minus.loc[:, 'Cumulative_expression'].tolist(), Genes_minus_minus.loc[:, 'Cumulative_expression'].tolist()]
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(8,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i<3:
            violins['bodies'][i].set_facecolor('#ff7762')
        else:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['MM', 'M-', '--', 'MM', 'M-', '--']
    set_axis_style(plt1, labels, pos1)
    
    yticknames1=np.arange(0.1, 47, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.1, 47)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f'    {Genes_plus_plus.shape[0]}',   xy=(0.5, 30), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_plus_minus.shape[0]}',  xy=(1.5, 27), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_minus_minus.shape[0]}', xy=(2.5, 27), xycoords='data', size=15, rotation=0)
    
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.30), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.30), xycoords='data', size=12, rotation=0)  
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(2.5, 0.30), xycoords='data', size=12, rotation=0)
        
    plt1.annotate(f'    {Genes_plus_plus.shape[0]}',   xy=(4.5, 9.7), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_plus_minus.shape[0]}',  xy=(5.5, 23), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_minus_minus.shape[0]}', xy=(6.5, 13), xycoords='data', size=15, rotation=0)
    
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.45), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(5.5, 0.45), xycoords='data', size=12, rotation=0)   
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(6.5, 0.4), xycoords='data', size=12, rotation=0)
    
    #Test EcTopoI FE difference between groups of IGRs.
    for i in range(3):
        for j in range(3):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j], dataset1[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i]),2)} Mean2={round(np.mean(dataset1[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    for i in range(3):
        for j in range(3):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j+3], dataset1[i+3], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i+3]),2)} Mean2={round(np.mean(dataset1[j+3]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')   
                
    
    #RNAP fold enrichment.
    #Draw violin-plots. Promoter factor.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['MM', 'M-', '--']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.023, 60, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.023, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f'     {Genes_plus_plus.shape[0]}', xy=(0.5, 30), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'    {Genes_plus_minus.shape[0]}', xy=(1.5, 30), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'   {Genes_minus_minus.shape[0]}', xy=(2.5, 30), xycoords='data', size=15, rotation=0)
    
    plt2.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "RNAP"].tolist()),2)}', xy=(0.5, 0.08), xycoords='data', size=12, rotation=0)
    plt2.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "RNAP"].tolist()),2)}', xy=(1.5, 0.06), xycoords='data', size=12, rotation=0)  
    plt2.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"RNAP"].tolist()),2)}', xy=(2.5, 0.055), xycoords='data', size=12, rotation=0)    
    
    #Test RNAP FE difference between groups of IGRs.
    for i in range(3):
        for j in range(3):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset2[j], dataset2[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[i]),2)} Mean2={round(np.mean(dataset2[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    #Expression level.
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['MM', 'M-', '--']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.0003, 8000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.0003, 8000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log') 
    
    plt3.annotate(f'     {Genes_plus_plus.shape[0]}', xy=(0.5, 700), xycoords='data', size=15, rotation=0)
    plt3.annotate(f'    {Genes_plus_minus.shape[0]}', xy=(1.5, 1300), xycoords='data', size=15, rotation=0)
    plt3.annotate(f'   {Genes_minus_minus.shape[0]}', xy=(2.5, 2100), xycoords='data', size=15, rotation=0)
    
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "Cumulative_expression"].tolist()),1)}', xy=(0.5, 0.007), xycoords='data', size=12, rotation=0)
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.0065), xycoords='data', size=12, rotation=0)  
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"Cumulative_expression"].tolist()),1)}', xy=(2.5, 0.0015), xycoords='data', size=12, rotation=0)    
    
    #Test Expression level difference between groups of IGRs.
    for i in range(3):
        for j in range(3):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset3[j], dataset3[i], equal_var=False)
                print(f'\nT-test Expression level means\nMean1={round(np.mean(dataset3[i]),2)} Mean2={round(np.mean(dataset3[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')     
    
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_Membrane_proteins_{DB_type}_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.png', dpi=400, figsize=(4.8, 6.5))       
    
    return

#read_test_membrane_orientation_GO_terms(PWD, 'PSORT')



#######
#Read final table, test genes-encode-membrane-targeting-protein factor effect on EcTopo signal.
#Data from PSORT database.
#######

def read_test_membrane_orientation_PSORT(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']
    input_data['G1_PSORT_mem']=(input_data['G1_PSORT_localization']=='CytoplasmicMembrane') | (input_data['G1_PSORT_localization']=='OuterMembrane')
    input_data['G2_PSORT_mem']=(input_data['G2_PSORT_localization']=='CytoplasmicMembrane') | (input_data['G2_PSORT_localization']=='OuterMembrane')
    

    #Classify gene pairs by membranness of the proteins they encode.
    Genes_plus_plus=input_data[(((input_data['G1_PSORT_mem']==True) | (input_data['G2_PSORT_mem']==True)) & (input_data['G2_strand']=='+') & (input_data['G1_strand']=='-'))]
    
    Genes_plus_minus=input_data[(~((input_data['G1_PSORT_mem']==True) | (input_data['G2_PSORT_mem']==True)) & (input_data['G2_strand']=='+') & (input_data['G1_strand']=='-'))]
    
    Genes_minus_minus=input_data[(input_data['G2_strand']=='-') | (input_data['G1_strand']=='+') ]
    
    print(Genes_plus_plus.shape, Genes_minus_minus.shape, Genes_plus_minus.shape)

    #Plot EcTopoI enrichment.
    pos1=[1, 2, 3, 5, 6, 7]
    dataset1=[Genes_plus_plus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_plus_minus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_minus_minus.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), 
              Genes_plus_plus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_plus_minus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_minus_minus.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    pos2=[1, 2, 3]
    dataset2=[Genes_plus_plus.loc[:, 'RNAP'].tolist(), Genes_plus_minus.loc[:, 'RNAP'].tolist(), Genes_minus_minus.loc[:, 'RNAP'].tolist()]    
    dataset3=[Genes_plus_plus.loc[:, 'Cumulative_expression'].tolist(), Genes_plus_minus.loc[:, 'Cumulative_expression'].tolist(), Genes_minus_minus.loc[:, 'Cumulative_expression'].tolist()]
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(8,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i<3:
            violins['bodies'][i].set_facecolor('#ff7762')
        else:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['M/<-->', '~M/<-->', '~<-->', 'M/<-->', '~M/<-->', '~<-->']
    set_axis_style(plt1, labels, pos1)
    
    yticknames1=np.arange(0.3, 35, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.3, 35)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f'      {Genes_plus_plus.shape[0]}', xy=(0.5, 8), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'     {Genes_plus_minus.shape[0]}', xy=(1.5, 7), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_minus_minus.shape[0]}', xy=(2.5, 9), xycoords='data', size=15, rotation=0)
    
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.55), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.48), xycoords='data', size=12, rotation=0)  
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(2.5, 0.5), xycoords='data', size=12, rotation=0)
        
    plt1.annotate(f'      {Genes_plus_plus.shape[0]}', xy=(4.5, 10), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'     {Genes_plus_minus.shape[0]}', xy=(5.5, 12), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_minus_minus.shape[0]}', xy=(6.5, 25), xycoords='data', size=15, rotation=0)
    
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.57), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(5.5, 0.5), xycoords='data', size=12, rotation=0)   
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(6.5, 0.45), xycoords='data', size=12, rotation=0)
    
    #Test EcTopoI FE difference between groups of IGRs.
    for i in range(3):
        for j in range(3):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j], dataset1[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i]),2)} Mean2={round(np.mean(dataset1[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    for i in range(3):
        for j in range(3):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j+3], dataset1[i+3], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i+3]),2)} Mean2={round(np.mean(dataset1[j+3]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')   
                
    
    #RNAP fold enrichment.
    #Draw violin-plots. Promoter factor.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['M/<-->', '~M/<-->', '~<-->']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.023, 60, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.023, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f'     {Genes_plus_plus.shape[0]}', xy=(0.5, 20), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'    {Genes_plus_minus.shape[0]}', xy=(1.5, 30), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'   {Genes_minus_minus.shape[0]}', xy=(2.5, 30), xycoords='data', size=15, rotation=0)
    
    plt2.annotate(r"     $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:, "RNAP"].tolist()),2)}', xy=(0.5, 0.08), xycoords='data', size=12, rotation=0)
    plt2.annotate(r"    $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "RNAP"].tolist()),2)}', xy=(1.5, 0.05), xycoords='data', size=12, rotation=0)  
    plt2.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:, "RNAP"].tolist()),2)}', xy=(2.5, 0.055), xycoords='data', size=12, rotation=0)    
    
    #Test RNAP FE difference between groups of IGRs.
    for i in range(3):
        for j in range(3):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset2[j], dataset2[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[i]),2)} Mean2={round(np.mean(dataset2[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    #Expression level.
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['M/<-->', '~M/<-->', '~<-->']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.0003, 8000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.0003, 8000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log') 
    
    plt3.annotate(f'     {Genes_plus_plus.shape[0]}', xy=(0.5, 700), xycoords='data', size=15, rotation=0)
    plt3.annotate(f'    {Genes_plus_minus.shape[0]}', xy=(1.5, 574), xycoords='data', size=15, rotation=0)
    plt3.annotate(f'   {Genes_minus_minus.shape[0]}', xy=(2.5, 2100), xycoords='data', size=15, rotation=0)
    
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_plus.loc[:,  "Cumulative_expression"].tolist()),1)}', xy=(0.5, 0.007), xycoords='data', size=12, rotation=0)
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_plus_minus.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.0013), xycoords='data', size=12, rotation=0)  
    plt3.annotate(r"   $\overline{X}$"+f'={round(np.mean(Genes_minus_minus.loc[:,"Cumulative_expression"].tolist()),1)}', xy=(2.5, 0.0017), xycoords='data', size=12, rotation=0)    
    
    #Test Expression level difference between groups of IGRs.
    for i in range(3):
        for j in range(3):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset3[j], dataset3[i], equal_var=False)
                print(f'\nT-test Expression level means\nMean1={round(np.mean(dataset3[i]),2)} Mean2={round(np.mean(dataset3[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')     
    
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_Membrane_proteins_and_EcTopoI_enrichment_in_IG_regions_no_rRNA_PSORT.png', dpi=400, figsize=(4.8, 6.5))       
    
    return

#read_test_membrane_orientation_PSORT(PWD)



#######
#Read final table, test IGRs complexity hypothesis (if IGR contain promoter or TF sites).
#######

def read_test_IGR_complexity(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']
    
    #Classify gene pairs by promoters and TFs.
    Genes_with_promoter=input_data[(input_data['Promoters_number']!=0)]
    Genes_no_promoter=input_data[(input_data['Promoters_number']==0)]
    Genes_with_TFs=input_data[(input_data['TF_number']!=0)]
    Genes_no_TFs=input_data[(input_data['TF_number']==0)]
    
    print(Genes_with_promoter.shape, Genes_no_promoter.shape, Genes_with_TFs.shape, Genes_no_TFs.shape)
    
    #Plot EcTopoI enrichment.
    pos=[1, 2, 4, 5]
    pos2=[1, 2]
    
    dataset1=[Genes_with_promoter.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_no_promoter.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),Genes_with_promoter.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_no_promoter.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    dataset2=[Genes_with_TFs.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_no_TFs.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_with_TFs.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_no_TFs.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    dataset3=[Genes_with_promoter.loc[:, 'RNAP'].tolist(), Genes_no_promoter.loc[:, 'RNAP'].tolist()]   
    dataset4=[Genes_with_promoter.loc[:, 'Cumulative_expression'].tolist(), Genes_no_promoter.loc[:, 'Cumulative_expression'].tolist()]    
    dataset5=[Genes_with_TFs.loc[:, 'RNAP'].tolist(), Genes_no_TFs.loc[:, 'RNAP'].tolist()]   
    dataset6=[Genes_with_TFs.loc[:, 'Cumulative_expression'].tolist(), Genes_no_TFs.loc[:, 'Cumulative_expression'].tolist()]        
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(5,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i<2:
            violins['bodies'][i].set_facecolor('#ff7762')
        else:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['+P', '-P', '+P', '-P']
    set_axis_style(plt1, labels, pos)    
    
    yticknames1=np.arange(0.2, 45, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.2, 45)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f' {Genes_with_promoter.shape[0]}',  xy=(0.5, 30), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_no_promoter.shape[0]}', xy=(1.5, 1.1), xycoords='data', size=15, rotation=0)

    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_with_promoter.loc[:,"EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.35), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_no_promoter.loc[:,  "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.53), xycoords='data', size=12, rotation=0)  
   
    plt1.annotate(f' {Genes_with_promoter.shape[0]}',  xy=(3.5, 23), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'    {Genes_no_promoter.shape[0]}', xy=(4.5, 1.3), xycoords='data', size=15, rotation=0)

    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_with_promoter.loc[:,"EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.45), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_no_promoter.loc[:,  "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.58), xycoords='data', size=12, rotation=0)   
    
    #Test EcTopoI FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[0], dataset1[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[0]),2)} Mean2={round(np.mean(dataset1[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[2], dataset1[3], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[2]),2)} Mean2={round(np.mean(dataset1[3]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    #RNAP fold enrichment.
    #Draw violin-plots.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['+P', '-P']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.02, 80, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.02, 80)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f'{Genes_with_promoter.shape[0]}',  xy=(0.5, 35), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'   {Genes_no_promoter.shape[0]}', xy=(1.5, 25), xycoords='data', size=15, rotation=0)

    plt2.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_with_promoter.loc[:,"RNAP"].tolist()),2)}', xy=(0.5, 0.055), xycoords='data', size=12, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_no_promoter.loc[:, "RNAP"].tolist()),2)}', xy=(1.5, 0.08), xycoords='data', size=12, rotation=0)   
    
    #Test RNAP FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset3[0], dataset3[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset3[0]),2)} Mean2={round(np.mean(dataset3[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    #Expression level.
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset4, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['+P', '-P']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.0001, 11000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.0001, 11000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log')
    
    plt3.annotate(f'{Genes_with_promoter.shape[0]}',  xy=(0.5, 3000), xycoords='data', size=15, rotation=0)
    plt3.annotate(f'   {Genes_no_promoter.shape[0]}', xy=(1.5, 20), xycoords='data', size=15, rotation=0)

    plt3.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_with_promoter.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(0.5, 0.0015), xycoords='data', size=12, rotation=0)
    plt3.annotate(r"  $\overline{X}$"+f'={round(np.mean(Genes_no_promoter.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.002), xycoords='data', size=12, rotation=0)   
    
    #Test Expression level difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset4[0], dataset4[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset4[0]),2)} Mean2={round(np.mean(dataset4[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')         
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_Promoter_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.png', dpi=400, figsize=(6.8, 6.5))       
    
    
    #EcTopoI fold enrichment.    
    #Draw violin-plots.
    fig1=plt.figure(figsize=(5,6.5), dpi=100)
    plt1=fig1.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset2, positions=pos, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i<2:
            violins['bodies'][i].set_facecolor('#ff7762')
        else:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['+TF', '-TF', '+TF', '-TF']
    set_axis_style(plt1, labels, pos)    
    
    yticknames1=np.arange(0.2, 45, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.2, 45)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')            
             
    plt1.annotate(f'  {Genes_with_TFs.shape[0]}', xy=(0.5, 30), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {Genes_no_TFs.shape[0]}',   xy=(1.5, 20), xycoords='data', size=15, rotation=0)             
    
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_with_TFs.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.35), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_no_TFs.loc[:,   "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.35), xycoords='data', size=12, rotation=0)               
      
    plt1.annotate(f'  {Genes_with_TFs.shape[0]}', xy=(3.5, 14), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {Genes_no_TFs.shape[0]}',   xy=(4.5, 25), xycoords='data', size=15, rotation=0)   
    
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_with_TFs.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.52), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_no_TFs.loc[:,   "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.40), xycoords='data', size=12, rotation=0)    
    
    #Test EcTopoI FE difference between groups of IGRs.   
    Intervals_stat=stats.ttest_ind(dataset2[0], dataset2[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[0]),2)} Mean2={round(np.mean(dataset2[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset2[2], dataset2[3], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[2]),2)} Mean2={round(np.mean(dataset2[3]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    #RNAP fold enrichment.
    #Draw violin-plots.
    plt2=fig1.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset5, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['+TF', '-TF']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.02, 80, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.02, 80)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f' {Genes_with_TFs.shape[0]}', xy=(0.5, 35), xycoords='data', size=15, rotation=0)
    plt2.annotate(f' {Genes_no_TFs.shape[0]}',   xy=(1.5, 36), xycoords='data', size=15, rotation=0)

    plt2.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_with_TFs.loc[:, "RNAP"].tolist()),2)}', xy=(0.5, 0.07), xycoords='data', size=12, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_no_TFs.loc[:,  "RNAP"].tolist()),2)}', xy=(1.5, 0.05), xycoords='data', size=12, rotation=0)   
    
    #Test RNAP FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset5[0], dataset5[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset5[0]),2)} Mean2={round(np.mean(dataset5[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    #Expression level.
    #Draw violin-plots. Promoter factor.
    plt3=fig1.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset6, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['+TF', '-TF']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.0001, 11000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.0001, 11000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log')
    
    plt3.annotate(f' {Genes_with_TFs.shape[0]}', xy=(0.5, 3000), xycoords='data', size=15, rotation=0)
    plt3.annotate(f' {Genes_no_TFs.shape[0]}',   xy=(1.5, 2000), xycoords='data', size=15, rotation=0)

    plt3.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_with_TFs.loc[:, "Cumulative_expression"].tolist()),0)}', xy=(0.5, 0.0013), xycoords='data', size=12, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_no_TFs.loc[:,  "Cumulative_expression"].tolist()),0)}', xy=(1.5, 0.003), xycoords='data', size=12, rotation=0)   
    
    #Test Expression level difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset6[0], dataset6[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset6[0]),2)} Mean2={round(np.mean(dataset6[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')        
           
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_TFs_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.png', dpi=400, figsize=(6.8, 6.5)) 
    
    return

#read_test_IGR_complexity(PWD)



#######
#Read final table, test effects of RNAP FE on EcTopoI signal.
#######

def read_test_RNAP_effect(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']    
    
    #Classify gene pairs by RNAP or Expression level.
    Genes_high_RNAP=input_data[input_data['RNAP']>2]
    Genes_low_RNAP=input_data[input_data['RNAP']<2]
    
    print(Genes_high_RNAP.shape, Genes_low_RNAP.shape)


    #Plot EcTopoI enrichment.
    pos1=[1, 2, 4, 5]
    dataset1=[Genes_high_RNAP.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_low_RNAP.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), 
              Genes_high_RNAP.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_low_RNAP.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(5,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i in [0,1]:
            violins['bodies'][i].set_facecolor('#ff7762')
        elif i in [2,3]:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
        
        #Taken from here: https://stackoverflow.com/questions/29776114/half-violin-plot/29781988#29781988
        #m=np.mean(violins['bodies'][i].get_paths()[0].vertices[:, 0])
        #violins['bodies'][i].get_paths()[0].vertices[:, 0]=np.clip(violins['bodies'][i].get_paths()[0].vertices[:, 0], -np.inf, m)
        #violins['bodies'][i].set_color('r')
        
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['RNAP>2', 'RNAP<2', 'RNAP>2', 'RNAP<2']
    set_axis_style(plt1, labels, pos1)    
    
    yticknames1=np.arange(0.2, 45, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.2, 45)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f'  {Genes_high_RNAP.shape[0]}', xy=(0.5, 31), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {Genes_low_RNAP.shape[0]}',  xy=(1.5, 29), xycoords='data', size=15, rotation=0)

    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_high_RNAP.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.35), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_low_RNAP.loc[:,  "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.32), xycoords='data', size=12, rotation=0)  
   
    plt1.annotate(f'  {Genes_high_RNAP.shape[0]}', xy=(3.5, 14), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {Genes_low_RNAP.shape[0]}',  xy=(4.5, 27), xycoords='data', size=15, rotation=0)

    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_high_RNAP.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.46), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_low_RNAP.loc[:,  "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.42), xycoords='data', size=12, rotation=0)   
    
    #Test EcTopoI FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[0], dataset1[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[0]),2)} Mean2={round(np.mean(dataset1[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[2], dataset1[3], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[2]),2)} Mean2={round(np.mean(dataset1[3]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    
    #RNAP fold enrichment.
    pos2=[1, 2]
    dataset2=[Genes_high_RNAP.loc[:, 'RNAP'].tolist(), Genes_low_RNAP.loc[:, 'RNAP'].tolist()]
    
    #Draw violin-plots. Promoter factor.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['RNAP>2', 'RNAP<2']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.05, 60, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.05, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f' {Genes_high_RNAP.shape[0]}', xy=(0.5, 34), xycoords='data', size=15, rotation=0)
    plt2.annotate(f' {Genes_low_RNAP.shape[0]}',  xy=(1.5, 2.42), xycoords='data', size=15, rotation=0)

    plt2.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_high_RNAP.loc[:, "RNAP"].tolist()),2)}', xy=(0.5, 1.12), xycoords='data', size=12, rotation=0)
    plt2.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_low_RNAP.loc[:,      "RNAP"].tolist()),2)}', xy=(1.5, 0.06), xycoords='data', size=12, rotation=0)   
    
    #Test RNAP FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset2[0], dataset2[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[0]),2)} Mean2={round(np.mean(dataset2[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    #Expression level.
    dataset3=[Genes_high_RNAP.loc[:, 'Cumulative_expression'].tolist(), Genes_low_RNAP.loc[:, 'Cumulative_expression'].tolist()]
    
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['RNAP>2', 'RNAP<2']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.001, 8000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.001, 8000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log')
    
    plt3.annotate(f' {Genes_high_RNAP.shape[0]}', xy=(0.5, 2586), xycoords='data', size=15, rotation=0)
    plt3.annotate(f' {Genes_low_RNAP.shape[0]}',  xy=(1.5, 1020), xycoords='data', size=15, rotation=0)

    plt3.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_high_RNAP.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(0.5, 0.007), xycoords='data', size=12, rotation=0)
    plt3.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_low_RNAP.loc[:,  "Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.0018), xycoords='data', size=12, rotation=0)   
    
    #Test Expression level difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset3[0], dataset3[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset3[0]),2)} Mean2={round(np.mean(dataset3[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')        
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_RNAP_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.png', dpi=400, figsize=(5, 6.5))       
    
    return

#read_test_RNAP_effect(PWD)



#######
#Read final table, test effects of expression level of adjacent genes on EcTopoI signal.
#######

def read_test_expression_level(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']    
    
    #Classify gene pairs by RNAP or Expression level.
    Genes_high_EL=input_data[input_data['Cumulative_expression']>5]
    Genes_low_EL=input_data[input_data['Cumulative_expression']<5]
    
    print(Genes_high_EL.shape, Genes_low_EL.shape)


    #Plot EcTopoI enrichment.
    pos1=[1, 2, 4, 5]
    dataset1=[Genes_high_EL.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_low_EL.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), 
              Genes_high_EL.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_low_EL.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(5,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i in [0,1]:
            violins['bodies'][i].set_facecolor('#ff7762')
        elif i in [2,3]:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
        
        #Taken from here: https://stackoverflow.com/questions/29776114/half-violin-plot/29781988#29781988
        #m=np.mean(violins['bodies'][i].get_paths()[0].vertices[:, 0])
        #violins['bodies'][i].get_paths()[0].vertices[:, 0]=np.clip(violins['bodies'][i].get_paths()[0].vertices[:, 0], -np.inf, m)
        #violins['bodies'][i].set_color('r')
        
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['EL>5', 'EL<5', 'EL>5', 'EL<5']
    set_axis_style(plt1, labels, pos1)    
    
    yticknames1=np.arange(0.1, 50, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.1, 50)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f'  {Genes_high_EL.shape[0]}', xy=(0.5, 32), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {Genes_low_EL.shape[0]}',  xy=(1.5, 13.5), xycoords='data', size=15, rotation=0)

    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_high_EL.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.37), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_low_EL.loc[:,  "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.27), xycoords='data', size=12, rotation=0)  
   
    plt1.annotate(f'  {Genes_high_EL.shape[0]}', xy=(3.5, 25), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {Genes_low_EL.shape[0]}',  xy=(4.5, 14), xycoords='data', size=15, rotation=0)

    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_high_EL.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.46), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_low_EL.loc[:,  "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.42), xycoords='data', size=12, rotation=0)   
    
    #Test EcTopoI FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[0], dataset1[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[0]),2)} Mean2={round(np.mean(dataset1[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[2], dataset1[3], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[2]),2)} Mean2={round(np.mean(dataset1[3]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    
    #RNAP fold enrichment.
    pos2=[1, 2]
    dataset2=[Genes_high_EL.loc[:, 'RNAP'].tolist(), Genes_low_EL.loc[:, 'RNAP'].tolist()]
    
    #Draw violin-plots. Promoter factor.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['EL>5', 'EL<5']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.05, 60, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.05, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f' {Genes_high_EL.shape[0]}', xy=(0.5, 34), xycoords='data', size=15, rotation=0)
    plt2.annotate(f' {Genes_low_EL.shape[0]}',  xy=(1.5, 32), xycoords='data', size=15, rotation=0)

    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_high_EL.loc[:, "RNAP"].tolist()),2)}', xy=(0.5, 0.12), xycoords='data', size=12, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_low_EL.loc[:,  "RNAP"].tolist()),2)}', xy=(1.5, 0.06), xycoords='data', size=12, rotation=0)   
    
    #Test RNAP FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset2[0], dataset2[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[0]),2)} Mean2={round(np.mean(dataset2[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    #Expression level.
    dataset3=[Genes_high_EL.loc[:, 'Cumulative_expression'].tolist(), Genes_low_EL.loc[:, 'Cumulative_expression'].tolist()]
    
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['EL>5', 'EL<5']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.001, 8000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.001, 8000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log')
    
    plt3.annotate(f' {Genes_high_EL.shape[0]}', xy=(0.5, 2586), xycoords='data', size=15, rotation=0)
    plt3.annotate(f' {Genes_low_EL.shape[0]}',  xy=(1.5, 7), xycoords='data', size=15, rotation=0)

    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_high_EL.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(0.5, 1.18), xycoords='data', size=12, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_low_EL.loc[:,  "Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.0018), xycoords='data', size=12, rotation=0)   
    
    #Test Expression level difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset3[0], dataset3[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset3[0]),2)} Mean2={round(np.mean(dataset3[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')        
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_Expression_level_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.png', dpi=400, figsize=(5, 6.5))       
    
    return

#read_test_expression_level(PWD)



#######
#Read final table, test expression factors combined together (RNAP and expression level),
#test together factors not related to expression directly (gene orientation, membrane-targeting proteins, transcription factors).
#######

def read_test_factors_combination(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']
    input_data['G1_PSORT_mem']=(input_data['G1_PSORT_localization']=='CytoplasmicMembrane') | (input_data['G1_PSORT_localization']=='OuterMembrane')
    input_data['G2_PSORT_mem']=(input_data['G2_PSORT_localization']=='CytoplasmicMembrane') | (input_data['G2_PSORT_localization']=='OuterMembrane')    
    
    #Classify gene pairs by non-expression factors(gene orientation, membrane-targeting proteins, transcription factors). 
    Genes_divergent_membrane_with_TF=input_data[(input_data['G1_strand']=='-') & (input_data['G2_strand']=='+') & (input_data['TF_number']>0) & ((input_data['G1_PSORT_mem']==True) | (input_data['G2_PSORT_mem']==True))]
    Genes_not_div_not_mem_no_TF=input_data[~((input_data['G1_strand']=='-') & (input_data['G2_strand']=='+')) & (input_data['TF_number']==0) & ~((input_data['G1_PSORT_mem']==True) | (input_data['G2_PSORT_mem']==True))]
    
    print(Genes_divergent_membrane_with_TF.shape, Genes_not_div_not_mem_no_TF.shape)

    #Plot EcTopoI enrichment.
    pos1=[1, 2, 4, 5]
    dataset1=[Genes_divergent_membrane_with_TF.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_not_div_not_mem_no_TF.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), 
              Genes_divergent_membrane_with_TF.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_not_div_not_mem_no_TF.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(5,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i in [0,1]:
            violins['bodies'][i].set_facecolor('#ff7762')
        elif i in [2,3]:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
        
        #Taken from here: https://stackoverflow.com/questions/29776114/half-violin-plot/29781988#29781988
        #m=np.mean(violins['bodies'][i].get_paths()[0].vertices[:, 0])
        #violins['bodies'][i].get_paths()[0].vertices[:, 0]=np.clip(violins['bodies'][i].get_paths()[0].vertices[:, 0], -np.inf, m)
        #violins['bodies'][i].set_color('r')
        
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['Div\n+TF\nMem', '-Div\n-TF\n-Mem', 'Div\n+TF\nMem', '-Div\n-TF\n-Mem']
    set_axis_style(plt1, labels, pos1)    
    
    yticknames1=np.arange(0.2, 45, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.2, 45)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f' {Genes_divergent_membrane_with_TF.shape[0]}', xy=(0.5, 30), xycoords='data', size=15, rotation=0)
    plt1.annotate(f' {Genes_not_div_not_mem_no_TF.shape[0]}',      xy=(1.5, 17), xycoords='data', size=15, rotation=0)

    plt1.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_divergent_membrane_with_TF.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.49), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_not_div_not_mem_no_TF.loc[:,     "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.3), xycoords='data', size=12, rotation=0)  
   
    plt1.annotate(f' {Genes_divergent_membrane_with_TF.shape[0]}', xy=(3.5, 12), xycoords='data', size=15, rotation=0)
    plt1.annotate(f' {Genes_not_div_not_mem_no_TF.shape[0]}',      xy=(4.5, 14.5), xycoords='data', size=15, rotation=0)

    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_divergent_membrane_with_TF.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.49), xycoords='data', size=12, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_not_div_not_mem_no_TF.loc[:,      "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.42), xycoords='data', size=12, rotation=0)   
    
    #Test EcTopoI FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[0], dataset1[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[0]),2)} Mean2={round(np.mean(dataset1[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[2], dataset1[3], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[2]),2)} Mean2={round(np.mean(dataset1[3]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    
    #RNAP fold enrichment.
    pos2=[1, 2]
    dataset2=[Genes_divergent_membrane_with_TF.loc[:, 'RNAP'].tolist(), Genes_not_div_not_mem_no_TF.loc[:, 'RNAP'].tolist()]
    
    #Draw violin-plots.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['Div\n+TF\nMem', '-Div\n-TF\n-Mem']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.05, 60, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.05, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f' {Genes_divergent_membrane_with_TF.shape[0]}', xy=(0.5, 18), xycoords='data', size=15, rotation=0)
    plt2.annotate(f' {Genes_not_div_not_mem_no_TF.shape[0]}',      xy=(1.5, 30.5), xycoords='data', size=15, rotation=0)

    plt2.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_divergent_membrane_with_TF.loc[:, "RNAP"].tolist()),2)}', xy=(0.5, 0.12), xycoords='data', size=12, rotation=0)
    plt2.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_not_div_not_mem_no_TF.loc[:,      "RNAP"].tolist()),2)}', xy=(1.5, 0.059), xycoords='data', size=12, rotation=0)   
    
    #Test RNAP FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset2[0], dataset2[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[0]),2)} Mean2={round(np.mean(dataset2[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    #Expression level.
    dataset3=[Genes_divergent_membrane_with_TF.loc[:, 'Cumulative_expression'].tolist(), Genes_not_div_not_mem_no_TF.loc[:, 'Cumulative_expression'].tolist()]
    
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['Div\n+TF\nMem', '-Div\n-TF\n-Mem']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.001, 8000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.001, 8000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log')
    
    plt3.annotate(f' {Genes_divergent_membrane_with_TF.shape[0]}', xy=(0.5, 800), xycoords='data', size=15, rotation=0)
    plt3.annotate(f' {Genes_not_div_not_mem_no_TF.shape[0]}',      xy=(1.5, 2000), xycoords='data', size=15, rotation=0)

    plt3.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_divergent_membrane_with_TF.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(0.5, 0.009), xycoords='data', size=12, rotation=0)
    plt3.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_not_div_not_mem_no_TF.loc[:,      "Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.004), xycoords='data', size=12, rotation=0)   
    
    #Test Expression level difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset3[0], dataset3[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset3[0]),2)} Mean2={round(np.mean(dataset3[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')        
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\Add_factors\\NEW_non_expression_factors_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.png', dpi=400, figsize=(5, 6.5))       
    
    
    
    #Classify gene pairs by expression-related factors (RNAP, expression level). 
    Genes_expression_high=input_data[(input_data['RNAP']>2) & (input_data['Cumulative_expression']>5)]
    Genes_expression_low=input_data[(input_data['RNAP']<2) & (input_data['Cumulative_expression']<5)]
    
    print(Genes_expression_high.shape, Genes_expression_low.shape)


    #Plot EcTopoI enrichment.
    pos1=[1, 2, 4, 5]
    dataset1=[Genes_expression_high.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_expression_low.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), 
              Genes_expression_high.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_expression_low.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(5,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i in [0,1]:
            violins['bodies'][i].set_facecolor('#ff7762')
        elif i in [2,3]:
            violins['bodies'][i].set_facecolor('#58ffb1')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
        
        #Taken from here: https://stackoverflow.com/questions/29776114/half-violin-plot/29781988#29781988
        #m=np.mean(violins['bodies'][i].get_paths()[0].vertices[:, 0])
        #violins['bodies'][i].get_paths()[0].vertices[:, 0]=np.clip(violins['bodies'][i].get_paths()[0].vertices[:, 0], -np.inf, m)
        #violins['bodies'][i].set_color('r')
        
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['EL>2\nRNAP>2', 'EL<2\nRNAP<2', 'EL>2\nRNAP>2', 'EL<2\nRNAP<2']
    set_axis_style(plt1, labels, pos1)    
    
    yticknames1=np.arange(0.2, 45, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.2, 45)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f' {Genes_expression_high.shape[0]}', xy=(0.5, 30), xycoords='data', size=15, rotation=0)
    plt1.annotate(f' {Genes_expression_low.shape[0]}',  xy=(1.5, 13), xycoords='data', size=15, rotation=0)

    plt1.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_expression_high.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.35), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_expression_low.loc[:,  "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.32), xycoords='data', size=12, rotation=0)  
   
    plt1.annotate(f' {Genes_expression_high.shape[0]}', xy=(3.5, 15), xycoords='data', size=15, rotation=0)
    plt1.annotate(f' {Genes_expression_low.shape[0]}',  xy=(4.5, 15), xycoords='data', size=15, rotation=0)

    plt1.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_expression_high.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.43), xycoords='data', size=12, rotation=0)
    plt1.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_expression_low.loc[:,  "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.42), xycoords='data', size=12, rotation=0)   
    
    #Test EcTopoI FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[0], dataset1[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[0]),2)} Mean2={round(np.mean(dataset1[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset1[2], dataset1[3], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[2]),2)} Mean2={round(np.mean(dataset1[3]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    
    #RNAP fold enrichment.
    pos2=[1, 2]
    dataset2=[Genes_expression_high.loc[:, 'RNAP'].tolist(), Genes_expression_low.loc[:, 'RNAP'].tolist()]
    
    #Draw violin-plots.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['EL>2\nRNAP>2', 'EL<2\nRNAP<2']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.05, 60, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.05, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f' {Genes_expression_high.shape[0]}', xy=(0.5, 35), xycoords='data', size=15, rotation=0)
    plt2.annotate(f' {Genes_expression_low.shape[0]}',  xy=(1.5, 2.5), xycoords='data', size=15, rotation=0)

    plt2.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_expression_high.loc[:, "RNAP"].tolist()),2)}', xy=(0.5, 1), xycoords='data', size=12, rotation=0)
    plt2.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_expression_low.loc[:,  "RNAP"].tolist()),2)}', xy=(1.5, 0.06), xycoords='data', size=12, rotation=0)   
    
    #Test RNAP FE difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset2[0], dataset2[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[0]),2)} Mean2={round(np.mean(dataset2[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    #Expression level.
    dataset3=[Genes_expression_high.loc[:, 'Cumulative_expression'].tolist(), Genes_expression_low.loc[:, 'Cumulative_expression'].tolist()]
    
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['EL>2\nRNAP>2', 'EL<2\nRNAP<2']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.001, 8000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.001, 8000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log')
    
    plt3.annotate(f' {Genes_expression_high.shape[0]}', xy=(0.5, 2100), xycoords='data', size=15, rotation=0)
    plt3.annotate(f' {Genes_expression_low.shape[0]}',  xy=(1.5, 8.63), xycoords='data', size=15, rotation=0)

    plt3.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_expression_high.loc[:, "Cumulative_expression"].tolist()),1)}', xy=(0.5, 1.09), xycoords='data', size=12, rotation=0)
    plt3.annotate(r"$\overline{X}$"+f'={round(np.mean(Genes_expression_low.loc[:,  "Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.0014), xycoords='data', size=12, rotation=0)   
    
    #Test Expression level difference between groups of IGRs.
    Intervals_stat=stats.ttest_ind(dataset3[0], dataset3[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset3[0]),2)} Mean2={round(np.mean(dataset3[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')        
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\Add_factors\\NEW_expression_factors_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.png', dpi=400, figsize=(5, 6.5))       
        
    return

#read_test_factors_combination(PWD)


#######
#Read final table, all factors together on one plot.
#######

def read_test_factors_combination_one_plot(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']
    input_data['G1_PSORT_mem']=(input_data['G1_PSORT_localization']=='CytoplasmicMembrane') | (input_data['G1_PSORT_localization']=='OuterMembrane')
    input_data['G2_PSORT_mem']=(input_data['G2_PSORT_localization']=='CytoplasmicMembrane') | (input_data['G2_PSORT_localization']=='OuterMembrane')    
    
    #Classify gene pairs by promoters and TFs.
    All=input_data
    EL_pos=input_data[(input_data['Cumulative_expression']>5)]
    RNAP_pos=input_data[(input_data['RNAP']>2)]
    EL_RNAP_pos=input_data[(input_data['RNAP']>2) & (input_data['Cumulative_expression']>5)]
    EL_RNAP_pos_div=input_data[(input_data['RNAP']>2) & (input_data['Cumulative_expression']>5) & (input_data['G1_strand']=='-') & (input_data['G2_strand']=='+')]
    input_data['EL_RNAP_high_div']=(input_data['RNAP']>2) & (input_data['Cumulative_expression']>5) & (input_data['G1_strand']=='-') & (input_data['G2_strand']=='+')
    EL_RNAP_pos_div_TF_pos=input_data[(input_data['RNAP']>2) & (input_data['Cumulative_expression']>5) & (input_data['G1_strand']=='-') & (input_data['G2_strand']=='+') & (input_data['TF_number']>0)]
    input_data['EL_RNAP_high_div_TF']=(input_data['RNAP']>2) & (input_data['Cumulative_expression']>5) & (input_data['G1_strand']=='-') & (input_data['G2_strand']=='+') & (input_data['TF_number']>0)
    EL_RNAP_pos_div_TF_pos_Mem=input_data[(input_data['RNAP']>2) & (input_data['Cumulative_expression']>5) & (input_data['G1_strand']=='-') & (input_data['G2_strand']=='+') & (input_data['TF_number']>0) & ((input_data['G1_PSORT_mem']==True) | (input_data['G2_PSORT_mem']==True))]
    input_data['EL_RNAP_high_div_TF_Mem']=(input_data['RNAP']>2) & (input_data['Cumulative_expression']>5) & (input_data['G1_strand']=='-') & (input_data['G2_strand']=='+') & (input_data['TF_number']>0) & ((input_data['G1_PSORT_mem']==True) | (input_data['G2_PSORT_mem']==True))
    
    EL_neg=input_data[(input_data['Cumulative_expression']<5)]
    RNAP_neg=input_data[(input_data['RNAP']<2)]
    EL_RNAP_neg=input_data[(input_data['RNAP']<2) & (input_data['Cumulative_expression']<5)]
    EL_RNAP_neg_not_div=input_data[(input_data['RNAP']<2) & (input_data['Cumulative_expression']<5) & ~((input_data['G1_strand']=='-') & (input_data['G2_strand']=='+'))]
    EL_RNAP_neg_not_div_TF_neg=input_data[(input_data['RNAP']<2) & (input_data['Cumulative_expression']<5) & ~((input_data['G1_strand']=='-') & (input_data['G2_strand']=='+')) & (input_data['TF_number']==0)]
    EL_RNAP_neg_not_div_TF_neg_not_Mem=input_data[(input_data['RNAP']<2) & (input_data['Cumulative_expression']<5) & ~((input_data['G1_strand']=='-') & (input_data['G2_strand']=='+')) & (input_data['TF_number']==0) & ~((input_data['G1_PSORT_mem']==True) | (input_data['G2_PSORT_mem']==True))]

    #Save new dataframe.
    input_data.to_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT_processed_1.xlsx', sheet_name='Intergenic_regions_info')
    
    
    #Plot EcTopoI enrichment. Positive factors.
    pos0=[1]
    dataset0=[All.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist()]
    
    #Draw violin-plots.
    fig=plt.figure(figsize=(7,4.5), dpi=100)
    plt0=plt.subplot2grid((2,4),(0,0), rowspan=2) 
    violins=plt0.violinplot(dataset0, positions=pos0, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff6863')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
        
        #Taken from here: https://stackoverflow.com/questions/29776114/half-violin-plot/29781988#29781988
        #m=np.mean(violins['bodies'][i].get_paths()[0].vertices[:, 0])
        #violins['bodies'][i].get_paths()[0].vertices[:, 0]=np.clip(violins['bodies'][i].get_paths()[0].vertices[:, 0], -np.inf, m)
        #violins['bodies'][i].set_color('r')
        
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=[0]
    set_axis_style(plt0, labels, pos0)    
    
    yticknames1=np.arange(0.25, 45, 5)
    plt0.set_yticks(yticknames1, minor=False)
    plt0.set_yticklabels(yticknames1)
    #plt0.set_ylabel('EcTopoI fold enrichment', size=15)
    plt0.set_ylim(0.25, 45)
    plt.setp(plt0.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt0.set_yscale('log')
    
    plt0.spines['top'].set_visible(False)
    plt0.spines['right'].set_visible(False)      
    
    plt0.annotate(f'   {All.shape[0]}', xy=(0.5, 30), xycoords='data', size=15, rotation=0)   

    plt0.annotate(r"    $\overline{X}$"+f'={round(np.mean(All.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.34), xycoords='data', size=12, rotation=0)     
    
    

    #Plot EcTopoI enrichment. Positive factors.
    pos1=[1, 2, 3, 4, 5, 6]
    dataset1=[EL_pos.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), RNAP_pos.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), EL_RNAP_pos.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), EL_RNAP_pos_div.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),
              EL_RNAP_pos_div_TF_pos.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), EL_RNAP_pos_div_TF_pos_Mem.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist()]
    
    #Draw violin-plots.
    plt1=plt.subplot2grid((2,4),(0,1), colspan=3) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9984')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
        
        #Taken from here: https://stackoverflow.com/questions/29776114/half-violin-plot/29781988#29781988
        #m=np.mean(violins['bodies'][i].get_paths()[0].vertices[:, 0])
        #violins['bodies'][i].get_paths()[0].vertices[:, 0]=np.clip(violins['bodies'][i].get_paths()[0].vertices[:, 0], -np.inf, m)
        #violins['bodies'][i].set_color('r')
        
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['1\n\n\n\n', 2, 3, 4, 5, 6]
    set_axis_style(plt1, labels, pos1)    
    
    yticknames1=np.arange(0.2, 65, 5)
    #plt1.set_yticks(yticknames1, minor=False)
    #plt1.set_yticklabels(yticknames1)
    #plt1.set_ylabel('EcTopoI FE', size=15)
    plt1.set_ylim(0.2, 65)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=0)   
    plt1.set_yscale('log')
    
    plt1.spines['top'].set_visible(False)
    plt1.spines['right'].set_visible(False)    
    
    plt1.annotate(f'  {EL_pos.shape[0]}',                     xy=(0.5, 31), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {RNAP_pos.shape[0]}',                   xy=(1.5, 31), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {EL_RNAP_pos.shape[0]}',                xy=(2.5, 31), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'  {EL_RNAP_pos_div.shape[0]}',            xy=(3.5, 31), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'   {EL_RNAP_pos_div_TF_pos.shape[0]}',     xy=(4.5, 31), xycoords='data', size=15, rotation=0)
    plt1.annotate(f'   {EL_RNAP_pos_div_TF_pos_Mem.shape[0]}', xy=(5.5, 31), xycoords='data', size=15, rotation=0)    

    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(EL_pos.loc[:,                     "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.25), xycoords='data', size=10, rotation=0)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(RNAP_pos.loc[:,                  "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.25), xycoords='data', size=10, rotation=0)  
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(EL_RNAP_pos.loc[:,                "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(2.5, 0.25), xycoords='data', size=10, rotation=0)
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(EL_RNAP_pos_div.loc[:,            "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.41), xycoords='data', size=10, rotation=0) 
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(EL_RNAP_pos_div_TF_pos.loc[:,     "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.44), xycoords='data', size=10, rotation=0)
    plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(EL_RNAP_pos_div_TF_pos_Mem.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(5.5, 0.60), xycoords='data', size=10, rotation=0)       
    
    #Test EcTopoI FE difference between groups of IGRs.
    for i in range(len(dataset1)):
        for j in range(len(dataset1)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[i], dataset1[j], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i]),2)} Mean2={round(np.mean(dataset1[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    
    
    #Plot EcTopoI enrichment. Negative factors.
    dataset2=[EL_neg.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), RNAP_neg.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), EL_RNAP_neg.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), EL_RNAP_neg_not_div.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),
              EL_RNAP_neg_not_div_TF_neg.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), EL_RNAP_neg_not_div_TF_neg_not_Mem.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist()]    
    
    #Draw violin-plots.
    plt2=plt.subplot2grid((2,4),(1,1), colspan=3) 
    violins=plt2.violinplot(dataset2, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#a3c0ff')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['1\n\n\n\n', 2, 3, 4, 5, 6]
    set_axis_style(plt2, labels, pos1)    
    
    yticknames1=np.arange(0.2, 60, 5)
    #plt2.set_yticks(yticknames1, minor=False)
    #plt2.set_yticklabels(yticknames1)
    #plt2.set_ylabel('EcTopoI FE', size=15)
    plt2.set_ylim(0.2, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=0)   
    plt2.set_yscale('log')
    
    plt2.spines['top'].set_visible(False)
    plt2.spines['right'].set_visible(False)    

    plt2.annotate(f'  {EL_neg.shape[0]}',                             xy=(0.5, 13), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'  {RNAP_neg.shape[0]}',                           xy=(1.5, 28), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'  {EL_RNAP_neg.shape[0]}',                        xy=(2.5, 13), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'  {EL_RNAP_neg_not_div.shape[0]}',                xy=(3.5, 13), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'  {EL_RNAP_neg_not_div_TF_neg.shape[0]}',         xy=(4.5, 13), xycoords='data', size=15, rotation=0)
    plt2.annotate(f'  {EL_RNAP_neg_not_div_TF_neg_not_Mem.shape[0]}', xy=(5.5, 13), xycoords='data', size=15, rotation=0)    

    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(EL_neg.loc[:,                             "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.25), xycoords='data', size=10, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(RNAP_neg.loc[:,                           "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.25), xycoords='data', size=10, rotation=0)  
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(EL_RNAP_neg.loc[:,                        "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(2.5, 0.25), xycoords='data', size=10, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(EL_RNAP_neg_not_div.loc[:,                "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.25), xycoords='data', size=10, rotation=0) 
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(EL_RNAP_neg_not_div_TF_neg.loc[:,         "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.25), xycoords='data', size=10, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(EL_RNAP_neg_not_div_TF_neg_not_Mem.loc[:, "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(5.5, 0.25), xycoords='data', size=10, rotation=0)   
    
    #Test EcTopoI FE difference between groups of IGRs.
    for i in range(len(dataset2)):
        for j in range(len(dataset2)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset2[i], dataset2[j], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[i]),2)} Mean2={round(np.mean(dataset2[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    Intervals_stat=stats.ttest_ind(dataset0[0], dataset1[0], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset0[0]),2)} Mean2={round(np.mean(dataset1[0]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')   
    
    Intervals_stat=stats.ttest_ind(dataset0[0], dataset1[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset0[0]),2)} Mean2={round(np.mean(dataset1[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')      
    
    Intervals_stat=stats.ttest_ind(dataset0[0], dataset2[0], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset0[0]),2)} Mean2={round(np.mean(dataset2[0]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')   
    
    Intervals_stat=stats.ttest_ind(dataset0[0], dataset2[1], equal_var=False)
    print(f'\nT-test FE means\nMean1={round(np.mean(dataset0[0]),2)} Mean2={round(np.mean(dataset2[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')      
    
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\Add_factors\\NEW_Factors_combinations_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.svg', dpi=400, figsize=(7, 4.5), transparent=True)       
    
    return

#read_test_factors_combination_one_plot(PWD)


#######
#Read final table, test different types of promoters by sigma-factor.
#######


def read_test_sigma_factor(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']
    
    #Classify intergenic regions by type of a promoter.
    Sigma70=[]
    Sigma54=[]
    Sigma32=[]
    Sigma24=[]
    Sigma28=[]
    Sigma38=[]
    unknown=[]
    unknown_only=[]
    mixed=[]
    for i in input_data.index:
        Sigma_list=input_data.loc[i, 'Promoters_sigmas']
        print(i, Sigma_list, isinstance(Sigma_list,float))
        if isinstance(Sigma_list,float)==False:
            Sigma_list=Sigma_list.split(';')
            Sigma70_check=0
            Sigma54_check=0
            Sigma32_check=0
            Sigma24_check=0
            Sigma28_check=0
            Sigma38_check=0
            Unknown_check=0
            #Detect different sigmas.
            for sigma_type in Sigma_list:
                if sigma_type=="Sigma70":
                    Sigma70_check=1
                if sigma_type=="Sigma54":
                    Sigma54_check=1
                if sigma_type=="Sigma32":
                    Sigma32_check=1
                if sigma_type=="Sigma24":
                    Sigma24_check=1
                if sigma_type=="Sigma28":
                    Sigma28_check=1
                if sigma_type=="Sigma38":
                    Sigma38_check=1
                if sigma_type=="unknown":
                    Unknown_check=1
            
            #If sigmas annotated.
            if Sigma70_check==1:
                Sigma70.append(1)
            if Sigma54_check==1:
                Sigma54.append(1)
            if Sigma32_check==1:
                Sigma32.append(1)
            if Sigma24_check==1:
                Sigma24.append(1)
            if Sigma28_check==1:
                Sigma28.append(1)  
            if Sigma38_check==1:
                Sigma38.append(1) 
            if Unknown_check==1:
                unknown.append(1)
            #If sigmas not annotated.    
            if Sigma70_check==0:
                Sigma70.append(0)
            if Sigma54_check==0:
                Sigma54.append(0)
            if Sigma32_check==0:
                Sigma32.append(0)
            if Sigma24_check==0:
                Sigma24.append(0)
            if Sigma28_check==0:
                Sigma28.append(0)  
            if Sigma38_check==0:
                Sigma38.append(0) 
            if Unknown_check==0:
                unknown.append(0)
            #Several different sigmas detected.    
            if Sigma70_check+Sigma54_check+Sigma32_check+Sigma24_check+Sigma28_check+Sigma38_check>1:
                mixed.append(1)
                unknown_only.append(0)
            #Just one type of sigma detected.
            if Sigma70_check+Sigma54_check+Sigma32_check+Sigma24_check+Sigma28_check+Sigma38_check==1:
                mixed.append(0)
                unknown_only.append(0)
            #No annotated sigmas detected.
            if Sigma70_check+Sigma54_check+Sigma32_check+Sigma24_check+Sigma28_check+Sigma38_check==0:
                mixed.append(0)
                unknown_only.append(1)
            
        else:
            Sigma70.append(0)
            Sigma54.append(0)
            Sigma32.append(0)
            Sigma24.append(0)
            Sigma28.append(0)
            Sigma38.append(0)
            unknown.append(0)
            unknown_only.append(0)
            mixed.append(0)          
     
    
    print(len(Sigma70), len(Sigma54), len(Sigma32), len(Sigma24), len(Sigma28), len(Sigma38), len(unknown), len(mixed))
    
    #Add data to dataframe.
    input_data['Sigma70']=Sigma70
    input_data['Sigma54']=Sigma54
    input_data['Sigma32']=Sigma32
    input_data['Sigma24']=Sigma24
    input_data['Sigma28']=Sigma28
    input_data['Sigma38']=Sigma38
    input_data['Sigma_unknown']=unknown
    input_data['Sigma_only_unknown']=unknown_only
    input_data['Sigma_mixed']=mixed
    
    #Create datasets.
    Genes_sigma70=input_data[(input_data['Sigma70']==1)]
    Genes_sigma54=input_data[(input_data['Sigma54']==1)]
    Genes_sigma32=input_data[(input_data['Sigma32']==1)]
    Genes_sigma24=input_data[(input_data['Sigma24']==1)]
    Genes_sigma28=input_data[(input_data['Sigma28']==1)]
    Genes_sigma38=input_data[(input_data['Sigma38']==1)]
    Genes_sigma_unknown=input_data[(input_data['Sigma_unknown']==1)]
    Genes_sigma_only_unknown=input_data[(input_data['Sigma_only_unknown']==1)]
    Genes_sigma_mixed=input_data[(input_data['Sigma_mixed']==1)]
    
    
    print(Genes_sigma70.shape, Genes_sigma54.shape, Genes_sigma32.shape, Genes_sigma24.shape, Genes_sigma28.shape, Genes_sigma38.shape, Genes_sigma_unknown.shape, Genes_sigma_mixed.shape)    
    
    #Plot EcTopoI enrichment.
    pos1=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]    
    dataset1=[input_data.loc[:, 'EcTopoI-Rif-CTD_FE'], Genes_sigma70.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_sigma38.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_sigma54.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_sigma32.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),  Genes_sigma24.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),
              Genes_sigma28.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),  Genes_sigma_unknown.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_sigma_only_unknown.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(), Genes_sigma_mixed.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist(),
              input_data.loc[:, 'EcTopoI+Rif-CTD_FE'], Genes_sigma70.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_sigma38.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_sigma54.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_sigma32.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(),  Genes_sigma24.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(),
              Genes_sigma28.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(),  Genes_sigma_unknown.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_sigma_only_unknown.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist(), Genes_sigma_mixed.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]    
    
    pos2=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10] 
    dataset2=[input_data.loc[:, 'RNAP'], Genes_sigma70.loc[:, 'RNAP'].tolist(), Genes_sigma38.loc[:, 'RNAP'].tolist(), Genes_sigma54.loc[:, 'RNAP'].tolist(), Genes_sigma32.loc[:, 'RNAP'].tolist(),  Genes_sigma24.loc[:, 'RNAP'].tolist(),
              Genes_sigma28.loc[:, 'RNAP'].tolist(),  Genes_sigma_unknown.loc[:, 'RNAP'].tolist(),  Genes_sigma_only_unknown.loc[:, 'RNAP'].tolist(), Genes_sigma_mixed.loc[:, 'RNAP'].tolist(),]
    dataset3=[input_data.loc[:, 'Cumulative_expression'], Genes_sigma70.loc[:, 'Cumulative_expression'].tolist(), Genes_sigma38.loc[:, 'Cumulative_expression'].tolist(), Genes_sigma54.loc[:, 'Cumulative_expression'].tolist(), Genes_sigma32.loc[:, 'Cumulative_expression'].tolist(),  Genes_sigma24.loc[:, 'Cumulative_expression'].tolist(),
              Genes_sigma28.loc[:, 'Cumulative_expression'].tolist(),  Genes_sigma_unknown.loc[:, 'Cumulative_expression'].tolist(),  Genes_sigma_only_unknown.loc[:, 'Cumulative_expression'].tolist(), Genes_sigma_mixed.loc[:, 'Cumulative_expression'].tolist(),]
    
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(15,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i<len(pos2):
            violins['bodies'][i].set_facecolor('#ff7762')
        else:
            violins['bodies'][i].set_facecolor('#58ffb1')        
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['All', 'Sigma70', 'Sigma38', 'Sigma54', 'Sigma32', 'Sigma24', 'Sigma28', 'Unknown', 'Only\nunknown', 'Mixed', 'All', 'Sigma70', 'Sigma38', 'Sigma54', 'Sigma32', 'Sigma24', 'Sigma28', 'Unknown', 'Only\nunknown', 'Mixed']
    set_axis_style(plt1, labels, pos1)
    
    yticknames1=np.arange(0.1, 45, 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.1, 45)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    plt1.annotate(f' {input_data.shape[0]}',               xy=(0.5, 30),    xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma70.shape[0]}',            xy=(1.5, 30), xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma38.shape[0]}',            xy=(2.5, 14), xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma54.shape[0]}',            xy=(3.5, 19), xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma32.shape[0]}',            xy=(4.5, 15),    xycoords='data', size=13, rotation=0) 
    plt1.annotate(f' {Genes_sigma24.shape[0]}',            xy=(5.5, 15), xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma28.shape[0]}',            xy=(6.5, 11), xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma_unknown.shape[0]}',      xy=(7.5, 30), xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma_only_unknown.shape[0]}', xy=(8.5, 17), xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma_mixed.shape[0]}',        xy=(9.5, 15),    xycoords='data', size=13, rotation=0)     
    
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(input_data.loc[:,              "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(0.5, 0.30), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma70.loc[:,           "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(1.5, 0.30), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma38.loc[:,           "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(2.5, 0.30), xycoords='data', size=9, rotation=0)   
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma54.loc[:,           "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(3.5, 0.30), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma32.loc[:,           "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(4.5, 0.30), xycoords='data', size=9, rotation=0) 
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma24.loc[:,           "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(5.5, 0.30), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma28.loc[:,           "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(6.5, 0.30), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_unknown.loc[:,     "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(7.5, 0.30), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_only_unknown.loc[:,"EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(8.5, 0.30), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_mixed.loc[:,       "EcTopoI-Rif-CTD_FE"].tolist()),2)}', xy=(9.5, 0.30), xycoords='data', size=9, rotation=0)      
    
    plt1.annotate(f' {input_data.shape[0]}',               xy=(11.5, 23.9),  xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma70.shape[0]}',            xy=(12.5, 13.6),  xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma38.shape[0]}',            xy=(13.5, 13.6),  xycoords='data', size=13, rotation=0)    
    plt1.annotate(f' {Genes_sigma54.shape[0]}',            xy=(14.5, 22.9),  xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma32.shape[0]}',            xy=(15.5, 6),     xycoords='data', size=13, rotation=0) 
    plt1.annotate(f' {Genes_sigma24.shape[0]}',            xy=(16.5, 12.7),  xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma28.shape[0]}',            xy=(17.5, 11.4),  xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma_unknown.shape[0]}',      xy=(18.5, 23.4),  xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma_only_unknown.shape[0]}', xy=(19.5, 9.03),  xycoords='data', size=13, rotation=0)
    plt1.annotate(f' {Genes_sigma_mixed.shape[0]}',        xy=(20.5, 13.04), xycoords='data', size=13, rotation=0) 
    
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(input_data.loc[:,               "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(11.5, 0.37), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma70.loc[:,            "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(12.5, 0.45), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma38.loc[:,            "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(13.5, 0.51), xycoords='data', size=9, rotation=0)    
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma54.loc[:,            "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(14.5, 0.52), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma32.loc[:,            "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(15.5, 0.45), xycoords='data', size=9, rotation=0) 
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma24.loc[:,            "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(16.5, 0.45), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma28.loc[:,            "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(17.5, 0.45), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_unknown.loc[:,      "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(18.5, 0.36), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_only_unknown.loc[:, "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(19.5, 0.45), xycoords='data', size=9, rotation=0)
    plt1.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_mixed.loc[:,        "EcTopoI+Rif-CTD_FE"].tolist()),2)}', xy=(20.5, 0.45), xycoords='data', size=9, rotation=0)   
    
    #Test EcTopoI FE difference between groups of IGRs.
    for i in range(len(pos2)):
        for j in range(len(pos2)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j], dataset1[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i]),2)} Mean2={round(np.mean(dataset1[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    for i in range(len(pos2)):
        for j in range(len(pos2)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j+len(pos2)], dataset1[i+len(pos2)], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i+len(pos2)]),2)} Mean2={round(np.mean(dataset1[j+len(pos2)]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')       
    
    
    
    #RNAP fold enrichment.
    #Draw violin-plots. Promoter factor.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['All', 'Sigma70', 'Sigma38', 'Sigma54', 'Sigma32', 'Sigma24', 'Sigma28', 'Unknown', 'Only\nunknown', 'Mixed']
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.01, 60, 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.01, 60)
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    plt2.annotate(f' {input_data.shape[0]}',               xy=(0.5, 31),    xycoords='data', size=13, rotation=0)
    plt2.annotate(f' {Genes_sigma70.shape[0]}',            xy=(1.5, 31),    xycoords='data', size=13, rotation=0)
    plt2.annotate(f' {Genes_sigma38.shape[0]}',            xy=(2.5, 31),    xycoords='data', size=13, rotation=0)
    plt2.annotate(f' {Genes_sigma54.shape[0]}',            xy=(3.5, 31),    xycoords='data', size=13, rotation=0)
    plt2.annotate(f' {Genes_sigma32.shape[0]}',            xy=(4.5, 24),    xycoords='data', size=13, rotation=0) 
    plt2.annotate(f' {Genes_sigma24.shape[0]}',            xy=(5.5, 31),    xycoords='data', size=13, rotation=0)
    plt2.annotate(f' {Genes_sigma28.shape[0]}',            xy=(6.5, 12),    xycoords='data', size=13, rotation=0)
    plt2.annotate(f' {Genes_sigma_unknown.shape[0]}',      xy=(7.5, 31),    xycoords='data', size=13, rotation=0)
    plt2.annotate(f' {Genes_sigma_only_unknown.shape[0]}', xy=(8.5, 31),    xycoords='data', size=13, rotation=0)
    plt2.annotate(f' {Genes_sigma_mixed.shape[0]}',        xy=(9.5, 30.37), xycoords='data', size=13, rotation=0)     

    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(input_data.loc[:,               "RNAP"].tolist()),2)}', xy=(0.5, 0.045), xycoords='data', size=9, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma70.loc[:,            "RNAP"].tolist()),2)}', xy=(1.5, 0.045), xycoords='data', size=9, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma38.loc[:,            "RNAP"].tolist()),2)}', xy=(2.5, 0.055), xycoords='data', size=9, rotation=0)    
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma54.loc[:,            "RNAP"].tolist()),2)}', xy=(3.5, 0.056), xycoords='data', size=9, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma32.loc[:,            "RNAP"].tolist()),2)}', xy=(4.5, 0.045), xycoords='data', size=9, rotation=0) 
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma24.loc[:,            "RNAP"].tolist()),2)}', xy=(5.5, 0.073), xycoords='data', size=9, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma28.loc[:,            "RNAP"].tolist()),2)}', xy=(6.5, 0.11),  xycoords='data', size=9, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_unknown.loc[:,      "RNAP"].tolist()),2)}', xy=(7.5, 0.073), xycoords='data', size=9, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_only_unknown.loc[:, "RNAP"].tolist()),2)}', xy=(8.5, 0.073), xycoords='data', size=9, rotation=0)
    plt2.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_mixed.loc[:,        "RNAP"].tolist()),2)}', xy=(9.5, 0.045), xycoords='data', size=9, rotation=0)         
    
    #Test RNAP FE difference between groups of IGRs.
    for i in range(len(pos2)):
        for j in range(len(pos2)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset2[j], dataset2[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[i]),2)} Mean2={round(np.mean(dataset2[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
        
    #Expression level.
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=['All', 'Sigma70', 'Sigma38', 'Sigma54', 'Sigma32', 'Sigma24', 'Sigma28', 'Unknown', 'Only\nunknown', 'Mixed']
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.0002, 8000, 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.0002, 8000)
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log') 
    
    plt3.annotate(f' {input_data.shape[0]}',               xy=(0.5, 2300), xycoords='data', size=13, rotation=0)
    plt3.annotate(f' {Genes_sigma70.shape[0]}',            xy=(1.5, 2300), xycoords='data', size=13, rotation=0)
    plt3.annotate(f' {Genes_sigma38.shape[0]}',            xy=(2.5, 200),  xycoords='data', size=13, rotation=0)
    plt3.annotate(f' {Genes_sigma54.shape[0]}',            xy=(3.5, 370),  xycoords='data', size=13, rotation=0)
    plt3.annotate(f' {Genes_sigma32.shape[0]}',            xy=(4.5, 506),  xycoords='data', size=13, rotation=0) 
    plt3.annotate(f' {Genes_sigma24.shape[0]}',            xy=(5.5, 788),  xycoords='data', size=13, rotation=0)
    plt3.annotate(f' {Genes_sigma28.shape[0]}',            xy=(6.5, 506),  xycoords='data', size=13, rotation=0)
    plt3.annotate(f' {Genes_sigma_unknown.shape[0]}',      xy=(7.5, 2300), xycoords='data', size=13, rotation=0)
    plt3.annotate(f' {Genes_sigma_only_unknown.shape[0]}', xy=(8.5, 2300), xycoords='data', size=13, rotation=0)
    plt3.annotate(f' {Genes_sigma_mixed.shape[0]}',        xy=(9.5, 419),  xycoords='data', size=13, rotation=0)     
    
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(input_data.loc[:,              "Cumulative_expression"].tolist()),1)}', xy=(0.5, 0.0012), xycoords='data', size=9, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma70.loc[:,           "Cumulative_expression"].tolist()),1)}', xy=(1.5, 0.0012), xycoords='data', size=9, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma38.loc[:,           "Cumulative_expression"].tolist()),1)}', xy=(2.5, 0.006),  xycoords='data', size=9, rotation=0)    
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma54.loc[:,           "Cumulative_expression"].tolist()),1)}', xy=(3.5, 0.026),  xycoords='data', size=9, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma32.loc[:,           "Cumulative_expression"].tolist()),1)}', xy=(4.5, 0.011),  xycoords='data', size=9, rotation=0) 
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma24.loc[:,           "Cumulative_expression"].tolist()),1)}', xy=(5.5, 0.011),  xycoords='data', size=9, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma28.loc[:,           "Cumulative_expression"].tolist()),1)}', xy=(6.5, 0.024),  xycoords='data', size=9, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_unknown.loc[:,     "Cumulative_expression"].tolist()),1)}', xy=(7.5, 0.011),  xycoords='data', size=9, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_only_unknown.loc[:,"Cumulative_expression"].tolist()),1)}', xy=(8.5, 0.011),  xycoords='data', size=9, rotation=0)
    plt3.annotate(r" $\overline{X}$"+f'={round(np.mean(Genes_sigma_mixed.loc[:,       "Cumulative_expression"].tolist()),1)}', xy=(9.5, 0.005),  xycoords='data', size=9, rotation=0) 
        
    
    #Test Expression level difference between groups of IGRs.
    for i in range(len(pos2)):
        for j in range(len(pos2)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset3[j], dataset3[i], equal_var=False)
                print(f'\nT-test Expression level means\nMean1={round(np.mean(dataset3[i]),2)} Mean2={round(np.mean(dataset3[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')                 
                
           
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_Sigma_factors_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.svg', dpi=400, figsize=(14, 6.5), transparent=True)         
    
    return


#read_test_sigma_factor(PWD)



#######
#Read final table, test different transcription factors.
#######


def read_test_transcription_factors(pwd):
    
    #Read input data table.
    input_data=pd.read_excel(pwd+'Intergenic_regions_NO_DPS_NO_rRNA_info_promoters_TF_sites_Membraness_EcTopoI_signal_new_RNAP_signal_PSORT.xlsx', header=0, index_col=0, sheet_name='Intergenic_regions_info')
    input_data['Cumulative_expression']=input_data['G1_expression']+input_data['G2_expression']
    
    #Classify intergenic regions by type of transcription factor.
    TF_dict={}
    no_TF_count=0
    first_TF=0
    for i in input_data.index:
        TF_list=input_data.loc[i, 'TF_names']
        print(i, TF_list, isinstance(TF_list, float))
        if isinstance(TF_list, float)==False:
            #Take TF names.
            TF_list=TF_list.split(';')
            
            #Unique TFs.
            TF_list=list(set(TF_list))
            print(TF_list)
            
            #Initiation step.
            if first_TF==0 and no_TF_count>0:
                Len_ar=no_TF_count   
                
            #Length of array on a previous step.
            elif first_TF!=0:
                Len_ar=len(TF_dict[list(TF_dict.keys())[0]])
            
            first_TF=1
               
            #Detect different TFs.
            for TF_type in TF_list:
                #If TF is already present, add new row.
                if TF_type in TF_dict:
                    TF_dict[TF_type].append(1)

                #If TF first time detected, add new list, add new row.
                else:
                    New_list=[0]*Len_ar
                    TF_dict[TF_type]=New_list
                    TF_dict[TF_type].append(1)
            
            #Add zeros.
            for TF, TF_ar in TF_dict.items():
                if len(TF_ar)==Len_ar:
                    TF_dict[TF].append(0)
                    
        else:
            no_TF_count+=1
            for TF, TF_ar in TF_dict.items():
                TF_dict[TF].append(0)            
                    
 
    #Add TF data to dataframe. 
    for TF, TF_ar in TF_dict.items():
        print(TF, len(TF_ar))
        input_data[TF]=TF_ar
        
    #Create datasets, position arrays, labels arrays.
    dataset_noRif_mean=[np.mean(input_data.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist())]
    dataset_noRif=[input_data.loc[:, 'EcTopoI-Rif-CTD_FE'].tolist()]
    dataset_Rif=[input_data.loc[:, 'EcTopoI+Rif-CTD_FE'].tolist()]
    dataset2=[input_data.loc[:, 'RNAP'].tolist()]
    dataset3=[input_data.loc[:, 'Cumulative_expression'].tolist()]
    pos2=[1]
    labels2=['All']
    i=2
    for TF, TF_ar in TF_dict.items():
        if input_data[(input_data[TF]==1)].shape[0]>15:
            pos2.append(i)
            labels2.append(TF)
            dataset_noRif.append(input_data[(input_data[TF]==1)].loc[:, 'EcTopoI-Rif-CTD_FE'].tolist())
            dataset_noRif_mean.append(np.mean(input_data[(input_data[TF]==1)].loc[:, 'EcTopoI-Rif-CTD_FE'].tolist()))
            dataset_Rif.append(input_data[(input_data[TF]==1)].loc[:, 'EcTopoI+Rif-CTD_FE'].tolist())
            dataset2.append(input_data[(input_data[TF]==1)].loc[:, 'RNAP'].tolist())
            dataset3.append(input_data[(input_data[TF]==1)].loc[:, 'Cumulative_expression'].tolist())
            i+=1
    
    #Sort all data by mean EcTopoI signal of IG sets.
    
    print(pos2)
    print(labels2)
    print(dataset_noRif_mean)
    
    dataset_noRif=[x for _,x in sorted(zip(dataset_noRif_mean,dataset_noRif))]
    dataset_Rif=[x for _,x in sorted(zip(dataset_noRif_mean,dataset_Rif))]
    dataset2=[x for _,x in sorted(zip(dataset_noRif_mean,dataset2))]
    dataset3=[x for _,x in sorted(zip(dataset_noRif_mean,dataset3))]
    labels2=[x for _,x in sorted(zip(dataset_noRif_mean,labels2))]
    pos2=np.arange(1,(len(labels2)+1),1).tolist()
    
    
    #Combine dataset, positions, labels for EcTopoI.        
    dataset1=dataset_noRif+dataset_Rif
    pos1=pos2+list(np.array(pos2)+len(pos2)+2)
    labels1=labels2+labels2
    
    print(pos1)
    print(labels1)
    
    #Report size of dataframes.
    for TF_set in dataset_noRif:
        print(len(TF_set))
    
    print(pos1, labels1, len(dataset1))
    
    #EcTopoI fold enrichment.
    #Draw violin-plots.
    fig=plt.figure(figsize=(14,6.5), dpi=100)
    plt1=fig.add_subplot(2,1,2) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        if i<len(pos2):
            violins['bodies'][i].set_facecolor('#ff7762')
        else:
            violins['bodies'][i].set_facecolor('#58ffb1')        
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=labels1
    set_axis_style(plt1, labels, pos1)
    
    yticknames1=np.arange(0.1, 5*max(max(dataset1)), 5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('EcTopoI fold enrichment', size=15)
    plt1.set_ylim(0.1, 5*max(max(dataset1)))
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt1.set_yscale('log')
    
    #Place set size annotation.
    for i in range(len(dataset1)):
        if i<len(dataset_noRif):
            plt1.annotate(f' {len(dataset1[i])}', xy=((i+0.5), (1.2*max(dataset1[i]))), xycoords='data', size=12, rotation=0)
            plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(dataset1[i]),2)}', xy=((i+0.75), (0.3*min(dataset1[i]))), xycoords='data', size=8, rotation=90)
        else:
            plt1.annotate(f' {len(dataset1[i])}', xy=((i+2.5), (1.2*max(dataset1[i]))), xycoords='data', size=12, rotation=0)
            plt1.annotate(r"  $\overline{X}$"+f'={round(np.mean(dataset1[i]),2)}', xy=((i+2.75), (0.3*min(dataset1[i]))), xycoords='data', size=8, rotation=90)           
    

    #Test EcTopoI FE difference between groups of IGRs.
    for i in range(len(dataset_noRif)):
        for j in range(len(dataset_noRif)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j], dataset1[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i]),2)} Mean2={round(np.mean(dataset1[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
    #Test EcTopoI Rif FE difference between groups of IGRs.
    for i in range(len(dataset_Rif)):
        for j in range(len(dataset_Rif)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[j+len(dataset_Rif)], dataset1[i+len(dataset_Rif)], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i+len(dataset_Rif)]),2)} Mean2={round(np.mean(dataset1[j+len(dataset_Rif)]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')       
    
    
    
    #RNAP fold enrichment.
    #Draw violin-plots. Promoter factor.
    plt2=fig.add_subplot(2,2,1) 
    violins=plt2.violinplot(dataset2, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff9df0')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)   
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=labels2
    set_axis_style(plt2, labels, pos2)    
    
    yticknames1=np.arange(0.01, 10*max(max(dataset2)), 5)
    plt2.set_yticks(yticknames1, minor=False)
    plt2.set_yticklabels(yticknames1)
    plt2.set_ylabel('RNAP fold enrichment', size=15)
    plt2.set_ylim(0.01, 10*max(max(dataset2)))
    plt.setp(plt2.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt2.set_yscale('log')
    
    #Place set size annotation.
    for i in range(len(dataset2)):
        plt2.annotate(f' {len(dataset2[i])}', xy=((i+0.5), (1.2*max(dataset2[i]))), xycoords='data', size=12, rotation=0)
        plt2.annotate(r"  $\overline{X}$"+f'={round(np.mean(dataset2[i]),2)}', xy=((i+0.75), (0.15*min(dataset2[i]))), xycoords='data', size=8, rotation=90)

    #Test RNAP FE difference between groups of IGRs.
    for i in range(len(dataset2)):
        for j in range(len(dataset2)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset2[j], dataset2[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset2[i]),2)} Mean2={round(np.mean(dataset2[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')

        
    #Expression level.
    #Draw violin-plots. Promoter factor.
    plt3=fig.add_subplot(2,2,2) 
    violins=plt3.violinplot(dataset3, positions=pos2, widths=0.77, showmeans=True, showmedians=True, points=500)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ffe294')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)  
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=labels2
    set_axis_style(plt3, labels, pos2)    
    
    yticknames1=np.arange(0.0002, 10*max(max(dataset3)), 5)
    plt3.set_yticks(yticknames1, minor=False)
    plt3.set_yticklabels(yticknames1)
    plt3.set_ylabel('Expression level', size=15)
    plt3.set_ylim(0.0002, 10*max(max(dataset3)))
    plt.setp(plt3.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    plt3.set_yscale('log') 
    
    #Place set size annotation.
    for i in range(len(dataset3)):
        plt3.annotate(f' {len(dataset3[i])}', xy=((i+0.5), (1.2*max(dataset3[i]))), xycoords='data', size=12, rotation=0)
        plt3.annotate(r"  $\overline{X}$"+f'={round(np.mean(dataset3[i]),1)}', xy=((i+0.75), (0.04*min(dataset3[i]))), xycoords='data', size=8, rotation=90)

    #Test Expression level difference between groups of IGRs.
    for i in range(len(dataset3)):
        for j in range(len(dataset3)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset3[j], dataset3[i], equal_var=False)
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset3[i]),2)} Mean2={round(np.mean(dataset3[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')
 
           
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}Factors_analysis\\NEW_Transcription_factors_and_EcTopoI_enrichment_in_IG_regions_no_rRNA.svg', dpi=400, figsize=(14, 6.5), transparent=True)         
    
    return


#read_test_transcription_factors(PWD)