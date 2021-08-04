###############################################
##Dmitry Sutormin, 2021##
##TopoI ChIP-Seq analysis##

#Script computes Fold Enrichment (FE) near TSS (-5000bp:500bp) and TES (-500bp:5000bp)

###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import random
from Bio import SeqIO


#Path to the directory with input files.
Path_to_input_files='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\\'
#Path to TUs groups file.
TUs_groups_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_transcripts\\Representative_transcripts\\"
#Path to the input annotation, type of annotation and name of TUs set.             

##1##                                   
Path_to_annotation_1=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_NO_DPS_EP_del_cor.txt'
Type_of_annot_1='broadPeak'             
Genes_set_name_1='All_TUs_no_dps_1660'    
##2##                                   
Path_to_annotation_2=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_NO_DPS_NO_RFA_no_tRNA_rRNA_EP_del_cor_HETU_200.txt'
Type_of_annot_2='broadPeak'             
Genes_set_name_2='HETU_no_dps_rfa_200'         
##3##                                   
Path_to_annotation_3=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_NO_DPS_no_tRNA_rRNA_EP_del_cor_LETU_200.txt'
Type_of_annot_3='broadPeak'             
Genes_set_name_3='LETU_no_dps_200'    

#Path to the file with regions to be omitted (e.g. deletions).
Deletions_inpath='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\Deletions_w3110_G_Mu_SGS.broadPeak'
#Width of US, DS regions.
US_DS_length=500
#Length of GB.
GB_length=500

#Dictionary of pathes to input data.
Dict_of_wigs_path={'RpoC_Borukhov' : Path_to_input_files+'Borukhov_RpoC_Pol_Sofi_LB_FE.wig',
                   'RpoB_Borukhov_wt' : Path_to_input_files + 'Borukhov_RpoB_BW28357_FE_av.wig',
                   'RpoB_Kahramanoglou' : Path_to_input_files+'Kahramanoglou_RpoB_IP_ME.wig',
                   'RpoS_Peano' : Path_to_input_files+'Peano_RpoS_FE_av.wig',
                   'RpoS_Seo' : Path_to_input_files+'Seo_RpoS_FE_Rep1.wig',
                   'RpoD_Myers' : Path_to_input_files + 'Myers_RpoD_FE_av.wig',
                   'HNS_Kahramanoglou' : Path_to_input_files+'Kahramanoglou_HNS_IP_ME.wig',
                   'Fis_Kahramanoglou' : Path_to_input_files+'Kahramanoglou_Fis_IP_ME.wig',
                   'MukB_Nolivos' : Path_to_input_files+'Nolivos_MukB_IP_av.wig',
                   'MatP_Nolivos' : Path_to_input_files+'Nolivos_MatP_IP_av.wig',
                   'TopA_CTD_minus_Rif_minus_av_1_2_3' : Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_123.wig',
                   'TopA_CTD_minus_Rif_plus_av_1_2_3' : Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_minus_Rif_plus_FE_av_123.wig',
                   'TopA_CTD_plus_Rif_minus_av_1_2_3' : Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_plus_Rif_minus_FE_av.wig',
                   'TopA_CTD_plus_Rif_plus_av_2_3' :   Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_plus_Rif_plus_FE_av.wig',                   
                   'Gyrase_Cfx' : Path_to_input_files+'Sutormin_Gyrase_Cfx_10mkM_FE_av.wig',
                   'Gyrase_Cfx_Rif' : Path_to_input_files+'Sutormin_Gyrase_RifCfx_122mkM_10mkM_FE_av.wig',
                   'TopoIV_Cfx' : Path_to_input_files+'Sutormin_TopoIV_Cfx_FE_av.wig',
                   'GC' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig',
                   'RNA_Seq' : Path_to_input_files+'Sutormin_RNA_Seq_Exponential_av.wig'
                   }

#######
#Checks if directory exists and if not creates it.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Path to the output directory.
Out_path='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Signal_over_TUs\Representative_transcripts\\'

#Output path.
def create_out_dirs(out_path, genes_set_name):
    Dir_check_create(out_path)
    Dir_check_create(out_path+'\Figures\Plots_TSS_TES\\'+genes_set_name)
    Dir_check_create(out_path+'\Figures\Histograms_TSS_TES\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TSS_TES_tab\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TSS_TES_wig\\'+genes_set_name)    
    return

create_out_dirs(Out_path, Genes_set_name_1)
create_out_dirs(Out_path, Genes_set_name_2)
create_out_dirs(Out_path, Genes_set_name_3)


#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    return NE_values


#######
#Reads annotation of particular set of genes .tab BroadPeak-like (determined on a basis of expression level).
#######

def parse_expression_annotation(annot_inpath):
    genes_annotation={}
    filein=open(annot_inpath, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['GeneID', 'OperonID', 'TU_ID']:
            TU_name=line[1].lstrip('"').rstrip(';"')
            TU_start=int(line[2])
            TU_end=int(line[3])
            TU_strand=line[4]
            TU_expression=float(line[5].replace(',','.'))
            genes_annotation[TU_name]=[TU_start, TU_end, TU_strand, TU_expression]
    filein.close()            
    return genes_annotation


#######
#Write .wig file.
#######

def write_wig(ar, fileout_path, name):
    fileout=open(fileout_path, 'w')
    fileout.write(f'track type=wiggle_0 name="{name}" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom=NC_007779.1_w3110_Mu start=1 step=1\n')
    for point in ar:
        fileout.write(f'{point}\n')
    fileout.close()
    return


#######
#Opens and reads BED file with deletions coordinates.
#Example:
#GenomeID\tStart\tEnd
#NC_007779.1_w3110_Mu\t274500\t372148
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar


#######
#Write .tab file with FE info for genes US, GB, and DS.
#######

def write_genes_FE(dict1, dict2, FE_track_name, path_out):
    fileout=open(path_out, 'w')
    fileout.write(f'Gene_name\tStart\tEnd\tStrand\t{FE_track_name}_FE_US\t{FE_track_name}_FE_DS\n')
    for k, v in dict1.items():
        fileout.write(f'{k}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[0]}\t{dict2[k][0]}\n')
    fileout.close()
    return

#######
#Convert dictionary to array, discard keys.
#######

def dict_to_ar(dictionary):
    ar=[]
    for k,v in dictionary.items():
        ar.append(v[0]) 
    return ar


#########
##Makes histogram for FE over TUs: US, GB, DS.
#########

def plot_FE_dist_UDB(ar0, name0, ar1, name1, pathout):
    #Plot distribution of FE values.
    
    mean_FE0=round(np.mean(ar0),2)
    print(f'Mean FE in {name0}={mean_FE0}')
    fig=plt.figure(figsize=(8, 3), dpi=100)
    bins0=np.arange(min(ar0+ar1), max(ar0+ar1), 0.25)
    plot0=plt.subplot2grid((1,2),(0,0), rowspan=1, colspan=1)
    plot0.hist(ar0, bins0, color='#ff878b', edgecolor='black', alpha=0.8, label=f'{name0}')
    plot0.annotate(f'Mean FE={mean_FE0}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot0.set_yscale('log')
    plot0.set_xlabel('Fold enrichment', size=17)
    plot0.set_ylabel('Number of TUs', size=17)
    plot0.set_title(name0, size=18)  
    #plot0.legend(fontsize=22)
       
    mean_FE1=round(np.mean(ar1),2)
    print(f'Mean FE in {name1}={mean_FE1}')
    bins1=np.arange(min(ar0+ar1), max(ar0+ar1), 0.25)
    plot1=plt.subplot2grid((1,2),(0,1), rowspan=1, colspan=1)     
    plot1.hist(ar1, bins1, color='#ffce91', edgecolor='black', alpha=0.5, label=f'{name1}')
    plot1.annotate(f'Mean FE={mean_FE1}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot1.set_yscale('log')
    plot1.set_xlabel('Fold enrichment', size=17)
    plot1.set_ylabel('Number of TUs', size=17)
    plot1.set_title(name1, size=18) 
    #plot1.legend(fontsize=22) 
    
    plt.tight_layout()
    plt.show()
    plt.savefig(pathout, dpi=300, figsize=(8, 3))
    plt.close() 
    return


#######
#Returns FE or Ded FE over the set of genes (US, GB, DS) - for each gene separately.
#######

def genes_and_FE(gene_annotation, genes_set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 

    #Calculate FE over genes.
    gene_US=np.array([0.0]*(win_width+length))
    gene_DS=np.array([0.0]*(win_width+length))
    gene_US_mean_dict={}
    gene_DS_mean_dict={}
    for gene_name, gene_info in gene_annotation.items():
        delited=0
        for deletion in deletions:
            if deletion[1]>=gene_info[0]>=deletion[0] or deletion[1]>=gene_info[1]>=deletion[0]:
                delited=1
        if delited==0:
            start=gene_info[0]
            end=gene_info[1]
            glen=len(FE_track)
            if gene_info[2]=='+':
                if start<win_width:
                    gene_US+=np.array(FE_track[glen-(win_width-start):] + FE_track[:start+length])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[glen-(win_width-start):] + FE_track[:start+length]), start, end, gene_info[2]]
                elif start+length>glen:
                    gene_US+=np.array(FE_track[start-win_width:] + FE_track[:start+length-glen])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:] + FE_track[:start+length-glen]), start, end, gene_info[2]]                    
                else:
                    gene_US+=np.array(FE_track[start-win_width:start+length])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start+length]), start, end, gene_info[2]]
                    
                if end+win_width>glen:
                    gene_DS+=np.array(FE_track[end-length:] + FE_track[:end+win_width-glen])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end-length:] + FE_track[:end+win_width-glen]), start, end, gene_info[2]]
                elif end-length<0:
                    gene_DS+=np.array(FE_track[glen-(length-end):]+FE_track[:end+win_width])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[glen-(length-end):]+FE_track[:end+win_width]), start, end, gene_info[2]]                    
                else:
                    gene_DS+=np.array(FE_track[end-length:end+win_width])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end-length:end+win_width]), start, end, gene_info[2]]

            elif gene_info[2]=='-':
                if start<win_width:
                    gene_DS+=np.array(FE_track[:start+length][::-1] + FE_track[glen-(win_width-start):][::-1])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[:start+length][::-1] + FE_track[glen-(win_width-start):][::-1]), start, end, gene_info[2]]
                elif start+length>glen:
                    gene_DS+=np.array(FE_track[:length-(glen-start)][::-1] + FE_track[start-win_width:][::-1])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[:length-(glen-start)][::-1] + FE_track[start-win_width:][::-1]), start, end, gene_info[2]]                    
                else:
                    gene_DS+=np.array(FE_track[start-win_width:start+length][::-1])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start+length][::-1]), start, end, gene_info[2]]
                if end+win_width>glen:
                    gene_US+=np.array(FE_track[:end+win_width-glen][::-1] + FE_track[end-length:][::-1])  
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[:end+win_width-glen][::-1] + FE_track[end-length:][::-1]), start, end, gene_info[2]]
                else:
                    gene_US+=np.array(FE_track[end-length:end+win_width][::-1]) 
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[end-length:end+win_width][::-1]), start, end, gene_info[2]]
                    
    Num_genes=len(gene_annotation)
    gene_US=gene_US/Num_genes
    gene_DS=gene_DS/Num_genes
    print(f'FE over TUs computed...')

    #Write wig-like file with FE over US, DS.
    print(f'Writing FE over TU, DS...')
    write_wig(np.concatenate((gene_US, gene_DS), axis=None), f'{out_path}\Signal_of_TSS_TES_wig\\{genes_set_name}\Signal_{FE_track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp.wig', f'{win_width}_{length}')

    #Plot FE over US, DS. 
    print(f'Plotting FE over TU, DS...')
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    positions_DS=np.arange(1+length,win_width+1+length+length,1)
    positions_US=np.arange(-1-win_width,-1+length,1)  
    plot1.plot(positions_US, gene_US, linestyle='-', color='#D17E7E', linewidth=1, label='Rep12') 
    plot1.plot(positions_DS, gene_DS, linestyle='-', color='#D17E7E', linewidth=1)
    ticks=np.arange(-win_width,win_width+length+length+1,length).tolist()
    ticks+=[0, 2*length]
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(2*length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)       
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{FE_track_name} fold enrichment', size=20)
    plot1.set_title(f'{FE_track_name} over the {genes_set_name}s', size=20)     
    plt.savefig(f'{out_path}\Figures\Plots_TSS_TES\\{genes_set_name}\\{FE_track_name}_over_{genes_set_name}_{win_width}bp.png', dpi=400, figsize=(10, 6))   
    plt.close()      

    #Make ar from dict.
    gene_US_mean=dict_to_ar(gene_US_mean_dict)
    gene_DS_mean=dict_to_ar(gene_DS_mean_dict)

    #Write table contains FE for US, GB, DS of TUs in a set.
    print(f'Writing FE for TUs\' TU, DS...')
    write_genes_FE(gene_US_mean_dict, gene_DS_mean_dict, FE_track_name, f'{out_path}\Signal_of_TSS_TES_tab\\{genes_set_name}\{FE_track_name}_over_{genes_set_name}_{win_width}bp.txt')

    #Plot distribution of mean TUs' FEs.
    print(f'Plotting FE distribution over TU, DS...')
    plot_FE_dist_UDB(gene_US_mean, f'{FE_track_name} US', gene_DS_mean, f'{FE_track_name} DS', f'{out_path}\Figures\Histograms_TSS_TES\\{genes_set_name}\Signal_distribution_{FE_track_name}_over_{genes_set_name}_{win_width}bp.png')
    print(len(gene_US_mean), len(gene_DS_mean))
    return gene_US, gene_DS, gene_US_mean, gene_DS_mean, gene_US_mean_dict, gene_DS_mean_dict

#######
#Wrapper: reads data and gene annotation, computes signal over US, GB, DS, plots and writes the output.
#######


def Wrapper_signal_over_TUs(dict_of_wigs_path, path_to_annotation, type_of_annot, genes_set_name, deletions_inpath, win_width, length, out_path):
    #Reads input data in wig files.
    dict_of_wigs={}
    for name, data_path in dict_of_wigs_path.items():
        dict_of_wigs[name]=wig_parsing(data_path)
    
    #Reads annotation.
    print(f'Now working with {path_to_annotation}')
    gene_annotation=parse_expression_annotation(path_to_annotation)
    
    #Calculate and plot signal over TUs.
    for FE_track_name, FE_track in dict_of_wigs.items():
        genes_and_FE(gene_annotation, genes_set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length)
    return

Wrapper_signal_over_TUs(Dict_of_wigs_path, Path_to_annotation_1, Type_of_annot_1, Genes_set_name_1, Deletions_inpath, US_DS_length, GB_length, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_path, Path_to_annotation_2, Type_of_annot_2, Genes_set_name_2, Deletions_inpath, US_DS_length, GB_length, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_path, Path_to_annotation_3, Type_of_annot_3, Genes_set_name_3, Deletions_inpath, US_DS_length, GB_length, Out_path)
