###############################################
##Dmitry Sutormin, 2020##
##TopoI ChIP-Seq analysis##

#Script extracts upstream regions of transcription units and alignes them by transcription start site.
#Then frequency of nucleotides is calculated and logo is constructed. Mean position-wise score of EcTopoI binding motif is
#calculated and aligned with the logo.
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
from Bio.Seq import Seq
from Bio import AlignIO
import weblogo
from weblogo import LogoOptions, LogoFormat, eps_formatter, read_seq_data, LogoData, png_formatter, pdf_formatter


#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\\"

#Pass to transcription units annotation, bed.
#Transcripts_annotation="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_transcripts\\Representative_transcripts\DY330_RNA-Seq_transcripts_representative_EP_del_cor.txt"
Transcripts_annotation=PWD + "Promoter_structure\DY330_RNA-Seq_transcripts_representative_EP_del_cor_forward_TUs.txt"
#Transcripts_annotation=PWD + "Promoter_structure\DY330_RNA-Seq_transcripts_representative_EP_del_cor_reverse_TUs.txt"
#Transcripts_annotation=PWD + "Promoter_structure\TopA_minus_CTD_plus_Rif_av123_TSS_TES_500_LS_0.88.txt"
#Transcripts_annotation=PWD + "Promoter_structure\TopA_minus_CTD_plus_Rif_av123_TSS_TES_500_HS_2.txt"


#DY330 genome scored by EcTopoI binding motif, wig.
#Scores_dict={'Minus' : PWD + 'Motif_scanning\SARUS_new\EcTopoI_noCTD_Rif_motif_1_w3110_Mu_scanned_minus.wig',
#             'Plus' :  PWD + 'Motif_scanning\SARUS_new\EcTopoI_noCTD_Rif_motif_1_w3110_Mu_scanned_plus.wig',
#             'Both' :  PWD + 'Motif_scanning\SARUS_new\EcTopoI_noCTD_Rif_motif_1_w3110_Mu_scanned_both.wig',}
Scores_dict={'GC'     : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig',
             'EcRpoC' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\Pol_Sofi_LB_w3110_for_Mu.wig'}

#Path to genome sequence, fasta.
Genome_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"

#Path to files with regions to be omitted, bed.
Deletions_inpath='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\Deletions_w3110_G_Mu_SGS.broadPeak'

#Output path.
Outpath=PWD + "Promoter_structure\\"


#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
def create_out_dirs(out_path, genes_set_name):
    Dir_check_create(out_path)
    Dir_check_create(out_path+'\Figures\Plots_TSS_TES\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TSS_TES_tab\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TSS_TES_wig\\'+genes_set_name)    
    Dir_check_create(out_path+'\Sequences_of_TSS_TES_fasta\\') 
    
    return


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
#Parses WIG file.
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
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
    return genomefa


#######
#Opens and reads BED file with deletions coordinates.
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
#Write .fasta file with sequences of TUs US and DS regions.
#######

def write_sequences_dict(seq_dict, win_width, length, outpath):
    fileout=open(outpath, 'w')
    for TU_name, sequence_info in seq_dict.items():
        fileout.write(f'>{TU_name}_{sequence_info[1]}_{sequence_info[2]}_{sequence_info[3]}_{win_width}_{length}\n{sequence_info[0]}\n')
    fileout.close()
    return


#######
#Reads SARUS output for mfa and classify data on "+" (plus) strand, "-" (minus) strand.
#######

def Read_sarus_plot(sarus_output, win_width, length, genes_set_name, out_path):
    
    #Read Sarus output.
    filein=open(sarus_output, 'r')
    
    init=0
    dict_plus_minus={}
    for line in filein:
        if line[0]=='>':

            if init==0:
                ar_plus=[]
                ar_minus=[]  
            elif init==1:
                seq_len=len(ar_plus)
                print(ar_plus)
                print(ar_minus)
                dict_plus_minus[seq_id]=[ar_plus, ar_minus]
                ar_plus=[]
                ar_minus=[]                
                
            line=line.lstrip('>').rstrip().split(' ')
            seq_id=line[0] 
            init=1
            
        else:
            line=line.rstrip().split('\t')
            score=float(line[0])
            position=int(line[1])
            strand=line[2]
            if strand=='+':
                ar_plus.append(score)
            elif strand=='-':
                ar_minus.append(score)

    filein.close()
    print(f'{sarus_output} was parsed succesfully')
    
    #Plot Sarus output as metagene plot.
    #Calculate signal over TSS/TES.
    gene_fw=np.array([0.0]*(seq_len))
    gene_rv=np.array([0.0]*(seq_len))   
    
    for gene_name, scan_info in dict_plus_minus.items():
        gene_fw+=np.array(scan_info[0])
        gene_rv+=np.array(scan_info[1])

    Num_genes=len(dict_plus_minus)
    gene_fw=gene_fw/Num_genes
    gene_rv=gene_rv/Num_genes
    print(f'Signal over anchor is computed...')

    #Write wig-like file with FE over US, DS.
    #print(f'Writing signal over anchor...')
    #write_wig(np.concatenate((gene_US, gene_DS), axis=None), f'{out_path}\Signal_of_TSS_TES_wig\\{genes_set_name}\Signal_{FE_track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp.wig', f'{win_width}_{length}')

    #Plot signal over anchor. 
    print(f'Plotting signal over anchor...')
    plt.figure(figsize=(6, 6), dpi=100)
    plot1=plt.subplot(111)  
    positions_US=np.arange(-1-win_width,-1-win_width+seq_len,1)  
    plot1.plot(positions_US, gene_fw, linestyle='-', color='#D17E7E', linewidth=2, label='fw') 
    plot1.plot(positions_US, gene_rv, linestyle='--', color='#D17E7E', linewidth=2, label='rv') 
    ticks=np.arange(-win_width,seq_len,10).tolist()
    ticks+=[0]
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)       
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{genes_set_name} score', size=20)     
    plt.savefig(f'{out_path}\\Try1\Sequences_of_TSS_TES__sarus_scan_noCTD_Rif_3_T_F\{genes_set_name}_{win_width}bp_{length}bp.png', dpi=400, figsize=(10, 6))   
    plt.close()         
    
    return


#######
#Reads promoters annotation from RegulonDB. Classify promoters by type of a sigma-factor.
#Aligns selected promoters and returns logo. Sends the promoters to SARUS and returns constructs metagene plot on SARUS output.
#######

def Select_promoters_run_sarus_plot(rdb_promoters_inpath, win_width, outpath):
    
    #Read regulonDB file.
    promoters_dict={'Strong' : {'Sigma70' : [], 'Sigma54' : [], 'Sigma32' : [], 'Sigma38' : [], 'Sigma28' : [], 'Sigma24' : [], 'mixed' : [], 'unknown' : []}, 'Weak' : {'Sigma70' : [], 'Sigma54' : [], 'Sigma32' : [], 'Sigma38' : [], 'Sigma28' : [], 'Sigma24' : [], 'mixed' : [], 'unknown' : []}}
    promoters_infile=open(rdb_promoters_inpath, 'r')
    for line in promoters_infile:
        line=line.rstrip().split('\t')
        if line[0] not in ['RegulonDB ID']:
            name=line[1]
            sigma=line[4]
            if ', ' in sigma:
                sigma='mixed'
            sequence=line[5]
            confidence=line[7]
                
            promoters_dict[confidence][sigma].append(sequence)
    
    #Write mfa files for SARUS.
    #Run SARUS.
    for conf_name, conf_set in promoters_dict.items():
        for sigma_name, sigma_set in conf_set.items():
            fileout=open(f'{outpath}\{conf_name}_{sigma_name}_promoters.fasta', 'w')
            for i in range(len(sigma_set)):
                fileout.write(f'>{i}_{conf_name}_{sigma_name}\n{sigma_set[i]}\n')
            fileout.close()
            #Tips are taken from https://datatofish.com/command-prompt-python/
            os.system(f'cmd /c "java -jar C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Dists\sarus-2.0.2.jar {outpath}\{conf_name}_{sigma_name}_promoters.fasta C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks_motifs_3.txt -100 > {outpath}\{conf_name}_{sigma_name}_sarus_scan_noCTD_Rif_3.txt"')
            
    
    #Read Sarus output.
    for conf_name, conf_set in promoters_dict.items():
        for sigma_name, sigma_set in conf_set.items():
            sarus_output=f'{outpath}\{conf_name}_{sigma_name}_sarus_scan_noCTD_Rif_3.txt'
            filein=open(sarus_output, 'r')
            
            init=0
            dict_plus_minus={}
            for line in filein:
                if line[0]=='>':
        
                    if init==0:
                        ar_plus=[]
                        ar_minus=[]  
                    elif init==1:
                        seq_len=len(ar_plus)
                        print(ar_plus)
                        print(ar_minus)
                        dict_plus_minus[seq_id]=[ar_plus, ar_minus]
                        ar_plus=[]
                        ar_minus=[]                
                        
                    line=line.lstrip('>').rstrip().split(' ')
                    seq_id=line[0] 
                    init=1
                    
                else:
                    line=line.rstrip().split('\t')
                    score=float(line[0])
                    position=int(line[1])
                    strand=line[2]
                    if strand=='+':
                        ar_plus.append(score)
                    elif strand=='-':
                        ar_minus.append(score)
        
            filein.close()
            print(f'{sarus_output} was parsed succesfully')
    
            #Plot Sarus output as metagene plot.
            #Calculate signal over TSS/TES.
            gene_fw=np.array([0.0]*(seq_len))
            gene_rv=np.array([0.0]*(seq_len))   
            
            for gene_name, scan_info in dict_plus_minus.items():
                gene_fw+=np.array(scan_info[0])
                gene_rv+=np.array(scan_info[1])
        
            Num_genes=len(dict_plus_minus)
            gene_fw=gene_fw/Num_genes
            gene_rv=gene_rv/Num_genes
            print(f'Signal over anchor is computed...')

            #Plot signal over anchor. 
            print(f'Plotting signal over anchor...')
            plt.figure(figsize=(13.1, 3), dpi=100)
            plot1=plt.subplot(111)  
            positions_US=np.arange(1,seq_len+1,1)  
            plot1.plot(positions_US, gene_fw, linestyle='-', color='black', linewidth=3, label='fw') 
            plot1.plot(positions_US, gene_rv, linestyle='--', color='black', linewidth=3, label='rv') 
            ticks=np.arange(0,seq_len+1,10).tolist()
            ticks+=[1]
            plot1.set_xticks(ticks)
            ticks_lables=ticks
            plt.axvline(x=win_width)
            plot1.legend(fontsize=12)    
            plot1.set_xlabel('Distance, bp', size=20)
            plot1.set_ylabel(f'{conf_name} {sigma_name} ' + '$\it{E. coli}$ promoters score' + f'({Num_genes})', size=10) 
            plt.xlim([0, seq_len+1])
            plt.savefig(f'{outpath}\{conf_name}_{sigma_name}_sarus_scan_noCTD_Rif_3.svg', dpi=400, figsize=(13.1, 3))   
            plt.close()         
            
            #Create LOGO for a set of promoters.
            Create_logo(f'{outpath}\{conf_name}_{sigma_name}_promoters.fasta', f'{outpath}\{conf_name}_{sigma_name}_promoters_LOGO.pdf')
    
    return


#######
#Returns anchored signal over the set of genes (US, DS) - for each gene separately.
#######

def TUs_and_signal(gene_annotation, genes_set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length):
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

    #Write table contains FE for US, GB, DS of TUs in a set.
    print(f'Writing FE for TUs\' TU, DS...')
    write_genes_FE(gene_US_mean_dict, gene_DS_mean_dict, FE_track_name, f'{out_path}\Signal_of_TSS_TES_tab\\{genes_set_name}\{FE_track_name}_over_{genes_set_name}_{win_width}bp.txt')
 
    return gene_US, gene_DS, gene_US_mean_dict, gene_DS_mean_dict


#######
#Returns anchored sequences over the set of genes (US, DS) - for each gene separately.
#######

def TUs_and_sequence(gene_annotation, genes_set_name, genome_path, out_path, deletions_inpath, win_width, length):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 
    
    #Reading the genome.
    genome_sequence=read_genome(genome_path)

    #Calculate FE over genes.
    gene_US_seq_dict={}
    gene_DS_seq_dict={}
    for gene_name, gene_info in gene_annotation.items():
        delited=0
        for deletion in deletions:
            if deletion[1]>=gene_info[0]>=deletion[0] or deletion[1]>=gene_info[1]>=deletion[0]:
                delited=1
        if delited==0:
            start=gene_info[0]
            end=gene_info[1]
            glen=len(genome_sequence)
            if gene_info[2]=='+':
                if start<win_width:
                    gene_US_seq=genome_sequence[glen-(win_width-start):] + genome_sequence[:start+length]
                    gene_US_seq_dict[gene_name]=[gene_US_seq, start, end, gene_info[2]]
                elif start+length>glen:
                    gene_US_seq=genome_sequence[start-win_width:] + genome_sequence[:start+length-glen]
                    gene_US_seq_dict[gene_name]=[gene_US_seq, start, end, gene_info[2]]                    
                else:
                    gene_US_seq=genome_sequence[start-win_width:start+length]
                    gene_US_seq_dict[gene_name]=[gene_US_seq, start, end, gene_info[2]]
                    
                if end+win_width>glen:
                    gene_DS_seq=genome_sequence[end-length:] + genome_sequence[:end+win_width-glen]
                    gene_DS_seq_dict[gene_name]=[gene_DS_seq, start, end, gene_info[2]]
                elif end-length<0:
                    gene_DS_seq=genome_sequence[glen-(length-end):]+genome_sequence[:end+win_width]
                    gene_DS_seq_dict[gene_name]=[gene_DS_seq, start, end, gene_info[2]]                    
                else:
                    gene_DS_seq=genome_sequence[end-length:end+win_width]
                    gene_DS_seq_dict[gene_name]=[gene_DS_seq, start, end, gene_info[2]]

            elif gene_info[2]=='-':
                if start<win_width:
                    gene_DS_seq=Seq(genome_sequence[:start+length]).reverse_complement() + Seq(genome_sequence[glen-(win_width-start):]).reverse_complement()
                    gene_DS_seq_dict[gene_name]=[gene_DS_seq, start, end, gene_info[2]]
                elif start+length>glen:
                    gene_DS_seq=Seq(genome_sequence[:length-(glen-start)]).reverse_complement() + Seq(genome_sequence[start-win_width:]).reverse_complement()
                    gene_DS_seq_dict[gene_name]=[gene_DS_seq, start, end, gene_info[2]]                    
                else:
                    gene_DS_seq=Seq(genome_sequence[start-win_width:start+length]).reverse_complement()
                    gene_DS_seq_dict[gene_name]=[gene_DS_seq, start, end, gene_info[2]]
                
                if end+win_width>glen:
                    gene_US_seq=Seq(genome_sequence[:end+win_width-glen]).reverse_complement() + Seq(genome_sequence[end-length:]).reverse_complement()  
                    gene_US_seq_dict[gene_name]=[gene_US_seq, start, end, gene_info[2]]
                else:
                    gene_US_seq=Seq(genome_sequence[end-length:end+win_width]).reverse_complement()
                    gene_US_seq_dict[gene_name]=[gene_US_seq, start, end, gene_info[2]]
                    
    Num_genes=len(gene_annotation)
    print(f'Sequences of TUs US and DS are extracting...')

    #Write extracted sequences in fasta format.
    print(f'Writing TUs US, DS sequences...')
    write_sequences_dict(gene_US_seq_dict, win_width, length, f'{out_path}\Sequences_of_TSS_TES_fasta\\{genes_set_name}_{win_width}bp_{length}bp_US.fasta')
    write_sequences_dict(gene_DS_seq_dict, win_width, length, f'{out_path}\Sequences_of_TSS_TES_fasta\\{genes_set_name}_{win_width}bp_{length}bp_DS.fasta')
 
    return gene_US_seq_dict, gene_DS_seq_dict


#######
#Creates logo for US and DS regions.
#######

def Create_logo(alig_inpath, out_path):
    MFA_data=open(alig_inpath)
    MFA_seqs=read_seq_data(MFA_data)
    logodata=LogoData.from_seqs(MFA_seqs)
    logooptions=LogoOptions(yaxis_scale=0.9, pad_right=True, stacks_per_line=90)
    logooptions.show_errorbars=False
    logoformat=LogoFormat(logodata, logooptions)
    pdf=weblogo.logo_formatter.pdf_formatter(logodata, logoformat)
    logout=open(out_path, 'wb')
    logout.write(pdf)
    logout.close()
    
    return


#######
#Wrapper: reads data and gene annotation, computes signal over US, DS, plots and writes the output; 
#extracts sequences at TSS TES, writes the output.
#######

def Wrapper_signal_sequence_at_TSS_TES(tracks_dict, transcripts_annotation, genes_set_name, genome_path, deletions_inpath, win_width, length, out_path):
    #Create output folders.
    create_out_dirs(out_path, genes_set_name)
    
    #Reads input data in wig files.
    dict_of_wigs={}
    for name, data_path in tracks_dict.items():
        dict_of_wigs[name]=wig_parsing(data_path)
    
    #Reads annotation.
    print(f'Now working with {transcripts_annotation}')
    TUs_annotation=parse_expression_annotation(transcripts_annotation)
    
    #Calculate and plot signal over TUs.
    for track_name, signal_track in dict_of_wigs.items():
        TUs_and_signal(TUs_annotation, genes_set_name, signal_track, track_name, out_path, deletions_inpath, win_width, length)
        
    #Extract sequences at TSS and TES.
    TUs_and_sequence(TUs_annotation, genes_set_name, genome_path, out_path, deletions_inpath, win_width, length)
    
    #Create logo.
    
    Create_logo(f'{out_path}\Sequences_of_TSS_TES_fasta\\{genes_set_name}_{win_width}bp_{length}bp_US.fasta', f'{out_path}\Sequences_of_TSS_TES_fasta\\{genes_set_name}_{win_width}bp_{length}bp_US.pdf')
    Create_logo(f'{out_path}\Sequences_of_TSS_TES_fasta\\{genes_set_name}_{win_width}bp_{length}bp_DS.fasta', f'{out_path}\Sequences_of_TSS_TES_fasta\\{genes_set_name}_{win_width}bp_{length}bp_DS.pdf')
    
    return

#Wrapper_signal_sequence_at_TSS_TES(Scores_dict, Transcripts_annotation, 'All_rep_TUs_fw', Genome_path, Deletions_inpath, 70, 30, Outpath)
#Read_sarus_plot("C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Promoter_structure\Try1\Sequences_of_TSS_TES__sarus_scan_noCTD_Rif_3_T_F\\TopoI_HS_2_TUs_70bp_30bp_US_sarus_scan_noCTD_Rif_3_T_F.txt", 70, 30, 'TopoI_HS_2_TUs_US_cooriented', Outpath)
Select_promoters_run_sarus_plot("C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Promoter_structure\RegulonDB_promoters.txt", 60, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Promoter_structure\Try3_regulonDB_promoters")