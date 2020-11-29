###############################################
##Dmitry Sutormin, 2020##
##TopoA ChIP-Seq analysis##

#Takes sets of transcription units (TUs) and calculate EcTopoI signal (wig) in upstreams, TU bodies and at TU starts.
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
import scipy
from scipy import stats


#Path to TUs data.
TUs_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_transcripts\\Representative_transcripts\\DY330_RNA-Seq_transcripts_representative_NO_DPS_NO_RFA_no_tRNA_rRNA_EP_del_cor_HETU_200.txt"

#Dict of path with signal data (wig).
Signal_wig_path_dict={'TopA_CTD_minus_Rif_minus_av_1_2_3': 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\\Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_123.wig',
                      'TopA_CTD_minus_Rif_plus_av_1_2_3':  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\\Sutormin_TopA_ChIP_CTD_minus_Rif_plus_FE_av_123.wig',
                      'TopA_CTD_plus_Rif_minus_av_1_2_3':  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\\Sutormin_TopA_ChIP_CTD_plus_Rif_minus_FE_av.wig',
                      'TopA_CTD_plus_Rif_plus_av_2_3':     'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\\Sutormin_TopA_ChIP_CTD_plus_Rif_plus_FE_av.wig',}

#Intervals: define US, TSS, GB borders in bp.
US_length=12000
TSS_halfwidth=700

#Output path.
Outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Signal_over_TUs\Representative_transcripts\Figures\Signal_US_TSS_GB\Rif_CTD_effects_4kb_half_distance_estimation\\"



#################################
#################################
######### Part 1.
######### Compare average FE of TUS, TSS, TUB.
#################################
#################################


#######
#Parses WIG file.
#######

def wig_parsing(wig_path_dict):

    wig_data_dict={}
    for name, wigfile in wig_path_dict.items():
        print('Now is processing: ' + str(wigfile))
        wigin=open(wigfile, 'r')
        NE_values=[]
        for line in wigin:
            line=line.rstrip().split(' ')
            if line[0] not in ['track', 'fixedStep']:
                NE_values.append(float(line[0]))
        wigin.close()
        
        wig_data_dict[name]=NE_values
    return wig_data_dict


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
#Genome binning with average peak width, compute distribution of GC% of bins.
#######

def genome_bin_mean_std_signal(genes_annotation, wig_data, TSS_halfwidth, US_width):
    
    #Mask deletions.
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]]
    mask=[0]*len(wig_data)
    for delition in deletions:
        mask[delition[0]:delition[1]]=[1]*(delition[1]-delition[0])  
    
    wig_data_masked=[]
    for i in range(len(mask)):
        if mask[i]==0:
            wig_data_masked.append(wig_data[i])
    
    #Calculate median length of TUBs (width of a bin).
    bins_width_TUB=[]
    for TU_name, TU_data in genes_annotation.items():
        TU_start=TU_data[0]
        TU_end=TU_data[1]
        TU_body=(TU_end-TU_start)-TSS_halfwidth
        if TU_body>0:
            bins_width_TUB.append(TU_body)
    
    median_len_TUB=int(np.median(bins_width_TUB))-1
    
    #Bin genome by TU body (TUB) length and keep signal of bins.
    number_of_genome_bin_TUB=int(len(wig_data_masked)/median_len_TUB)-1
    genome_signal_values_TUB=[]
    for i in range(number_of_genome_bin_TUB):
        genome_signal_values_TUB.append(np.mean(wig_data_masked[median_len_TUB*i:median_len_TUB*(i+1)]))
    
    Signal_mean_genome=np.mean(wig_data_masked)
    Signal_std_genome=np.std(wig_data_masked)
    
    Signal_mean_genome_binned_TUB=np.mean(genome_signal_values_TUB)
    Signal_std_genome_binned_TUB=np.std(genome_signal_values_TUB)    
    
    print(f'Median width of TUB bins: {median_len_TUB}')
    print(f'Global mean: {Signal_mean_genome}; Global std: {Signal_std_genome}') 
    print(f'TUB binned mean: {Signal_mean_genome_binned_TUB}; TUB binned std: {Signal_std_genome_binned_TUB}') 
    
    TUB_estimators=[Signal_mean_genome, Signal_std_genome, Signal_mean_genome_binned_TUB, Signal_std_genome_binned_TUB]
    
    #Bin genome by transcription start site (TSS) length and keep signal of bins.
    number_of_genome_bin_TSS=int(len(wig_data_masked)/(TSS_halfwidth*2))-1
    genome_signal_values_TSS=[]
    for i in range(number_of_genome_bin_TSS):
        genome_signal_values_TSS.append(np.mean(wig_data_masked[(TSS_halfwidth*2)*i:(TSS_halfwidth*2)*(i+1)]))
    
    Signal_mean_genome_binned_TSS=np.mean(genome_signal_values_TSS)
    Signal_std_genome_binned_TSS=np.std(genome_signal_values_TSS)    
    
    print(f'Median width of TSS bins: {2*TSS_halfwidth}')
    print(f'TSS binned mean: {Signal_mean_genome_binned_TSS}; TSS binned std: {Signal_std_genome_binned_TSS}')  
    
    TSS_estimators=[Signal_mean_genome, Signal_std_genome, Signal_mean_genome_binned_TSS, Signal_std_genome_binned_TSS]    
    
    #Bin genome by TU upstream (TUU) length and keep signal of bins.
    number_of_genome_bin_TUU=int(len(wig_data_masked)/(US_width))-1
    genome_signal_values_TUU=[]
    for i in range(number_of_genome_bin_TUU):
        genome_signal_values_TUU.append(np.mean(wig_data_masked[(US_width)*i:(US_width)*(i+1)]))
    
    Signal_mean_genome_binned_TUU=np.mean(genome_signal_values_TUU)
    Signal_std_genome_binned_TUU=np.std(genome_signal_values_TUU)    
    
    print(f'Median width of TUU bins: {US_width}')
    print(f'TUU binned mean: {Signal_mean_genome_binned_TUU}; TUU binned std: {Signal_std_genome_binned_TUU}')  
    
    TUU_estimators=[Signal_mean_genome, Signal_std_genome, Signal_mean_genome_binned_TUU, Signal_std_genome_binned_TUU]        
    
    return TUB_estimators, TSS_estimators, TUU_estimators


#######
#Fetch relative signal from US, TSS, GB.
#Relative signal X-mean(X).
#######

def fetch_signal_relative(wig_data_dict, genes_annotation, US_width, TSS_halfwidth):
    
    Signal_dict={}
    
    for wig_name, wig_data in wig_data_dict.items():
        genome_len=4647454
        
        #Try using estimators calculated for binned genome (like for GC in Return_sequences_under_peaks...)
        TUB_estimators, TSS_estimators, TUU_estimators=genome_bin_mean_std_signal(genes_annotation, wig_data, TSS_halfwidth, US_width)
        
        WIG_mean_TUB=TUB_estimators[2]
        WIG_std_TUB=TUB_estimators[3]
        WIG_mean_TSS=TSS_estimators[2]
        WIG_std_TSS=TSS_estimators[3]
        WIG_mean_TUU=TUU_estimators[2]
        WIG_std_TUU=TUU_estimators[3]        
    
        US_signal_ar=[]
        TSS_signal_ar=[]
        GB_signal_ar=[]
        
        for TU_name, TU_data in genes_annotation.items():
            TU_start=TU_data[0]
            TU_end=TU_data[1]
            TU_strand=TU_data[2]
            
            if TU_strand=="+":
                #print(TU_start-TSS_halfwidth-US_width, TU_start-TSS_halfwidth, TU_start+TSS_halfwidth, TU_end)
                
                if (TU_start-TSS_halfwidth-US_width>0) & (TU_start-TSS_halfwidth>0):
                    US_signal=(np.mean(wig_data[TU_start-TSS_halfwidth-US_width:TU_start-TSS_halfwidth]) - WIG_mean_TUU)
                    TSS_signal=(np.mean(wig_data[TU_start-TSS_halfwidth:TU_start+TSS_halfwidth]) - WIG_mean_TSS)
                elif (TU_start-TSS_halfwidth-US_width<0) & (TU_start-TSS_halfwidth>0):
                    US_signal=(np.mean(wig_data[TU_start-TSS_halfwidth-US_width:]+wig_data[:TU_start-TSS_halfwidth]) - WIG_mean_TUU)
                    TSS_signal=(np.mean(wig_data[TU_start-TSS_halfwidth:TU_start+TSS_halfwidth]) - WIG_mean_TSS)
                elif (TU_start-TSS_halfwidth-US_width<0) & (TU_start-TSS_halfwidth<0):
                    US_signal=(np.mean(wig_data[TU_start-TSS_halfwidth-US_width:TU_start-TSS_halfwidth]) - WIG_mean_TUU)
                    TSS_signal=(np.mean(wig_data[TU_start-TSS_halfwidth:]+wig_data[:TU_start]) - WIG_mean_TSS)
                        
                if TU_start+TSS_halfwidth<TU_end:
                    GB_signal=(np.mean(wig_data[TU_start+TSS_halfwidth:TU_end]) - WIG_mean_TUB)
                elif TU_start+TSS_halfwidth>=TU_end:
                    GB_signal=-1
                    
            elif TU_strand=="-":
                #print(TU_start, TU_end-TSS_halfwidth, TU_end+TSS_halfwidth, TU_end+TSS_halfwidth+US_width)
                
                US_signal=(np.mean(wig_data[TU_end+TSS_halfwidth:TU_end+TSS_halfwidth+US_width]) - WIG_mean_TUU)
                TSS_signal=(np.mean(wig_data[TU_end-TSS_halfwidth:TU_end+TSS_halfwidth]) - WIG_mean_TSS)
                if TU_end-TSS_halfwidth>TU_start:
                    GB_signal=(np.mean(wig_data[TU_start:TU_end-TSS_halfwidth]) - WIG_mean_TUB)
                else:
                    GB_signal=-1
                    
            US_signal_ar.append(US_signal)
            TSS_signal_ar.append(TSS_signal)
            if GB_signal!=-1:
                GB_signal_ar.append(GB_signal)
        
        Signal_dict[wig_name]=[US_signal_ar, TSS_signal_ar, GB_signal_ar]
    
    return Signal_dict


#######
#Fetch relative signal from US, TSS, GB.
#Normalized signal (X-mean(X))/std(X).
#######

def fetch_signal_normalized(wig_data_dict, genes_annotation, US_width, TSS_halfwidth):
    
    Signal_dict={}
    
    for wig_name, wig_data in wig_data_dict.items():
        genome_len=4647454
        
        #Try using estimators calculated for binned genome (like for GC in Return_sequences_under_peaks...)
        TUB_estimators, TSS_estimators, TUU_estimators=genome_bin_mean_std_signal(genes_annotation, wig_data, TSS_halfwidth, US_width)
        
        WIG_mean_TUB=TUB_estimators[2]
        WIG_std_TUB=TUB_estimators[3]
        WIG_mean_TSS=TSS_estimators[2]
        WIG_std_TSS=TSS_estimators[3]
        WIG_mean_TUU=TUU_estimators[2]
        WIG_std_TUU=TUU_estimators[3]        
    
        US_signal_ar=[]
        TSS_signal_ar=[]
        GB_signal_ar=[]
        
        for TU_name, TU_data in genes_annotation.items():
            TU_start=TU_data[0]
            TU_end=TU_data[1]
            TU_strand=TU_data[2]
            
            if TU_strand=="+":
                #print(TU_start-TSS_halfwidth-US_width, TU_start-TSS_halfwidth, TU_start+TSS_halfwidth, TU_end)
                
                if (TU_start-TSS_halfwidth-US_width>0) & (TU_start-TSS_halfwidth>0):
                    US_signal=(np.mean(wig_data[TU_start-TSS_halfwidth-US_width:TU_start-TSS_halfwidth]) - WIG_mean_TUU)/WIG_std_TUU
                    TSS_signal=(np.mean(wig_data[TU_start-TSS_halfwidth:TU_start+TSS_halfwidth]) - WIG_mean_TSS)/WIG_std_TSS
                elif (TU_start-TSS_halfwidth-US_width<0) & (TU_start-TSS_halfwidth>0):
                    US_signal=(np.mean(wig_data[TU_start-TSS_halfwidth-US_width:]+wig_data[:TU_start-TSS_halfwidth]) - WIG_mean_TUU)/WIG_std_TUU
                    TSS_signal=(np.mean(wig_data[TU_start-TSS_halfwidth:TU_start+TSS_halfwidth]) - WIG_mean_TSS)/WIG_std_TSS
                elif (TU_start-TSS_halfwidth-US_width<0) & (TU_start-TSS_halfwidth<0):
                    US_signal=(np.mean(wig_data[TU_start-TSS_halfwidth-US_width:TU_start-TSS_halfwidth]) - WIG_mean_TUU)/WIG_std_TUU
                    TSS_signal=(np.mean(wig_data[TU_start-TSS_halfwidth:]+wig_data[:TU_start]) - WIG_mean_TSS)/WIG_std_TSS
                        
                if TU_start+TSS_halfwidth<TU_end:
                    GB_signal=(np.mean(wig_data[TU_start+TSS_halfwidth:TU_end]) - WIG_mean_TUB)/WIG_std_TUB
                elif TU_start+TSS_halfwidth>=TU_end:
                    GB_signal=-1
                    
            elif TU_strand=="-":
                US_signal=(np.mean(wig_data[TU_end+TSS_halfwidth:TU_end+TSS_halfwidth+US_width]) - WIG_mean_TUU)/WIG_std_TUU
                TSS_signal=(np.mean(wig_data[TU_end-TSS_halfwidth:TU_end+TSS_halfwidth]) - WIG_mean_TSS)/WIG_std_TSS
                if TU_end-TSS_halfwidth>TU_start:
                    GB_signal=(np.mean(wig_data[TU_start:TU_end-TSS_halfwidth]) - WIG_mean_TUB)/WIG_std_TUB
                else:
                    GB_signal=-1
                    
            US_signal_ar.append(US_signal)
            TSS_signal_ar.append(TSS_signal)
            if GB_signal!=-1:
                GB_signal_ar.append(GB_signal)
        
        Signal_dict[wig_name]=[US_signal_ar, TSS_signal_ar, GB_signal_ar]
    
    return Signal_dict


#######
#Fetch absolute signal from TSS, GB.
#Absolute signal.
#######

def fetch_signal_abs(wig_data_dict, genes_annotation, US_width, TSS_halfwidth):
    
    Signal_dict={}
    
    for wig_name, wig_data in wig_data_dict.items():
        genome_len=4647454
    
        US_signal_ar=[]
        TSS_signal_ar=[]
        GB_signal_ar=[]
        
        for TU_name, TU_data in genes_annotation.items():
            TU_start=TU_data[0]
            TU_end=TU_data[1]
            TU_strand=TU_data[2]
            
            if TU_strand=="+":
                #print(TU_start-TSS_halfwidth-US_width, TU_start-TSS_halfwidth, TU_start+TSS_halfwidth, TU_end)
                
                if (TU_start-TSS_halfwidth-US_width>0) & (TU_start-TSS_halfwidth>0):
                    US_signal=np.mean(wig_data[TU_start-TSS_halfwidth-US_width:TU_start-TSS_halfwidth])
                    TSS_signal=np.mean(wig_data[TU_start-TSS_halfwidth:TU_start+TSS_halfwidth])
                elif (TU_start-TSS_halfwidth-US_width<0) & (TU_start-TSS_halfwidth>0):
                    US_signal=np.mean(wig_data[TU_start-TSS_halfwidth-US_width:]+wig_data[:TU_start-TSS_halfwidth])
                    TSS_signal=np.mean(wig_data[TU_start-TSS_halfwidth:TU_start+TSS_halfwidth])
                elif (TU_start-TSS_halfwidth-US_width<0) & (TU_start-TSS_halfwidth<0):
                    US_signal=np.mean(wig_data[TU_start-TSS_halfwidth-US_width:TU_start-TSS_halfwidth])
                    TSS_signal=np.mean(wig_data[TU_start-TSS_halfwidth:]+wig_data[:TU_start])
                        
                if TU_start+TSS_halfwidth<TU_end:
                    GB_signal=np.mean(wig_data[TU_start+TSS_halfwidth:TU_end])
                elif TU_start+TSS_halfwidth>=TU_end:
                    GB_signal=-1
                    
            elif TU_strand=="-":
                #print(TU_start, TU_end-TSS_halfwidth, TU_end+TSS_halfwidth, TU_end+TSS_halfwidth+US_width)
                
                US_signal=np.mean(wig_data[TU_end+TSS_halfwidth:TU_end+TSS_halfwidth+US_width])
                TSS_signal=np.mean(wig_data[TU_end-TSS_halfwidth:TU_end+TSS_halfwidth])
                if TU_end-TSS_halfwidth>TU_start:
                    GB_signal=np.mean(wig_data[TU_start:TU_end-TSS_halfwidth])
                else:
                    GB_signal=-1
                    
            US_signal_ar.append(US_signal)
            TSS_signal_ar.append(TSS_signal)
            if GB_signal!=-1:
                GB_signal_ar.append(GB_signal)
        
        Signal_dict[wig_name]=[US_signal_ar, TSS_signal_ar, GB_signal_ar]

    return Signal_dict


#######
#Plot data.
#######

#For relative data X-mean(X).
def plot_proportions_US_TSS_GB_rel(Signal_dict, outpath):
    
    #Compare US, TSS, TUB.
    
    fig, plot_av=plt.subplots(1,1,figsize=(8,3), dpi=100)

    #EcTopoI vs EcTopoI Rif+
    Ec_TopoI=[np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][0]), np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif=[np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][0]), np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][2])]
    Ec_TopoI_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][0]), np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][0]), np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][2])]

    
    Conditions=['US', 'TSS', 'TUB', 'US', 'TSS', 'TUB', 'US', 'TSS', 'TUB', 'US', 'TSS', 'TUB']
    
    X_coords_1=[1,2,3]
    X_coords_2=[5,6,7]
    X_coords_3=[9,10,11]
    X_coords_4=[13,14,15]
    X_coords=[1,2,3, 5,6,7, 9,10,11, 13,14,15]
      
    plot_av.bar(X_coords_1, Ec_TopoI, width=0.9, color='#ff6f8e', edgecolor='k', linewidth=0.6, label='EcTopoI')
    plot_av.bar(X_coords_2, Ec_TopoI_Rif, width=0.9, color='#75595e', edgecolor='k', linewidth=0.6, label='EcTopoI Rif')
    plot_av.bar(X_coords_3, Ec_TopoI_CTD, width=0.9, color='#6be1ff', edgecolor='k', linewidth=0.6, label='EcTopoI CTD')
    plot_av.bar(X_coords_4, Ec_TopoI_Rif_CTD, width=0.9, color='#ffcab8', edgecolor='k', linewidth=0.6, label='EcTopoI Rif CTD')
    plot_av.set_ylabel('Relative FE (X-' + r"$\overline{X}$" + ')', size=16)
    plot_av.set_xticks(X_coords)
    plot_av.set_xticklabels(Conditions, rotation=0, size=7)     
    plt.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.5, 1))
    plt.tight_layout(rect=[0,0,1,1])   
    plt.show()
    plt.savefig(outpath + 'HETU_EcTopoI_Rif_CTD_effects_relative_barplot_US_12000_TSS_700_GB.png', dpi=300, size=(8,3))    
    
    
    #Compare TSS, TUB.
    fig, plot_av=plt.subplots(1,1,figsize=(8,3), dpi=100)
    
    Ec_TopoI=[np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif=[np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][2])]
    Ec_TopoI_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][2])]    
    
    Conditions=['TSS', 'TUB', 'TSS', 'TUB', 'TSS', 'TUB', 'TSS', 'TUB']
    
    X_coords_1=[1,2]
    X_coords_2=[4,5]
    X_coords_3=[7,8]
    X_coords_4=[10,11]
    X_coords=[1,2, 4,5, 7,8, 10,11]
      
    plot_av.bar(X_coords_1, Ec_TopoI, width=0.9, color='#ff6f8e', edgecolor='k', linewidth=0.6, label='EcTopoI')
    plot_av.bar(X_coords_2, Ec_TopoI_Rif, width=0.9, color='#75595e', edgecolor='k', linewidth=0.6, label='EcTopoI Rif')
    plot_av.bar(X_coords_3, Ec_TopoI_CTD, width=0.9, color='#6be1ff', edgecolor='k', linewidth=0.6, label='EcTopoI CTD')
    plot_av.bar(X_coords_4, Ec_TopoI_Rif_CTD, width=0.9, color='#ffcab8', edgecolor='k', linewidth=0.6, label='EcTopoI Rif CTD')
    plot_av.set_ylabel('Relative FE (X-' + r"$\overline{X}$" + ')', size=16)
    plot_av.set_xticks(X_coords)
    plot_av.set_xticklabels(Conditions, rotation=0, size=7)     
    plt.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.5, 1))
    plt.tight_layout(rect=[0,0,1,1]) 
    plt.show()
    plt.savefig(outpath + 'HETU_EcTopoI_Rif_CTD_effects_relative_barplot_TSS_700_GB.png', dpi=300, size=(8,3))  
    
    #Test whether differences between dataset are significant.
    
    #Welch t-test. CTD-Rif- vs CTD-Rif+
    #TSS
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1], Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][1], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1])}, Sample size CTD-/Rif+: {len(Signal_dict["TopA_CTD_minus_Rif_plus_av_1_2_3"][1])}')
    print(f'\nT-test TSS RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1]),2)} CTD-/Rif+={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_plus_av_1_2_3"][1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    #TUB
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2], Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][2], equal_var=False)
    print(f'\nT-test TUB- RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][2]),2)} CTD-/Rif+={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_plus_av_1_2_3"][2]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    #Welch t-test. CTD-Rif- vs CTD+Rif-
    #TSS
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1], Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][1], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1])}, Sample size CTD+/Rif-: {len(Signal_dict["TopA_CTD_plus_Rif_minus_av_1_2_3"][1])}')
    print(f'\nT-test TSS RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1]),2)} CTD+/Rif-={round(np.mean(Signal_dict["TopA_CTD_plus_Rif_minus_av_1_2_3"][1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    #TUB
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2], Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][2], equal_var=False)
    print(f'\nT-test TUB RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][2]),2)} CTD+/Rif-={round(np.mean(Signal_dict["TopA_CTD_plus_Rif_minus_av_1_2_3"][2]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    #Welch t-test. CTD-Rif- vs CTD+Rif+
    #TSS
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1], Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][1], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1])}, Sample size CTD+/Rif+: {len(Signal_dict["TopA_CTD_plus_Rif_plus_av_2_3"][1])}')
    print(f'\nT-test TSS RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1]),2)} CTD+/Rif+={round(np.mean(Signal_dict["TopA_CTD_plus_Rif_plus_av_2_3"][1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    #TUB
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2], Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][2], equal_var=False)
    print(f'\nT-test TUB RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][2]),2)} CTD+/Rif+={round(np.mean(Signal_dict["TopA_CTD_plus_Rif_plus_av_2_3"][2]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    return

#For normalized data (X-mean(X))/std(X).
def plot_proportions_US_TSS_GB_norm(Signal_dict, outpath):
    
    #Compare US, TSS, TUB.
    
    fig, plot_av=plt.subplots(1,1,figsize=(8,3), dpi=100)

    #EcTopoI vs EcTopoI Rif+
    Ec_TopoI=[np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][0]), np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif=[np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][0]), np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][2])]
    Ec_TopoI_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][0]), np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][0]), np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][2])]

    
    Conditions=['US', 'TSS', 'TUB', 'US', 'TSS', 'TUB', 'US', 'TSS', 'TUB', 'US', 'TSS', 'TUB']
    
    X_coords_1=[1,2,3]
    X_coords_2=[5,6,7]
    X_coords_3=[9,10,11]
    X_coords_4=[13,14,15]
    X_coords=[1,2,3, 5,6,7, 9,10,11, 13,14,15]
      
    plot_av.bar(X_coords_1, Ec_TopoI, width=0.9, color='#ff6f8e', edgecolor='k', linewidth=0.6, label='EcTopoI')
    plot_av.bar(X_coords_2, Ec_TopoI_Rif, width=0.9, color='#75595e', edgecolor='k', linewidth=0.6, label='EcTopoI Rif')
    plot_av.bar(X_coords_3, Ec_TopoI_CTD, width=0.9, color='#6be1ff', edgecolor='k', linewidth=0.6, label='EcTopoI CTD')
    plot_av.bar(X_coords_4, Ec_TopoI_Rif_CTD, width=0.9, color='#ffcab8', edgecolor='k', linewidth=0.6, label='EcTopoI Rif CTD')
    plot_av.set_ylabel('Normalized FE (X-' + r"$\overline{X}$" + ')/'  + r"$\sigma$", size=16)
    plot_av.set_xticks(X_coords)
    plot_av.set_xticklabels(Conditions, rotation=0, size=7)     
    plt.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.5, 1))
    plt.tight_layout(rect=[0,0,1,1]) 
    plt.show()
    plt.savefig(outpath + 'HETU_EcTopoI_Rif_CTD_effects_normalized_barplot_US_12000_TSS_700_GB.png', dpi=300, size=(8,3))    
    
    
    #Compare TSS, TUB.
    fig, plot_av=plt.subplots(1,1,figsize=(8,3), dpi=100)
    
    Ec_TopoI=[np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif=[np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][2])]
    Ec_TopoI_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][2])]    
    
    Conditions=['TSS', 'TUB', 'TSS', 'TUB', 'TSS', 'TUB', 'TSS', 'TUB']
    
    X_coords_1=[1,2]
    X_coords_2=[4,5]
    X_coords_3=[7,8]
    X_coords_4=[10,11]
    X_coords=[1,2, 4,5, 7,8, 10,11]
      
    plot_av.bar(X_coords_1, Ec_TopoI, width=0.9, color='#ff6f8e', edgecolor='k', linewidth=0.6, label='EcTopoI')
    plot_av.bar(X_coords_2, Ec_TopoI_Rif, width=0.9, color='#75595e', edgecolor='k', linewidth=0.6, label='EcTopoI Rif')
    plot_av.bar(X_coords_3, Ec_TopoI_CTD, width=0.9, color='#6be1ff', edgecolor='k', linewidth=0.6, label='EcTopoI CTD')
    plot_av.bar(X_coords_4, Ec_TopoI_Rif_CTD, width=0.9, color='#ffcab8', edgecolor='k', linewidth=0.6, label='EcTopoI Rif CTD')
    plot_av.set_ylabel('Normalized FE (X-' + r"$\overline{X}$" + ')/'  + r"$\sigma$", size=16)
    plot_av.set_xticks(X_coords)
    plot_av.set_xticklabels(Conditions, rotation=0, size=7)     
    plt.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.5, 1))
    plt.tight_layout(rect=[0,0,1,1]) 
    plt.show()
    plt.savefig(outpath + 'HETU_EcTopoI_Rif_CTD_effects_normalized_barplot_TSS_700_GB.png', dpi=300, size=(8,3))  
    
    #Test whether differences between dataset are significant.
    
    #Welch t-test. CTD-Rif- vs CTD-Rif+
    #TSS
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1], Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][1], equal_var=True)
    print(f'Sample size CTD-/Rif-: {len(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1])}, Sample size CTD-/Rif+: {len(Signal_dict["TopA_CTD_minus_Rif_plus_av_1_2_3"][1])}')
    print(f'\nT-test TSS RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1]),2)} CTD-/Rif+={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_plus_av_1_2_3"][1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    #TUB
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2], Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][2], equal_var=True)
    print(f'\nT-test TUB RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][2]),2)} CTD-/Rif+={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_plus_av_1_2_3"][2]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    #Welch t-test. CTD-Rif- vs CTD+Rif-
    #TSS
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1], Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][1], equal_var=True)
    print(f'Sample size CTD-/Rif-: {len(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1])}, Sample size CTD+/Rif-: {len(Signal_dict["TopA_CTD_plus_Rif_minus_av_1_2_3"][1])}')
    print(f'\nT-test TSS RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1]),2)} CTD+/Rif-={round(np.mean(Signal_dict["TopA_CTD_plus_Rif_minus_av_1_2_3"][1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    #TUB
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2], Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][2], equal_var=True)
    print(f'\nT-test TUB RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][2]),2)} CTD+/Rif-={round(np.mean(Signal_dict["TopA_CTD_plus_Rif_minus_av_1_2_3"][2]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    #Welch t-test. CTD-Rif- vs CTD+Rif+
    #TSS
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1], Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][1], equal_var=True)
    print(f'Sample size CTD-/Rif-: {len(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1])}, Sample size CTD+/Rif+: {len(Signal_dict["TopA_CTD_plus_Rif_plus_av_2_3"][1])}')
    print(f'\nT-test TSS RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][1]),2)} CTD+/Rif+={round(np.mean(Signal_dict["TopA_CTD_plus_Rif_plus_av_2_3"][1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    #TUB
    Intervals_stat=stats.ttest_ind(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2], Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][2], equal_var=True)
    print(f'\nT-test TUB RFE means\nCTD-/Rif-={round(np.mean(Signal_dict["TopA_CTD_minus_Rif_minus_av_1_2_3"][2]),2)} CTD+/Rif+={round(np.mean(Signal_dict["TopA_CTD_plus_Rif_plus_av_2_3"][2]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    
    return


#For absolute data X.
def plot_proportions_TSS_GB(Signal_dict, outpath):
    
    fig, plot_av=plt.subplots(1,1,figsize=(8,3), dpi=100)

    #EcTopoI vs EcTopoI Rif+
    Ec_TopoI=[np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif=[np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'][2])]
    Ec_TopoI_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'][2])]
    Ec_TopoI_Rif_CTD=[np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][1]), np.mean(Signal_dict['TopA_CTD_plus_Rif_plus_av_2_3'][2])]

    
    Conditions=['TSS', 'TUB', 'TSS', 'TUB', 'TSS', 'TUB', 'TSS', 'TUB']
    
    X_coords_1=[1,2,]
    X_coords_2=[4,5]
    X_coords_3=[7,8]
    X_coords_4=[10,11]
    X_coords=[1,2, 4,5, 7,8, 10,11]
      
    plot_av.bar(X_coords_1, Ec_TopoI, width=0.9, color='#ff6f8e', edgecolor='k', linewidth=0.6, label='EcTopoI')
    plot_av.bar(X_coords_2, Ec_TopoI_Rif, width=0.9, color='#75595e', edgecolor='k', linewidth=0.6, label='EcTopoI Rif')
    plot_av.bar(X_coords_3, Ec_TopoI_CTD, width=0.9, color='#6be1ff', edgecolor='k', linewidth=0.6, label='EcTopoI CTD')
    plot_av.bar(X_coords_4, Ec_TopoI_Rif_CTD, width=0.9, color='#ffcab8', edgecolor='k', linewidth=0.6, label='EcTopoI Rif CTD')
    plot_av.set_ylabel('EcTopoI FE', size=16)
    plot_av.set_xticks(X_coords)
    plot_av.set_xticklabels(Conditions, rotation=0, size=7)     
    plt.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.5, 1))
    plt.tight_layout(rect=[0,0,1,1]) 
    plt.show()
    plt.savefig(outpath + 'HETU_EcTopoI_Rif_CTD_effects_absolute_barplot_TSS_700_GB.png', dpi=300, size=(8,3))    
    
    #Make pie charts.
    fig, plot_av=plt.subplots(1,4,figsize=(9,3), dpi=100)
    
    labels=f'TSS', f'TUB'
    colors=['lightcoral', 'lightskyblue']
    explode=(0.1, 0)  # explode 1st slice
    plot_av[0].pie(Ec_TopoI, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=-70)  
    plot_av[1].pie(Ec_TopoI_Rif, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=-70)
    plot_av[2].pie(Ec_TopoI_CTD, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=-70)
    plot_av[3].pie(Ec_TopoI_Rif_CTD, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=-70)
    plt.show()
    plt.savefig(outpath + 'HETU_EcTopoI_Rif_CTD_effects_absolute_piechart_TSS_700_GB.png', dpi=300, size=(9,3))       
    
    return

WIG_data_dict=wig_parsing(Signal_wig_path_dict)
TUs_dictionary=parse_expression_annotation(TUs_data_path)
#Relative_Signal_dictionary=fetch_signal_relative(WIG_data_dict, TUs_dictionary, US_length, TSS_halfwidth)
#Normalized_Signal_dictionary=fetch_signal_normalized(WIG_data_dict, TUs_dictionary, US_length, TSS_halfwidth)
#Absolute_Signal_dictionary=fetch_signal_abs(WIG_data_dict, TUs_dictionary, US_length, TSS_halfwidth)


#plot_proportions_US_TSS_GB_rel(Relative_Signal_dictionary, Outpath)
#plot_proportions_US_TSS_GB_norm(Normalized_Signal_dictionary, Outpath)
#plot_proportions_TSS_GB(Absolute_Signal_dictionary, Outpath)



#################################
#################################
######### Part 2.
######### Compute dependence of FE from the distance from TSS. Do statistics.
#################################
#################################


#######
#Compute half-distance (distance from TSS, where signal reaches one-half of the initial signal (at the TSS)).
#######

def compute_half_distance(wig_data_dict, genes_annotation):
    
    Half_data_dict={}
    
    for wig_name, wig_data in wig_data_dict.items():
        
        print(f'Now on {wig_name}')
        
        Half_signal_ar=[]
        Half_distance_ar=[]
        
        for TU_name, TU_data in genes_annotation.items():  
            TU_start=TU_data[0]
            TU_end=TU_data[1]
            TU_strand=TU_data[2]
            
            if TU_strand=="+":
                        
                GB=wig_data[TU_start:TU_end]
                    
            elif TU_strand=="-":
                
                GB=wig_data[TU_start:TU_end][::-1]            
    
            Half_signal=((GB[0]-1)/2)+1
            Half_distance=0
            for position in range(len(GB)):
                if GB[position]<Half_signal:
                    Half_distance=position
                    break  
            
            if Half_distance!=0:
                Half_signal_ar.append(Half_signal)
                Half_distance_ar.append(Half_distance)  
        
        Half_data_dict[wig_name]=[Half_signal_ar, Half_distance_ar]
        
    HD_mm=Half_data_dict['TopA_CTD_minus_Rif_minus_av_1_2_3'] 
    HD_mp=Half_data_dict['TopA_CTD_minus_Rif_plus_av_1_2_3'] 
    HD_pm=Half_data_dict['TopA_CTD_plus_Rif_minus_av_1_2_3'] 
    HD_pp=Half_data_dict['TopA_CTD_plus_Rif_plus_av_2_3'] 
    
    #Welch t-test. CTD-Rif- vs CTD-Rif+
    #Half signal.
    Intervals_stat=stats.ttest_ind(HD_mm[0], HD_mp[0], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(HD_mm[0])}, Sample size CTD-/Rif+: {len(HD_mp[0])}')
    print(f'\nT-test half signal means\nCTD-/Rif-={round(np.mean(HD_mm[0]),2)} CTD-/Rif+={round(np.mean(HD_mp[0]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')  
    #Half distance.
    Intervals_stat=stats.ttest_ind(HD_mm[1], HD_mp[1], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(HD_mm[1])}, Sample size CTD-/Rif+: {len(HD_mp[1])}')
    print(f'\nT-test half signal means\nCTD-/Rif-={round(np.mean(HD_mm[1]),2)} CTD-/Rif+={round(np.mean(HD_mp[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    
    
    #Welch t-test. CTD-Rif- vs CTD+Rif-
    #Half signal.
    Intervals_stat=stats.ttest_ind(HD_mm[0], HD_pm[0], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(HD_mm[0])}, Sample size CTD+/Rif-: {len(HD_pm[0])}')
    print(f'\nT-test half signal means\nCTD-/Rif-={round(np.mean(HD_mm[0]),2)} CTD+/Rif-={round(np.mean(HD_pm[0]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')  
    #Half distance.
    Intervals_stat=stats.ttest_ind(HD_mm[1], HD_pm[1], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(HD_mm[1])}, Sample size CTD+/Rif-: {len(HD_pm[1])}')
    print(f'\nT-test half signal means\nCTD-/Rif-={round(np.mean(HD_mm[1]),2)} CTD+/Rif-={round(np.mean(HD_pm[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
    
    #Welch t-test. CTD-Rif- vs CTD+Rif+
    #Half signal.
    Intervals_stat=stats.ttest_ind(HD_mm[0], HD_pp[0], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(HD_mm[0])}, Sample size CTD+/Rif+: {len(HD_pp[0])}')
    print(f'\nT-test half signal means\nCTD-/Rif-={round(np.mean(HD_mm[0]),2)} CTD+/Rif+={round(np.mean(HD_pp[0]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')  
    #Half distance.
    Intervals_stat=stats.ttest_ind(HD_mm[1], HD_pp[1], equal_var=False)
    print(f'Sample size CTD-/Rif-: {len(HD_mm[1])}, Sample size CTD+/Rif+: {len(HD_pp[1])}')
    print(f'\nT-test half signal means\nCTD-/Rif-={round(np.mean(HD_mm[1]),2)} CTD+/Rif+={round(np.mean(HD_pp[1]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')     
    
    
    return Half_data_dict


#######
#Take FE values from some distance from TSS going into TU body.
#Compute mean and std, plot for several conditions (e.g. -CTD vs +CTD).
#######

def TSS_distance_FE(wig_data_dict, genes_annotation):
    
    Signal_dict={}
    
    for wig_name, wig_data in wig_data_dict.items():
        
        print(f'Now on {wig_name}')
        
        GB_mean_ar=[]
        GB_std_ar=[]
        Distance_ar=[]
        Sample_size=[]
        
        for TSS_distance in range(4000):

            GB_signal_ar=[]
            
            for TU_name, TU_data in genes_annotation.items():
                TU_start=TU_data[0]
                TU_end=TU_data[1]
                TU_strand=TU_data[2]
                
                if TU_strand=="+":
                            
                    if TU_start+TSS_distance<TU_end:
                        GB_signal=wig_data[TU_start+TSS_distance]
                        GB=wig_data[TU_start:TU_start+TSS_distance]
                    elif TU_start+TSS_distance>=TU_end:
                        GB_signal=-1
                        GB=wig_data[TU_start:TU_end]
                        
                elif TU_strand=="-":
                    
                    if TU_end-TSS_distance>TU_start:
                        GB_signal=wig_data[TU_end-TSS_distance]
                        GB=wig_data[TU_end-TSS_distance:TU_end][::-1]
                    else:
                        GB_signal=-1
                        GB=wig_data[TU_start:TU_end][::-1]
                        
                if GB_signal!=-1:
                    GB_signal_ar.append(GB_signal)
                
                
                    
            GB_mean=np.mean(GB_signal_ar)
            GB_std=np.std(GB_signal_ar)
            
            GB_mean_ar.append(GB_mean)
            GB_std_ar.append(GB_std)
            Distance_ar.append(TSS_distance)
            Sample_size.append(len(GB_signal_ar))
            
        Signal_dict[wig_name]=[GB_mean_ar, GB_std_ar, Distance_ar, Sample_size]
     
    
    return  Signal_dict


#######
#Plot FE TSS distance dependence.
#######

def plot_TSS_distance_FE(Signal_dict, Half_data_dict, output_path):
    
    #Plot FE over genes.
    plt.figure(figsize=(6, 5), dpi=100)
    plot1=plt.subplot(111)

    #Standard order of colors: #B6B8BD, #333738, #b08642, #757d8b.
    
    
    
    #Dataset 1
    dataset1='TopA_CTD_minus_Rif_minus_av_1_2_3'
    plot1.plot(Signal_dict[dataset1][2], Signal_dict[dataset1][0], linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=r"HETU $\overline{FE}\pm$SE", zorder=6)
    
    Upper_conf_interval=np.array(Signal_dict[dataset1][0])+(np.array(Signal_dict[dataset1][1])/np.sqrt(Signal_dict[dataset1][3]))
    Lower_conf_interval=np.array(Signal_dict[dataset1][0])-(np.array(Signal_dict[dataset1][1])/np.sqrt(Signal_dict[dataset1][3]))
    
    plot1.plot(Signal_dict[dataset1][2], Upper_conf_interval, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, zorder=5)       
    plot1.plot(Signal_dict[dataset1][2], Lower_conf_interval, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, zorder=5) 
    plot1.fill_between(Signal_dict[dataset1][2], Lower_conf_interval, Upper_conf_interval, facecolor='#43c287', alpha=0.4, interpolate=True, zorder=4)   
    
    Half_signal_1=np.mean(Half_data_dict[dataset1][0])
    Half_distance_1=np.mean(Half_data_dict[dataset1][1])
    Half_signal_SE_1=np.std(Half_data_dict[dataset1][0])/np.sqrt(len(Half_data_dict[dataset1][0]))
    Half_distance_SE_1=np.std(Half_data_dict[dataset1][1])/np.sqrt(len(Half_data_dict[dataset1][1]))  
    
          
    #Dataset 2
    dataset2='TopA_CTD_plus_Rif_plus_av_2_3'
    plot1.plot(Signal_dict[dataset2][2], Signal_dict[dataset2][0], linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=r"HETU Rif CTD $\overline{FE}\pm$ SE", zorder=3)
    
    Upper_conf_interval=np.array(Signal_dict[dataset2][0])+(np.array(Signal_dict[dataset2][1])/np.sqrt(Signal_dict[dataset2][3]))
    Lower_conf_interval=np.array(Signal_dict[dataset2][0])-(np.array(Signal_dict[dataset2][1])/np.sqrt(Signal_dict[dataset2][3]))
    
    plot1.plot(Signal_dict[dataset2][2], Upper_conf_interval, linestyle='-', color='#333738', linewidth=1, alpha=1, zorder=2)       
    plot1.plot(Signal_dict[dataset2][2], Lower_conf_interval, linestyle='-', color='#333738', linewidth=1, alpha=1, zorder=2)     
    plot1.fill_between(Signal_dict[dataset2][2], Lower_conf_interval, Upper_conf_interval, facecolor='#7ce0ff', alpha=0.4, interpolate=True, zorder=1) 
    
    Half_signal_2=np.mean(Half_data_dict[dataset2][0])
    Half_distance_2=np.mean(Half_data_dict[dataset2][1])
    Half_signal_SE_2=np.std(Half_data_dict[dataset2][0])/np.sqrt(len(Half_data_dict[dataset2][0]))
    Half_distance_SE_2=np.std(Half_data_dict[dataset2][1])/np.sqrt(len(Half_data_dict[dataset2][1]))    
    
    #Sample size.
    plot2=plot1.twinx() 
    plot2.plot(Signal_dict[dataset1][2], Signal_dict[dataset1][3], linestyle='--', color='red', linewidth=1.5, alpha=1, label='Sample size')
    
                    
    ticks=sorted(np.arange(0, len(Signal_dict[dataset2][2])-100, 1000).tolist() + [len(Signal_dict[dataset2][2])])
    plot1.set_xticks(ticks, minor=False)
    ticks_lables=sorted(np.arange(0, len(Signal_dict[dataset2][2])-100, 1000).tolist() + [len(Signal_dict[dataset2][2])])
    TSS_index=ticks_lables.index(0)
    TES_index=ticks_lables.index(len(Signal_dict[dataset2][2]))
    ticks_lables[TSS_index]='TS'
    ticks_lables[TES_index]='TE'
    print(ticks)
    print(ticks_lables)
    plot1.set_xticklabels(ticks_lables, minor=False) 
    
    yticks1=[Half_signal_1-Half_signal_SE_1, Half_signal_1, Half_signal_1+Half_signal_SE_1, Half_signal_2-Half_signal_SE_2, Half_signal_2, Half_signal_2+Half_signal_SE_2]
    xticks1=[Half_distance_1-Half_distance_SE_1, Half_distance_1, Half_distance_1+Half_distance_SE_1, Half_distance_2-Half_distance_SE_2, Half_distance_2, Half_distance_2+Half_distance_SE_2]
    print(yticks1)
    print(xticks1)
    
    plot1.set_yticks(yticks1, minor=True)
    plot1.set_xticks(xticks1, minor=True)
    
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(len(Signal_dict[dataset2][2]), color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axhline(1, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_xlim([-100, 4100])  
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI fold enrichment', size=20) 
    plot2.set_ylabel('Number of TUs', size=20)
    plot1.legend(fontsize=12, frameon=False, loc='upper right', bbox_to_anchor=(0.6, 1.3))
    plot2.legend(fontsize=12, frameon=False, loc='upper left', bbox_to_anchor=(0.6, 1.3))
    plt.tight_layout(rect=[0,0,1,0.9])    
    plt.savefig(f'{output_path}EcTopoI_FE_over_HETU_Rif_CTD_plus_minus.png', dpi=400, figsize=(6, 5), transparent=False)
    plt.show()
    plt.close()   
    
    return


Signal_form_TSS_dict=TSS_distance_FE(WIG_data_dict, TUs_dictionary)
Half_data_dictionary=compute_half_distance(WIG_data_dict, TUs_dictionary)
plot_TSS_distance_FE(Signal_form_TSS_dict, Half_data_dictionary, Outpath)



#################################
#################################
######### Part 3.
######### Normalize a signal.
######### Compute dependence of FE from the distance from TSS. Do statistics.
#################################
#################################

#######
#Normalize signal data.
#######

def signal_normalize(wig_data_dict, width):
    
    Transformed_data_dict={}
    
    #Estimate mean and standard deviation of datasets.
    for wig_name, wig_data in wig_data_dict.items():
    
        #Mask deletions.
        deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]]
        mask=[0]*len(wig_data)
        for delition in deletions:
            mask[delition[0]:delition[1]]=[1]*(delition[1]-delition[0])  
        
        wig_data_masked=[]
        for i in range(len(mask)):
            if mask[i]==0:
                wig_data_masked.append(wig_data[i])
        
        #Bin genome with bins of length width.
        number_of_genome_bin=int(len(wig_data_masked)/width)-1
        genome_signal_values_binned=[]
        for i in range(number_of_genome_bin):
            genome_signal_values_binned.append(np.mean(wig_data_masked[width*i:width*(i+1)]))
        
        Signal_mean_genome=np.mean(wig_data_masked)
        Signal_std_genome=np.std(wig_data_masked)
        
        Signal_mean_genome_binned=np.mean(genome_signal_values_binned)
        Signal_std_genome_binned=np.std(genome_signal_values_binned)    
        
        print(f'Global mean: {Signal_mean_genome}; Global std: {Signal_std_genome}') 
        print(f'Signal binned mean: {Signal_mean_genome_binned}; Signal binned std: {Signal_std_genome_binned}') 
        
        Signal_estimators=[Signal_mean_genome, Signal_std_genome, Signal_mean_genome_binned, Signal_std_genome_binned]
        
        #Compute relative signal X-mean(X).
        wig_relative_global=np.array(wig_data)-Signal_mean_genome
        wig_normalized_global=(np.array(wig_data)-Signal_mean_genome)/Signal_std_genome
        
        wig_relative_binned=np.array(wig_data)-Signal_mean_genome_binned
        wig_normalized_binned=(np.array(wig_data)-Signal_mean_genome_binned)/Signal_std_genome_binned
        
        wig_transformed=[wig_relative_global, wig_normalized_global, wig_relative_binned, wig_normalized_binned]
        
        Transformed_data_dict[wig_name]={'Estimators' : Signal_estimators, 'Transformed data' : wig_transformed}
        
    
    return Transformed_data_dict


#######
#Take FE values from some distance from TSS going into TU body. 
#If we reach TU end, discard data.
#Compute mean and std, plot for several conditions (e.g. -CTD vs +CTD).
#######

def TSS_distance_FE_normalized(Transformed_data_dict, genes_annotation):
    
    for wig_name, data in Transformed_data_dict.items():
        
        print(f'Now on a {wig_name}.')
        
        Prepared_data_ar=[]
    
        for wig_data in data['Transformed data']:
            
            GB_mean_ar=[]
            GB_std_ar=[]
            Distance_ar=[]
            Sample_size=[]
            
            for TSS_distance in range(15000):
        
                GB_signal_ar=[]
                
                for TU_name, TU_data in genes_annotation.items():
                    TU_start=TU_data[0]
                    TU_end=TU_data[1]
                    TU_strand=TU_data[2]
                    
                    if TU_strand=="+":
                                
                        if TU_start+TSS_distance<TU_end:
                            GB_signal=wig_data[TU_start+TSS_distance]
                        elif TU_start+TSS_distance>=TU_end:
                            GB_signal=-1
                            
                    elif TU_strand=="-":
                        
                        if TU_end-TSS_distance>TU_start:
                            GB_signal=wig_data[TU_end-TSS_distance]
                        else:
                            GB_signal=-1
                            
                    if GB_signal!=-1:
                        GB_signal_ar.append(GB_signal)
                        
                GB_mean=np.mean(GB_signal_ar)
                GB_std=np.std(GB_signal_ar)
                
                GB_mean_ar.append(GB_mean)
                GB_std_ar.append(GB_std)
                Distance_ar.append(TSS_distance)
                Sample_size.append(len(GB_signal_ar))
                
            Prepared_data_ar.append([GB_mean_ar, GB_std_ar, Distance_ar, Sample_size])
                
        Transformed_data_dict[wig_name]['Prepared FE']=Prepared_data_ar
    
    return  Transformed_data_dict


#######
#Take FE values from some distance from TSS going into TU body.
#Compute mean and std, plot for several conditions (e.g. -CTD vs +CTD).
#######

def TSS_distance_FE_normalized_no_TU_end(Transformed_data_dict, genes_annotation):
    
    for wig_name, data in Transformed_data_dict.items():
        
        print(f'Now on a {wig_name}.')
        
        Prepared_data_ar=[]
    
        for wig_data in data['Transformed data']:
            
            GB_mean_ar=[]
            GB_std_ar=[]
            Distance_ar=[]
            Sample_size=[]
            
            for TSS_distance in range(4000):
        
                GB_signal_ar=[]
                
                for TU_name, TU_data in genes_annotation.items():
                    TU_start=TU_data[0]
                    TU_end=TU_data[1]
                    TU_strand=TU_data[2]
                    
                    if TU_strand=="+":
                                
                        if TU_start+TSS_distance<len(wig_data):
                            GB_signal=wig_data[TU_start+TSS_distance]
                        elif TU_start+TSS_distance>=len(wig_data):
                            GB_signal=wig_data[TU_start+TSS_distance-len(wig_data)]
                            
                    elif TU_strand=="-":
                        
                        GB_signal=wig_data[TU_end-TSS_distance]
                            
                    GB_signal_ar.append(GB_signal)
                        
                GB_mean=np.mean(GB_signal_ar)
                GB_std=np.std(GB_signal_ar)
                
                GB_mean_ar.append(GB_mean)
                GB_std_ar.append(GB_std)
                Distance_ar.append(TSS_distance)
                Sample_size.append(len(GB_signal_ar))
                
            Prepared_data_ar.append([GB_mean_ar, GB_std_ar, Distance_ar, Sample_size])
                
        Transformed_data_dict[wig_name]['Prepared FE']=Prepared_data_ar
    
    return  Transformed_data_dict


#######
#Plot normalized FE TSS distance dependence.
#######

def plot_TSS_distance_FE_normalized(Signal_dict, output_path):
    
    
    ####################################################
    #Retrive data. Relative FE (X-mean(X)), global mean.
    ####################################################
    Distance_ar_rgm=Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][0][2]
    Mean_ar_rgm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][0][0])
    Std_ar_rgm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][0][1])
    Sample_size_ar_rgm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][0][3])
    SE_ar_rgm=Std_ar_rgm/np.sqrt(Sample_size_ar_rgm)
    Upper_conf_interval_rgm=Mean_ar_rgm+SE_ar_rgm
    Lower_conf_interval_rgm=Mean_ar_rgm-SE_ar_rgm
    
    Distance_ar_rgp=Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][0][2]
    Mean_ar_rgp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][0][0])
    Std_ar_rgp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][0][1])
    Sample_size_ar_rgp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][0][3])
    SE_ar_rgp=Std_ar_rgp/np.sqrt(Sample_size_ar_rgp)
    Upper_conf_interval_rgp=Mean_ar_rgp+SE_ar_rgp
    Lower_conf_interval_rgp=Mean_ar_rgp-SE_ar_rgp    
    
    #Plot FE over genes.
    plt.figure(figsize=(6, 5), dpi=100)
    plot1=plt.subplot(111)

    #Standard order of colors: #B6B8BD, #333738, #b08642, #757d8b.
    
    #CTD-/Rif-
    plot1.plot(Distance_ar_rgm, Mean_ar_rgm, linestyle='-', color='#B6B8BD', linewidth=1.5, alpha=1, label=r"HETU $\overline{relFE}\pm$SE", zorder=6)
    plot1.plot(Distance_ar_rgm, Upper_conf_interval_rgm, linestyle='-', color='#B6B8BD', linewidth=0.5, alpha=0.3, zorder=5)       
    plot1.plot(Distance_ar_rgm, Lower_conf_interval_rgm, linestyle='-', color='#B6B8BD', linewidth=0.5, alpha=0.3, zorder=5) 
    plot1.fill_between(Distance_ar_rgm, Lower_conf_interval_rgm, Upper_conf_interval_rgm, facecolor='#43c287', alpha=0.3, interpolate=True, zorder=4)        
          
    #CTD+/Rif-
    plot1.plot(Distance_ar_rgp, Mean_ar_rgp, linestyle='-', color='#333738', linewidth=1.5, alpha=1, label=r"HETU Rif $\overline{relFE}\pm$ SE", zorder=3)
    plot1.plot(Distance_ar_rgp, Upper_conf_interval_rgp, linestyle='-', color='#333738', linewidth=0.5, alpha=0.3, zorder=2)       
    plot1.plot(Distance_ar_rgp, Lower_conf_interval_rgp, linestyle='-', color='#333738', linewidth=0.5, alpha=0.3, zorder=2)     
    plot1.fill_between(Distance_ar_rgp, Lower_conf_interval_rgp, Upper_conf_interval_rgp, facecolor='#7ce0ff', alpha=0.3, interpolate=True, zorder=1) 
    
    #Sample size.
    plot2=plot1.twinx() 
    plot2.plot(Distance_ar_rgm, Sample_size_ar_rgm, linestyle='--', color='red', linewidth=1.2, alpha=1, label='Sample size')
                 
    ticks=sorted(np.arange(0, len(Distance_ar_rgm)+100, 1000).tolist() + [len(Distance_ar_rgm)])
    plot1.set_xticks(ticks, minor=False)
    ticks_lables=sorted(np.arange(0, len(Distance_ar_rgm)+100, 1000).tolist() + [len(Distance_ar_rgm)])
    TSS_index=ticks_lables.index(0)
    TES_index=ticks_lables.index(len(Distance_ar_rgm))
    ticks_lables[TSS_index]='TS'
    ticks_lables[TES_index]='TE'
    plot1.set_xticklabels(ticks_lables, minor=False) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(len(Distance_ar_rgm), color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axhline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_xlim([-100, 4100])    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI relative FE', size=20) 
    plot2.set_ylabel('Number of TUs', size=20)
    plot1.legend(fontsize=12, frameon=False, loc='upper right', bbox_to_anchor=(0.6, 1.3))
    plot2.legend(fontsize=12, frameon=False, loc='upper left', bbox_to_anchor=(0.6, 1.3))
    plt.tight_layout(rect=[0,0,1,0.9])     
    plt.savefig(f'{output_path}EcTopoI_global_relative_FE_over_HETU_Rif_plus_minus.png', dpi=400, figsize=(6, 5), transparent=False)
    plt.show()
    plt.close()  
    
    
    ####################################################
    #Retrive data. Normalized FE (X-mean(X))/sigma(X), global mean, std.
    ####################################################
    Distance_ar_ngm=Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][1][2]
    Mean_ar_ngm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][1][0])
    Std_ar_ngm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][1][1])
    Sample_size_ar_ngm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][1][3])
    SE_ar_ngm=Std_ar_ngm/np.sqrt(Sample_size_ar_ngm)
    Upper_conf_interval_ngm=Mean_ar_ngm+SE_ar_ngm
    Lower_conf_interval_ngm=Mean_ar_ngm-SE_ar_ngm
    
    Distance_ar_ngp=Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][1][2]
    Mean_ar_ngp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][1][0])
    Std_ar_ngp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][1][1])
    Sample_size_ar_ngp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][1][3])
    SE_ar_ngp=Std_ar_ngp/np.sqrt(Sample_size_ar_ngp)
    Upper_conf_interval_ngp=Mean_ar_ngp+SE_ar_ngp
    Lower_conf_interval_ngp=Mean_ar_ngp-SE_ar_ngp    
    
    #Plot FE over genes.
    plt.figure(figsize=(6, 5), dpi=100)
    plot1=plt.subplot(111)

    #Standard order of colors: #B6B8BD, #333738, #b08642, #757d8b.
    
    #CTD-/Rif-
    plot1.plot(Distance_ar_ngm, Mean_ar_ngm, linestyle='-', color='#B6B8BD', linewidth=1.5, alpha=1, label=r"HETU $\overline{normFE}\pm$SE", zorder=6)
    plot1.plot(Distance_ar_ngm, Upper_conf_interval_ngm, linestyle='-', color='#B6B8BD', linewidth=0.5, alpha=0.3, zorder=5)       
    plot1.plot(Distance_ar_ngm, Lower_conf_interval_ngm, linestyle='-', color='#B6B8BD', linewidth=0.5, alpha=0.3, zorder=5) 
    plot1.fill_between(Distance_ar_ngm, Lower_conf_interval_ngm, Upper_conf_interval_ngm, facecolor='#43c287', alpha=0.3, interpolate=True, zorder=4)        
          
    #CTD+/Rif-
    plot1.plot(Distance_ar_ngp, Mean_ar_ngp, linestyle='-', color='#333738', linewidth=1.5, alpha=1, label=r"HETU Rif $\overline{normFE}\pm$ SE", zorder=3)
    plot1.plot(Distance_ar_ngp, Upper_conf_interval_ngp, linestyle='-', color='#333738', linewidth=0.5, alpha=0.3, zorder=2)       
    plot1.plot(Distance_ar_ngp, Lower_conf_interval_ngp, linestyle='-', color='#333738', linewidth=0.5, alpha=0.3, zorder=2)     
    plot1.fill_between(Distance_ar_ngp, Lower_conf_interval_ngp, Upper_conf_interval_ngp, facecolor='#7ce0ff', alpha=0.3, interpolate=True, zorder=1) 
    
    #Sample size.
    plot2=plot1.twinx() 
    plot2.plot(Distance_ar_ngm, Sample_size_ar_ngm, linestyle='--', color='red', linewidth=1.2, alpha=1, label='Sample size')
                 
    ticks=sorted(np.arange(0, len(Distance_ar_ngm)+100, 1000).tolist() + [len(Distance_ar_ngm)])
    plot1.set_xticks(ticks, minor=False)
    ticks_lables=sorted(np.arange(0, len(Distance_ar_ngm)+100, 1000).tolist() + [len(Distance_ar_ngm)])
    TSS_index=ticks_lables.index(0)
    TES_index=ticks_lables.index(len(Distance_ar_ngm))
    ticks_lables[TSS_index]='TS'
    ticks_lables[TES_index]='TE'
    plot1.set_xticklabels(ticks_lables, minor=False) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(len(Distance_ar_ngm), color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axhline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_xlim([-100, 4100])   
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI normalized FE', size=20) 
    plot2.set_ylabel('Number of TUs', size=20)
    plot1.legend(fontsize=12, frameon=False, loc='upper right', bbox_to_anchor=(0.6, 1.3))
    plot2.legend(fontsize=12, frameon=False, loc='upper left', bbox_to_anchor=(0.6, 1.3))
    plt.tight_layout(rect=[0,0,1,0.9])       
    plt.savefig(f'{output_path}EcTopoI_global_normalized_FE_over_HETU_Rif_plus_minus.png', dpi=400, figsize=(6, 5), transparent=False)
    plt.show()
    plt.close()       
    
    
    ####################################################
    #Retrive data. Relative FE (X-mean(X)), binned mean.
    ####################################################
    Distance_ar_rbm=Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][2][2]
    Mean_ar_rbm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][2][0])
    Std_ar_rbm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][2][1])
    Sample_size_ar_rbm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][2][3])
    SE_ar_rbm=Std_ar_rbm/np.sqrt(Sample_size_ar_rbm)
    Upper_conf_interval_rbm=Mean_ar_rbm+SE_ar_rbm
    Lower_conf_interval_rbm=Mean_ar_rbm-SE_ar_rbm
    
    Distance_ar_rbp=Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][2][2]
    Mean_ar_rbp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][2][0])
    Std_ar_rbp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][2][1])
    Sample_size_ar_rbp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][2][3])
    SE_ar_rbp=Std_ar_rbp/np.sqrt(Sample_size_ar_rbp)
    Upper_conf_interval_rbp=Mean_ar_rbp+SE_ar_rbp
    Lower_conf_interval_rbp=Mean_ar_rbp-SE_ar_rbp    
    
    #Plot FE over genes.
    plt.figure(figsize=(6, 5), dpi=100)
    plot1=plt.subplot(111)

    #Standard order of colors: #B6B8BD, #333738, #b08642, #757d8b.
    
    #CTD-/Rif-
    plot1.plot(Distance_ar_rbm, Mean_ar_rbm, linestyle='-', color='#B6B8BD', linewidth=1.5, alpha=1, label=r"HETU $\overline{relFE}\pm$SE", zorder=6)
    plot1.plot(Distance_ar_rbm, Upper_conf_interval_rbm, linestyle='-', color='#B6B8BD', linewidth=0.5, alpha=0.3, zorder=5)       
    plot1.plot(Distance_ar_rbm, Lower_conf_interval_rbm, linestyle='-', color='#B6B8BD', linewidth=0.5, alpha=0.3, zorder=5) 
    plot1.fill_between(Distance_ar_rbm, Lower_conf_interval_rbm, Upper_conf_interval_rbm, facecolor='#43c287', alpha=0.3, interpolate=True, zorder=4)        
          
    #CTD+/Rif-
    plot1.plot(Distance_ar_rbp, Mean_ar_rbp, linestyle='-', color='#333738', linewidth=1.5, alpha=1, label=r"HETU Rif $\overline{relFE}\pm$ SE", zorder=3)
    plot1.plot(Distance_ar_rbp, Upper_conf_interval_rbp, linestyle='-', color='#333738', linewidth=0.5, alpha=0.3, zorder=2)       
    plot1.plot(Distance_ar_rbp, Lower_conf_interval_rbp, linestyle='-', color='#333738', linewidth=0.5, alpha=0.3, zorder=2)     
    plot1.fill_between(Distance_ar_rbp, Lower_conf_interval_rbp, Upper_conf_interval_rbp, facecolor='#7ce0ff', alpha=0.3, interpolate=True, zorder=1) 
    
    #Sample size.
    plot2=plot1.twinx() 
    plot2.plot(Distance_ar_rbm, Sample_size_ar_rbm, linestyle='--', color='red', linewidth=1.2, alpha=1, label='Sample size')
                 
    ticks=sorted(np.arange(0, len(Distance_ar_rbm)+100, 1000).tolist() + [len(Distance_ar_rbm)])
    plot1.set_xticks(ticks, minor=False)
    ticks_lables=sorted(np.arange(0, len(Distance_ar_rbm)+100, 1000).tolist() + [len(Distance_ar_rbm)])
    TSS_index=ticks_lables.index(0)
    TES_index=ticks_lables.index(len(Distance_ar_rbm))
    ticks_lables[TSS_index]='TS'
    ticks_lables[TES_index]='TE'
    plot1.set_xticklabels(ticks_lables, minor=False) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(len(Distance_ar_rbm), color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axhline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_xlim([-100, 4100])   
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI relative FE', size=20) 
    plot2.set_ylabel('Number of TUs', size=20)
    plot1.legend(fontsize=12, frameon=False, loc='upper right', bbox_to_anchor=(0.6, 1.3))
    plot2.legend(fontsize=12, frameon=False, loc='upper left', bbox_to_anchor=(0.6, 1.3))
    plt.tight_layout(rect=[0,0,1,0.9])       
    plt.savefig(f'{output_path}EcTopoI_binned_relative_FE_over_HETU_Rif_plus_minus.png', dpi=400, figsize=(6, 5), transparent=False)
    plt.show()
    plt.close()  
    
    
    ####################################################
    #Retrive data. Normalized FE (X-mean(X))/sigma(X), binned mean, sigma.
    ####################################################
    Distance_ar_nbm=Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][3][2]
    Mean_ar_nbm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][3][0])
    Std_ar_nbm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][3][1])
    Sample_size_ar_nbm=np.array(Signal_dict['TopA_CTD_minus_Rif_minus_av_1_2_3']['Prepared FE'][3][3])
    SE_ar_nbm=Std_ar_nbm/np.sqrt(Sample_size_ar_nbm)
    Upper_conf_interval_nbm=Mean_ar_nbm+SE_ar_nbm
    Lower_conf_interval_nbm=Mean_ar_nbm-SE_ar_nbm
    
    Distance_ar_nbp=Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][3][2]
    Mean_ar_nbp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][3][0])
    Std_ar_nbp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][3][1])
    Sample_size_ar_nbp=np.array(Signal_dict['TopA_CTD_minus_Rif_plus_av_1_2_3']['Prepared FE'][3][3])
    SE_ar_nbp=Std_ar_nbp/np.sqrt(Sample_size_ar_nbp)
    Upper_conf_interval_nbp=Mean_ar_nbp+SE_ar_nbp
    Lower_conf_interval_nbp=Mean_ar_nbp-SE_ar_nbp    
    
    #Plot FE over genes.
    plt.figure(figsize=(6, 5), dpi=100)
    plot1=plt.subplot(111)

    #Standard order of colors: #B6B8BD, #333738, #b08642, #757d8b.
    
    #CTD-/Rif-
    plot1.plot(Distance_ar_nbm, Mean_ar_nbm, linestyle='-', color='#B6B8BD', linewidth=1.5, alpha=1, label=r"HETU $\overline{normFE}\pm$SE", zorder=6)
    plot1.plot(Distance_ar_nbm, Upper_conf_interval_nbm, linestyle='-', color='#B6B8BD', linewidth=0.5, alpha=0.3, zorder=5)       
    plot1.plot(Distance_ar_nbm, Lower_conf_interval_nbm, linestyle='-', color='#B6B8BD', linewidth=0.5, alpha=0.3, zorder=5) 
    plot1.fill_between(Distance_ar_nbm, Lower_conf_interval_nbm, Upper_conf_interval_nbm, facecolor='#43c287', alpha=0.3, interpolate=True, zorder=4)        
          
    #CTD+/Rif-
    plot1.plot(Distance_ar_nbp, Mean_ar_nbp, linestyle='-', color='#333738', linewidth=1.5, alpha=1, label=r"HETU Rif $\overline{normFE}\pm$ SE", zorder=3)
    plot1.plot(Distance_ar_nbp, Upper_conf_interval_nbp, linestyle='-', color='#333738', linewidth=0.5, alpha=0.3, zorder=2)       
    plot1.plot(Distance_ar_nbp, Lower_conf_interval_nbp, linestyle='-', color='#333738', linewidth=0.5, alpha=0.3, zorder=2)     
    plot1.fill_between(Distance_ar_nbp, Lower_conf_interval_nbp, Upper_conf_interval_nbp, facecolor='#7ce0ff', alpha=0.3, interpolate=True, zorder=1) 
    
    #Sample size.
    plot2=plot1.twinx() 
    plot2.plot(Distance_ar_nbm, Sample_size_ar_nbm, linestyle='--', color='red', linewidth=1.2, alpha=1, label='Sample size')
                 
    ticks=sorted(np.arange(0, len(Distance_ar_nbm)+100, 1000).tolist() + [len(Distance_ar_nbm)])
    plot1.set_xticks(ticks, minor=False)
    ticks_lables=sorted(np.arange(0, len(Distance_ar_nbm)+100, 1000).tolist() + [len(Distance_ar_nbm)])
    TSS_index=ticks_lables.index(0)
    TES_index=ticks_lables.index(len(Distance_ar_nbm))
    ticks_lables[TSS_index]='TS'
    ticks_lables[TES_index]='TE'
    plot1.set_xticklabels(ticks_lables, minor=False) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(len(Distance_ar_nbm), color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axhline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_xlim([-100, 4100])    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI normalized FE', size=20) 
    plot2.set_ylabel('Number of TUs', size=20)
    plot1.legend(fontsize=12, frameon=False, loc='upper right', bbox_to_anchor=(0.6, 1.3))
    plot2.legend(fontsize=12, frameon=False, loc='upper left', bbox_to_anchor=(0.6, 1.3))
    plt.tight_layout(rect=[0,0,1,0.9])       
    plt.savefig(f'{output_path}EcTopoI_binned_normalized_FE_over_HETU_Rif_plus_minus.png', dpi=400, figsize=(6, 5), transparent=False)
    plt.show()
    plt.close()      
    
    return



#Transformed_data_dictionary=signal_normalize(WIG_data_dict, 200)
#Signal_normalized_from_TSS_dict=TSS_distance_FE_normalized(Transformed_data_dictionary, TUs_dictionary)
#Signal_normalized_from_TSS_TU_no_end_dict=TSS_distance_FE_normalized_no_TU_end(Transformed_data_dictionary, TUs_dictionary)
#plot_TSS_distance_FE_normalized(Signal_normalized_from_TSS_dict, Outpath)


