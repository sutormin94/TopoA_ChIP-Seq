###############################################
##Dmitry Sutormin, 2020##
##ChIP-Seq analysis##

####
#Filtering of transcription units, identification of representative ones by RNA-Seq data.
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import scipy


#Import RNA-Seq data, wig.
RNAseq_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Sutormin_RNA_Seq_Exponential_av.wig"

#Transcription units.
TUs_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_transcripts\DY330_RNA-Seq_transcripts_EP_del_cor.txt"

#Output path for filtered data.
TUs_filtered_outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_transcripts\DY330_RNA-Seq_transcripts_representative_EP_del_cor.txt"

#######
#Read wig data.
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
    print(len(NE_values))
    return NE_values

RNASeq_data=wig_parsing(RNAseq_data_path)

#######
#Read TUs data, filter.
#######

def read_TUs(TU_inpath, wig_data, outfile_path):
    
    #Read TUs.
    TUs_dict={}
    filein=open(TU_inpath, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['TU_ID']:
            TU_ID=line[0]
            TU_name=line[1]
            TU_start=int(line[2])
            TU_end=int(line[3])
            TU_strand=line[4]
            TU_expression=float(line[5])
            Description=line[6]
            if TU_name in TUs_dict:
                TUs_dict[TU_name][1].append([TU_start, TU_end, TU_expression])
            else:
                TUs_dict[TU_name]=[[TU_ID, TU_name, TU_strand, Description], [[TU_start, TU_end, TU_expression]]]
    filein.close()
    
    
    #Write filetered data.
    fileout=open(outfile_path, 'w')
    fileout.write('TU_ID\tGene_name\tStart\tEnd\tStrand\tExpression_E\tDescription\n')
    
    #Sort TUs by Start site.
    for TU_type, TU_data in TUs_dict.items():
        if len(TU_data[1])>1:
            #print(f'Before sorting {TU_data[0][2]} {TU_data[1]}')
            #Transcription start refinment.
            if TU_data[0][2]=="+":
                TU_data[1].sort(key=lambda x: int(x[0]))
                
                #Analysis of 3' UTRs coverage.
                TU_start_adjusted=TU_data[1][-1][0]
                for i in range(len(TU_data[1])-1):
                    #print(TU_data[1][-(i+2)][0], TU_data[1][-(i+1)][0])
                    reg_s=TU_data[1][-(i+2)][0]
                    reg_e=TU_data[1][-(i+1)][0]
                    reg_signal=np.mean(wig_data[reg_s:reg_e])
                    body_signal=np.mean(wig_data[reg_e:TU_data[1][-(i+1)][1]])
                    #print(f' UTR: {reg_signal} Body: {body_signal}')
                    if reg_signal>0.2*body_signal:
                        TU_start_adjusted=reg_s
                    
            elif TU_data[0][2]=="-":
                TU_data[1].sort(key=lambda x: int(x[1]))
                
                #Analysis of 3' UTRs coverage.
                TU_start_adjusted=TU_data[1][0][1]
                for i in range(len(TU_data[1])-1):
                    #print(TU_data[1][i][1], TU_data[1][i+1][1])
                    reg_s=TU_data[1][i][1]
                    reg_e=TU_data[1][i+1][1]
                    reg_signal=np.mean(wig_data[reg_s:reg_e])
                    body_signal=np.mean(wig_data[TU_data[1][i][0]:reg_s])
                    #print(f' UTR: {reg_signal} Body: {body_signal}') 
                    if reg_signal>0.2*body_signal:
                        TU_start_adjusted=reg_e                 
                 
            print(f'After sorting {TU_data[0][2]} {TU_data[1]}')
            print(f'Refined TU_start: {TU_start_adjusted}')   
            
            #Filter possible TUs by adjusted TU_start.
            TU_data_start_adjusted=[]
            if TU_data[0][2]=="+":
                for i in range(len(TU_data[1])):
                    if TU_start_adjusted==TU_data[1][i][0]:
                        TU_data_start_adjusted.append(TU_data[1][i])
            elif TU_data[0][2]=="-":
                for i in range(len(TU_data[1])):
                    if TU_start_adjusted==TU_data[1][i][1]:
                        TU_data_start_adjusted.append(TU_data[1][i])
                
            if len(TU_data_start_adjusted)>1:
                print(f'TU_end refinment is needed: {TU_data_start_adjusted}')
            
            #Transcription end refinment.
            if TU_data[0][2]=="+":
                TU_data_start_adjusted.sort(key=lambda x: int(x[1])) 
                
                #Analysis of 5' UTRs coverage.
                TU_end_adjusted=TU_data_start_adjusted[0][1]
                for i in range(len(TU_data_start_adjusted)-1):
                    #print(TU_data[1][-(i+2)][0], TU_data[1][-(i+1)][0])
                    reg_s=TU_data_start_adjusted[i][1]
                    reg_e=TU_data_start_adjusted[i+1][1]
                    reg_signal=np.mean(wig_data[reg_s:reg_e])
                    body_signal=np.mean(wig_data[TU_data_start_adjusted[i][0]:reg_s])
                    print(f' UTR: {reg_signal} Body: {body_signal}')
                    if reg_signal>0.2*body_signal:
                        TU_end_adjusted=reg_e
                        
                
            elif TU_data[0][2]=="-":
                TU_data_start_adjusted.sort(key=lambda x: int(x[0])) 
                
                #Analysis of 5' UTRs coverage.
                TU_end_adjusted=TU_data_start_adjusted[-1][0]
                for i in range(len(TU_data_start_adjusted)-1):
                    #print(TU_data[1][i][1], TU_data[1][i+1][1])
                    reg_s=TU_data_start_adjusted[-(i+2)][0]
                    reg_e=TU_data_start_adjusted[-(i+1)][0]
                    reg_signal=np.mean(wig_data[reg_s:reg_e])
                    body_signal=np.mean(wig_data[reg_e:TU_data_start_adjusted[-(i+1)][1]])
                    print(f' UTR: {reg_signal} Body: {body_signal}') 
                    if reg_signal>0.2*body_signal:
                        TU_end_adjusted=reg_s                  
            
            print(f'Refined TU_end: {TU_end_adjusted}') 
            
            #Write data.
            if TU_data[0][2]=="+":
                fileout.write(f'{TU_data[0][0]}\t{TU_data[0][1]}\t{TU_start_adjusted}\t{TU_end_adjusted}\t{TU_data[0][2]}\t{np.mean(wig_data[TU_start_adjusted:TU_end_adjusted])}\t{TU_data[0][3]}\n')
            elif TU_data[0][2]=="-":
                fileout.write(f'{TU_data[0][0]}\t{TU_data[0][1]}\t{TU_end_adjusted}\t{TU_start_adjusted}\t{TU_data[0][2]}\t{np.mean(wig_data[TU_end_adjusted:TU_start_adjusted])}\t{TU_data[0][3]}\n')
        
        #Write TU if only one transcript is reported.
        else:
            print(TU_data)
            fileout.write(f'{TU_data[0][0]}\t{TU_data[0][1]}\t{TU_data[1][0][0]}\t{TU_data[1][0][1]}\t{TU_data[0][2]}\t{TU_data[1][0][2]}\t{TU_data[0][3]}\n')
            
    fileout.close()    
    
    return


read_TUs(TUs_path, RNASeq_data, TUs_filtered_outpath)