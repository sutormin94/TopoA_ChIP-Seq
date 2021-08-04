###############################################
##Dmitry Sutormin, 2021##
##Translate coordinates of MG1655 chip probes to DY330 genome##

####
#Reads chip results.
#Calculates chip signal: signal/background.
#BLAST probes sequences against target genome (DY330).
#Return new coordinates.
#Integrate chip data and return WIG track with coverage depth.
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import re
import matplotlib.pyplot as plt
import os

#######
#Data to be used.
#######

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\Mooney_RNAP_NusAG_Rho\\"

#Dictionary of IP path.
Dict_IP_path={'NusA_1'   : PWD + 'Raw_data\GSM350983_Cy5_pair_NusA_IP_1.txt',
              'NusA_2'   : PWD + 'Raw_data\GSM351012_Cy3_pair_NusA_IP_2.txt',
              'NusA_3'   : PWD + 'Raw_data\GSM351013_Cy3_pair_NusA_IP_3.txt',
              'NusG_1'   : PWD + 'Raw_data\GSM351002_Cy5_pair_NusG_IP_1.txt',
              'NusG_2'   : PWD + 'Raw_data\GSM351014_Cy3_pair_NusG_IP_2.txt',
              'NusG_3'   : PWD + 'Raw_data\GSM351015_Cy3_pair_NusG_IP_3.txt',
              'RpoC_1'   : PWD + 'Raw_data\GSM350984_Cy5_pair_RpoC_IP_1.txt',
              'RpoC_2'   : PWD + 'Raw_data\GSM351006_Cy5_pair_RpoC_IP_2.txt',
              'RpoC_3'   : PWD + 'Raw_data\GSM351009_Cy3_pair_RpoC_IP_3.txt',              
              'RpoD_1'   : PWD + 'Raw_data\GSM350985_Cy5_pair_RpoD_IP_1.txt',
              'RpoD_2'   : PWD + 'Raw_data\GSM351010_Cy3_pair_RpoD_IP_2.txt',
              'RpoD_3'   : PWD + 'Raw_data\GSM351011_Cy3_pair_RpoD_IP_3.txt', 
              'Rho_1'    : PWD + 'Raw_data\GSM351001_Cy5_pair_Rho_IP_1.txt',
              'Rho_2'    : PWD + 'Raw_data\GSM351007_Cy5_pair_Rho_IP_2.txt', 
              'RpoC_Rif' : PWD + 'Raw_data\GSM351003_Cy5_pair_RpoC_Rif_IP_1.txt',
              'RpoD_Rif' : PWD + 'Raw_data\GSM351004_Cy5_pair_RpoD_Rif_IP_1.txt',
              'RpoC_BCM' : PWD + 'Raw_data\GSM351005_Cy5_pair_RpoC_BCM_IP_1.txt',
              }

#Dictionary of Mock path.
Dict_Mock_path={'NusA_1'   : PWD + 'Raw_data\GSM350983_Cy3_pair_NusA_Mock_1.txt',
                'NusA_2'   : PWD + 'Raw_data\GSM351012_Cy5_pair_NusA_Mock_2.txt',
                'NusA_3'   : PWD + 'Raw_data\GSM351013_Cy5_pair_NusA_Mock_3.txt',
                'NusG_1'   : PWD + 'Raw_data\GSM351002_Cy3_pair_NusG_Mock_1.txt',
                'NusG_2'   : PWD + 'Raw_data\GSM351014_Cy5_pair_NusG_Mock_2.txt',
                'NusG_3'   : PWD + 'Raw_data\GSM351015_Cy5_pair_NusG_Mock_3.txt',
                'RpoC_1'   : PWD + 'Raw_data\GSM350984_Cy3_pair_RpoC_Mock_1.txt',
                'RpoC_2'   : PWD + 'Raw_data\GSM351006_Cy3_pair_RpoC_Mock_2.txt',
                'RpoC_3'   : PWD + 'Raw_data\GSM351009_Cy5_pair_RpoC_Mock_3.txt',              
                'RpoD_1'   : PWD + 'Raw_data\GSM350985_Cy3_pair_RpoD_Mock_1.txt',
                'RpoD_2'   : PWD + 'Raw_data\GSM351010_Cy5_pair_RpoD_Mock_2.txt',
                'RpoD_3'   : PWD + 'Raw_data\GSM351011_Cy5_pair_RpoD_Mock_3.txt', 
                'Rho_1'    : PWD + 'Raw_data\GSM351001_Cy3_pair_Rho_Mock_1.txt',
                'Rho_2'    : PWD + 'Raw_data\GSM351007_Cy3_pair_Rho_Mock_2.txt', 
                'RpoC_Rif' : PWD + 'Raw_data\GSM351003_Cy3_pair_RpoC_Rif_Mock_1.txt',
                'RpoD_Rif' : PWD + 'Raw_data\GSM351004_Cy3_pair_RpoD_Rif_Mock_1.txt',
                'RpoC_BCM' : PWD + 'Raw_data\GSM351005_Cy3_pair_RpoC_BCM_Mock_1.txt',}

#Chip array annotation.
Chip_array_annotation_path=PWD + "Array_description\GPL7790.ndf"

#Outpath for fasta file with sequences of chip probes.
Chip_annot_fasta_outpath=PWD + "Array_description\GPL7790_probes.fasta"

#Outpath for chip probes BLAST results.
Probes_BLAST_res_path=PWD + "Array_description\GPL7790_probes_DY330_blast_res.txt"

#Path to the genome sequence.
Genome=PWD + "Genome\E_coli_w3110_G_Mu.fasta"


#######
#Read ChIP-chip data.
#######

def read_chip_data(pathin):
    #Read input file with chip data.
    filein=open(pathin, 'r')
    Chip_data_dict={}
    for line in filein:
        if line[0]!='#':
            line=line.rstrip().split('\t')
            block=line[1]
            ref_seq_id=line[2]
            prob_id=line[3]
            if ref_seq_id=='U00096':
                signal=float(line[9])
                Chip_data_dict[f'{ref_seq_id}_{block}_{prob_id}']=signal
    filein.close()  
    return Chip_data_dict


#######
#Calculate fold enrichment.
#######

def calc_FU_for_probes(dict_IP, dict_Mock):
    dict_FE={}
    for probe_name, probe_signal_IP in dict_IP.items():
        if probe_name in dict_Mock:
            if dict_Mock[probe_name]>0:
                probe_signal_FE=probe_signal_IP/dict_Mock[probe_name]
                dict_FE[probe_name]=probe_signal_FE
    return dict_FE


#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Create BLAST DB.
    make_blast_db=NcbimakeblastdbCommandline(dbtype="nucl", input_file=genome_path)
    print('Prepare index for target genome: ' + str(make_blast_db))
    make_blast_db()
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
        genomeid=record.id
    genome.close()
    return genomeid, genomefa


#######
#Read chip array annotation, get sequences of probes.
#######

def read_chip_annot_write_fasta(chip_annot_path, fasta_pathout):
    #Read input file with chip annotation.
    filein=open(chip_annot_path, 'r')
    #Write fasta file with oligos sequences (and chip signal).
    fileout=open(fasta_pathout, 'w')
    
    Chip_annot_dict={}
    for line in filein:
        line=line.rstrip().split('\t')
        block=line[1]
        ref_seq_id=line[4]
        prob_id=line[12]
        sequence=line[5]
        if ref_seq_id=='U00096':
            Chip_annot_dict[f'{ref_seq_id}_{block}_{prob_id}']=sequence
            fileout.write(f'>{ref_seq_id}_{block}_{prob_id}\n{sequence}\n')
    
    filein.close()    
    fileout.close()
    return Chip_annot_dict


#######
#Blast flank sequences over reference genome, return coodinates of flanks.
#######

def blast_flanks_return_abs_coords(flanks_pathout, genome_db, flank_blast_res_path):
    #Blast flanks.
    flank_blast=NcbiblastnCommandline(query=flanks_pathout, db=genome_db, out=flank_blast_res_path, outfmt=6, dust='no', task='blastn-short', max_hsps=20, num_threads=6)  
    print('Run blast of flanking sequences: ' + str(flank_blast))
    flank_blast()
    return


#######
#Reads blast result file.
#######

def read_blast_results(blast_data_path, chip_probes_dict):
    
    #Read and take probes coordinates from the blast results.
    filein=open(blast_data_path, 'r')
    
    #Probes coordinates dict.
    Probes_coords_dict={}
    
    line_counter=0
    for line in filein:
        #Count lines.
        line_counter+=1
        if line_counter%10000==0:
            print('Now on ' + str(line_counter) + 'th line')
        #Parse blast results.
        line=line.rstrip().split('\t')
        probe_name=line[0]
        coord_start=min([int(line[8]), int(line[9])])-1
        coord_end=max([int(line[8]), int(line[9])])
        #Get original probe length.
        probe_len=len(chip_probes_dict[probe_name])
        if (float(line[2])==100.0) & ((coord_end-coord_start)>=probe_len-1):
            if probe_name not in Probes_coords_dict:
                Probes_coords_dict[probe_name]=[[coord_start, coord_end]]
            else:
                Probes_coords_dict[probe_name].append([coord_start, coord_end])
        
    filein.close()

    return Probes_coords_dict


#######
#Convert chip signal to WIG for a new reference genome using BLAST data.
#######

def convert_chip_to_wig(set_FE_data, probes_coords_dict, genome_len):
    
    #Create wig array to keep signal FE and coverage data.
    signal_ar=np.array([float(0)]*genome_len)
    coverage_depth_ar=np.array([float(0)]*genome_len)
    
    #Take probe coordinates from BLAST results and FE from chip data.
    for probe_name, probe_FE in set_FE_data.items():
        
        signal=probe_FE
        if probe_name in probes_coords_dict:
            for coord_pair in probes_coords_dict[probe_name]:
                coord_start=coord_pair[0]
                coord_end=coord_pair[1]
                coverage_depth_ar[coord_start:coord_end]+=1
                signal_ar[coord_start:coord_end]+=signal
    
    #Normalize signal by coverage depth.
    signal_ar_norm=[float(0)]*genome_len
    for i in range(len(signal_ar)):
        if coverage_depth_ar[i]!=0:
            signal_ar_norm[i]=signal_ar[i]/coverage_depth_ar[i]
        else:
            signal_ar_norm[i]=0
    return signal_ar, coverage_depth_ar, signal_ar_norm


#######
#Write WIG.
#######

def write_wig(ar, pathout, name, Chromosome_name):
    #Write wig file.
    fileout=open(pathout, 'w')
    fileout.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
    for i in range(len(ar)):
        fileout.write(str(ar[i])+'\n')   
    fileout.close()
    return

#######
#Wrapper function.
#######

def wrapper(pwd, dict_IP_path, dict_Mock_path, chip_array_annotation_path, chip_annot_fasta_outpath, probes_BLAST_res_path, genome):
    
    #Read ChIP-chip data.
    Chip_data_dict={}
    for set_name, set_data_path in dict_IP_path.items():
        print(f'Now reading chip data of the {set_name} dataset')
        Data_IP=read_chip_data(set_data_path)
        Data_Mock=read_chip_data(dict_Mock_path[set_name])
        Data_FE=calc_FU_for_probes(Data_IP, Data_Mock)
        Chip_data_dict[set_name]=Data_FE
    
    #Read chip array annotation, get probe sequence, write to fasta file. 
    Chip_probes_dict=read_chip_annot_write_fasta(chip_array_annotation_path, chip_annot_fasta_outpath)
    
    #Read target genome, create BLAST database.
    genome_id, genome_fasta=read_genome(genome)
    
    #BLAST oligos against target genome.
    blast_flanks_return_abs_coords(chip_annot_fasta_outpath, genome, probes_BLAST_res_path)
    
    #Read probes blast results.
    Probes_coords_dict=read_blast_results(probes_BLAST_res_path, Chip_probes_dict)
    
    #Convert chip-chip signal to WIG format, write wig.
    for set_name, set_FE_data in Chip_data_dict.items():
        print(f'Now converting signal of the {set_name} dataset')
        signal_ar, coverage_depth_ar, signal_ar_norm=convert_chip_to_wig(set_FE_data, Probes_coords_dict, len(genome_fasta))
    
        #Write wig file.
        write_wig(signal_ar,         f'{pwd}FE\{set_name}_raw_signal_FE.wig',    f'{set_name}_raw_signal_FE',    genome_id)
        write_wig(coverage_depth_ar, f'{pwd}FE\{set_name}_probes_cov_depth.wig', f'{set_name}_probes_cov_depth', genome_id)
        write_wig(signal_ar_norm,    f'{pwd}FE\{set_name}_norm_signal_FE.wig',   f'{set_name}_norm_signal_FE',   genome_id)
        
    return

wrapper(PWD, Dict_IP_path, Dict_Mock_path, Chip_array_annotation_path, Chip_annot_fasta_outpath, Probes_BLAST_res_path, Genome)