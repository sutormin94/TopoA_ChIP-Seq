###############################################
##Dmitry Sutormin, 2020##
##ChIP-Seq analysis##

####
#Converts wig file with some genomic feature to Circos-compatible format (bed-like).
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import scipy

#Input files.
Input_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig"

#Bin width, nt.
Bin_width=100

#Chromosome ID.
Chromosome_id='chr1'

#Output file.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Manuals\circos-tutorials-0.67\E_coli_test\E_coli_w3110_Mu_GC_133bp_bin_100.txt"


#######
##Parses WIG file.
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

#######
##Read wig, bin, convert to Circus format.
#######

def read_bin_write(inpath, bwidth, chrid, outpath):
    #Read wig data.
    wig_data=wig_parsing(inpath)
    
    #Bin wig data.
    wig_binned=[]
    start=0
    while start+bwidth<len(wig_data):
        wig_binned.append(np.mean(wig_data[start:(start+bwidth)]))
        start+=bwidth
    wig_binned.append(np.mean(wig_data[start:]))
    
    #Write binned data in Circos format.
    outfile=open(outpath, 'w')
    for i in range(len(wig_binned)):
        if ((i+1)*bwidth)<len(wig_data):
            outfile.write(f'{chrid} {i*bwidth} {(i+1)*bwidth} {wig_binned[i]}\n')
        else:
            outfile.write(f'{chrid} {i*bwidth} {len(wig_data)} {wig_binned[i]}\n')
    
    outfile.close()
    return

read_bin_write(Input_path, Bin_width, Chromosome_id, Output_path)