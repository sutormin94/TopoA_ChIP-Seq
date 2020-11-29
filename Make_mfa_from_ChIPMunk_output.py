###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes output from ChIP-Munk (table with seqs)
#and prepares multy-fasta file ready for logo visualization (e.g. with WebLogo online service).
###############################################

#######
#Packages to be imported.
#######

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio import AlignIO, motifs
from Bio.Align import AlignInfo
import weblogo
from weblogo import LogoOptions, LogoFormat, eps_formatter, read_seq_data, LogoData, png_formatter, pdf_formatter

#######
#Variables to be defined.
#######

#Filename.
Filename="EcTopoI_CTD_Rif_rep12_nm_0.001_peaks_motifs_3"

#Path to ChIP-Munk results and multi-fasta output.
Data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\\"

def read_convert(datapath, filename):
    input_path=datapath + filename + '.txt'
    output_path=datapath + filename + ".fasta"
    output_path_RC=datapath + filename + "_RC.fasta"
    #Read ChIP-Munk output, extract sequences and write mfa file.
    filein=open(input_path, 'r')
    fileout=open(output_path, 'w')
    fileout_RC=open(output_path_RC, 'w')
    for line in filein:
        line=line.rstrip().split('|')
        if line[0]=='WORD':
            data=line[1].split('\t')
            seq=data[2]
            seq_bio=Seq(seq, generic_dna)
            seq_bio_rc=seq_bio.reverse_complement()
            fileout.write(f'>{data[0]}\n{seq}\n')
            fileout_RC.write(f'>{data[0]}\n{str(seq_bio_rc)}\n')
    filein.close()
    fileout.close()
    fileout_RC.close()
    print('ChIP-Munk output converted to mfa succesfully.')
    
    #Read mfa file, identify consensus sequences.
    alignment=AlignIO.read(output_path, "fasta")
    alignment_summary=AlignInfo.SummaryInfo(alignment)
    consensus=alignment_summary.dumb_consensus(threshold=0.35,  ambiguous='X')
    print('Consensus sequence:\n' + consensus)
    consensus_rc=consensus.reverse_complement()
    print('Reverse-complement consensus sequence:\n' + consensus_rc)
    print('Done!')
    
    #Read mfa file, draw motif. + strand.
    MFA_data=open(output_path)
    MFA_seqs=read_seq_data(MFA_data)
    logodata=LogoData.from_seqs(MFA_seqs)
    logooptions=LogoOptions()
    logooptions.show_errorbars=False
    logoformat=LogoFormat(logodata, logooptions)
    pdf=weblogo.logo_formatter.pdf_formatter(logodata, logoformat)
    logout=open(datapath + filename + ".pdf", 'wb')
    logout.write(pdf)
    logout.close()
    
    #Read mfa file, draw motif. - strand.
    MFA_data_RC=open(output_path_RC)
    MFA_seqs_RC=read_seq_data(MFA_data_RC)
    logodata_RC=LogoData.from_seqs(MFA_seqs_RC)
    logooptions=LogoOptions()
    logooptions.show_errorbars=False
    logoformat=LogoFormat(logodata_RC, logooptions)
    pdf_RC=weblogo.logo_formatter.pdf_formatter(logodata_RC, logoformat)
    logout=open(datapath + filename + "_RC.pdf", 'wb')
    logout.write(pdf_RC)
    logout.close()    
    
    
    return

read_convert(Data_path, Filename)