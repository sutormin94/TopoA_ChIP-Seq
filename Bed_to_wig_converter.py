###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\BW25113_topA_mutants\Set3\\"
#Path to the input file
filein_path_dict={'1' :  PWD + "Cov_depth\\D10_S1_L001.bed",
                  '2' :  PWD + "Cov_depth\\D11_S2_L001.bed",
                  '3' :  PWD + "Cov_depth\\D12_S3_L001.bed",
                  '4' :  PWD + "Cov_depth\\D13_S4_L001.bed",
                  '5' :  PWD + "Cov_depth\\D14_S5_L001.bed",
                  '6' :  PWD + "Cov_depth\\D15_S6_L001.bed",
                  '7' :  PWD + "Cov_depth\\D16_S7_L001.bed",
                  '8' :  PWD + "Cov_depth\\D17_S8_L001.bed",  
                  '9' :  PWD + "Cov_depth\\D18_S9_L001.bed",
                  '10' :  PWD + "Cov_depth\\D19_S10_L001.bed",                
                  }

#Path to the output file.
fileout_path_dict={'1' :  PWD + "Cov_depth\\D10_S1_L001.wig",
                   '2' :  PWD + "Cov_depth\\D11_S2_L001.wig",
                   '3' :  PWD + "Cov_depth\\D12_S3_L001.wig",
                   '4' :  PWD + "Cov_depth\\D13_S4_L001.wig",
                   '5' :  PWD + "Cov_depth\\D14_S5_L001.wig",
                   '6' :  PWD + "Cov_depth\\D15_S6_L001.wig",
                   '7' :  PWD + "Cov_depth\\D16_S7_L001.wig",
                   '8' :  PWD + "Cov_depth\\D17_S8_L001.wig",  
                   '9' :  PWD + "Cov_depth\\D18_S9_L001.wig",
                   '10' :  PWD + "Cov_depth\\D19_S10_L001.wig",                 
                    }

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  'D10_S1',
           '2' :  'D11_S2',
           '3' :  'D12_S3',
           '4' :  'D13_S4',
           '5' :  'D14_S5',
           '6' :  'D15_S6',
           '7' :  'D16_S7',
           '8' :  'D17_S8',  
           '9' :  'D18_S9',
          '10' :  'D19_S10'         
           }      

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)


def read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual):
    for sample_name, sample_path in filein_path_dict.items():
        print(f'Now is processing: {sample_path}')
        print(f'Progress: {sample_name}/{len(filein_path_dict)}')
        
        filein=open(filein_path_dict[sample_name], 'r')
        fileout=open(fileout_path_dict[sample_name], 'w')
        
        Ar_of_Cromosome_names=[]
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in Ar_of_Cromosome_names:
                if Auto_or_manual==0:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+str(line[0])+' start=1 step=1\n')
                elif Auto_or_manual==1:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
                Ar_of_Cromosome_names.append(line[0])
            else:
                fileout.write(line[2]+'\n')
            
        filein.close()
        fileout.close()    
    return


read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual)