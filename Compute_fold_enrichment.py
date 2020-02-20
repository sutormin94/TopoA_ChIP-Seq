###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to compute fold enrichment of IP over the Mock DNA control (FE).
####

###############################################

import numpy as np

#Path to the working directory.
PWD="C:\Sutor\Science\E_coli_ChIP-Seqs\\"
#Path to the file with IP data
IP_path_dict={'1' :  PWD +        "Singh_CRP\WIG\Singh_CRP_EE_IP_Rep1_SRR5099114.wig",
              '2' :  PWD +        "Singh_CRP\WIG\Singh_CRP_EE_IP_Rep2_SRR5099115.wig",  
              '3' :  PWD +        "Singh_CRP\WIG\Singh_CRP_ME_IP_Rep1_SRR5099118.wig",
              '4' :  PWD +        "Singh_CRP\WIG\Singh_CRP_ME_IP_Rep2_SRR5099119.wig", 
              '5' :  PWD +        "Singh_CRP\WIG_nodup\Singh_CRP_EE_IP_Rep1_SRR5099114_nodup.wig",
              '6' :  PWD +        "Singh_CRP\WIG_nodup\Singh_CRP_EE_IP_Rep2_SRR5099115_nodup.wig",
              '7' :  PWD +        "Singh_CRP\WIG_nodup\Singh_CRP_ME_IP_Rep1_SRR5099118_nodup.wig",
              '8' :  PWD +        "Singh_CRP\WIG_nodup\Singh_CRP_ME_IP_Rep2_SRR5099119_nodup.wig",              
              }

#Path to the file Mock control data
Mock_path_dict={'1' :  PWD +           "Singh_CRP\WIG\Singh_CRP_EE_Mock_Rep1_SRR5099116.wig",
                '2' :  PWD +           "Singh_CRP\WIG\Singh_CRP_EE_Mock_Rep2_SRR5099117.wig",
                '3' :  PWD +           "Singh_CRP\WIG\Singh_CRP_ME_Mock_Rep1_SRR5099120.wig",
                '4' :  PWD +           "Singh_CRP\WIG\Singh_CRP_ME_Mock_Rep2_SRR5099121.wig",   
                '5' :  PWD +           "Singh_CRP\WIG_nodup\Singh_CRP_EE_Mock_Rep1_SRR5099116_nodup.wig",
                '6' :  PWD +           "Singh_CRP\WIG_nodup\Singh_CRP_EE_Mock_Rep2_SRR5099117_nodup.wig",
                '7' :  PWD +           "Singh_CRP\WIG_nodup\Singh_CRP_ME_Mock_Rep1_SRR5099120_nodup.wig",
                '8' :  PWD +           "Singh_CRP\WIG_nodup\Singh_CRP_ME_Mock_Rep2_SRR5099121_nodup.wig",                
                }


#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "Singh_CRP_EE_IP_Rep1/Singh_CRP_EE_Mock_Rep1",
           '2' :  "Singh_CRP_EE_IP_Rep2/Singh_CRP_EE_Mock_Rep2", 
           '3' :  "Singh_CRP_ME_IP_Rep1/Singh_CRP_ME_Mock_Rep1",
           '4' :  "Singh_CRP_ME_IP_Rep2/Singh_CRP_ME_Mock_Rep2", 
           '5' :  "Singh_CRP_EE_IP_Rep1_nodup/Singh_CRP_EE_Mock_Rep1_nodup",
           '6' :  "Singh_CRP_EE_IP_Rep2_nodup/Singh_CRP_EE_Mock_Rep2_nodup",
           '7' :  "Singh_CRP_ME_IP_Rep1_nodup/Singh_CRP_ME_Mock_Rep1_nodup",
           '8' :  "Singh_CRP_ME_IP_Rep2_nodup/Singh_CRP_ME_Mock_Rep2_nodup",            
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name_manual=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)
#Output path to the final file (fold enrichment).
FE_file_path_dict={'1' :  PWD +           "Singh_CRP\FE\Singh_CRP_EE_FE_Rep1.wig",
                   '2' :  PWD +           "Singh_CRP\FE\Singh_CRP_EE_FE_Rep2.wig",
                   '3' :  PWD +           "Singh_CRP\FE\Singh_CRP_ME_FE_Rep1.wig",
                   '4' :  PWD +           "Singh_CRP\FE\Singh_CRP_ME_FE_Rep2.wig",
                   '5' :  PWD +           "Singh_CRP\FE_nodup\Singh_CRP_EE_FE_Rep1_nodup.wig",
                   '6' :  PWD +           "Singh_CRP\FE_nodup\Singh_CRP_EE_FE_Rep2_nodup.wig",
                   '7' :  PWD +           "Singh_CRP\FE_nodup\Singh_CRP_ME_FE_Rep1_nodup.wig",
                   '8' :  PWD +           "Singh_CRP\FE_nodup\Singh_CRP_ME_FE_Rep2_nodup.wig",                   
                   }


#######
#Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    Dict_of_chromosomes_data={}
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0]=='fixedStep':
            chrom_name=line[1].split('=')[1]
            Dict_of_chromosomes_data[chrom_name]=[]
        if line[0] not in ['track', 'fixedStep']:
            Dict_of_chromosomes_data[chrom_name].append(float(line[0]))
    wigin.close()
    
    for Chromosome_name, data in Dict_of_chromosomes_data.items():
        data_array=np.array(data)
        data_mean=np.mean(data_array)
        print(f'Mean coverage of {Chromosome_name}: {data_mean}')
        data_array_scaled=data_array/data_mean
        Dict_of_chromosomes_data[Chromosome_name]=data_array_scaled
    return Dict_of_chromosomes_data


def read_files(input_dict):
    Data_dict={}
    for name, path in input_dict.items():
        Data_dict[name]=wig_parsing(path)
        print(f'Progress: {name}/{len(input_dict)}')
    return Data_dict

IP_dict=read_files(IP_path_dict)
Mock_dict=read_files(Mock_path_dict)


def divide_write(IP_dict, Mock_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict):
    for sample_name, sample_data in IP_dict.items():
        print(f'Now is processing: {sample_name}')
        print(f'Progress: {sample_name}/{len(IP_dict)}')
        FE_out=open(FE_file_path_dict[sample_name], 'w')
        #Write file with fold enrichment data.
        for Chromosome_name, data in sample_data.items():
            print(f'Average normalized covarage of IP: {np.mean(data)}')
            print(f'Average normalized covarage of Mock: {np.mean(Mock_dict[sample_name][Chromosome_name])}')
            if Auto_or_manual==0:
                FE_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
            elif Auto_or_manual==1:
                FE_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name_manual+' start=1 step=1\n')
            for i in range(len(data)):
                if Mock_dict[sample_name][Chromosome_name][i]!=0:
                    FE_data_position=IP_dict[sample_name][Chromosome_name][i]/Mock_dict[sample_name][Chromosome_name][i]
                    FE_out.write(str(FE_data_position)+'\n')
                else:
                    FE_out.write(str(0)+'\n')
            
        FE_out.close()        
    return

divide_write(IP_dict, Mock_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict)
