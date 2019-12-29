###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to compute fold enrichment of IP over the Mock DNA control (FE).
####

###############################################

import numpy as np

#Path to the file with IP data
IP_path_dict={'1' : "C:\Sutor\Science\BREX\WIG_nodup\EMit2_S16_nodup.wig",
              '2' : "C:\Sutor\Science\BREX\WIG_nodup\EMit4_S18_nodup.wig",
              '3' : "C:\Sutor\Science\BREX\WIG_nodup\EMit6_S20_nodup.wig",
              '4' : "C:\Sutor\Science\BREX\WIG_nodup\EMit8_S22_nodup.wig",
              '5' : "C:\Sutor\Science\BREX\WIG_nodup\EMit9_S97_nodup.wig",
              '6' : "C:\Sutor\Science\BREX\WIG_nodup\EMit10_S98_nodup.wig",
              '7' : "C:\Sutor\Science\BREX\WIG_nodup\EMit11_S99_nodup.wig",
              '8' : "C:\Sutor\Science\BREX\WIG_nodup\EMit12_S100_nodup.wig",}

#Path to the file Mock control data
Mock_path_dict={'1' : "C:\Sutor\Science\BREX\WIG_nodup\EMit1_S15_nodup.wig",
                '2' : "C:\Sutor\Science\BREX\WIG_nodup\EMit3_S17_nodup.wig",
                '3' : "C:\Sutor\Science\BREX\WIG_nodup\EMit5_S19_nodup.wig",
                '4' : "C:\Sutor\Science\BREX\WIG_nodup\EMit7_S21_nodup.wig",
                '5' : "C:\Sutor\Science\BREX\WIG_nodup\EMit13_S101_nodup.wig",
                '6' : "C:\Sutor\Science\BREX\WIG_nodup\EMit14_S102_nodup.wig",
                '7' : "C:\Sutor\Science\BREX\WIG_nodup\EMit15_S103_nodup.wig",
                '8' : "C:\Sutor\Science\BREX\WIG_nodup\EMit16_S104_nodup.wig",}


#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' : "EMit2_S16_nodup/EMit1_S15_nodup",
           '2' : "EMit4_S18_nodup/EMit3_S17_nodup",
           '3' : "EMit6_S20_nodup/EMit5_S19_nodup",
           '4' : "EMit8_S22_nodup/EMit7_S21_nodup",
           '5' : "EMit9_S97_nodup/EMit13_S101_nodup",
           '6' : "EMit10_S98_nodup/EMit14_S102_nodup",  
           '7' : "EMit11_S99_nodup/EMit15_S103_nodup",
           '8' : "EMit12_S100_nodup/EMit16_S104_nodup",}

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name_manual=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)
#Output path to the final file (fold enrichment).
FE_file_path_dict={'1' : "C:\Sutor\Science\BREX\Fold_enrichment_nodup\EMit2_1_nodup.wig",
                   '2' : "C:\Sutor\Science\BREX\Fold_enrichment_nodup\EMit4_3_nodup.wig",
                   '3' : "C:\Sutor\Science\BREX\Fold_enrichment_nodup\EMit6_5_nodup.wig",
                   '4' : "C:\Sutor\Science\BREX\Fold_enrichment_nodup\EMit8_7_nodup.wig",
                   '5' : "C:\Sutor\Science\BREX\Fold_enrichment_nodup\EMit9_13_nodup.wig",
                   '6' : "C:\Sutor\Science\BREX\Fold_enrichment_nodup\EMit10_14_nodup.wig",  
                   '7' : "C:\Sutor\Science\BREX\Fold_enrichment_nodup\EMit11_15_nodup.wig",
                   '8' : "C:\Sutor\Science\BREX\Fold_enrichment_nodup\EMit12_16_nodup.wig",}


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
