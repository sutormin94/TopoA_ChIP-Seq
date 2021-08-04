###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to compute fold enrichment of IP over the Mock DNA control (FE).
####

###############################################

import numpy as np

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\BW25113_EcTopoI_mutants_ChIP-Seq\\"
#Path to the file with IP data
IP_path_dict={'1' :  PWD + "WIG\\EcTopoI_BW25113_delta11_IP_1_Jan_4_S28.wig",
              '2' :  PWD + "WIG\\EcTopoI_BW25113_delta11_IP_2_Jan_5_S29.wig",
              '3' :  PWD + "WIG\\EcTopoI_BW25113_delta11_IP_3_Jan_6_S30.wig",  
              '4' :  PWD + "WIG\\EcTopoI_BW25113_delta14_IP_1_Jan_10_S34.wig",
              '5' :  PWD + "WIG\\EcTopoI_BW25113_delta14_IP_2_Jan_11_S35.wig",
              '6' :  PWD + "WIG\\EcTopoI_BW25113_delta14_IP_3_Jan_23_S47.wig", 
              '7' :  PWD + "WIG\\EcTopoI_BW25113_delta30_IP_1_Jan_16_S40.wig",
              '8' :  PWD + "WIG\\EcTopoI_BW25113_delta30_IP_2_Jan_17_S41.wig",
              '9' :  PWD + "WIG\\EcTopoI_BW25113_delta30_IP_3_Jan_18_S42.wig",    
              '10' : PWD + "WIG\\EcTopoI_DY330_IP_1_Jan_22_S46.wig",
              '11' : PWD + "WIG\\EcTopoI_DY330_IP_2_Jan_12_S36.wig",
              '12' : PWD + "WIG\\EcTopoI_DY330_IP_3_Jan_24_S48.wig", 
              }

#Path to the file Mock control data
Mock_path_dict={'1' :  PWD + "WIG\\EcTopoI_BW25113_delta11_Mock_1_Jan_1_S25.wig",
                '2' :  PWD + "WIG\\EcTopoI_BW25113_delta11_Mock_2_Jan_2_S26.wig",
                '3' :  PWD + "WIG\\EcTopoI_BW25113_delta11_Mock_3_Jan_3_S27.wig",  
                '4' :  PWD + "WIG\\EcTopoI_BW25113_delta14_Mock_1_Jan_7_S31.wig",
                '5' :  PWD + "WIG\\EcTopoI_BW25113_delta14_Mock_2_Jan_8_S32.wig",
                '6' :  PWD + "WIG\\EcTopoI_BW25113_delta14_Mock_3_Jan_9_S33.wig", 
                '7' : PWD + "WIG\\EcTopoI_BW25113_delta30_Mock_1_Jan_13_S37.wig",
                '8' : PWD + "WIG\\EcTopoI_BW25113_delta30_Mock_2_Jan_14_S38.wig",
                '9' : PWD + "WIG\\EcTopoI_BW25113_delta30_Mock_3_Jan_15_S39.wig",    
                '10' : PWD + "WIG\\EcTopoI_DY330_Mock_1_Jan_19_S43.wig", 
                '11' : PWD + "WIG\\EcTopoI_DY330_Mock_2_Jan_20_S44.wig",
                '12' : PWD + "WIG\\EcTopoI_DY330_Mock_3_Jan_21_S45.wig", 
                }


#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "EcTopoI_BW25113_delta11_IP_1_Jan_4_S28/EcTopoI_BW25113_delta11_Mock_1_Jan_1_S25",
           '2' :  "EcTopoI_BW25113_delta11_IP_2_Jan_5_S29/EcTopoI_BW25113_delta11_Mock_2_Jan_2_S26", 
           '3' :  "EcTopoI_BW25113_delta11_IP_3_Jan_6_S30/EcTopoI_BW25113_delta11_Mock_3_Jan_3_S27",
           '4' :  "EcTopoI_BW25113_delta14_IP_1_Jan_10_S34/EcTopoI_BW25113_delta14_Mock_1_Jan_7_S31",
           '5' :  "EcTopoI_BW25113_delta14_IP_2_Jan_11_S35/EcTopoI_BW25113_delta14_Mock_2_Jan_8_S32", 
           '6' :  "EcTopoI_BW25113_delta14_IP_3_Jan_23_S47/EcTopoI_BW25113_delta14_Mock_3_Jan_9_S33",
           '7' :  "EcTopoI_BW25113_delta30_IP_1_Jan_16_S40/EcTopoI_BW25113_delta30_Mock_1_Jan_13_S37",
           '8' :  "EcTopoI_BW25113_delta30_IP_2_Jan_17_S41/EcTopoI_BW25113_delta30_Mock_2_Jan_14_S38", 
           '9' :  "EcTopoI_BW25113_delta30_IP_3_Jan_18_S42/EcTopoI_BW25113_delta30_Mock_3_Jan_15_S39",
           '10' : "EcTopoI_DY330_IP_1_Jan_22_S46/EcTopoI_DY330_Mock_1_Jan_19_S43",
           '11' : "EcTopoI_DY330_IP_2_Jan_12_S36/EcTopoI_DY330_Mock_2_Jan_20_S44", 
           '12' : "EcTopoI_DY330_IP_3_Jan_24_S48/EcTopoI_DY330_Mock_3_Jan_21_S45",                                 
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name_manual=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)
#Output path to the final file (fold enrichment).
FE_file_path_dict={'1' :  PWD + "WIG\\EcTopoI_BW25113_delta11_FE_1.wig",
                   '2' :  PWD + "WIG\\EcTopoI_BW25113_delta11_FE_2.wig",
                   '3' :  PWD + "WIG\\EcTopoI_BW25113_delta11_FE_3.wig", 
                   '4' :  PWD + "WIG\\EcTopoI_BW25113_delta14_FE_1.wig",
                   '5' :  PWD + "WIG\\EcTopoI_BW25113_delta14_FE_2.wig",
                   '6' :  PWD + "WIG\\EcTopoI_BW25113_delta14_FE_3.wig",
                   '7' :  PWD + "WIG\\EcTopoI_BW25113_delta30_FE_1.wig",
                   '8' :  PWD + "WIG\\EcTopoI_BW25113_delta30_FE_2.wig",
                   '9' :  PWD + "WIG\\EcTopoI_BW25113_delta30_FE_3.wig",    
                   '10' : PWD + "WIG\\EcTopoI_DY330_FE_1.wig", 
                   '11' : PWD + "WIG\\EcTopoI_DY330_FE_2.wig",
                   '12' : PWD + "WIG\\EcTopoI_DY330_FE_3.wig", 
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
