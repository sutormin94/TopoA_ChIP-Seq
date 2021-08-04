###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to compute fold enrichment of IP over the Mock DNA control (FE).
####

###############################################

import numpy as np

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\\"
#Path to the file with IP data
IP_path_dict={'1' :  PWD + "WIG\\DSu_14_S14_edt_N3E.wig",
              '2' :  PWD + "WIG\\DSu_16_S16_edt_N3E.wig",
              '3' :  PWD + "WIG\\DSu_18_S18_edt_N3E.wig",  
              '4' :  PWD + "WIG\\DSu_20_S20_edt_N3E.wig",
              '5' :  PWD + "WIG\\DSu_22_S22_edt_N3E.wig",
              '6' :  PWD + "WIG\\DSu_24_S24_edt_N3E.wig", 
              '7' :  PWD + "WIG\\DSu_14_S14_edt_N5E.wig",
              '8' :  PWD + "WIG\\DSu_16_S16_edt_N5E.wig",
              '9' :  PWD + "WIG\\DSu_18_S18_edt_N5E.wig",    
              '10' : PWD + "WIG\\DSu_20_S20_edt_N5E.wig",
              '11' : PWD + "WIG\\DSu_22_S22_edt_N5E.wig",
              '12' : PWD + "WIG\\DSu_24_S24_edt_N5E.wig", 
              }

#Path to the file Mock control data
Mock_path_dict={'1' :  PWD + "WIG\\DSu_13_S13_edt_N3E.wig",
                '2' :  PWD + "WIG\\DSu_15_S15_edt_N3E.wig",
                '3' :  PWD + "WIG\\DSu_17_S17_edt_N3E.wig",  
                '4' :  PWD + "WIG\\DSu_19_S19_edt_N3E.wig",
                '5' :  PWD + "WIG\\DSu_21_S21_edt_N3E.wig",
                '6' :  PWD + "WIG\\DSu_23_S23_edt_N3E.wig", 
                '7' :  PWD + "WIG\\DSu_13_S13_edt_N5E.wig",
                '8' :  PWD + "WIG\\DSu_15_S15_edt_N5E.wig",
                '9' :  PWD + "WIG\\DSu_17_S17_edt_N5E.wig",    
                '10' : PWD + "WIG\\DSu_19_S19_edt_N5E.wig", 
                '11' : PWD + "WIG\\DSu_21_S21_edt_N5E.wig",
                '12' : PWD + "WIG\\DSu_23_S23_edt_N5E.wig", 
                }


#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "DSu_14_S14_edt_N3E/DSu_13_S13_edt_N3E",
           '2' :  "DSu_16_S16_edt_N3E/DSu_15_S15_edt_N3E", 
           '3' :  "DSu_18_S18_edt_N3E/DSu_17_S17_edt_N3E",
           '4' :  "DSu_20_S20_edt_N3E/DSu_19_S19_edt_N3E",
           '5' :  "DSu_22_S22_edt_N3E/DSu_21_S21_edt_N3E", 
           '6' :  "DSu_24_S24_edt_N3E/DSu_23_S23_edt_N3E",
           '7' :  "DSu_14_S14_edt_N5E/DSu_13_S13_edt_N5E",
           '8' :  "DSu_16_S16_edt_N5E/DSu_15_S15_edt_N5E", 
           '9' :  "DSu_18_S18_edt_N5E/DSu_17_S17_edt_N5E",
           '10' : "DSu_20_S20_edt_N5E/DSu_19_S19_edt_N5E",
           '11' : "DSu_22_S22_edt_N5E/DSu_21_S21_edt_N5E", 
           '12' : "DSu_24_S24_edt_N5E/DSu_23_S23_edt_N5E",                                 
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name_manual=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)
#Output path to the final file (fold enrichment).
FE_file_path_dict={'1' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_Ara_N3E_FE_1.wig",
                   '2' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_N3E_FE_1.wig",
                   '3' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_Ara_N3E_FE_2.wig", 
                   '4' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_N3E_FE_2.wig",
                   '5' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_Ara_N3E_FE_3.wig",
                   '6' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_N3E_FE_3.wig",
                   '7' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_Ara_N5E_FE_1.wig",
                   '8' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_N5E_FE_1.wig",
                   '9' :  PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_Ara_N5E_FE_2.wig",    
                   '10' : PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_N5E_FE_2.wig", 
                   '11' : PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_Ara_N5E_FE_3.wig",
                   '12' : PWD + "FE_masked\\EcTopoI_G116S_M320V_Topo_Seq_N5E_FE_3.wig", 
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
        print(f'Estimated pseudocount value: {data_mean/2}')
        data_array=data_array+(data_mean/2)
        data_array_scaled=data_array/(data_mean+(data_mean/2))
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
