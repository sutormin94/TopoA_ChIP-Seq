###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\BW25113_EcTopoI_mutants_ChIP-Seq\\"
#Path to the input file
filein_path_dict={'1' :  PWD + "Cov_depth\\Jan_1_S25.bed",
                  '2' :  PWD + "Cov_depth\\Jan_2_S26.bed",
                  '3' :  PWD + "Cov_depth\\Jan_3_S27.bed",
                  '4' :  PWD + "Cov_depth\\Jan_4_S28.bed",
                  '5' :  PWD + "Cov_depth\\Jan_5_S29.bed",
                  '6' :  PWD + "Cov_depth\\Jan_6_S30.bed",  
                  '7' :  PWD + "Cov_depth\\Jan_7_S31.bed",
                  '8' :  PWD + "Cov_depth\\Jan_8_S32.bed",
                  '9' :  PWD + "Cov_depth\\Jan_9_S33.bed",
                  '10' : PWD + "Cov_depth\\Jan_10_S34.bed",
                  '11' : PWD + "Cov_depth\\Jan_11_S35.bed",
                  '12' : PWD + "Cov_depth\\Jan_12_S36.bed", 
                  '13' : PWD + "Cov_depth\\Jan_13_S37.bed",
                  '14' : PWD + "Cov_depth\\Jan_14_S38.bed",
                  '15' : PWD + "Cov_depth\\Jan_15_S39.bed",
                  '16' : PWD + "Cov_depth\\Jan_16_S40.bed",
                  '17' : PWD + "Cov_depth\\Jan_17_S41.bed",
                  '18' : PWD + "Cov_depth\\Jan_18_S42.bed",    
                  '19' : PWD + "Cov_depth\\Jan_19_S43.bed", 
                  '20' : PWD + "Cov_depth\\Jan_20_S44.bed",
                  '21' : PWD + "Cov_depth\\Jan_21_S45.bed",
                  '22' : PWD + "Cov_depth\\Jan_22_S46.bed",
                  '23' : PWD + "Cov_depth\\Jan_23_S47.bed",
                  '24' : PWD + "Cov_depth\\Jan_24_S48.bed",                   
                  }

#Path to the output file.
fileout_path_dict={'1' :  PWD + "WIG\\EcTopoI_BW25113_delta11_Mock_1_Jan_1_S25.wig",
                   '2' :  PWD + "WIG\\EcTopoI_BW25113_delta11_Mock_2_Jan_2_S26.wig",
                   '3' :  PWD + "WIG\\EcTopoI_BW25113_delta11_Mock_3_Jan_3_S27.wig",
                   '4' :  PWD + "WIG\\EcTopoI_BW25113_delta11_IP_1_Jan_4_S28.wig",
                   '5' :  PWD + "WIG\\EcTopoI_BW25113_delta11_IP_2_Jan_5_S29.wig",
                   '6' :  PWD + "WIG\\EcTopoI_BW25113_delta11_IP_3_Jan_6_S30.wig",  
                   '7' :  PWD + "WIG\\EcTopoI_BW25113_delta14_Mock_1_Jan_7_S31.wig",
                   '8' :  PWD + "WIG\\EcTopoI_BW25113_delta14_Mock_2_Jan_8_S32.wig",
                   '9' :  PWD + "WIG\\EcTopoI_BW25113_delta14_Mock_3_Jan_9_S33.wig",
                   '10' : PWD + "WIG\\EcTopoI_BW25113_delta14_IP_1_Jan_10_S34.wig",
                   '11' : PWD + "WIG\\EcTopoI_BW25113_delta14_IP_2_Jan_11_S35.wig",
                   '12' : PWD + "WIG\\EcTopoI_BW25113_delta14_IP_3_Jan_12_S36.wig", 
                   '13' : PWD + "WIG\\EcTopoI_BW25113_delta30_Mock_1_Jan_13_S37.wig",
                   '14' : PWD + "WIG\\EcTopoI_BW25113_delta30_Mock_2_Jan_14_S38.wig",
                   '15' : PWD + "WIG\\EcTopoI_BW25113_delta30_Mock_3_Jan_15_S39.wig",
                   '16' : PWD + "WIG\\EcTopoI_BW25113_delta30_IP_1_Jan_16_S40.wig",
                   '17' : PWD + "WIG\\EcTopoI_BW25113_delta30_IP_2_Jan_17_S41.wig",
                   '18' : PWD + "WIG\\EcTopoI_BW25113_delta30_IP_3_Jan_18_S42.wig",    
                   '19' : PWD + "WIG\\EcTopoI_DY330_Mock_1_Jan_19_S43.wig", 
                   '20' : PWD + "WIG\\EcTopoI_DY330_Mock_2_Jan_20_S44.wig",
                   '21' : PWD + "WIG\\EcTopoI_DY330_Mock_3_Jan_21_S45.wig",
                   '22' : PWD + "WIG\\EcTopoI_DY330_IP_1_Jan_22_S46.wig",
                   '23' : PWD + "WIG\\EcTopoI_DY330_IP_2_Jan_23_S47.wig",
                   '24' : PWD + "WIG\\EcTopoI_DY330_IP_3_Jan_24_S48.wig",                
                    }

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "EcTopoI_BW25113_delta11_Mock_1_Jan_1_S25",
           '2' :  "EcTopoI_BW25113_delta11_Mock_2_Jan_2_S26",
           '3' :  "EcTopoI_BW25113_delta11_Mock_3_Jan_3_S27",
           '4' :  "EcTopoI_BW25113_delta11_IP_1_Jan_4_S28",
           '5' :  "EcTopoI_BW25113_delta11_IP_2_Jan_5_S29",
           '6' :  "EcTopoI_BW25113_delta11_IP_3_Jan_6_S30",  
           '7' :  "EcTopoI_BW25113_delta14_Mock_1_Jan_7_S31",
           '8' :  "EcTopoI_BW25113_delta14_Mock_2_Jan_8_S32",
           '9' :  "EcTopoI_BW25113_delta14_Mock_3_Jan_9_S33",
           '10' : "EcTopoI_BW25113_delta14_IP_1_Jan_10_S34",
           '11' : "EcTopoI_BW25113_delta14_IP_2_Jan_11_S35",
           '12' : "EcTopoI_BW25113_delta14_IP_3_Jan_12_S36", 
           '13' : "EcTopoI_BW25113_delta30_Mock_1_Jan_13_S37",
           '14' : "EcTopoI_BW25113_delta30_Mock_2_Jan_14_S38",
           '15' : "EcTopoI_BW25113_delta30_Mock_3_Jan_15_S39",
           '16' : "EcTopoI_BW25113_delta30_IP_1_Jan_16_S40",
           '17' : "EcTopoI_BW25113_delta30_IP_2_Jan_17_S41",
           '18' : "EcTopoI_BW25113_delta30_IP_3_Jan_18_S42",    
           '19' : "EcTopoI_DY330_Mock_1_Jan_19_S43", 
           '20' : "EcTopoI_DY330_Mock_2_Jan_20_S44",
           '21' : "EcTopoI_DY330_Mock_3_Jan_21_S45",
           '22' : "EcTopoI_DY330_IP_1_Jan_22_S46",
           '23' : "EcTopoI_DY330_IP_2_Jan_23_S47",
           '24' : "EcTopoI_DY330_IP_3_Jan_24_S48",       
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