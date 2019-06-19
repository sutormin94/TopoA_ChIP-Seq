###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the input file
filein_path_dict={'1' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970849.bed",
                  '2' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970850.bed",
                  '3' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970851.bed",
                  '4' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970852.bed",
                  '5' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970853.bed",
                  '6' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970854.bed",
                  '7' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970855.bed",
                  '8' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970856.bed",
                  '9' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970857.bed",
                  '10' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970858.bed",
                  '11' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970859.bed",
                  '12' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970860.bed",  
                  '13' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970861.bed", 
                  '14' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970862.bed", 
                  '15' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970863.bed", 
                  '16' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970864.bed", 
                  '17' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970865.bed", 
                  '18' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970866.bed",
                  '19' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970867.bed", 
                  '20' : "F:\Sayyed_TopoIV_data\Cov_depth\SRR2970868.bed"}

#Path to the output file.
fileout_path_dict={'1' : "F:\Sayyed_TopoIV_data\Cov_depth\ParC_IP_1.wig",
                   '2' : "F:\Sayyed_TopoIV_data\Cov_depth\ParC_Input_1.wig",
                   '3' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_IP_1.wig",
                   '4' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_Input_1.wig",
                   '5' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_IP_2.wig",
                   '6' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_Input_2.wig",
                   '7' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_IP_G1.wig",
                   '8' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_Input_G1.wig",
                   '9' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_IP_S20min.wig",
                   '10' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_Input_S20min.wig",
                   '11' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_IP_S40min.wig",
                   '12' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_Input_S40min.wig",  
                   '13' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_IP_G2.wig", 
                   '14' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_Input_G2.wig", 
                   '15' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParC_IP_1.wig", 
                   '16' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParC_Input_1.wig", 
                   '17' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParE_IP_1.wig", 
                   '18' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParE_Input_1.wig",
                   '19' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParE_IP_2.wig", 
                   '20' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParE_Input_2.wig"}
#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' : "ParC_IP_1",
           '2' : "ParC_Input_1",
           '3' : "ParE_IP_1",
           '4' : "ParE_Input_1",
           '5' : "ParE_IP_2",
           '6' : "ParE_Input_2",
           '7' : "ParE_IP_G1",
           '8' : "ParE_Input_G1",
           '9' : "ParE_IP_S20min",
           '10' : "ParE_Input_S20min",
           '11' : "ParE_IP_S40min",
           '12' : "ParE_Input_S40min",  
           '13' : "ParE_IP_G2", 
           '14' : "ParE_Input_G2", 
           '15' : "NorflIP_ParC_IP_1", 
           '16' : "NorflIP_ParC_Input_1", 
           '17' : "NorflIP_ParE_IP_1", 
           '18' : "NorflIP_ParE_Input_1",
           '19' : "NorflIP_ParE_IP_2", 
           '20' : "NorflIP_ParE_Input_2"}
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