###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the input file
filein_path_dict={'1' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit1_S15_ns_fm_ps_nodup.bed",
                  '2' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit2_S16_ns_fm_ps_nodup.bed",
                  '3' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit3_S17_ns_fm_ps_nodup.bed",
                  '4' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit4_S18_ns_fm_ps_nodup.bed",
                  '5' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit5_S19_ns_fm_ps_nodup.bed",
                  '6' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit6_S20_ns_fm_ps_nodup.bed",
                  '7' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit7_S21_ns_fm_ps_nodup.bed",
                  '8' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit8_S22_ns_fm_ps_nodup.bed",
                  '9' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit_9_S97_ns_fm_ps_nodup.bed",
                  '10' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit_10_S98_ns_fm_ps_nodup.bed",
                  '11' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit_11_S99_ns_fm_ps_nodup.bed",
                  '12' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit_12_S100_ns_fm_ps_nodup.bed",  
                  '13' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit_13_S101_ns_fm_ps_nodup.bed", 
                  '14' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit_14_S102_ns_fm_ps_nodup.bed", 
                  '15' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit_15_S103_ns_fm_ps_nodup.bed",
                  '16' : "C:\Sutor\Science\BREX\Cov_depth_nodup\EMit_16_S104_ns_fm_ps_nodup.bed",}

#Path to the output file.
fileout_path_dict={'1' : "C:\Sutor\Science\BREX\WIG_nodup\EMit1_S15_nodup.wig",
                  '2' : "C:\Sutor\Science\BREX\WIG_nodup\EMit2_S16_nodup.wig",
                  '3' : "C:\Sutor\Science\BREX\WIG_nodup\EMit3_S17_nodup.wig",
                  '4' : "C:\Sutor\Science\BREX\WIG_nodup\EMit4_S18_nodup.wig",
                  '5' : "C:\Sutor\Science\BREX\WIG_nodup\EMit5_S19_nodup.wig",
                  '6' : "C:\Sutor\Science\BREX\WIG_nodup\EMit6_S20_nodup.wig",
                  '7' : "C:\Sutor\Science\BREX\WIG_nodup\EMit7_S21_nodup.wig",
                  '8' : "C:\Sutor\Science\BREX\WIG_nodup\EMit8_S22_nodup.wig",
                  '9' : "C:\Sutor\Science\BREX\WIG_nodup\EMit9_S97_nodup.wig",
                  '10' : "C:\Sutor\Science\BREX\WIG_nodup\EMit10_S98_nodup.wig",
                  '11' : "C:\Sutor\Science\BREX\WIG_nodup\EMit11_S99_nodup.wig",
                  '12' : "C:\Sutor\Science\BREX\WIG_nodup\EMit12_S100_nodup.wig",  
                  '13' : "C:\Sutor\Science\BREX\WIG_nodup\EMit13_S101_nodup.wig", 
                  '14' : "C:\Sutor\Science\BREX\WIG_nodup\EMit14_S102_nodup.wig", 
                  '15' : "C:\Sutor\Science\BREX\WIG_nodup\EMit15_S103_nodup.wig",
                  '16' : "C:\Sutor\Science\BREX\WIG_nodup\EMit16_S104_nodup.wig",}

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' : "EMit1_S15_nodup",
           '2' : "EMit2_S16_nodup",
           '3' : "EMit3_S17_nodup",
           '4' : "EMit4_S18_nodup",
           '5' : "EMit5_S19_nodup",
           '6' : "EMit6_S20_nodup",
           '7' : "EMit7_S21_nodup",
           '8' : "EMit8_S22_nodup",
           '9' : "EMit9_S97_nodup",
           '10' : "EMit10_S98_nodup",
           '11' : "EMit11_S99_nodup",
           '12' : "EMit12_S100_nodup",  
           '13' : "EMit13_S101_nodup", 
           '14' : "EMit14_S102_nodup", 
           '15' : "EMit15_S103_nodup",
           '16' : "EMit16_S104_nodup",}

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