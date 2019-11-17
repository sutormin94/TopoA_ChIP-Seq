###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the input file
filein_path_dict={'1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\ChiP_CTDplusRif-1_S55.bed",
                  '2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\ChiP_CTDplusRif-2_S53.bed",
                  '3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\ChiP_CTDplusRif-3_S51.bed",
                  '4' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\ChiP_CTDplusRifplus2_S59.bed",
                  '5' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\ChiP_CTDplusRifplus3_S57.bed",
                  '6' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\ChiP_CTD-Rif-3_S47.bed",
                  '7' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\ChiP_CTD-Rifplus3_S49.bed",
                  '8' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\Mock_CTDplusRif-1_S54.bed",
                  '9' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\Mock_CTDplusRif-2_S52.bed",
                  '10' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\Mock_CTDplusRif-3_S50.bed",
                  '11' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\Mock_CTDplusRifplus1_S60.bed",
                  '12' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\Mock_CTDplusRifplus2_S58.bed",  
                  '13' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\Mock_CTDplusRifplus3_S56.bed", 
                  '14' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\Mock_CTD-Rif-3_S46.bed", 
                  '15' : "C:\Sutor\Science\TopoI-ChIP-Seq\Cov_depth\Mock_CTD-Rifplus3_S48.bed", }

#Path to the output file.
fileout_path_dict={'1' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\ChiP_CTDplusRif-1_S55.wig",
                  '2' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\ChiP_CTDplusRif-2_S53.wig",
                  '3' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\ChiP_CTDplusRif-3_S51.wig",
                  '4' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\ChiP_CTDplusRifplus2_S59.wig",
                  '5' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\ChiP_CTDplusRifplus3_S57.wig",
                  '6' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\ChiP_CTD-Rif-3_S47.wig",
                  '7' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\ChiP_CTD-Rifplus3_S49.wig",
                  '8' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\Mock_CTDplusRif-1_S54.wig",
                  '9' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\Mock_CTDplusRif-2_S52.wig",
                  '10' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\Mock_CTDplusRif-3_S50.wig",
                  '11' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\Mock_CTDplusRifplus1_S60.wig",
                  '12' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\Mock_CTDplusRifplus2_S58.wig",  
                  '13' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\Mock_CTDplusRifplus3_S56.wig", 
                  '14' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\Mock_CTD-Rif-3_S46.wig", 
                  '15' : "C:\Sutor\Science\TopoI-ChIP-Seq\WIG\Mock_CTD-Rifplus3_S48.wig", }
#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' : "ChiP_CTDplusRif-1_S55",
           '2' : "ChiP_CTDplusRif-2_S53",
           '3' : "ChiP_CTDplusRif-3_S51",
           '4' : "ChiP_CTDplusRifplus2_S59",
           '5' : "ChiP_CTDplusRifplus3_S57",
           '6' : "ChiP_CTD-Rif-3_S47",
           '7' : "ChiP_CTD-Rifplus3_S49",
           '8' : "Mock_CTDplusRif-1_S54",
           '9' : "Mock_CTDplusRif-2_S52",
           '10' : "Mock_CTDplusRif-3_S50",
           '11' : "Mock_CTDplusRifplus1_S60",
           '12' : "Mock_CTDplusRifplus2_S58",  
           '13' : "Mock_CTDplusRifplus3_S56", 
           '14' : "Mock_CTD-Rif-3_S46", 
           '15' : "Mock_CTD-Rifplus3_S48", }
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