###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the input file
filein_path="F:\V_cholerae_ToxR\Cov_depth_nodup\SRR2188532_1_name_sorted_fm_ps_nd.bed"
#Path to the output file.
fileout_path="F:\V_cholerae_ToxR\Cov_depth_nodup\ToxR_ChIP_control_nodup.wig"
#ID or short description of the track (will be the name of a track in IGV).
name='ToxR_ChIP_control_nodup'
#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)

filein=open(filein_path, 'r')
fileout=open(fileout_path, 'w')


Ar_of_Cromosome_names=[]
for line in filein:
    line=line.rstrip().split('\t')
    if line[0] not in Ar_of_Cromosome_names:
        if Auto_or_manual==0:
            fileout.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+str(line[0])+' start=1 step=1\n')
        elif Auto_or_manual==1:
            fileout.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
        Ar_of_Cromosome_names.append(line[0])
    else:
        fileout.write(line[2]+'\n')
    
filein.close()
fileout.close()