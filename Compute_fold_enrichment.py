###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to compute fold enrichment of IP over the Mock DNA control (FE).
####

###############################################


#Path to the file with IP data
IP_path="F:\Gyrase_time-course_experiment\Reads_eq\WIG_files\Raw_wig\Coverage\+IP+Cfx_R2\DSu_48_S120_eq_edt_for_rev_depth.wig"
#Path to the file Mock control data
Mock_path="F:\Gyrase_time-course_experiment\Reads_eq\WIG_files\Raw_wig\Coverage\-IP+Cfx_R2\DSu_40_S112_eq_edt_for_rev_depth.wig"
#ID or short description of the track (will be the name of a track in IGV).
name='Gyrase_tc_Topo-Seq_30min_2_FE'
#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''
#Output path to the final file (fold enrichment).
FE_file_path="F:\Gyrase_time-course_experiment\Reads_eq\WIG_files\Fold_enrichment\Gyrase_tc_Topo-Seq_30min_2_FE.wig"


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
    return Dict_of_chromosomes_data
    
IP=wig_parsing(IP_path)
Mock=wig_parsing(Mock_path)


FE_out=open(FE_file_path, 'w')
#Write file with fold enrichment data.
for Chromosome_name, data in IP.items():
    FE_out.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
    for i in range(len(data)):
        if Mock[Chromosome_name][i]!=0:
            FE_data_position=IP[Chromosome_name][i]/Mock[Chromosome_name][i]
            FE_out.write(str(FE_data_position)+'\n')
        else:
            FE_out.write(str(0)+'\n')
    
FE_out.close()