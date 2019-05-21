###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to compute fold enrichment of IP over the Mock DNA control (FE).
####

###############################################


#Path to the file with IP data
IP_path="Path\IP_data.wig"
#Path to the file Mock control data
Mock_path="Path\Mock_data.wig"
#ID or short description of the track (will be the name of a track in IGV).
name=''
#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''
#Output path to the final file (fold enrichment).
FE_file_path="Path\Fold_enrichment.wig"


#######
#Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    return NE_values
    
IP=wig_parsing(IP_path)
Mock=wig_parsing(Mock_path)


FE_out=open(FE_file_path, 'w')
#Write file with fold enrichment data.
FE_out.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\r\nfixedStep chrom='+Chromosome_name+' start=1 step=1\r\n')

for i in range(len(IP)):
    if Mock[i]!=0:
        FE_data_position=IP[i]/Mock[i]
        FE_out.write(str(FE_data_position)+'\r\n')
    else:
        FE_out.write(str(0)+'\r\n')
    
FE_out.close()