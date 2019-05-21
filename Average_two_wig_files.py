###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to compute by-position average of a set of wig files.
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np

#Dictionary of replicas 
#'Replica name' : 'Path to wig file'
Dict_of_replicas={'Replic 1' : "Path\Replica_1.wig",
                  'Replic 2' : "Path\Replica_2.wig",}

#ID or short description of the track (will be the name of a track in IGV).
name=''
#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''
#Output path two the final file.
average_file_path="Path\Average.wig"


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

#Contains data of all replicas in separate arrays.
ar_of_replicas=[]
for replica_name, replica_path in Dict_of_replicas.items():
    ar_of_replicas.append(wig_parsing(replica_path))


#Write file with avaraged data.
average_out=open(average_file_path, 'w')
average_out.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\r\nfixedStep chrom='+Chromosome_name+' start=1 step=1\r\n')

for i in range(len(ar_of_replicas[0])):
    av_data_position=[]
    for ar in ar_of_replicas:
        av_data_position.append(ar[i])
    average_out.write(str(np.mean(av_data_position))+'\r\n')

average_out.close()