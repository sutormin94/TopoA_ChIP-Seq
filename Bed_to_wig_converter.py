###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the input file
filein_path="Path\file.bed"
#Path to the output file.
fileout_path="Path\file.wig"
#ID or short description of the track (will be the name of a track in IGV).
name=''
#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''

filein=open(filein_path, 'r')
fileout=open(fileout_path, 'w')

fileout.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\r\nfixedStep chrom='+Chromosome_name+' start=1 step=1\r\n')

for line in filein:
    line=line.rstrip().split('\t')
    fileout.write(line[2]+'\r\n')
    
filein.close()
fileout.close()