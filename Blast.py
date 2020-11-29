###############################################
##Dmitry Sutormin, 2020##
##ChIP-Seq analysis##

####
#Run blast.
####

###############################################

from os import listdir
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\BW25113_topA_mutants\Assembly\\"

DB_to_create=PWD + "Final_assemblies\\"

Query_file=PWD + "topA_search\\topA_full.fasta"

Output_file=PWD + "topA_search\\"


for file in listdir(DB_to_create):
    #Create blast database.
    Make_IG_sequences_db=NcbimakeblastdbCommandline(dbtype="nucl", input_file=DB_to_create + file)    
    print('Making blast database: ' + str(Make_IG_sequences_db))
    Make_IG_sequences_db()
    
    #Blast query.
    Query_blast=NcbiblastnCommandline(query=Query_file, db=DB_to_create + file, out=Output_file + file, outfmt=6)  
    print('Run blast of query sequences: ' + str(Query_blast))
    Query_blast()
       