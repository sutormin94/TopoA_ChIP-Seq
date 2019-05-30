pathin="path\to\the\file"
pathout="path\to\the\file"

filein=open(pathin, 'r')

Dict_of_spacers={}
for line in filein:
    line=line.rstrip().split('\t')
    coordinate=line[0]
    sequence=line[1]
    strand=line[2]
    if sequence not in Dict_of_spacers:
        Dict_of_spacers[sequence]=[1, [coordinate], [strand]]
    else:
        Dict_of_spacers[sequence][0]+=1
        if coordinate not in Dict_of_spacers[sequence][1]:
            Dict_of_spacers[sequence][1].append(coordinate)
        if strand not in Dict_of_spacers[sequence][2]:    
            Dict_of_spacers[sequence][2].append(strand)
filein.close()
     
fileout=open(pathout, 'w')   
Dict_of_spacers_unique_mapping={}  
count_non_unique_mapping=0
count_non_unique_strand=0
for sequence, data in Dict_of_spacers.items():
    check_quality=0
    if len(data[1])>1:
        count_non_unique_mapping+=1
        check_quality+=1
    if len(data[2])>1:
        count_non_unique_strand+=1
        check_quality+=2
    if check_quality==0:
        Dict_of_spacers_unique_mapping[sequence]=data
        fileout.write(f'{data[1][0]}\t{data[0]}\t{sequence}\t{data[2][0]}\n')

fileout.close()

print(f'Number of raw spacers: {len(Dict_of_spacers)}')
print(f'Number of non-unique mapped spacers: {count_non_unique_mapping}')
print(f'Number of non-unique stranded spacers: {count_non_unique_strand}')
print(f'Number of filtered spacers: {len(Dict_of_spacers_unique_mapping)}')
    

