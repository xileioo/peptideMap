#!/usr/bin/python
# 2022-22-07 by lixue
import re
import pickle

proteinID_dict = dict()
ID_convert_FILE = 'uniprot_proteinEntryName_to_GeneNames.tsv'
ID_convert_file = open(ID_convert_FILE, 'r')
next(ID_convert_file)

ID_output_FILE = 'proteinID_geneID_table.csv'
ID_output_file = open(ID_output_FILE, 'w')
ID_output_file.write('proteinID,geneID\n')

count = 0
for line in ID_convert_file:
    currentline = line.split('\t')
    if gene := re.match('\w+', currentline[4]):
        geneID = gene.group(0)
        proteinID = currentline[2]
        ID_output_file.write(proteinID + ',' + geneID + '\n')
        proteinID_dict[proteinID] = geneID
        count += 1
print(count)


with open('proteinID_geneID_dict.pickle', 'wb') as handle:
    pickle.dump(proteinID_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('proteinID_geneID_dict.pickle', 'rb') as handle:
    proteinID_dict2 = pickle.load(handle)

print(proteinID_dict == proteinID_dict2)
print(proteinID_dict2['WDR70_HUMAN'])

