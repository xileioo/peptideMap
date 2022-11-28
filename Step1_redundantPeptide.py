#!/usr/bin/python
# 2022-22-07 by lixue
# example： python Step1_redundantPeptide.py KVLEH***V 4

# need install biopython
# conda install -c conda-forge biopython
import sys
from Bio import SeqIO
import re
from collections import defaultdict

def window(fseq, window_size=9):
    """
    Using sliding window on protein sequence to get kmer
    example: 
    for seq,start in window('CGAGDGHCAAGDTHNUKHGFCAAGD',9): print(seq + ';' + str(start) + '\n')
    """
    for i in range(len(fseq) - window_size + 1):
        yield fseq[i:i+window_size], i

def peptideCompare(kmer, qury_str):
    """
    Calculate the match count of amino acid for peptide fragment compare each kmer
    X or * represent any amino acid
    example: 
    print( peptideCompare('CGAGDGHCAAGD','**AGAG**AA*D')) # output should be 6
    """
    matchCount = 0
    for i in range(len(qury_str)):
        if qury_str[i] == '*' or qury_str[i] == 'X':
            continue
        elif qury_str[i] == kmer[i]:
            matchCount = matchCount + 1
        else:
            continue
    return matchCount


#-------------- define input and parameter -------------------------------------------------------
ref_file = 'ReferenceProteomes_EMBL_EBI/UP000005640_9606.fasta'  # huamn reference - no need to change
#query_peptide_fragment = 'KVLEH***V'   # query mapping peptide fragment
#atleast_matchAA_count = 4                     # query_peptide compare 9-kmer that match at least three/four AA
query_peptide_fragment = sys.argv[1]
atleast_matchAA_count  = sys.argv[2]
output_file = sys.argv[3]
peptide_length = len(query_peptide_fragment)      # peptide length

#output_file = 'kmerSearchRes_XXLEHVVXX_9kmer_atleastMatch4.csv'
#output_file =  'Step1_redundantPeptide_' + re.sub('\\*', 'X', query_peptide_fragment) + '_' + str(peptide_length) + 'kmer' + '_atleastMatch' + str(atleast_matchAA_count) + '.csv'
print('Your query peptide is: ' + query_peptide_fragment)
print('At least match amino acid: ' + str(atleast_matchAA_count))
print('Save output to FILE: ' + output_file)


# ----------- main ------------------------------------------------------------- 
myout = open(output_file, 'w')
myout.write('Protein,matchCount,Peptide_start,Peptide_sequence,Protein_sequence\n')
fasta_seq = SeqIO.parse(open(ref_file), 'fasta')

#i = 0
for fasta in fasta_seq:
    #i = i + 1
    #if(i == 300):
    #    break
    name, sequence = fasta.id, str(fasta.seq)
    for seq,start_loci in window(sequence,peptide_length):
        matchCount = peptideCompare(seq, query_peptide_fragment)
        if matchCount >= int(atleast_matchAA_count):
            myout.write(name + ',' + str(matchCount) + ',' + str(start_loci) + ',' + seq + ',' + sequence + '\n') # output_file

myout.close()

