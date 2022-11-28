#!/usr/bin/python
# 2022-10-26 by lixue
# demo: python Step3_Protein_expression.py step1Out expFile outFile

import sys
import re
import pandas as pd
import numpy as np


# define input and parameter
step1Out = sys.argv[1]
expFile = sys.argv[2]
outFile = sys.argv[3]
#step1Out = 'D:/peptideMap/KVLEHXXXV_atleastMatch4/Step1_matchedRedundantPeptides.csv'
#expFile = 'D:/peptideMap/exp_HPA/My_Merged_normal_and_rna_tissue_consensus.tsv'
#outFile = "D:/peptideMap/KVLEHXXXV_atleastMatch4/Step3_out.csv"

#----------------------------------------------------------
# read IHC-exp-file and Step1-peptide-list
#----------------------------------------------------------
EXP = pd.read_table(expFile, sep='\t')
PEP = pd.read_csv(step1Out)
#print(EXP.shape)
#print(EXP)

tmp = PEP['Protein'].str.split("|", expand = True)
PEP['proteinID'] = tmp[2]
PEP = PEP.iloc[:,[5,0,1,2,3,4]]

result = pd.merge(PEP, EXP, how="left", on="proteinID")
result.to_csv(outFile, index=False)
