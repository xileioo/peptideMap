#!/usr/bin/python
# 2022-10-26 by lixue
# demo: python Step4_finalMerge.py Step2Out step3Out outFile

#import sys
import sys
from cmath import nan
import re
import pandas as pd
import numpy as np


# define input and parameter
Step2Out = sys.argv[1]
Step3Out = sys.argv[2]
outFile = sys.argv[3]
#Step2Out = 'D:/peptideMap/KVLEHXXXV_atleastMatch4/Step2_matchedRedundantPeptides.netMHCpan.out'
#Step3Out = 'D:/peptideMap/KVLEHXXXV_atleastMatch4/Step3_out.csv'
#outFile = "D:/peptideMap/KVLEHXXXV_atleastMatch4/Step4_out.csv"

#----------------------------------------------------------
# Wrangling Step2 Out
#----------------------------------------------------------
df_list = []
MHC_res = open(Step2Out,'r')
for line in MHC_res:
    if re.match('.+<=', line):
        myarr = line.split()
        df_list.append(myarr)
df = pd.DataFrame(df_list)

df = df[[1,2,11,12,14]]
df.columns = ['HLA','Peptide','Score_EL','Rank_EL','BindLevel']
df = df.drop_duplicates()
#print(df)

#----------------------------------------------------------
# Wrangling Step3 Out
#----------------------------------------------------------
PEP = pd.read_csv(Step3Out)
#print(PEP['Level'].value_counts())
# Medium    28991
# Low       17199
# High      12533

# KEEP main tissue
#print(PEP['Tissue'].value_counts())
keep_tissues = ['colon','breast','lung','kidney','pancreas','lymph node','spleen',
'liver','rectum','duodenum','small intestine','retina','adrenal gland','thymus',
'adipose tissue','smooth muscle','heart muscle','skin','thyroid gland','skeletal muscle','bone marrow']
PEP = PEP[PEP['Tissue'].isin(keep_tissues)]

# IHC Level convert to digital
PEP.loc[PEP['Level'] == 'Low','Level'] = 1
PEP.loc[PEP['Level'] == 'Medium','Level'] = 2
PEP.loc[PEP['Level'] == 'High','Level'] = 3
PEP['Level'] = pd.to_numeric(PEP['Level'])
PEP = PEP.rename({'Peptide_sequence': 'Peptide'}, axis='columns')
# print(PEP.dtypes)
# print(PEP['Level'].value_counts())


#----------------------------------------------------------
# Final Wrangling Step3 Out expression data
#----------------------------------------------------------
# Level = max of tissues
# Level = sum of tissues
# nTPM = mean of tissues
group_by = ['proteinID','GeneName','Peptide','matchCount']
myres = PEP[[*group_by,'Tissue','Level','nTPM']].groupby([*group_by,'Tissue'],as_index=False).agg({'Level':'max','nTPM':'max'})

myres = myres.groupby(group_by,as_index=False).agg({'Level':['max','sum'],'nTPM':'mean'})
myres.columns = ['proteinID','GeneName','Peptide','matchCount','IHC_Level_max','IHC_Level_sum','nTPM_mean']
#print(myres)

# format Concatenate Tissue strings
mytissue = PEP.loc[np.isnan(PEP['Level']) == 0,[*group_by,'Tissue']]
mytissue = mytissue.drop_duplicates()
mytissue['Tissue'] = mytissue.groupby(group_by)['Tissue'].transform(lambda x: ','.join(x))
mytissue = mytissue.drop_duplicates()
#print(mytissue)

# merge Tissue, Level and nTPM
myres = pd.merge(mytissue, myres, how="right", on=group_by)
#print(myres)

#----------------------------------------------------------
# Merge MHCpan (df) and Expression (myres)
#----------------------------------------------------------
#myres_filter = myres[myres['Peptide'].isin(df['Peptide'])]
final_res = pd.merge(myres, df, how="left", on='Peptide')
final_res.to_csv(outFile)

## check
#a = set(myres['Peptide']).difference(set(df['Peptide']))
#a = set(set(df['Peptide'])).difference(myres_filter['Peptide'])
#outTest1 = "D:/peptideMap/KVLEHXXXV_atleastMatch4/Step4_test1.csv"
#outTest2 = "D:/peptideMap/KVLEHXXXV_atleastMatch4/Step4_test2.csv"
#myres.to_csv(outTest1)
#df.to_csv(outTest2)

