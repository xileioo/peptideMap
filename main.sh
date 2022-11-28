#!/bin/sh

# conda create -n BIO python=3.8
# conda activate BIO
# conda install -c conda-forge biopython
# install netMHCpan-4.1
# cd ~/pipetideMap 

#---------------------------------------------#
#  Required parameters                        #
#---------------------------------------------#
query="KVLEHXXXV"
atleast_matchAA=4
MHCallele="HLA-A02:01"

#---------------------------------------------#
#  prepare files                              #
#---------------------------------------------#
outDir="${query}_atleastMatch${atleast_matchAA}"
mkdir -p $outDir
output1="$outDir/Step1_matchedRedundantPeptides.csv"
output2="$outDir/Step2_matchedRedundantPeptides.netMHCpan.out"
output3="$outDir/Step3_matchedRedundantPeptides.exp.csv"
finalOut="$outDir/Step4_final_result.csv"

processing_temp="$outDir/Step2_netMHCpan_input.tmp"

exp_HPA="exp_HPA/My_Merged_normal_and_rna_tissue_consensus.tsv"

#---------------------------------------------#
#   Step1: sliding window                     #
#---------------------------------------------#
python Step1_redundantPeptide.py $query $atleast_matchAA $output1


#---------------------------------------------#
#   Step2: netMHCpan                          #
#---------------------------------------------#
awk -F "," '{print $4}' $output1 > $processing_temp
/home/innovent/public/netMHCpan-4.1/netMHCpan -p $processing_temp -a $MHCallele > $output2 
rm $processing_temp


#---------------------------------------------#
#   Step3: IHC and TPM                        #
#---------------------------------------------#
python Step3_Protein_expression.py $output1 $exp_HPA $output3


#---------------------------------------------#
#   Step4: Final merge ICH/TPM and netMHCpan--#
#---------------------------------------------#
python Step4_finalMerge.py $output2 $output3 $finalOut

