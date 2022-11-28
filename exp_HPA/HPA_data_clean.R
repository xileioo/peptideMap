
#exp_ihc[exp_ihc$Gene == "ENSG00000006611" & exp_ihc$Tissue == "kidney",]
library(tidyverse)
mainPath <- 'D:/work/task_TS/project_2022/20220527_TCR_peptide/peptideMap/exp_HPA/'
setwd(mainPath)

mydata <- read.csv("../ReferenceProteomes_EMBL_EBI/proteinID_geneID_table.csv",header = T)
exp_ihc <- read_tsv("normal_tissue.tsv") # [1] 1193218       6
exp_con <- read_tsv("rna_tissue_consensus.tsv") # [1] 1084208       4
colnames(exp_ihc)[2] = "GeneName"
colnames(exp_ihc)[4] = "CellType"
colnames(exp_con)[2] = "GeneName"
head(exp_ihc)
#              Gene  GeneName         Tissue            CellType        Level Reliability
# 1 ENSG00000000003    TSPAN6 adipose tissue          adipocytes Not detected    Approved
# 2 ENSG00000000003    TSPAN6  adrenal gland     glandular cells Not detected    Approved
# 3 ENSG00000000003    TSPAN6       appendix     glandular cells       Medium    Approved
# 4 ENSG00000000003    TSPAN6       appendix     lymphoid tissue Not detected    Approved
# 5 ENSG00000000003    TSPAN6    bone marrow hematopoietic cells Not detected    Approved
# 6 ENSG00000000003    TSPAN6         breast          adipocytes Not detected    Approved
head(exp_con)
#              Gene  GeneName         Tissue nTPM
# 1 ENSG00000000003    TSPAN6 adipose tissue 29.6
# 2 ENSG00000000003    TSPAN6  adrenal gland 18.2
# 3 ENSG00000000003    TSPAN6       amygdala  7.2
# 4 ENSG00000000003    TSPAN6       appendix  4.9
# 5 ENSG00000000003    TSPAN6  basal ganglia  7.5
# 6 ENSG00000000003    TSPAN6    bone marrow  0.6

table(exp_ihc$Level)

# Ascending         Descending               High                Low 
#       172                 73             139347             183718 
# Medium                N/A       Not detected Not representative 
# 302624               1867             565394                 23 

exp_ihc <- exp_ihc[exp_ihc$Level %in% c("High","Low","Medium"),]

ihc_tissue <- table(exp_ihc$Tissue)
con_tissue <- table(exp_con$Tissue)
keep_tissue <- intersect(names(ihc_tissue), names(con_tissue))
# [1]  "adipose tissue"    "adrenal gland"     "appendix"          "bone marrow"       "breast"           
# [6]  "cerebellum"        "cerebral cortex"   "cervix"            "choroid plexus"    "colon"            
# [11] "duodenum"          "epididymis"        "esophagus"         "fallopian tube"    "gallbladder"      
# [16] "heart muscle"      "hypothalamus"      "kidney"            "liver"             "lung"             
# [21] "lymph node"        "ovary"             "pancreas"          "parathyroid gland" "pituitary gland"  
# [26] "placenta"          "prostate"          "rectum"            "retina"            "salivary gland"   
# [31] "seminal vesicle"   "skeletal muscle"   "skin"              "small intestine"   "smooth muscle"    
# [36] "spleen"            "testis"            "thymus"            "thyroid gland"     "tonsil"           
# [41] "urinary bladder"   "vagina"  

setdiff(names(ihc_tissue), names(con_tissue))
# [1]  "bronchus"         "cartilage"        "caudate"          "dorsal raphe"     "endometrium 1"   
# [6]  "endometrium 2"    "eye"              "hair"             "hippocampus"      "lactating breast"
# [11] "N/A"              "nasopharynx"      "oral mucosa"      "skin 1"           "skin 2"          
# [16] "soft tissue 1"    "soft tissue 2"    "stomach 1"        "stomach 2"   

setdiff(names(con_tissue), names(ihc_tissue))
# [1] "amygdala"              "basal ganglia"         "endometrium"           "hippocampal formation"
# [5] "medulla oblongata"     "midbrain"              "pons"                  "spinal cord"          
# [9] "stomach"               "thalamus"              "tongue"                "white matter"  

exp_ihc <- exp_ihc[exp_ihc$Tissue %in% keep_tissue,]
exp_con <- exp_con[exp_con$Tissue %in% keep_tissue,]
dim(exp_ihc)
#[1] 458537      6
dim(exp_con)
#[1] 843128      4

exp_total <- merge(exp_ihc, exp_con, 
                   by = c("Gene","GeneName","Tissue"),
                   all = T)
dim(exp_total)
# [1] 989461      7


### merge protein ID and gene ID data
keep_protein <- intersect(mydata$geneID, exp_total$GeneName)
length(keep_protein)
#[1] 18686

exp_total2 <- merge(exp_total, mydata,
                    by.x = "GeneName", by.y = "geneID",
                    all.x = T)
exp_total2 <- exp_total2[,c(8,1:7)]

#write_tsv(exp_total,"My_Merged_normal_and_rna_tissue_consensus.tsv")
write_tsv(exp_total2,"My_Merged_normal_and_rna_tissue_consensus.tsv")

mean(exp_total2$nTPM,na.rm = T)
#[1] 66.84245



