# Gather cognitive scores and confounds into one file.
# Output:
# pred_all.csv Cognitive scores and confounds.
# dum_all.csv Cognitive scores and confounds with dummy variables generated.
# dum_sel.csv Cognitive scores and confounds with dummy variables generated and reference column removed.

#Load packages.
library(tidyverse)

#Set outpath.
subgroup <- 'full'
outpath <- paste0('../outputs/c_cognition/',subgroup,'/')
dir.create(outpath,recursive=T)

#Read in subjects.
subjects <- scan('r_full_submain.txt',character(),quote='') 

#Read in behavioral subjects.
behave_subjects <- scan('ourfa_subs.txt',character(),quote='') 

#Read in HCP cognitive variables and demographics.
inpath <- '../inputs/data/hcp/'
infile <- paste0(inpath,'HCP1200_Data_Dictionary.csv')
hcpopen <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject))
infile <- paste0(inpath,'RESTRICTED_HCP1200_Data_Dictionary.csv')
hcpclose <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject))
infile <- paste0(inpath,'unrestricted_hcp_freesurfer.csv')
hcpfs <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject))
hcpfull <- hcpopen %>%
  full_join(hcpclose) %>%
  full_join(hcpfs)
predall <- hcpfull %>%
  select(Subject,
         CogTotalComp_Unadj,
         CogFluidComp_Unadj,
         CogCrystalComp_Unadj,
         PMAT24_A_CR, 
         PicVocab_Unadj,
         ProcSpeed_Unadj,
         CardSort_Unadj,
         Flanker_Unadj,
         ListSort_Unadj,
         PicSeq_Unadj,
         ReadEng_Unadj,
         IWRD_TOT,
         VSPLOT_TC,
         Gender,
         Age_in_Yrs,
         Race) %>%
  rename('NF'=CogFluidComp_Unadj,
         'NC'=CogCrystalComp_Unadj,
         'P24_CR'=PMAT24_A_CR, 
         'PV'=PicVocab_Unadj,
         'PS'=ProcSpeed_Unadj,
         'CardSort'=CardSort_Unadj,
         'Flanker'=Flanker_Unadj,
         'ListSort'=ListSort_Unadj,
         'PicSeq'=PicSeq_Unadj,
         'ReadEng'=ReadEng_Unadj,
         'Age'=Age_in_Yrs) %>%
  mutate_if(is.character,as.factor) %>%
  mutate(Gender=recode_factor(Gender,
                              'F'='F',
                              'M'='M'),
         Race=recode_factor(Race,
                            'Unknown or Not Reported'='Unknown',
                            'More than one'='Multiple',
                            'Asian/Nat. Hawaiian/Othr Pacific Is.'='Asian',
                            'Black or African Am.'='Black',
                            'White'='White',
                            'Am. Indian/Alaskan Nat.'='Indian'))

#Add a g from PCA column.
infile <- paste0('../outputs/c_cognition/ourfa_subs/cogscores.csv')
gPCA <- read.csv(infile) %>%
  mutate(Subject=as.character(Subject)) %>%
  select(c('Subject','gPCA'))    
colnames(gPCA) <- c('Subject','gPCA')
predall <- predall %>%
  left_join(gPCA)

#Add NF and NC.
infile <- '../inputs/data/hcp/HCP1200_Behavioral.csv'
hcpbehave <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject)) %>%
  select(Subject,
        CogFluidComp_Unadj,
        CogCrystalComp_Unadj) %>%
  rename('NF'=CogFluidComp_Unadj,
         'NC'=CogCrystalComp_Unadj)
hcpbehave <- data.frame(hcpbehave)
rownames(hcpbehave) <- hcpbehave$Subject
hcpbehave <- hcpbehave[behave_subjects,]
rownames(hcpbehave) <- NULL

#Save the matrix.
write.csv(predall,paste0(outpath,'pred_all.csv'),row.names=F)

#Isolate categorical variables.
dummat <- predall %>%
  select(Gender,Race,Ethnicity,Recon) %>%
  mutate_if(is.character,as.factor) 

#Produce full dummy variables.
fulldum <- as.data.frame(model.matrix(~.+0,data=dummat,
                                      contrasts.arg=lapply(dummat,contrasts,contrasts=FALSE)))
rownames(fulldum) <- predall[,'Subject']

#Save the matrix.
write.csv(fulldum,paste0(outpath,'dum_all.csv'),row.names=T)

#Produce reference column subtracted dummy variables which all of them have.
seldum <- fulldum %>%
  select(-GenderF,-RaceWhite,-EthnicityNonhispanic,-Reconr227)

#Save the matrix.
write.csv(seldum,paste0(outpath,'dum_sel.csv'),row.names=T)
