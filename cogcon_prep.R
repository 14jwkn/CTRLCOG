#Set code as working directory.

#Load packages.
library(tidyverse)

# Functions ---------------------------------------------------------------

#Residualizer function.
resid_maker <- function(cvar,cont_list,cat_list,cmat,resflag) {
  
  #If built-in R function without intercept (demeaning).
  if (resflag == 'demean') {
    
    #Do test with built-in R function.
    conform <- paste0(c(cont_list,cat_list),collapse='+')
    inform <- as.formula(paste0(cvar,'~',conform))
    valresid <- data.frame(resid(lm(inform,data=cmat))) 
    colnames(valresid) <- cvar
  
  #If raw value with intercept kept (no demeaning).
  } else if (resflag == 'raw') {
    
    #Isolate target variable.
    valtrue <- cmat %>%
      select(all_of(cvar))
    selmat <- valtrue
    
    #Isolate categorical variables, produce dummy variables, remove reference column.
    if (length(cat_list) != 0) {
      catmat <- cmat %>%
        select(all_of(cat_list)) %>%
        mutate_if(is.character,as.factor) 
      # fullcat <- as.data.frame(model.matrix(~.+0,data=catmat,
      #                                       contrasts.arg=lapply(catmat,contrasts,contrasts=FALSE)))
      selcat <- as.data.frame(model.matrix(~.,data=catmat))
      selcat[,'(Intercept)'] <- NULL
      selmat <- cbind(selmat,selcat)
    }
    
    #Isolate continuous variables.
    if (length(cont_list) != 0) {
      selcont <- cmat %>%
        select(all_of(cont_list))
      selmat <- cbind(selmat,selcont)
    }

    #Fit.
    resid_out <- colnames(selmat)[colnames(selmat)!=cvar]
    resid_form <- paste0(resid_out,collapse='+')
    inform <- as.formula(paste0(cvar,'~',resid_form))
    cfit <- lm(inform,data=selmat)
    
    #Use confound coefficients to build predicted value.
    valpred <- 0
    for (cconf in resid_out) {
      valpred <- valpred + cfit$coefficients[cconf]*selmat[,cconf]
    }
    # valpred <- valpred + cfit$coefficients['(Intercept)']
    
    #Produce true - predicted.
    valresid <- data.frame(valtrue-valpred)
  }

  #Return the true - predicted value without labels.
  return(valresid)
}

# Driver ------------------------------------------------------------------

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
  rename('NT'=CogTotalComp_Unadj,
         'NF'=CogFluidComp_Unadj,
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

#Add a P24 correct and no gPCA column.
infile <- paste0('../outputs/c_cognition/ourfa_subs/cogscores.csv')
hcpbehave <- read.csv(infile) %>%
  mutate(Subject=as.character(Subject)) 
cval <- 'PMAT24'
cont_list <- c('gPCA')
cat_list <- c()
cmat <- hcpbehave
gFngPCA2 <- data.frame(cmat$Subject,resid_maker(cval,cont_list,cat_list,cmat,'raw'))
colnames(gFngPCA) <- c('Subject','gFngPCA')
predall <- predall %>%
  left_join(gFngPCA)

#Add a PicVoc correct and no gPCA column.
infile <- paste0('../outputs/c_cognition/ourfa_subs/cogscores.csv')
hcpbehave <- read.csv(infile) %>%
  mutate(Subject=as.character(Subject)) 
cval <- 'PicVoc'
cont_list <- c('gPCA')
cat_list <- c()
cmat <- hcpbehave
gCngPCA <- data.frame(cmat$Subject,resid_maker(cval,cont_list,cat_list,cmat,'raw'))
colnames(gCngPCA) <- c('Subject','gCngPCA')
predall <- predall %>%
  left_join(gCngPCA2)

#Add a P24 correct and no gCFA column.
infile <- paste0('../outputs/c_cognition/ourfa_subs/cogscores.csv')
hcpbehave <- read.csv(infile) %>%
  mutate(Subject=as.character(Subject))  %>%
  inner_join(du_CFA)
cval <- 'PMAT24'
cont_list <- c('gCFA')
cat_list <- c()
cmat <- hcpbehave
gFngCFA <- data.frame(cmat$Subject,resid_maker(cval,cont_list,cat_list,cmat,'raw'))
colnames(gFngCFA) <- c('Subject','gFngCFA')
predall <- predall %>%
  left_join(gFngCFA)

#Add a PicVoc and no gCFA column.
infile <- paste0('../outputs/c_cognition/ourfa_subs/cogscores.csv')
hcpbehave <- read.csv(infile) %>%
  mutate(Subject=as.character(Subject))  %>%
  inner_join(du_CFA)
cval <- 'PicVoc'
cont_list <- c('gCFA')
cat_list <- c()
cmat <- hcpbehave
gCngCFA <- data.frame(cmat$Subject,resid_maker(cval,cont_list,cat_list,cmat,'raw'))
colnames(gCngCFA) <- c('Subject','gCngCFA')
predall <- predall %>%
  left_join(gCngCFA)

#Add a P24 correct and no gEFA column.
infile <- paste0('../outputs/c_cognition/ourfa_subs/cogscores.csv')
hcpbehave <- read.csv(infile) %>%
  mutate(Subject=as.character(Subject))  %>%
  inner_join(du_EFA)
cval <- 'PMAT24'
cont_list <- c('gEFA')
cat_list <- c()
cmat <- hcpbehave
gFngEFA <- data.frame(cmat$Subject,resid_maker(cval,cont_list,cat_list,cmat,'raw'))
colnames(gFngEFA) <- c('Subject','gFngEFA')
predall <- predall %>%
  left_join(gFngEFA)

#Add a PicVoc and no gEFA column.
infile <- paste0('../outputs/c_cognition/ourfa_subs/cogscores.csv')
hcpbehave <- read.csv(infile) %>%
  mutate(Subject=as.character(Subject))  %>%
  inner_join(du_EFA)
cval <- 'PicVoc'
cont_list <- c('gEFA')
cat_list <- c()
cmat <- hcpbehave
gCngEFA <- data.frame(cmat$Subject,resid_maker(cval,cont_list,cat_list,cmat,'raw'))
colnames(gCngEFA) <- c('Subject','gCngEFA')
predall <- predall %>%
  left_join(gCngEFA)

#Generate full NT, NF, and NC.
infile <- '../inputs/data/hcp/HCP1200_Behavioral.csv'
hcpbehave <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject)) %>%
  select(Subject,
        CogTotalComp_Unadj,
        CogFluidComp_Unadj,
        CogCrystalComp_Unadj) %>%
  rename('NT'=CogTotalComp_Unadj,
         'NF'=CogFluidComp_Unadj,
         'NC'=CogCrystalComp_Unadj)
hcpbehave <- data.frame(hcpbehave)
rownames(hcpbehave) <- hcpbehave$Subject
hcpbehave <- hcpbehave[behave_subjects,]
rownames(hcpbehave) <- NULL

#Add a NF and no NT column.
cval <- 'NF'
cont_list <- c('NT')
cat_list <- c()
cmat <- hcpbehave
outmat <- data.frame(cmat$Subject,resid_maker(cval,cont_list,cat_list,cmat,'raw'))
colnames(outmat) <- c('Subject','NFnNT')
predall <- predall %>%
  left_join(outmat)

#Add a NC and no NT column.
cval <- 'NC'
cont_list <- c('NT')
cat_list <- c()
cmat <- hcpbehave
outmat <- data.frame(cmat$Subject,resid_maker(cval,cont_list,cat_list,cmat,'raw'))
colnames(outmat) <- c('Subject','NCnNT')
predall <- predall %>%
  left_join(outmat)

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
