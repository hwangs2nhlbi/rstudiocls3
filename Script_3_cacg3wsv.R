library(coxme)
library(kinship2)


#load surrogate variables
#this is the data set we obtained using R Script 2
 load("/restricted/projectnb/cvdmrx/cacscnt2/gen3/Sample_SV_gen3cac.RData")
 SV=svcacg3[[1]]
#load phenotype data set
# phe_AN=phe[,c("ID","age8","femalesex", "bmi8","LIP8","CD4T","CD8T","NK","Bcell","Mono","Gran","sbp8","dbp8", "HRx8","curr_diab8","tc8","hdl8","calc_ldl8","dmRX8","currsmk8","statinuse")]
#obtain data matrix with non-missing phenotype and covariates

phe=read.csv("/restricted/projectnb/cvdmrx/cacscnt2/g3cacscan2a2.csv", header=T, as.is=T, na.strings="")

phe=phe[,c("framid","lgcac2ad1","agescan2","femalesex","lipidrxscan2","bmiscan2","sbpscan2","hrxscan2","totalscan2","hdlscan2","diabscan2","logcac1","prgbr","yrdif2ct","chgpyr","currsmk","frmsmk","excessdrk","CD8T", "CD4T","NK", "Bcell","Mono","Gran")]

#load data set that contains the ID number of subjects who have non-missing methylation, phenotype and covariates
#this is the data sets we obtained using R script 1
load("/restricted/projectnb/cvdmrx/cacscnt2/gen3/Sample_Design_Matrix_cacg3.RData")

phe=phe[phe$framid%in% all_id, ]
phe=phe[order(phe$framid),]

## Get SVs that are associated with phenotype (CALC_LDL) at the level of 0.1
 p_val=c()
 for (i in 1:ncol(SV))
 {
 p_val[i]=cor.test(phe$lgcac2ad1, SV[,i])$p.value
 }

 SV=SV[, which(p_val<0.1)]
 SV=data.frame(SV) 
 
#load technical variables (For obtaining technical variables data, please see: /restricted/projectnb/fhs-methylation/cleaned/scripts/shuo/Script_Create_Technical_Variables.R)
tech=read.csv("/restricted/projectnb/fhs-methylation/cleaned/scripts/shuo/sample_tech.csv", header=T, as.is=T, na.strings="")

#load pedigree matrix
pedigree=read.csv("/restricted/projectnb/fhs_data/go550k/SubjectInfo/fam_unrel_comb_08122014.csv", header=T)
kmat= makekinship(famid=pedigree$famid, id=pedigree$id, mother.id=pedigree$mo, father.id=pedigree$fa)

dat<-phe
# dat=merge(phe, tech[, c(1,6)], by=c("frammid"), all.x=T)
# dat$femalesex=as.factor(dat$femalesex)
# dat$currsmk=as.factor(dat$currsmk)	
# dat$frmsmk=as.factor(dat$frmsmk)	
# dat$ppoor2=as.factor(dat$ppoor2)
# dat$ppoor1=as.factor(dat$ppoor1)
# dat$excessdrk=as.factor(dat$excessdrk)
# dat$schlg1=as.factor(dat$schlg1)

#dat$femalesex=as.factor(dat$femalesex)
#dat$lipidrxscan2=as.factor(dat$lipidrxscan2)	
#dat$hrxscan2=as.factor(dat$hrxscan2)	
#dat$diabscan2=as.factor(dat$diabscan2)
#dat$prgbr=as.factor(dat$prgbr)
#dat$currsmk=as.factor(dat$currsmk)
#dat$frmsmk=as.factor(dat$frmsmk)
#dat$excessdrk=as.factor(dat$excessdrk)
 dat=cbind(dat, SV)

# a is a formula, that is X1+бн+X50 (depend on how many SV to use)
 a="X1"
 for (i in 2:ncol(SV))
 {
  a=paste(a, " +X",i, sep="")
 }

#remove data set and variables that will not be used for lmekin model
# rm(i, pedigree, phe,  p_val, svlprxany, mod, mod0)

#save the data frames used for lme model fitting
#the data set can be found in script directory: /restricted/projectnb/fhs-methylation/cleaned/scripts/shuo/
save.image(file="/restricted/projectnb/cvdmrx/cacscnt2/gen3/Sample_Data_cacg3wsv.RData")
