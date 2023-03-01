#load one part of the methylation data
dat=read.csv("/restricted/projectnb/cvdmrx/cacscnt2/gen2/data_7.csv", header=T, as.is=T, na.strings="")

#obtain the ID of subject, who has methylation data
g_ID=dat[, 1]
# how to make sure framid is the first column
#load phenotype data set (the phenotype data is user-specified)
#the colunm name for subject ID in this sample phenotype data is "ID"
phe=read.csv("/restricted/projectnb/cvdmrx/cacscnt2/g2cacscan2a2.csv", header=T, as.is=T, na.strings="")
	
phe$femalesex=as.factor(phe$femalesex)
phe$lipidrxscan2=as.factor(phe$lipidrxscan2)	
phe$hrxscan2=as.factor(phe$hrxscan2)	
phe$diabscan2=as.factor(phe$diabscan2)
phe$prgbr=as.factor(phe$prgbr)
phe$currsmk=as.factor(phe$currsmk)
phe$frmsmk=as.factor(phe$frmsmk)
phe$excessdrk=as.factor(phe$excessdrk)

#select only phenotype and covariates to adjust, drop other variables
## 
##select only phenotype and covariates to adjust, drop other variables
## phe_AN=phe[,c("framid","age2", "sex", "calc_ldl2", "lprxany","site")]
## framid,IDTYPE,ID,dmany,hrxany,"lprxany,CD8T,CD4T,NK,Bcell,Mono,Gran,medcvd,sex,age,femalesex,diab,hrx,sbp,dbp,bmi,tc,hdl,calc_ldl,currsmk,shareid,offcoh
phe_AN=phe[,c("ID","lgcac2ad1","agescan2","femalesex","lipidrxscan2","bmiscan2","sbpscan2","hrxscan2","totalscan2","hdlscan2","diabscan2","logcac1","prgbr","yrdif2ct","chgpyr","currsmk","frmsmk","excessdrk","CD8T", "CD4T","NK", "Bcell","Mono","Gran")]
#obtain data matrix with non-missing phenotype and covariates
ind=which(rowSums(is.na(phe_AN))==0)
phe_AN=phe_AN[ind, ]

#obtain data matrix with non-missing methylation, phenotype and covariates
all_id=intersect(g_ID, phe_AN$ID)
all_id=sort(all_id)
phe_AN=phe_AN[phe_AN$ID %in% all_id,]
phe_AN=phe_AN[with(phe_AN, order(ID)), ]

#obtain design matrix for calculating surrogate variables 
mod=model.matrix(~ CD4T + CD8T+ NK +Bcell + Mono + Gran + lgcac2ad1+agescan2+femalesex+lipidrxscan2+bmiscan2+sbpscan2+hrxscan2+totalscan2+hdlscan2+diabscan2+logcac1 +prgbr + yrdif2ct +chgpyr + currsmk + frmsmk+excessdrk, data=phe_AN)
mod0=model.matrix(~ CD4T + CD8T+ NK +Bcell + Mono + Gran +         agescan2+femalesex+lipidrxscan2+bmiscan2+sbpscan2+hrxscan2+totalscan2+hdlscan2+diabscan2+logcac1 +prgbr+ yrdif2ct +chgpyr+currsmk + frmsmk +excessdrk , data=phe_AN)

#save data 
#the data set can be found in script directory: /restricted/projectnb/fhs-methylation/cleaned/scripts/shuo/
save(mod, mod0, all_id, file="/restricted/projectnb/cvdmrx/cacscnt2/gen2/Sample_Design_Matrix_cacg2.RData")

