# read in which part to run
# args <- commandArgs(TRUE) 
# part <- args[1]
rm(list=ls())
#the script takes about 100 mins to finish for each single data set with 1000 probes
setwd("/restricted/projectnb/cvdmrx/cacscnt2/gen2/")
library(kinship2) 
library(coxme)
library(quadprog)
library(bdsmatrix)
library(nlme)
# define working directories
datadir<-c("/restricted/projectnb/cvdmrx/cacscnt2/gen2/")
outputdir<-c("/restricted/projectnb/cvdmrx/cacscnt2/gen2/cacasuse2/")
 num = as.numeric(Sys.getenv("SGE_TASK_ID"))

#############################################################################
# 1. Load in the data and remove & record probe extreme outliers (+/- 5SD) 
# read in which part to run
args <- commandArgs(TRUE) 
part <- args[1]

#the script takes about 100 mins to finish for each single data set with 1000 probes

#############################################################################
# 1. Load in the data and remove & record probe extreme outliers (+/- 5SD)  #
#############################################################################


#load data sets (phenotype, SVs and pedigree matrix) this is the data sets we obtained from R Script 3
load(paste(datadir,"Sample_Data_cacg2wsv.RData",sep=""))
##  add code  ##

#load one part of the DNA methylation data (1000 probes)
res=read.csv(paste("/restricted/projectnb/fhs-methylation/cleaned/data/data_",num,".csv", sep=""), header=T, as.is=T, na.strings="")
#extract probe names
nam=colnames(res)[-1]

dat=merge(dat, res, by=c("ID"), all.x=T)

 #
#############################################################################

load(paste(datadir,"Sample_Data_cacg2wsv.RData",sep=""))

res=read.csv(paste("/restricted/projectnb/fhs-methylation/cleaned/data/data_",num,".csv", sep=""), header=T, as.is=T, na.strings="")
#extract probe names'
nam=colnames(res)[-1]

dat=merge(dat, res, by=c("ID"), all.x=T)
out_cacg2=data.frame(matrix(NA, (ncol(res)-1), 5))
out_cacg2[,1]=nam
# "age2","sex", "BMI2","lprxany","CALC_LDL2","CD4T","CD8T","NK","Bcell","Mono","Gran","TC2","HDL2","TRIG2","SBP2","DBP2", "hrxany","dmany","CIURRSMK2","curr_diab2","f_BG2"
# +lipidrxscan2+ bmiscan2+ sbpscan2+ hrxscan2+ totalscan2+ hdlscan2+ diabscan2

#looping probes for model fitting
for(i in 1:length(nam))
{
  
f1=formula(paste(nam[i], "~ lgcac2ad1 + agescan2 + femalesex   +  CD4T + CD8T + NK + Bcell + Mono + Gran +  ",a,"  +   (1|ID)", sep="" ))
L1=try(lmekin(f1, data=dat, varlist=list(kmat), na.action=na.omit))
if (class(L1)!="try-error")
   {
    N1=L1$n
    Beta1=fixef(L1)[2]
    SE1=sqrt(L1$var[2,2])
    p1=pnorm(abs(Beta1/SE1), lower.tail = F)*2
    out_cacg2[i, 2:5]=c(Beta1, SE1, p1, N1)
  }
  else {}
  
  if(i %%100 ==0)
  {
    cat("Now, it is the ", i,"-th iteration\n", sep=" ")
  }
  else{}
}

#output results
write.table(out_cacg2, file=paste('/restricted/projectnb/cvdmrx/cacscnt2/gen2/cacasuse2/Final_cacg2as_',num,'.csv', sep=""), col.names=F, row.names=F, sep=",", na="", quote=F)
}
