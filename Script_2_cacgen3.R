library(sva)
#load design matrix for calculating surrogate variables 
#this is the data sets we obtained using R Script 1
load("/restricted/projectnb/cvdmrx/cacscnt2/gen3/Sample_Design_Matrix_cacg3.RData")
#load original full DNA methylation matrix (443252 probes, which include data_1.csv to data_444.csv in the directory: /restricted/projectnb/fhs-methylation/cleaned/data/; data_445.csv and data_446.csv are not included due to the excluded probes for JHU or UMN subjects)
load("/restricted/projectnb/cvdmrx/forgen3/allg3may17.RData")

#obtain DNA methylation matrix with subjects who have non-missing methylation, phenotype and covariates
forallg3=forallg3[forallg3$framid %in% all_id, ]
forallg3=forallg3[with(forallg3, order(framid)), ]
forallg3=t(forallg3[, -1])
colnames(forallg3)=all_id
dim(forallg3)

#two-steps for calculating surrogate variables (set permutation times=10 and only include the 100,000 most variable probes)
n.sv=num.sv(forallg3, mod, method="be", vfilter=100000, B=10)
n.sv
svcacg3 = sva(forallg3, mod, mod0, n.sv=n.sv, method="irw", vfilter=100000, B=10)

#save calculated surrogate variables
#the data set can be found script directory: /restricted/projectnb/fhs-methylation/cleaned/scripts/shuo/
save(svcacg3, file="/restricted/projectnb/cvdmrx/cacscnt2/gen3/Sample_SV_gen3cac.RData")
system(paste("qsub -P ",project," -l mem_free=10G -l h_rt=720:00:00 -cwd)) 
