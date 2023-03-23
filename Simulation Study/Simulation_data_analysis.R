###This is a demo to illustrate the application of CARMA to the simulated data, which can be generated through the code "Genotype_to_Z-score_generation.R"

###Create a directory
#mkdir CARMA
#cd CARMA
##### Download and save the demo data in folder `CARMA'
#wget -O Sample_data.tar.gz https://osf.io/5gqz8/download
##### or download the file from https://osf.io/4t2bz/
#tar -zxf Sample_data.tar.gz
#unzip gz file
#gzip -d Sample_data.tar.gz
##Load libraries

library(Matrix )
library(MASS )
library(CARMA )
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(glmnet)
library(dplyr)
library(data.table)
rm(list=ls())



Sys.setenv("PKG_CXXFLAGS"="-std=c++11")#compile functions that use C++11 in R


setwd('CARMA')#Set the working directory

######The example data are in the folder CARMA, which could be downloaded through the code at the beginning of this document.
###### load the GWAS summary statistics
sumstat<- fread(file = "Sample_data/sumstats_chr1_200937832_201937832.txt.gz",
                           sep = "\t", header = T, check.names = F, data.table = F,
                           stringsAsFactors = F)
###### load the pair-wise LD matrix (assuming the same order of variants
###### as the variants in sumstat file)
ld =  fread(file = "Sample_data/sumstats_chr1_200937832_201937832_ld.txt.gz",
                       sep = "\t", header = F, check.names = F, data.table = F,
                       stringsAsFactors = F)

##### set up the input files for CARMA
z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1 ####setting eta=1

##### run CARMA
##### We are using in-sample LD here, therefore, the outlier detection is off
CARMA.results<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=F)
###### Posterior inclusion probability (PIP) and credible sets (CS)
sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}
###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = "Sample_data/sumstats_chr1_200937832_201937832_carma.txt.gz",
       sep = "\t", quote = F, na = "NA", row.names = F, col.names = T,
       compress = "gzip")
