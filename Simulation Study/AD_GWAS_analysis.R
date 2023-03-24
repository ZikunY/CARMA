###This is a demo to illustrate the application of CARMA to two loci, ADAMTS4 and CR1, on chromosome 1 from an Alzheimer's disease GWAS (https://ctg.cncr.nl/software/summary_statistics).

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

###### load the GWAS summary statistics (part of AD GWAS sumstats from Jansen et al., 2019)
sumstat.1 =  fread(file = "Sample_data/ADAMTS4_sumstats.txt.gz",
                   sep = "\t", header = T, check.names = F, data.table = F,
                   stringsAsFactors = F)
sumstat.2 =  fread(file = "Sample_data/CR1_sumstats.txt.gz",
                   sep = "\t", header = T, check.names = F, data.table = F,
                   stringsAsFactors = F)
 
###### load the functional annotations for the variants included in
###### GWAS summary statistics (assuming the same order of variants
###### as the variants in sumstat file)
annot.1 =  fread(file = "Sample_data/ADAMTS4_annotations.txt.gz",
                 sep = "\t", header = T, check.names = F, data.table = F,
                 stringsAsFactors = F)
annot.2 =  fread(file = "Sample_data/CR1_annotations.txt.gz",
                 sep = "\t", header = T, check.names = F, data.table = F,
                 stringsAsFactors = F)
 
###### load the pair-wise LD matrix (assuming the same order of variants
###### as in sumstat file)
ld.1 =  fread(file = "Sample_data/ADAMTS4_ld.txt.gz",
              sep = "\t", header = F, check.names = F, data.table = F,
              stringsAsFactors = F)
ld.2 =  fread(file = "Sample_data/CR1_ld.txt.gz",
              sep = "\t", header = F, check.names = F, data.table = F,
              stringsAsFactors = F)
              

###The format of input for CARMA is list. For each element of the list, e.g., the list of Z scores or the list of LD, it contains the corresponding input data for each locus.
z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat.1$Z #Z scores for ADAMTS4 region
z.list[[2]]<-sumstat.2$Z #Z scores for CR1 region
ld.list[[1]]<-as.matrix(ld.1) #LD for ADAMTS4 region
ld.list[[2]]<-as.matrix(ld.2) #LD for CR1 region
lambda.list[[1]]<-1 #setting eta=1 for CARMA for ADAMTS4 region
lambda.list[[2]]<-1 #setting eta=1 for CARMA for CR1 region

###### Without annotations
###### The outlier detection is turned on. 
##### Here we run CARMA main function without annotations.
CARMA.results_no_annot<-CARMA_fixed_sigma(z.list,ld.list,lambda.list =  lambda.list,
                                          outlier.switch=T)
 
###### With annotations
###### Exclude the variant information (such as positions and REF/ALT alleles) columns in the annotation file
annot.list<-list()
annot.list[[1]]<-as.matrix(cbind(1, annot.1 %>% select(-(uniqID.a1a2:SNP)))) #Functional annotations for ADAMTS4 region
annot.list[[2]]<-as.matrix(cbind(1, annot.2 %>% select(-(uniqID.a1a2:SNP)))) #Functional annotations for CR1 region


####Here we run CARMA main function with annotations in annot.list.
CARMA.results_annot<-CARMA_fixed_sigma(z.list,ld.list,w.list=annot.list,
                                       lambda.list =  lambda.list,
                                       input.alpha=0, outlier.switch=T)
 
###### Posterior inclusion probability (PIP) and credible sets (CS) for ADAMTS4 region
sumstat.1 = sumstat.1 %>% mutate(PIP = CARMA.results_annot[[1]]$PIPs, CS = 0)
sumstat.1$CS[CARMA.results_annot[[1]]$`Credible set`[[2]][[1]]] = 1
###### Posterior inclusion probability (PIP) and credible sets (CS) for CR1 region
sumstat.2 = sumstat.2 %>% mutate(PIP = CARMA.results_annot[[2]]$PIPs, CS = 0)
sumstat.2$CS[CARMA.results_annot[[2]]$`Credible set`[[2]][[1]]] = 1
 
###### write the GWAS summary statistics with PIPs and CS
fwrite(x = sumstat.1, file = "Sample_data/ADAMTS4_carma.txt.gz",
        sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
fwrite(x = sumstat.2, file = "Sample_data/CR1_carma.txt.gz",
        sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
