library(PRROC )
library(pROC )
library(pryr )
library(MASS )
library(Matrix )
library(CARMA )
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(RcppProgress)
library(glmnet)
rm(list=ls())


Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
.libPaths("/ifs/scratch/msph/eigen/zy2412/R_libs/")
args = commandArgs(TRUE)

index<-1

ref.table<-read.csv('AD_data_loci.csv')

z.list<-list()
ld.list<-list()
full.w.list<-list()
lambda.list<-c()
label.list<-list()
input.prior.list<-list()

for(g in 1:length(which(ref.table$chr==index))){
      print(ref.table$locus.name[(which(ref.table$chr==index))[g]])
      pro.data<-read.table(paste0(ref.table$locus.name[(which(ref.table$chr==index))[g]],'.processed'),header = T)
      pro.ld<-read.table(paste0(ref.table$locus.name[(which(ref.table$chr==index))[g]],'.ld'))
      anno<-read.table(paste0(ref.table$locus.name[(which(ref.table$chr==index))[g]],'_polyfun.annot'),header = T)
      full.w.list[[g]]<-as.matrix(data.frame(1,anno[,-(1:6)]))
      
      ld.matrix<-as.matrix(pro.ld)
      ld.list[[g]]<-ld.matrix
      
      z<-pro.data[,7,drop=F]
      z.list[[g]]<-as.matrix(z)
      
      colnames(ld.list[[g]])<-rownames(ld.list[[g]])<-1:nrow(ld.list[[g]])
      lambda.list[[g]]<-1
      label.list[[g]]<-ref.table$locus.name[(which(ref.table$chr==index))[g]]
      full.w.list[[g]]<-as.matrix(data.frame(1,anno[,-(1:6)]))
    }




n=max(pro.data$Nsum);phi=0.01
tau.value=n*phi/qchisq(0.95,1)
tau.value=1/tau.value

    

a<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,w.list=full.w.list,outlier.switch=T,
effect.size.prior='Ind.Normal',outlier.BF.index = 1/3.2,tau=tau.value)

saveRDS(a,paste0('/chr',index,'.RDS') )


for(s in 1:length(a)){
    save.list<-list()
    save.list[[1]]<-a[[s]]
saveRDS(save.list,paste0(ref.table$locus.name[(which(ref.table$chr==index))[s]],'.RDS') )
}
