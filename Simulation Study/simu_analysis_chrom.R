
library(PRROC )
library(pROC )
library(Matrix )
library(pryr )
library(MASS )
library(CARMA )
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(RcppProgress)
library(glmnet)
rm(list=ls())



Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

args = commandArgs(TRUE)

annot<-(args[4])
cor.index<-(args[3])
lambda.power<-as.numeric(args[2])
chr<-as.numeric(args[1])
print(args[1])
print(args)
setwd('Set your working direction')
   ref.table<-read.csv('Simulation loci.csv')
labels='large'
i=1

chor.index<-which(ref.table$chr==paste0('chr',chr))
z.list<-list()
ld.list<-list()  
full.w.list<-list()
lambda.list<-c()
label.list<-list()
g=1
for(s in chor.index){
pro.data<-as.matrix(read.table(paste0('data/',labels[i],'/', ref.table$chr[s],'_',
                                as.character(ref.table$region_start [s]),'_',ref.table$region_end[s]),header=T))

                                
pro.ld<-as.matrix(read.table(paste0(ref.table$chr[s],'_',as.character(ref.table$region_start [s]),'_',ref.table$region_end[s],'.ld'),header = F))



anno<-read.table(paste0('data/',labels[i],'/',ref.table$chr[s],'_',
         as.character(ref.table$region_start [s]),'_',ref.table$region_end[s],'.deepsea'),header=T)

 
      print(ref.table$gene.names[s])
       print(nrow(pro.data))
       print(nrow(pro.ld))

      z<-pro.data[,1,drop=F]
      z.list[[g]]<-as.matrix(z)
      ld.matrix<-as.matrix(pro.ld)
      ld.list[[g]]<-ld.matrix
      colnames(ld.list[[g]])<-rownames(ld.list[[g]])<-1:nrow(ld.list[[g]])

       lambda.list[[g]]<-lambda.power
      label.list[[g]]<-paste0(ref.table$chr[s],'_',as.character(ref.table$region_start [s]),'_',ref.table$region_end[s],'_',s)
        full.w.list[[g]]<-as.matrix(anno)
   
      g=g+1
}
print(annot)
if(annot=='no'){
  full.w.list<-NULL
}
  




a<-CARMA_fixed_sigma(z.list,ld.list,lambda.list=lambda.list,w.list=full.w.list,
effect.size.prior='Ind.Normal',label.list=label.list,
output.labels=paste0('Output/',annot,'_anno_chrom_',lambda.power),
outlier.switch=F)
