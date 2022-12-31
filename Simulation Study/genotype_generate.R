 library(hapsim,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/" )
library(stringr,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/")
library(cli,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/")
library(glue,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/")
library(crayon,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/")
library(readr,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/")
library(sim1000G,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/")
library(Matrix,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/")
library(pryr,lib="/ifs/scratch/msph/eigen/zy2412/R_libs/")
args = commandArgs(TRUE)
i=as.numeric(args[1])
#rm(list=ls())
set.seed(i)
num.snps=100000;n=10000

print(i)
 setwd('Set your working direction')
    


ref.table<-read.csv('Simulation loci.csv')

print(ref.table[i,])

pop.names<-c('EUR')#Population index, can be changed to EAS, AFR, etc..
for(s in 1:length(pop.names)){

    a.vcf<-readVCF(paste0('reduced_vcf_file_complete_locus/',
                          ref.table$chr[i],'_',
                          as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.vcf' ),
                   maxNumberOfVariants = num.snps,min_maf=0.01,max_maf = 0.99)#####Pre process the 1000G vcf files by using PLINK/1.9.
map.file<-read.table('integrated_call_samples_v3.20130502.ALL.panel',header=T)
pop.index<-which(map.file$super_pop==pop.names[s])
a.vcf$gt1<-a.vcf$gt1[,pop.index]
a.vcf$gt2<-a.vcf$gt2[,pop.index]
a.vcf$individual_ids<-a.vcf$individual_ids[pop.index]
maf<-rowSums(a.vcf$gt1+a.vcf$gt2)/(2*dim(a.vcf$gt1)[1])
maf.index<-which(maf>0.01&maf<0.99)
a.vcf$gt1<-a.vcf$gt1[maf.index,]
a.vcf$gt2<-a.vcf$gt2[maf.index,]
a.vcf$varid<-a.vcf$varid[maf.index]

readGeneticMapFromFile(paste0("map_file/genetic_map_GRCh37_",ref.table$chr[i],".txt.gz"))

startSimulation(a.vcf, totalNumberOfIndividuals =n)

ids = generateUnrelatedIndividuals(n)
genotype = retrieveGenotypes(ids)
genotype<-as(genotype,"dgCMatrix")
print(object_size(genotype))
writeMM(genotype,file=paste0('genotype_matrix/',pop.names[s],'/',
                             ref.table$chr[i],'_',
                             as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.mtx'))
saveRDS(a.vcf,file=paste0('vcf_files/',pop.names[s],'/',
                          ref.table$chr[i],'_',
                        as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.RData'))
#}
}
