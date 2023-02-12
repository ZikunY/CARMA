 library(hapsim )
library(stringr)
library(cli)
library(glue)
library(crayon)
library(readr)
library(sim1000G)
library(Matrix)
library(pryr)
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

    a.vcf<-readVCF(paste0(ref.table$chr[i],'_',
                          as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.vcf' ),
                   maxNumberOfVariants = num.snps,min_maf=0.01,max_maf = 0.99)#####Pre process the 1000G vcf files by using PLINK/1.9.
map.file<-read.table('integrated_call_samples_v3.20130502.ALL.panel',header=T)
pop.index<-which(map.file$super_pop==pop.names[s])


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
                        as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.RDS'))
#}
}
