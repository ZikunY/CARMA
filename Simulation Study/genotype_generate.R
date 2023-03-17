#Loading libraries
library(hapsim )
library(stringr)
library(cli)
library(glue)
library(crayon)
library(readr)
library(sim1000G)
library(Matrix)
library(pryr)
#setting up jobs from cluster
args = commandArgs(TRUE)
i=as.numeric(args[1])
set.seed(i)


num.snps=100000;n=10000##Upper limit of extracted SNPs; Total number of simulated individuals

setwd('Set your working direction')
    
ref.table<-read.csv('Simulation loci.csv') ## The file contains the position information of loci from breast cancer GWAS, used in CARMA paper


pop.names<-c('EUR')#Population index, can be changed to EAS, AFR, etc..


a.vcf<-readVCF(paste0(ref.table$chr[i],'_',
                          as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.vcf' ),
                   maxNumberOfVariants = num.snps,min_maf=0.01,max_maf = 0.99)#####reading the loucs from 1000G vcf file, which is pre-processed  by using PLINK/1.9.
map.file<-read.table('integrated_call_samples_v3.20130502.ALL.panel',header=T)###1000G samples file
pop.index<-which(map.file$super_pop==pop.names) ### match the population in 1000G


readGeneticMapFromFile(paste0("map_file/genetic_map_GRCh37_",ref.table$chr[i],".txt.gz")) ###1000G map file

startSimulation(a.vcf, totalNumberOfIndividuals =n) ###simulating n individuals

ids = generateUnrelatedIndividuals(n)
genotype = retrieveGenotypes(ids)
genotype<-as(genotype,"dgCMatrix")

#Save the simulated genotype matrix and variants info

writeMM(genotype,file=paste0('genotype_matrix/',pop.names,'/',
                             ref.table$chr[i],'_',
                             as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.mtx'))
saveRDS(a.vcf,file=paste0('vcf_files/',pop.names,'/',
                          ref.table$chr[i],'_',
                        as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.RDS'))
