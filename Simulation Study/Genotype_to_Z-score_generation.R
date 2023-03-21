###Making the dirctory
#git clone https://github.com/ZikunY/CARMA.git
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

##Upper limit of extracted SNPs; Total number of simulated individuals
num.snps=100000;n=100

setwd('CARMA/Simulation Study')
## The file contains the position information of loci from breast cancer GWAS, used in CARMA paper
ref.table<-read.csv('Simulation loci.csv')

#Population index, can be changed to EAS, AFR, etc..
pop.names<-c('EUR')

#####1000G vcf file can be found in the link /vol1/ftp/release/20130502
#####reading the loucs from 1000G vcf file, which is pre-processed  by using PLINK/1.9.
#####Here we use the first locus in the csv file to demonstrate the simulation.

a.vcf<-readVCF(paste0(ref.table$chr[1],'_',
                          as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.vcf' ),
                   maxNumberOfVariants = num.snps,min_maf=0.01,max_maf = 0.99)
###1000G samples file
map.file<-read.table('integrated_call_samples_v3.20130502.ALL.panel',header=T)
### match the population in 1000G
pop.index<-which(map.file$super_pop==pop.names)
###The vcf file is pre-processed that only contains 503 EUR individuals; The map file is the complete list of 2504 individuals in 1000G;

###1000G map file
readGeneticMapFromFile(paste0("genetic_map_GRCh37_",ref.table$chr[1],".txt.gz"))

###simulating n individuals
startSimulation(a.vcf, totalNumberOfIndividuals =n)

ids = generateUnrelatedIndividuals(n)
genotype = retrieveGenotypes(ids)
genotype<-as(genotype,"dgCMatrix")



#####Reading DeepSEA
#####DeepSEA file can be generated in http://deepsea.princeton.edu/job/analysis/create/
whole.deepsea<-read.table(paste0(ref.table$chr[1],'_',
                                 as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.infile.vcf.out.evalue'),header=T)
whole.deepsea<-whole.deepsea[order(whole.deepsea$pos),] # order


### The folder name of scenarios regarding number of causal variants per locus
batch.index<-c('first_batch_genome_wise_noisy',
               'sec_batch_genome_wise_noisy',
               'third_batch_genome_wise_noisy')
### Here we set the number of causal variant equal to 1 as an example. 
num.causal=1
## learning chr of the locus
chr<-chor.index<-as.numeric(substr(ref.table[1,1],4,5))


#explained variance
phi_all=0.0075


s=1
ld.list<-list()
p.list<-list()
deepsea.list<-list()
true.x.list<-list()
w.list<-list()
deepsea.list<-list()
maf.list<-list()


true.x<-genotype

true.x<-as.matrix(true.x)
p.list[[s]]<-ncol(true.x)
maf<-apply(true.x,2,function(x){sum(x)/(length(x)*2)})
maf.index<-which(maf==0|maf==1)### Prevention measure, if n is very small, there would not be any alt allele.
print(maf.index)
if(length(maf.index)!=0){
  maf.prob<-c(length(which(true.x[,-maf.index]==0)),length(which(true.x[,-maf.index]==1)),length(which(true.x[-maf.index]==2)))/
    (prod(dim(true.x[,-maf.index])))
  for(h in 1:length(maf.index)){
    true.x[,maf.index[h]]<-sample(c(0,1,2),n,prob = maf.prob,replace = T)
  }
}
maf.list[[s]]<-maf
rm(maf)
true.x.list[[s]]<-true.x
x<-scale(true.x)
cor.x<-cor(true.x) # compute ld
print(ref.table[1,])
print(paste0('This is check; ',sum(is.na(cor.x))))
colnames(cor.x)<-row.names(cor.x)<-1:(p.list[[s]])
ld.list[[s]]<-cor.x

# Extract annotations for the simulated variants
deepsea<-whole.deepsea[,(ncol(whole.deepsea)-918):ncol(whole.deepsea)]
deepsea<-log(deepsea)
deepsea<-abs(deepsea)
deepsea<-scale(deepsea)
deepsea.list[[s]]<-deepsea

# theta vector that was previously sampled, so it is genome-wise true theta, all loci use same file
theta<-read.csv("genome_wise_theta.csv")
theta<-as.matrix(theta)
##########genome-wise
# Same true 100 functional annotations.
deepsea.index<-round(seq(1,919,length.out=100))
prior.prob.list<-list()
true.index.list<-list()


#save annotaitons
w.list[[s]]<-as.matrix(cbind(rep(1,p.list[[s]]),as.data.frame(deepsea.list[[s]][,deepsea.index])))
#save prior probability
prior.prob.list[[s]]<-(w.list[[s]]%*%theta)
top.index<-1:10
### Sample the true causal for the locus, subjected to the LD restriction and prior probabilities
chor.index<-1 #Chr 1
repeat{
  for(s in 1:length(chor.index)){
    if(num.causal>1){
      
      temp.index<-apply(as.matrix(order(prior.prob.list[[s]],decreasing = T)[top.index]),1,function(x){(sum(abs(ld.list[[s]][x,])>0.90)>=2)& (sum(abs(ld.list[[s]][x,])>0.90)<20)})
      if(sum(temp.index)>=num.causal){
        temp.ld<-ld.list[[s]][order(prior.prob.list[[s]],decreasing = T)[c(top.index)[temp.index]],
                              order(prior.prob.list[[s]],decreasing = T)[c(top.index)[temp.index]]]
        
        posi.ch<-which(apply( combn(1:nrow(temp.ld),num.causal),2,function(x){max(abs((temp.ld[x,x])[lower.tri(matrix(1,num.causal,num.causal))]))<0.1}))
        if(any(posi.ch)){
          
          temp.rank<-apply( (combn(1:nrow(temp.ld),num.causal)[,posi.ch,drop=F]),2,function(x){  match(as.numeric(colnames(temp.ld))[x],order(prior.prob.list[[s]],decreasing = T)) })
          posi.index<-which.min(colSums(temp.rank))
          true.index.list[[s]]<-as.numeric(colnames(temp.ld)[combn(1:nrow(temp.ld),num.causal)[,posi.ch[posi.index]]])
        }
      }
    }else{
      temp.index<-apply(as.matrix(order(prior.prob.list[[s]],decreasing = T)[top.index]),1,function(x){(sum(abs(ld.list[[s]][x,])>0.90)>=2)& (sum(abs(ld.list[[s]][x,])>0.90)<20)})
      if(sum(temp.index)>=num.causal){
        temp.index<-which(temp.index)
        true.index.list[[s]]<-order(prior.prob.list[[s]],decreasing = T)[c(top.index)[temp.index[1]]]
      }
    }
    
    
  }
  
  if(length(true.index.list)==length(chor.index)){
    break
  }
  top.index<-1:(max(top.index)+1)
  print(top.index)
}


top.rank<-match(true.index.list[[s]],order(prior.prob.list[[s]],decreasing = T))


if(!dir.exists(paste0(batch.index[num.causal],'/',
                      'data/data_plot/'))){
  dir.create(paste0(batch.index[num.causal],'/',
                    'data/data_plot/'),recursive = T)
}
setwd(paste0(batch.index[num.causal],'/'))
# making figure for Z score and the corresponding prior probability generated by annotations
pdf(paste0('data/data_plot/', ref.table$chr[1],'_',
           as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.pdf') ,width = 10,height=15)
par(mfrow=c(2,length(phi_all)))

beta.index<-rep(0,p.list[[s]])
beta.index[true.index.list[[s]]]<-1


plot(prior.prob.list[[s]],main=paste0( paste0(ref.table[1,],collapse = ';')),ylab='Prior Probs.')
points(which(beta.index==1),(prior.prob.list[[s]])[which(beta.index==1)],col=2,pch=19)


beta.value<-rnorm(sum(beta.index),0,sd=0.5)
true.beta<-beta.index
true.beta[which(beta.index==1)]=beta.value
true.mu<-as.matrix(true.x.list[[s]]%*%as.matrix(true.beta))
phi=phi_all[s]
###true regression model
y<-rnorm(n,mean=true.mu,sd=sqrt(var(true.mu)*(1-phi)/phi))
print(paste0('Proportion; ', (var(true.mu))/
               (t(y-true.x.list[[s]]%*%true.beta)%*%(y-true.x.list[[s]]%*%true.beta)/n+var(true.mu))))
y<-scale(y)     
stand.x<-scale(true.x.list[[s]])
lm.b.se<-apply(as.matrix(stand.x),2,function(i) {l<-lm(y~as.matrix(i)-1);s<-summary(l);return(c(coefficients(s)[1:2]))})

#compute Z scores
lm.z<-lm.b.se[1,]/lm.b.se[2,]
plot(lm.z,main='Summary statistics')
points(which(beta.index==1),(lm.z)[which(beta.index==1)],col=2,pch=19)

z.score<-data.frame(Zscore=lm.z) 
if(dir.exists(paste0('data/'))==F){
  dir.create(paste0('data/'))
}

############Annotations############
full_anno<-cbind(1,deepsea.list[[s]][,c(deepsea.index,(1:919)[-deepsea.index][1:100])])# Adding 100 noisy annotations.
colnames(full_anno)[1]<-'Coding'
###############Save CARMA data
write.table(z.score, paste0('data/',
                            ref.table$chr[1],'_',
                            as.character(ref.table$region_start [1]),'_',ref.table$region_end[1]),row.names = F,quote = F)

write.table(y, paste0('data/',
                      ref.table$chr[1],'_',
                      as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.y'),row.names = F,quote = F,col.names = F)    
write.table(full_anno,paste0('data/',
                             ref.table$chr[1],'_',
                             as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.deepsea'),
            row.names = F,quote = F,col.names = T)

write.table(top.rank,paste0('data/',
                            ref.table$chr[1],'_',
                            as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'_rank.txt'),
            row.names = F,quote = F,col.names = T)

#################Save SuSiE data
if(dir.exists(paste0('data/SuSiE/' ))==F){
  dir.create(paste0('data/SuSiE/' ),recursive = T)
}
write.table(lm.b.se[1,],paste0('data/SuSiE/',
                               ref.table$chr[1],'_',
                               as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.beta')
            ,row.names = F,quote = F,col.names = F)
write.table(lm.b.se[2,],paste0('data/SuSiE/',
                               ref.table$chr[1],'_',
                               as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.se')
            ,row.names = F,quote = F,col.names = F)
#####Save the true causal effect size
write.table(true.beta,file=paste0('data/',
                                  ref.table$chr[1],'_',
                                  as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.truebeta'),row.names = F,
            col.names = F,quote = F)

dev.off()  



