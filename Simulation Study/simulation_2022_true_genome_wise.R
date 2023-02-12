library(hapsim )
library(stringr)
library(cli)
library(crayon)
library(readr)
library(sim1000G)
library(Matrix)
library(pryr)
library(susieR)
library(CARMA )
rm(list=ls())
args = commandArgs(TRUE)
locus=as.numeric(args[1])
num.causal=as.numeric(args[2])
set.seed(args+50)


setwd('Set your working direction')
   
ref.table<-read.csv('Simulation loci.csv')

n=10000
top.index<-1:10
total.num=nrow(ref.table)
batch.index<-c('first_batch_genome_wise_noisy',
               'sec_batch_genome_wise_noisy',
               'third_batch_genome_wise_noisy')
chor.index<-locus
chr<-ref.table[chor.index,1]
chr<-as.numeric(substr(chr,4,5))



whole.deepsea<-read.table(paste0('deepsea/',paste0('chr',chr),'/infile.vcf.out.evalue'),sep = ',',header=T)
whole.deepsea<-whole.deepsea[order(whole.deepsea$pos),]

phi_all=c(0.3)
phi_all=phi_all/(40)
labels=c('large')


s=1
ld.list<-list()
p.list<-list()
deepsea.list<-list()
true.x.list<-list()
w.list<-list()
deepsea.list<-list()
maf.list<-list()
for(i in chor.index){

true.x<-readMM(paste0('genotype_matrix/',
                      ref.table$chr[i],'_',
                      as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.mtx'))

true.x<-as.matrix(true.x)
p.list[[s]]<-ncol(true.x)
maf<-apply(true.x,2,function(x){sum(x)/(length(x)*2)})
maf.index<-which(maf==0)### Prevention measure.
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
cor.x<-cor(true.x)
print(ref.table[i,])
print(paste0('This is check; ',sum(is.na(cor.x))))
colnames(cor.x)<-row.names(cor.x)<-1:(p.list[[s]])
ld.list[[s]]<-cor.x
s=s+1
}
s=1


for(i in chor.index){
pos.index<-which((whole.deepsea$pos<=ref.table$region_end[i] & whole.deepsea$pos>=ref.table$region_start[i]))
if(length(pos.index)<=p.list[[s]]){
  if(length(pos.index)!=0){
    sup.index<-sample((1:nrow(whole.deepsea))[-pos.index],p.list[[s]]-length(pos.index))
    pos.index<-c(pos.index,sup.index)
  }else{
    pos.index<-sample(1:nrow(whole.deepsea),p.list[[s]])
  }
}else{
  pos.index<-sample(pos.index,p.list[[s]],replace = F)
}

deepsea<-whole.deepsea[pos.index,(ncol(whole.deepsea)-918):ncol(whole.deepsea)]
deepsea<-log(deepsea)
deepsea<-abs(deepsea)
deepsea<-scale(deepsea)
deepsea.list[[s]]<-deepsea


s=s+1
}

s=1
theta<-read.csv("genome_wise_theta.csv")# theta vector that was previously sampled, so it is genome-wise true theta.
theta<-as.matrix(theta)
##########genome-wise
deepsea.index<-round(seq(1,919,length.out=100))# Same true 100 functional annotations.
prior.prob.list<-list()
true.index.list<-list()
w.list[[s]]<-as.matrix(cbind(rep(1,p.list[[s]]),as.data.frame(deepsea.list[[s]][,deepsea.index])))
prior.prob.list[[s]]<-(w.list[[s]]%*%theta)

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


ss=s=1
top.rank<-match(true.index.list[[s]],order(prior.prob.list[[s]],decreasing = T))

for(i in chor.index){

    
    setwd(paste0(batch.index[num.causal]))

pdf(paste0('data/data_plot/', ref.table$chr[i],'_',
           as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.pdf') ,width = 10,height=15)
par(mfrow=c(2,length(phi_all)))

  beta.index<-rep(0,p.list[[ss]])
  beta.index[true.index.list[[ss]]]<-1
  
  
  plot(prior.prob.list[[ss]],main=paste0( paste0(ref.table[i,],collapse = ';'),'; ',labels[s]),ylab='Prior Probs.')
  points(which(beta.index==1),(prior.prob.list[[ss]])[which(beta.index==1)],col=2,pch=19)
  

  beta.value<-rnorm(sum(beta.index),0,sd=0.5)
  true.beta<-beta.index
  true.beta[which(beta.index==1)]=beta.value
  true.mu<-as.matrix(true.x.list[[ss]]%*%as.matrix(true.beta))
  phi=phi_all[s]
  y<-rnorm(n,mean=true.mu,sd=sqrt(var(true.mu)*(1-phi)/phi))
  print(paste0('Proportion; ', (var(true.mu))/
                 (t(y-true.x.list[[ss]]%*%true.beta)%*%(y-true.x.list[[ss]]%*%true.beta)/n+var(true.mu))))
  y<-scale(y)     
  stand.x<-scale(true.x.list[[ss]])
  lm.b.se<-apply(as.matrix(stand.x),2,function(i) {l<-lm(y~as.matrix(i)-1);s<-summary(l);return(c(coefficients(s)[1:2]))})
   
  lm.z<-lm.b.se[1,]/lm.b.se[2,]
  plot(lm.z,main='Summary statistics')
  points(which(beta.index==1),(lm.z)[which(beta.index==1)],col=2,pch=19)
  
  z.score<-data.frame(Zscore=lm.z) 
  if(dir.exists(paste0('data/',labels[s]))==F){
    dir.create(paste0('data/',labels[s]))
  }
  
  ############Annotations############
  full_anno<-cbind(1,deepsea.list[[ss]][,c(deepsea.index,(1:919)[-deepsea.index][1:100])])# Adding 100 noisy annotations.
  colnames(mid_anno)[1]<-'Coding'
  colnames(full_anno)[1]<-'Coding'
  ############### CS data
  write.table(z.score, paste0('data/',labels[s],'/',
                              ref.table$chr[i],'_',
                              as.character(ref.table$region_start [i]),'_',ref.table$region_end[i]),row.names = F,quote = F)
  
  write.table(y, paste0('data/',labels[s],'/',
                        ref.table$chr[i],'_',
                        as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.y'),row.names = F,quote = F,col.names = F)    
  write.table(full_anno,paste0('data/',labels[s],'/',
                               ref.table$chr[i],'_',
                               as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.deepsea'),
              row.names = F,quote = F,col.names = T)
 
  write.table(top.rank,paste0('data/',labels[s],'/',
  ref.table$chr[i],'_',
  as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'_rank.txt'),
row.names = F,quote = F,col.names = T)
  
  #################SuSiE data
  if(dir.exists(paste0('data/SuSiE/',labels[s] ))==F){
    dir.create(paste0('data/SuSiE/',labels[s] ),recursive = T)
  }
  write.table(lm.b.se[1,],paste0('data/SuSiE/',labels[s],'/',
                                 ref.table$chr[i],'_',
                                 as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.beta')
              ,row.names = F,quote = F,col.names = F)
  write.table(lm.b.se[2,],paste0('data/SuSiE/',labels[s],'/',
                                 ref.table$chr[i],'_',
                                 as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.se')
              ,row.names = F,quote = F,col.names = F)
  ################PAINTOR data
  write.table(true.beta,file=paste0('data/',labels[s],'/',
                                    ref.table$chr[i],'_',
                                    as.character(ref.table$region_start [i]),'_',ref.table$region_end[i],'.truebeta'),row.names = F,
              col.names = F,quote = F)

dev.off()  
ss=ss+1
}
