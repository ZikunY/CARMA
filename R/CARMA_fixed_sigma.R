#' CARMA (fixed variance)
#' 
#' Performs a Bayesian fine-mapping model in order to identify putative causal variants at GWAS loci. The model requires the summary statistics
#' of the SNPs in the testing loci, the corresponding LD matrices for fine-mapping, and an estimated variance of traits. Functional annotations can be included as the prior 
#' information of the causality of the testing SNPs. The model also provides a procedure of outlier detection, which resolves the discrepancies
#' between the summary statistics and the LD matrix extracted from reference panels. The model can be executed chromosome-wise to increase power. 
#'@param z.list Input list of the summary statistics of the testing loci, and each element of the list is the summary statistics of each individual locus.
#'@param ld.list Input list of the LD correlation matrix of the testing loci, and each element of the list is the LD matrix of each individual locus.
#'@param w.list Input list of the functional annotations of the testing loci, and each element of the list is the functional annotation matrix of each individual locus.
#'@param lambda.list The default input is NULL as \eqn{\eta} is chosen through an adaptive procedure.
#'@param label.list Input list of the names of the testing loci. Default is NULL. 
#'@param effect.size.prior The prior of the effect size. The choice are 'Cauchy' and 'Hyper-g' priors, where the Cauchy prior is the default prior.
#'@param input.alpha The elastic net mixing parameter, where \eqn{0\le}\eqn{\alpha}\eqn{\le 1}.
#'@param y.var The input variance of the summary statistics, the default value is 1 as the summary statistics are standardized. 
#'@param rho.index A number between 0 and 1 specifying \eqn{\rho} of the estimated credible sets.
#'@param BF.index  A number greater than 1 to specifying the threshold of the Bayes factor of the estimated credible models.
#'@param outlier.switch The indicator variable of whether turn on the outlier detection. We suggest that the detection should always turn on. 
#'@param outlier.threshold The Bayes threshold of the hypothesis testing of determining outliers.
#'@param outlier.cor.range.threshold The correlation threshold of defining the highly correlated group, within which the outlier detection run. The default is 0.9.
#'@param num.causal The maximum number of causal variants assumed per locus.
#'@param Max.Model.Dim Maximum number of the top candidate models based on the unnormalized posterior probability. 
#'@param all.inner.iter Maximum iterations for Shotgun algorithm to run per iteration within EM algorithm.
#'@param all.iter Maximum iterations for EM algorithm to run.
#'@param output.labels Output directory where output will be written while CARMA is running. Default is NULL.
#'@param epsilon.threshold Convergence threshold measured by average of Bayes factors.
#'@return The form of the return is a list, for each list:
#'\itemize{
#'\item pip - The posterior inclusion probability of each individual locus.
#'\item Credibleset - The information regarding the credible set given a threshold \eqn{\rho}.
#'\item Credible model - The information regarding the credible model given a threshold  of the Bayes factor.
#'\item Outliers - The information regarding the detected outliers and the corresponding testing statistics of each detected outliers.
#'}
#'@details The function performs a Bayesian fine-mapping method. 
#'@examples 
#'# Example 
#'set.seed(1)
#'n = 400
#'p = 500
#'beta = rep(0,p)
#'beta[1] = 1
#'X = matrix(rnorm(n*p),nrow = n,ncol = p)
#'X = scale(X,center = TRUE,scale = TRUE)
#'y = drop(X %*% beta + rnorm(n))
#'SS=compute_summary_statistics(X,y)
#'z.list<-list()
#'z.list[[1]]<-(SS$betahat/SS$sebetahat)
#'ld.list<-list()
#'ld.list[[1]]<-cov(X)
#'CARMA.result<-CARMA_fixed_sigma(z.list,ld.list=ld.list,effect.size.prior='Hyper-g')
CARMA_fixed_sigma<-function(z.list,ld.list,w.list=NULL,lambda.list=NULL,output.labels=NULL,label.list=NULL,
                                 effect.size.prior='Cauchy',rho.index=0.99,BF.index=10,
                                Max.Model.Dim=1e+4,all.iter=10,all.inner.iter=10,input.alpha=0.5,epsilon.threshold=1e-3,
                                 num.causal=10,y.var=1,outlier.switch=T,outlier.threshold=1e-3,outlier.cor.range.threshold=0.9){
  EM.M.step.func<-function(Model.space=NULL,w=w,input.alpha=0.5){
      count.index<-Model.space
        print('this is running!!!')
        cv.poisson<-cv.glmnet(w,count.index,family = 'poisson',alpha=input.alpha,type.measure='deviance' )
        print('this is running!!!')
        cv.index<-which(cv.poisson$lambda==cv.poisson$lambda.min)
        print('this is running!!!')
        glm.beta<-as.matrix(c(cv.poisson$glmnet.fit$a0[cv.index],cv.poisson$glmnet.fit$beta[-1,cv.index]))
     
        #print(cv.index)   
    
      print(paste0('This is starting theta: ',glm.beta[1:ifelse(length(glm.beta)>10,10,length(glm.beta)),]))
  
    return(glm.beta=glm.beta)
  }
  credible.set.fun.improved<-function(pip,ld,true.beta=NULL,rho=0.99){
    candidate.r<-c((seq(from=.5,0.95,by=0.05)),seq(0.96,0.99,0.01))
    #pip=cs.results[[1]]$PIPs
    snp.list<-list()
    colnames(ld)<-rownames(ld)<-1:nrow(ld)
    for(r in 1:length(candidate.r)){
      working.ld<-ld
      cor.threshold<-candidate.r[r]
      pip.order<-order(pip,decreasing = T)
      snp.list[[r]]<-list()
      group.index<-1
      s<-1
      while(sum(pip[pip.order[s:length(pip.order)]])>rho){
        
        #print(s)
        cor.group<-as.numeric(names(which(abs(working.ld[which(pip.order[s]==as.numeric(colnames(working.ld))),])>cor.threshold)))
        #print(sum(as.numeric(colnames(working.ld))==sort(pip.order)))
        if(sum(pip[cor.group])>rho){
          group.pip<- pip[cor.group]
          snp.index<- cor.group[order(group.pip,decreasing = T)][1: min(which(cumsum( sort(group.pip,decreasing = T))>rho))]
          snp.list[[r]][[group.index]]<-snp.index
          group.index<-group.index+1
          pip.order<-pip.order[-match(snp.index,pip.order)]
          working.ld<-working.ld[-match(snp.index,as.numeric(colnames(working.ld))),-match(snp.index,as.numeric(colnames(working.ld)))]
          #dim(working.ld)
          #print(length(snp.index))
        }else{
          s=s+1
        }
      }
      
    }
    if(sum(sapply(snp.list,length))!=0){
      group.index<-max(which(sapply(snp.list,length)==max(sapply(snp.list,length))))
      credible.set.list<-snp.list[[group.index]]
      if(!is.null(true.beta)){
        purity<-c()
        for(s in 1:length(credible.set.list)){
          purity<-c(purity, (mean(ld[ credible.set.list[[s]], credible.set.list[[s]]]^2)))
        }
        causal.snp<-ifelse(length(na.omit(match(unlist(credible.set.list),true.beta)))!=0,
                           length(na.omit(match(unlist(credible.set.list),true.beta))),0)
        length.credible<-length(credible.set.list)
        return(list(c(causal.snp,
                      ifelse(length.credible==0,NA,length.credible),
                      mean(sapply(credible.set.list,length)),
                      mean(purity)),credible.set.list))
        
      }else{
        purity<-c()
        for(s in 1:length(credible.set.list)){
          purity<-c(purity, (mean(ld[ credible.set.list[[s]], credible.set.list[[s]]]^2)))
        }
        length.credible<-length(credible.set.list)
        return(list(c(ifelse(length.credible==0,NA,length.credible),
                      mean(sapply(credible.set.list,length)),
                      mean(purity)),
                    credible.set.list))
      }
    }else{
      return(list(rep(0,4),list()))
    }
  }
  credible.model.fun<-function(likelihood,model.space,bayes.threshold=10){
    post.like.temp<-likelihood-likelihood[1]
    post.prob<-exp(post.like.temp)/(sum(exp(post.like.temp)))
    bayes.f<-post.prob[1]/post.prob
    candidate.model<-1
    
    credible.model.list<-list()
    credible.model.list[[1]]<-list()
    input.rs<-c()
    for(ss in 1:length(which(bayes.f<bayes.threshold))){
      credible.model.list[[1]][[ss]]<-which(model.space[ss,]==1)
      input.rs<-c(input.rs,which(model.space[ss,]==1))
    }
    credible.model.list[[2]]<-data.frame(Posterior.Prob=post.prob[which(bayes.f<bayes.threshold)])
    credible.model.list[[3]]<-unique(input.rs)
    return(credible.model.list )
    
    
  }
  eta.adapt.function<-function(input.z,input.ld,quan=c(0.99)){
    
    eta.setting.func<-function(p.power){
      dim.model<-1
      result<-dim.model*log(p^(-p.power))-p^(-p.power)+lfactorial(p-dim.model)-lfactorial(p)
      return(result)
    }
    
    
    
    tau.sample<-rgamma(1e+5,.5,.5)
    input.z<-c(input.z)
    p<-length(input.z)
   
    eta.quan<-c()
    for(s in 1:length(quan)){
      z.index<-which(abs(input.z)==min(abs(input.z)[which(abs(input.z)>=quantile(abs(input.z),quan[s]))]))[1]
      
      all.index<-which(eta.setting.func(seq(0,10,by=0.01))+marginal_likelihood(z.index,input.ld,input.z,tau=tau.sample,p_S=1,1)>0)
      if(length(all.index)!=0){
        eta.quan<-c(eta.quan,seq(0,10,by=0.01)[max(all.index)])
      }else{
        eta.quan<-c(eta.quan,0)
      }
    }
    result.lambda.power<-c(eta.quan)
    return(result.lambda.power)
    
  }
  ########## Input data#########
  {
    log.2pi<-log(2*pi)
    L<-length(z.list)
    p.list<-list()
    for(i in 1:L){
      print(nrow(z.list[[i]]))
      z.list[[i]]<-as.matrix(z.list[[i]])
      p.list[[i]]<-nrow(z.list[[i]])
    }
    B<-Max.Model.Dim
    all.B.list<-list()
    for(i in 1:L){
    all.B.list[[i]]<-list()
    all.B.list[[i]][[1]]<-integer(0)
    all.B.list[[i]][[2]]<-Matrix(nrow = 0,ncol=p.list[[i]],data=0,sparse = T)
    }
    q.list<-list()
    if(!is.null(w.list)){
      for(i in 1:L){
      q.list[[i]]<-ncol(w.list[[i]])
      w.list[[i]]<-as.matrix(cbind(1,scale(w.list[[i]][,-1])))
      }
    }
    if(is.null(label.list)){
      for(i in 1:L){
        label.list[[i]]=paste0('locus_',i)
      }
    }
    Sigma.list<-list()
    for(i in 1:L){
    Sigma.list[[i]]<-as.matrix(ld.list[[i]])
    }
    S.list<-list()
    for(i in 1:L){
    S.list[[i]]<-integer(0)
    }
    
    all.C.list<-list()
    for(i in 1:L){
      all.C.list[[i]]<-list()
      all.C.list[[i]][[1]]<-integer(0)
      all.C.list[[i]][[2]]<-Matrix(nrow = 0,ncol=p.list[[i]],data=0,sparse = T)
    }
   
    all.epsilon.threshold<-0
    epsilon.list<-list()
    for(i in 1:L){
      epsilon.list[[i]]<-epsilon.threshold*p.list[[i]]
      all.epsilon.threshold<-all.epsilon.threshold+  epsilon.threshold*p.list[[i]]
    }
    model.prior='Poisson'
    standardize.model.space=T
    
    if(is.null(lambda.list)){
    if(effect.size.prior=='Cauchy'){
      marginal_likelihood=Cauchy_fixed_sigma_marginal
    }
    if(effect.size.prior=='Hyper-g'){
      marginal_likelihood=hyper_g_fixed_sigma_marginal
    }
     for(i in 1:L){
       print(eta.adapt.function(z.list[[i]],input.ld = Sigma.list[[i]]))
       lambda.list[[i]]<-1/( p.list[[i]]^eta.adapt.function(z.list[[i]],input.ld = Sigma.list[[i]]))
       }
    }
    print(lambda.list)
  }

#########Module model##########
  
  Module.Cauchy.Shotgun<-function(z,ld.matrix,Max.Model.Dim=1e+4,input.S=NULL,lambda,label,
                                              num.causal=10,output.labels,y.var=1,effect.size.prior=effect.size.prior,model.prior=model.prior,
                                              outlier.switch,outlier.cor.range.threshold,outlier.threshold=1e-3,input.conditional.S.list=NULL,
                                              C.list=NULL,prior.prob=NULL,epsilon=1e-3,inner.all.iter=10){
       {
       if(model.prior=='input.prob'){
         p<-nrow(z);
         poisson.prior<-function(dim.model){
           result<-dim.model*log(lambda)-lambda+lfactorial(p-dim.model)-lfactorial(p)
           return(result)
         }
         prior.prob.list<-list()
         posi.log.pro.list<-list()
         nega.log.pro.list<-list()
         for(pp in 1:100){
           top.index<-order(prior.prob,decreasing = T)[1:pp]
           if((poisson.prior(pp)-(sum(log(1-prior.prob[-top.index]))+sum(log(prior.prob[top.index]))))>0){
             posi.log.pro.list[[pp]]<-log(prior.prob)
             nega.log.pro.list[[pp]]<-log(1-prior.prob)
           }else{
             balance.value<-(poisson.prior(pp)-(sum(log(1-prior.prob[-top.index]))+sum(log(prior.prob[top.index]))))/p
             balance.value<-exp(balance.value)
             prior.prob.list[[pp]]<-1-balance.value+balance.value* prior.prob
             prior.prob.list[[pp]][top.index]<-balance.value* prior.prob[top.index]
             posi.log.pro.list[[pp]]<-log(prior.prob.list[[pp]])
             nega.log.pro.list[[pp]]<-log(1-prior.prob.list[[pp]])
           }
           
         }
         
         input.prior.dist<-function(x){
           variable.index<-which(x==1)
           pp<-length(variable.index)
           if(any(variable.index)){
             return( sum((nega.log.pro.list[[pp]][-variable.index]))+sum(( posi.log.pro.list[[pp]][variable.index])))
           }else{
             return(poisson.prior(0))
           }
         }
         prior.dist<-input.prior.dist
       }
       
      
       if(model.prior=='Poisson'){
        Poisson.prior.dist<-function(t){
          dim.model<-sum(t)
          result<-dim.model*log(lambda)-lambda+lfactorial(p-dim.model)-lfactorial(p)
          return(result)
        }
        prior.dist<-Poisson.prior.dist
      }
      
      ###### Define marginal likelihood#########  
      
         if(effect.size.prior=='Cauchy'){
           marginal_likelihood=Cauchy_fixed_sigma_marginal
         }
         if(effect.size.prior=='Hyper-g'){
           marginal_likelihood=hyper_g_fixed_sigma_marginal
         }
      
      
      p<-nrow(z);
      log.2pi<-log(2*pi)
      tau.sample<-rgamma(1e+5,0.5,rate=0.5)
      B<-Max.Model.Dim
      stored.result.prob<-rep(0,p)
      stored.bf<-0
      Sigma<-as.matrix(ld.matrix)
      #########
      if(!is.null(input.S)){
        S<-input.S
      }else{
        S<-integer(0)
      }
      conditional.S=NULL
      null.model<-Matrix(nrow = 1,ncol=p,data=0,sparse = T)
      null.margin<-prior.dist(null.model)
      if(is.null(C.list)){
        C.list<-list()
        C.list[[2]]<-list()
        C.list[[1]]<-list()
      }
           B.list<-list()
           B.list[[1]]<-prior.dist(null.model)
           B.list[[2]]<-Matrix(nrow = 1,ncol=p,data=0,sparse = T)
         if(length(input.conditional.S.list)==0){
           conditional.S.list<-list()
           conditional.S=NULL
         }else{
           conditional.S.list<-input.conditional.S.list
           conditional.S<-c()
           for(t in 1:length(conditional.S.list)){
             conditional.S<-c(conditional.S,conditional.S.list[[t]]$Index)
           }
           conditional.S<-unique(conditional.S)
           S<-conditional.S
         }
    }
    
  
    
    set.gamma.func<-function(input.S,condition.index=NULL){
       set.gamma.func.base<-function(S){
         add.function<-function(y){results<-(apply(as.matrix(S_sub),1,function(x){return(sort(c(x,y)))}))
         return(t(results))
         }
         
         set.gamma<-list()
         for(i in 1:3){
           set.gamma[[i]]<-c()
         }
         
         #set of gamma-
         if(length(S)==0){
           S_sub<-(1:p)
           set.gamma[[1]]<-c()
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(c(x,S))})
           set.gamma[[2]]<-as.matrix(set.gamma[[2]])
           set.gamma[[3]]<-c()
         }
         if(length(S)>1){
           S_sub<-(1:p)[-S]
           set.gamma[[1]]<-t(combn(S,length(S)-1))
           if(length(S)>2){
             set.gamma[[1]]<-t(apply(as.matrix(set.gamma[[1]]),1,sort))
           }
           #set of gamma+
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
           set.gamma[[2]]<-t(set.gamma[[2]])
           #set of gamma=
           set.gamma[[3]]<-add.function(set.gamma[[1]][1,])
           for(i in 2:nrow(set.gamma[[1]])){
             set.gamma[[3]]<-rbind(set.gamma[[3]],add.function(set.gamma[[1]][i,]))  
           }
         }
         if(length(S)==1){
           S_sub<-(1:p)[-S]
           set.gamma[[1]]<-t(combn(S,length(S)-1))
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
           set.gamma[[2]]<-t(set.gamma[[2]])
           set.gamma[[3]]<-t(add.function(set.gamma[[1]][1,]))
         }
         return(set.gamma)
       }
       set.gamma.func.conditional<-function(input.S,condition.index){
         add.function<-function(y){results<-(apply(as.matrix(S_sub),1,function(x){return(sort(c(x,y)))}))
         return(t(results))
         }
         
         set.gamma<-list()
         for(i in 1:3){
           set.gamma[[i]]<-c()
         }
         S=input.S[-match(condition.index,input.S)]
         
         #set of gamma-
         if(length(S)==0){
           S=integer(0)
           S_sub<-(1:p)[-condition.index]
           set.gamma[[1]]<-c()
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(c(x,S))})
           set.gamma[[2]]<-as.matrix(set.gamma[[2]])
           set.gamma[[3]]<-c()
         }
         if(length(S)==1){
           S_sub<-(1:p)[-input.S]
           set.gamma[[1]]<-t(combn(S,length(S)-1))
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
           set.gamma[[2]]<-as.matrix((t(set.gamma[[2]])))
           set.gamma[[3]]<-t(add.function(set.gamma[[1]][1,]))
         }
         if(length(S)>1){
           S_sub<-(1:p)[-input.S]
           if(length(S)>2){
             set.gamma[[1]]<-t(combn(S,length(S)-1))
             set.gamma[[1]]<-t(apply(as.matrix(set.gamma[[1]]),1,sort))
           }else{
             set.gamma[[1]]<-t(combn(S,length(S)-1))
           }
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
           set.gamma[[2]]<-as.matrix(t(set.gamma[[2]]))
           set.gamma[[3]]<-add.function(set.gamma[[1]][1,])
           factorial(3)
           for(i in 2:nrow(set.gamma[[1]])){
             set.gamma[[3]]<-rbind(set.gamma[[3]],add.function(set.gamma[[1]][i,]))  
           }
         }
         
         return(set.gamma)
       }
       
      if(is.null(condition.index)){
        results<-set.gamma.func.base(input.S)
      }else{
        results<-set.gamma.func.conditional(input.S,condition.index)
      }
      return(results)
    }
    duplicated.dgCMatrix <- function (dgCMat, MARGIN) {
      MARGIN <- as.integer(MARGIN)
      n <- nrow(dgCMat)
      p <- ncol(dgCMat)
      J <- rep(1:p, diff(dgCMat@p))
      I <- dgCMat@i + 1
      x <- dgCMat@x
      if (MARGIN == 1L) {
        ## check duplicated rows
        names(x) <- J
        RowLst <- split(x, I)
        is_empty <- setdiff(1:n, I)
        result <- duplicated.default(RowLst)
      } else if (MARGIN == 2L) {
        ## check duplicated columns
        names(x) <- I
        ColLst <- split(x, J)
        is_empty <- setdiff(1:p, J)
        result <- duplicated.default(ColLst)
      } else {
        warning("invalid MARGIN; return NULL")
        result <- NULL
      }
      
      if(any(is_empty)){
        out <- logical(if(MARGIN == 1L) n else p)
        out[-is_empty] <- result
        if(length(is_empty) > 1)
          out[is_empty[-1]] <- TRUE
        result <- out
      }
      
      result
    }
    match.dgCMatrix <- function (dgCMat1,dgCMat2) {
      #  dgCMat1=B.list[[2]]
      # dgCMat2=add.B[[2]]
      n1 <- nrow(dgCMat1)
      p1 <- ncol(dgCMat1)
      J1 <- rep(1:p1, diff(dgCMat1@p))
      I1 <- dgCMat1@i + 1
      x1 <- dgCMat1@x
      n2 <- nrow(dgCMat2)
      p2 <- ncol(dgCMat2)
      J2 <- rep(1:p2, diff(dgCMat2@p))
      I2 <- dgCMat2@i + 1
      x2 <- dgCMat2@x
      ## check duplicated rows
      names(x1) <- J1
      RowLst1 <- split(J1, I1)
      is_empty1<- setdiff(1:n1, I1)
      ## check duplicated rows
      names(x2) <- J2
      RowLst2 <- split(J2, I2)
      is_empty2 <- setdiff(1:n2, I2)
      result<-(match(RowLst2,RowLst1))
      if(any(which(result>is_empty1))){
        result[which(result>=is_empty1)]=result[which(result>=is_empty1)]+1
      }
      if(any(is_empty1)){
        if(any(is_empty2)){
          result<-c(is_empty1,result)  
        }
      }else{
        if(any(is_empty2)){
          result<-c(NA,result)  
        }
      }
      return(result)
      
    }
    PIP.func<-function(likeli,model.space){
      infi.index<-which(is.infinite(likeli))
      if(length(infi.index)!=0){
        likeli<-likeli[-infi.index]
        model.space<-model.space[-infi.index,]
      }
      aa<-likeli-max(likeli,na.rm=T)
      prob.sum<-sum(exp(aa))
      result.prob<-rep(NA,p)
      for(i in 1:p){
        result.prob[i]<-sum(exp(aa[ which(model.space[,i]==1)]))/prob.sum
      }
      
      # if(sum(result.prob==max(result.prob))==1){
      #   max.index<-which.max(result.prob)
      #   if(sum(Sigma[max.index,]==1)>1){
      #     perfect.index<-which(Sigma[max.index,]==1)
      #     result.prob[perfect.index]=mean( result.prob[perfect.index])
      #   }
      # }
      return(result.prob)
    }
    index.fun.inner<-function(x){
   #     print(nrow(x));print(p)
      m<-as(matrix(0,nrow=nrow(x),ncol=p),'dgTMatrix')
      m@i<-as.integer(rep(1:nrow(x)-1,each=ncol(x)))
      m@j<-as.integer(c(t(x))-1)
      m@x=rep(1,nrow(x)*ncol(x))
      m<-as(m,"dgCMatrix")
      return(m)
    }
    index.fun<-function(outer.x,Max.Model.Dimins=10){
     #   print(dim(outer.x))
     if(nrow(outer.x)>1000){
    index.bins<-which((1:nrow(outer.x))%%floor(nrow(outer.x)/Max.Model.Dimins)==0)
    #print(index.bins)
    result.m<-index.fun.inner(outer.x[1:index.bins[1],,drop=F])
    for(b in 1:(length(index.bins)-1)){
      result.m<-rbind(result.m,index.fun.inner(outer.x[(index.bins[b]+1):index.bins[b+1],,drop=F]))
    }
    if(index.bins[length(index.bins)]!=nrow(outer.x)){
      result.m<-rbind(result.m,index.fun.inner(outer.x[(index.bins[length(index.bins)]+1):nrow(outer.x),,drop=F]))
    }
     }else{
         result.m<-index.fun.inner(outer.x)
     }
    return(result.m)
     }
    ######################################################
    for(l in 1:inner.all.iter){
      for(h in 1:10){
        ##############COMPUTATION ############
        { 
          set.gamma<-set.gamma.func(S,conditional.S)  
          if(is.null(conditional.S)){
            working.S=S
            base.model<-null.model
            base.model.margin<-null.margin
          }else{
            working.S=S[-match(conditional.S,S)]
            if(length(working.S)!=0){
            base.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'dgCMatrix')
            base.model[,working.S]<-1
            p_S=length(working.S);
            base.model.margin<-marginal_likelihood(working.S,Sigma,z,tau=tau.sample,p_S=p_S,y.var)+prior.dist(base.model)
            }else{
              base.model<-null.model
              base.model.margin<-null.margin
            }
            }
          set.gamma.margin<-list()
          set.gamma.prior<-list()
          matrix.gamma<-list()
          if(length(working.S)!=0){
            S.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'dgCMatrix')
            S.model[,working.S]<-1
            p_S=length(working.S);
            current.log.margin<-marginal_likelihood(working.S,Sigma,z,tau=tau.sample,p_S=p_S,y.var)+prior.dist(S.model)
          }else{
            current.log.margin<-prior.dist(null.model)
          }
          
          if(length(working.S)>1){
            for(i in  1:length(set.gamma)){
              t0=Sys.time()
              matrix.gamma[[i]]<-index.fun(set.gamma[[i]])
              
              if(length(C.list[[2]])<ncol(set.gamma[[i]])){
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'dgCMatrix')
                C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
                computed.index<-integer(0)
              }else{
                computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
              }
              
              
              p_S=dim(set.gamma[[i]])[2]
              if(length(na.omit(computed.index))==0){
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                }
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
                set.gamma.margin[[i]]<- set.gamma.margin[[i]]+apply(matrix.gamma[[i]],1,prior.dist)
              }
              t1=Sys.time()-t0
              #print(t1)
              #print(paste0(sum(!is.na(computed.index))/length(computed.index),'%'))
            }
            
            add.B<-list()
            add.B[[1]]<-c(set.gamma.margin[[1]],
                          set.gamma.margin[[2]],
                          set.gamma.margin[[3]])
            add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'dgCMatrix')
            for(i in 1:3){
              add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[i]])
            }
           
          }
          
          if(length(working.S)==1){
            set.gamma.margin[[1]]<-null.margin
            matrix.gamma[[1]]<-null.model
            for(i in 2:3){
              
              matrix.gamma[[i]]<-index.fun(set.gamma[[i]])
              
              if(length(C.list[[2]])<ncol(set.gamma[[i]])){
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'dgCMatrix')
                C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
                computed.index<-integer(0)
              }else{
                computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
              }
              
              
              
              p_S=dim(set.gamma[[i]])[2]
              if(length(na.omit(computed.index))==0){
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                }
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
                set.gamma.margin[[i]]<- set.gamma.margin[[i]]+apply(matrix.gamma[[i]],1,prior.dist)
              }
            }
            
            
            add.B<-list()
            add.B[[1]]<-c(set.gamma.margin[[1]],
                          set.gamma.margin[[2]],
                          set.gamma.margin[[3]])
            add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'dgCMatrix')
            for(i in 1:3){
              add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[i]])
            }
            #   add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[3]])
            
            
            #print(paste0('length(S)=',length(S)))
          }
          if(length(working.S)==0){
            
            for(i in 2){
              matrix.gamma[[i]]<-index.fun(set.gamma[[i]])
              
              if(length(C.list[[2]])<ncol(set.gamma[[i]])){
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'dgCMatrix')
                C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
                computed.index<-integer(0)
              }else{
                computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
              }
              
              
              p_S=dim(set.gamma[[i]])[2]
              if(length(na.omit(computed.index))==0){
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                }
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
                set.gamma.margin[[i]]<- set.gamma.margin[[i]]+apply(matrix.gamma[[i]],1,prior.dist)
              }
            }
            
            add.B<-list()
            add.B[[1]]<-c(set.gamma.margin[[2]])
            add.B[[2]]<-matrix.gamma[[2]]
            #print(paste0('length(S)=',length(S)))
          }
          ########## add computed model into the storage of model###############
          ##########
          
          
          ## Add B.list
          
          add.index<-match.dgCMatrix(B.list[[2]],add.B[[2]])
          if(length(which(!is.na(add.index)))>10){
            check.index<-sample(which(!is.na(add.index)),10)
            for(g in 1:4){
              #print(which(add.B[[2]][check.index[g],]==1))
              #print(which(B.list[[2]][(add.index[check.index][g]),]==1))
            }
          }
          if(length(na.omit(add.index))!=0){
            B.list[[1]]<-c((B.list[[1]]),(add.B[[1]][is.na(add.index)]))
            B.list[[2]]<-rbind(B.list[[2]],add.B[[2]][is.na(add.index),,drop=F])
          }else{
            B.list[[1]]<-c((B.list[[1]]),(add.B[[1]]))
            B.list[[2]]<-rbind(B.list[[2]],add.B[[2]])
          }
          B.list[[2]]<-B.list[[2]][order(B.list[[1]],decreasing = T),]
          B.list[[1]]<-B.list[[1]][order(B.list[[1]],decreasing = T)]
          ###################Select Next S###############
          if(length(working.S)!=0){
            set.star<-data.frame(set.index=1:3,gamma.set.index=rep(NA,3),margin=rep(NA,3))
            for(i in 1){
              aa<-set.gamma.margin[[i]]-current.log.margin
              aa<-aa-aa[which.max(aa)]
              if(length(which(is.nan(aa)))!=0){
              aa[which(is.nan(aa))]<-min(aa)
              }
              set.star$gamma.set.index[i] <-c(sample(1:length(set.gamma.margin[[i]]),1,prob=exp(aa)))
              set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
              rm(aa)
            }
            for(i in 2:length(set.gamma)){
              aa<-set.gamma.margin[[i]]-current.log.margin
              aa<-aa-aa[which.max(aa)]
              if(length(which(is.nan(aa)))!=0){
                aa[which(is.nan(aa))]<-min(aa)
              }
              set.star$gamma.set.index[i]<-c(sample((1:length(set.gamma.margin[[i]]))[order(exp(aa),decreasing = T)[1:(floor(length(set.gamma.margin[[i]])/2))]],
                                                    1,prob=exp(aa)[order(exp(aa),decreasing = T)[1:(floor(length(set.gamma.margin[[i]])/2))]]))
              set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
              rm(aa)
            }
            
            
           print(set.star)
            ##second sample
            aa<-set.star$margin-current.log.margin-max(set.star$margin-current.log.margin)
            if(length(working.S)==num.causal){
              aa<-aa-max(aa[c(1,3)])
              sec.sample<-sample(c(1,3),1,prob=exp(aa)[-2] )
              S<-set.gamma[[sec.sample]][set.star$gamma.set.index[[sec.sample]] ,]
            }else{
            sec.sample<-sample(1:3,1,prob=exp(aa) )
            S<-set.gamma[[sec.sample]][set.star$gamma.set.index[[sec.sample]] ,]
            }
            
            }else{
                set.star<-data.frame(set.index=rep(1,3),gamma.set.index=rep(NA,3),margin=rep(NA,3))
                         aa<-set.gamma.margin[[2]]-current.log.margin
                         aa<-aa-aa[which.max(aa)]
                         if(length(which(is.nan(aa)))!=0){
                           aa[which(is.nan(aa))]<-min(aa)
                         }
                         
                         set.star$gamma.set.index[2] <-c(sample((1:length(set.gamma.margin[[2]]))[order(exp(aa),decreasing = T)[1:(floor(p/2))]],
                                                             1,prob=exp(aa)[order(exp(aa),decreasing = T)[1:(floor(p/2))]]))
                         set.star$margin[2]<-set.gamma.margin[[2]][  set.star$gamma.set.index[2]]
                         
                         S<-set.gamma[[2]][set.star$gamma.set.index[2],]
                         print(set.star)
          }
          print(paste0('this is running S: ',paste0(S,collapse = ',')))
          S<-c(S,conditional.S)
          print(paste0('this is all S: ',paste0(S,collapse = ',')))
          #print(h)
          }
         if(outlier.switch){
           inner.cor.threhold=0.99999
           test.S<-set.gamma[[which.max(set.star$margin)]][set.star$gamma.set.index[[which.max(set.star$margin)]] ,]
           
            outlier.table<-data.frame(matrix(NA,nrow = 0,ncol=5))
          
          for(t in test.S){
              
              ss.index<-which(abs(Sigma[t,])>outlier.cor.range.threshold)
              ss.index<-ss.index[is.na(match(ss.index,conditional.S))]
              if(length(ss.index)>1){
              test.SS<-ss.index
               repeat{
              
                   f.ld<-as.matrix(Sigma[test.SS,test.SS])
                   off.diag.index<-which(f.ld==1)[is.na(match(which(f.ld==1),seq(1,length(test.SS)^2,by=length(test.SS)+1)))]
                   if(length(off.diag.index)!=0){
                     f.ld[off.diag.index]=inner.cor.threhold
                   }
                  test.ld<-f.ld
                 alt.test.ld<-matrix(inner.cor.threhold,length(test.SS),length(test.SS))
                 diag(alt.test.ld)<-1
                 test.z<-as.matrix(z[test.SS,])
                 test.size<-nrow(test.z)-1
                 vector1<-as.matrix(rep(1,nrow(test.z)-1))
                 result.table<-data.frame(matrix(NA,nrow=length(test.SS),ncol=5))
                 colnames(result.table)<-c('Index','Z','a','m1','B1')
                 for(s.index in 1:nrow(test.z)){
                   ss_inverse<-ginv(test.ld[-s.index,-s.index])
                   a=test.z[s.index]-test.ld[s.index,-s.index,drop=F]%*%ss_inverse%*%test.z[-s.index,,drop=F]
                   b=1-test.ld[s.index,-s.index,drop=F]%*%ss_inverse%*%vector1
                   if(b<=0){
                     b=1-alt.test.ld[s.index,-s.index,drop=F]%*%solve(alt.test.ld[-s.index,-s.index])%*%vector1
                   }
                   m1=t(test.z[-s.index,,drop=F])%*%ss_inverse%*%vector1/(test.size+1/test.size)
                   sigma1=1/(test.size+1/test.size)
                   sigma0=1-test.ld[s.index,-s.index,drop=F]%*%ss_inverse%*%test.ld[-s.index,s.index,drop=F]
                   if(sigma0<=0){
                     sigma0=1-alt.test.ld[s.index,-s.index,drop=F]%*%solve(alt.test.ld[-s.index,-s.index])%*%alt.test.ld[-s.index,s.index,drop=F]
                   }
                   if(sigma0<0){
                          sigma0=0
                        }
                   B1=sqrt(exp(1))*sqrt(
                     (a-m1*b)^2/
                       (b^2*sigma1+sigma0)
                   )*exp(
                     -abs(a-m1*b)/(2*sqrt(b^2*sigma1+sigma0))
                   )
                   result.table[s.index,]<-c(test.SS[s.index],test.z[s.index,],a,m1,B1)
                 }
                 print(result.table[which.min(result.table$B1),])
                 if((!any(result.table$B1<outlier.threshold))){
                   break
                 }
                 if(nrow(result.table)==2){
                   outlier.table<-rbind(outlier.table,result.table[which.min(abs(result.table$Z)),])
                   break
                 }else{
                         equal.index<-which(result.table$B1<outlier.threshold)
                         drop.index<-equal.index[which.min(abs(result.table$Z[equal.index]))]
                         print(result.table[  drop.index,])
                         outlier.table<-rbind(outlier.table,result.table[drop.index,])
                         test.SS<-test.SS[-drop.index]
                 
                 }
                 }
               
             }
           }
            conditional.S<-c(conditional.S,outlier.table$Index)
            if(length(outlier.table$Index)!=0){
            conditional.S.list[[length(conditional.S.list)+1]]<-outlier.table
            }
            conditional.S<-unique(conditional.S)
            S<-unique(c(S,conditional.S))
            
        
          }
          
        
        }
      ######Result
    
      
      result.B.list<-list()
      if(!is.null(conditional.S)){
      all.c.index<-c()
      for(t in 1:length(conditional.S.list)){
        for(tt in 1:length(conditional.S.list[[t]]$Index)){
        c.index<-(B.list[[2]]@i[(B.list[[2]]@p[conditional.S.list[[t]]$Index[tt]]+1):B.list[[2]]@p[conditional.S.list[[t]]$Index[tt]+1]])+1
       all.c.index<-c(all.c.index,c.index)
        }
        }
      all.c.index<-unique(all.c.index)
      temp.B.list<-list()
      temp.B.list[[1]]<-B.list[[1]][-all.c.index]
      temp.B.list[[2]]<-B.list[[2]][-all.c.index,]
      }else{
        temp.B.list<-list()
        temp.B.list[[1]]<-B.list[[1]]
        temp.B.list[[2]]<-B.list[[2]]
      }
      result.B.list<-list()
      result.B.list[[1]]<-temp.B.list[[1]][(1:min(B,nrow(temp.B.list[[2]])))]
      result.B.list[[2]]<-temp.B.list[[2]][(1:min(B,nrow(temp.B.list[[2]]))),]
      
      result.prob<-PIP.func(result.B.list[[1]],result.B.list[[2]])
      if(!is.null(output.labels)){
        if(dir.exists(output.labels)==F ){
          dir.create(output.labels,recursive = T)
        }
      write.table(result.B.list[[1]],file=paste0(output.labels,'/post_',label,'_poi_likeli','.txt'),row.names = F,col.names = F)
      writeMM(result.B.list[[2]],file=paste0(output.labels,'/post_',label,'_poi_gamma','.mtx'))
      write.table((result.prob),file=paste0(output.labels,'/post_', label,'.txt'),row.names = F,append = F,col.names = F)
      if(outlier.switch){
        saveRDS(conditional.S.list,file=paste0(output.labels,'/post_', label,'_','outliers.RData'))
      }
      }

      
      difference<-abs(mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])-stored.bf)
       print(difference)
       if(difference<epsilon){
         break
       }else{
         stored.bf<-mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])
       }
    }
    return(list(result.B.list,C.list,result.prob,conditional.S.list))
  }
  ######## All burning###########
  previous.result<-list()
  
  for(i in 1:L){
    t0=Sys.time()
   
   all.C.list[[i]]<-Module.Cauchy.Shotgun(z.list[[i]],ld.list[[i]],epsilon=epsilon.list[[i]],
                                           Max.Model.Dim=Max.Model.Dim,lambda = lambda.list[[i]],
                                          outlier.switch=outlier.switch,outlier.threshold=outlier.threshold,outlier.cor.range.threshold=outlier.cor.range.threshold,
                                           num.causal = num.causal,y.var=y.var,
                                           label = label.list[[i]],output.labels = output.labels,
                                           effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)

   
   t1=Sys.time()-t0
   print(paste0('This is locus ',i,' burning time'))
   print((t1))
  }
  for(g in 1:all.iter){ 
    if(outlier.switch){
      delete.list<-list()
      for(i in 1:L){
        delete.list[[i]]<-integer(0)
        if(length(all.C.list[[i]][[4]])!=0){
        temp.delete.list<-c()
        for(t in 1:length(all.C.list[[i]][[4]])){
          temp.delete.list<-c(temp.delete.list,all.C.list[[i]][[4]][[t]]$Index) 
        }
        delete.list[[i]]<-temp.delete.list
        }
      }
    }else{
        delete.list<-list()
        for(i in 1:L){
          delete.list[[i]]<-all.C.list[[i]][[4]]
        }
    }
    if(!is.null(w.list)){
      w<-matrix(NA,nrow=0,ncol=ncol(w.list[[1]]))
      colnames(w)<-colnames(w.list[[1]])
      for(i in 1:L){
        if(length(delete.list[[i]])!=0){
          w<-rbind(w,w.list[[i]][-delete.list[[i]],])
        }else{
          w<-rbind(w,w.list[[i]])
        }
      }
    }
  
    for(i in 1:L){
      previous.result[[i]]<-mean(all.C.list[[i]][[1]][[1]][1:round(quantile(1:length(all.C.list[[i]][[1]][[1]]),probs = 0.25))])
    }
    ################################
    if(!is.null(w.list)){
    if(!standardize.model.space){
      model.space.count<-c()
      for(i in 1:L){
        if(length(delete.list[[i]])!=0){
        model.space.count<-c(model.space.count,colSums(all.C.list[[i]][[1]][[2]][,-delete.list[[i]]]))
        }else{
          model.space.count<-c(model.space.count,colSums(all.C.list[[i]][[1]][[2]]))
        }
        
      }
    }else{
      model.space.count<-c()
      for(i in 1:L){
        if(length(delete.list[[i]])!=0){
        indi.count<-colSums(all.C.list[[i]][[1]][[2]][,-delete.list[[i]]])
        indi.count<-floor(colSums(all.C.list[[i]][[1]][[2]][,-delete.list[[i]]])/nrow(all.C.list[[i]][[1]][[2]][,-delete.list[[i]]])*Max.Model.Dim)
        }else{
          indi.count<-colSums(all.C.list[[i]][[1]][[2]])
          indi.count<-floor(colSums(all.C.list[[i]][[1]][[2]])/nrow(all.C.list[[i]][[1]][[2]])*Max.Model.Dim)
        }
        model.space.count<-c(model.space.count,indi.count)
      }
    }
    
    try.index<-try(glm.beta<-EM.M.step.func(Model.space = model.space.count ,w=w,input.alpha=input.alpha))
    prior.prob.list<-list()
    
    if(class(try.index)[1]!='try-error'){
      for(i in 1:L){
        glm.beta[1]=log((min(Max.Model.Dim,nrow(all.C.list[[i]][[1]][[2]])))*lambda.list[[i]]/(lambda.list[[i]]+p.list[[i]]))
        prior.prob.list[[i]]<-(exp(w.list[[i]]%*%glm.beta)/(max(1+max(exp(w.list[[i]]%*%glm.beta)),min(Max.Model.Dim,nrow(all.C.list[[i]][[1]][[2]])))))
        if(!is.null(output.labels)){
        write.table((glm.beta),file=paste0(output.labels,'/post_', label.list[[i]],'_theta.txt'),row.names = F,append = F,col.names = F)
        }
        plot(prior.prob.list[[i]])
      }
      model.prior='input.prob'  
    }else{
      model.prior='Poisson'
      prior.prob.list<-list()
      for(i in 1:L){
        prior.prob.list[[i]]<-list(NULL)
      }
      
    }
    
    }else{
      prior.prob.list<-list()
      for(i in 1:L){
        prior.prob.list[[i]]<-list(NULL)
      }
    }
    ########    
    for(i in 1:L){
      t0=Sys.time()
      all.C.list[[i]]<-Module.Cauchy.Shotgun(z=z.list[[i]],ld.list[[i]],input.conditional.S.list = all.C.list[[i]][[4]],
                                             Max.Model.Dim=Max.Model.Dim,y.var=y.var,num.causal = num.causal,epsilon=epsilon.list[[i]],
                                             C.list = all.C.list[[i]][[2]],prior.prob = prior.prob.list[[i]],
                                             outlier.switch=outlier.switch,outlier.threshold=outlier.threshold,outlier.cor.range.threshold=outlier.cor.range.threshold,
                                             lambda = lambda.list[[i]],
                                             label = label.list[[i]],output.labels = output.labels,
                                             effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
      t1=Sys.time()-t0
      print(paste0('This is locus ',i,' computing time'))
      print((t1))
    }
   #  plot(previous.result[[i]])
   #  plot(all.C.list[[i]][[3]])
    difference<-0
     for(i in 1:L){
     #  difference<-difference+sum(abs(previous.result[[i]]-all.C.list[[i]][[3]]))
       difference<-difference+abs(previous.result[[i]]-mean(all.C.list[[i]][[1]][[1]][1:round(quantile(1:length(all.C.list[[i]][[1]][[1]]),probs = 0.25))]))
     }
    print(paste0('This is difference; ',difference))
     if(difference<all.epsilon.threshold){
       break
     }
  }
  results.list<-list()
  for(i in 1:L){
    results.list[[i]]<-list()
    pip=all.C.list[[i]][[3]]
    credible.set<-credible.set.fun.improved(pip,ld.list[[i]],rho=rho.index)
    credible.model<-credible.model.fun(all.C.list[[i]][[1]][[1]],all.C.list[[i]][[1]][[2]],bayes.threshold = BF.index)
    results.list[[i]][[1]]<-pip
    results.list[[i]][[2]]<-credible.set
    results.list[[i]][[3]]<-credible.model
    if(outlier.switch ){
      results.list[[i]][[4]]<-all.C.list[[i]][[4]]
      names(results.list[[i]])<-c('PIPs','Credible set','Credible model','Outliers')
    }else{
      names(results.list[[i]])<-c('PIPs','Credible set','Credible model')
    }
    }
  return(results.list)
}
