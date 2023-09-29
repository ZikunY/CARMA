#' CARMA
#' 
#' Performs a Bayesian fine-mapping model in order to identify putative causal variants at GWAS loci. The model requires the summary statistics
#' for the SNPs at the testing loci, the corresponding LD matrices for fine-mapping, and an estimated trait variance. Functional annotations can be included as prior 
#' information on the causality of the testing SNPs. The model also provides a procedure of outlier detection, which aims to identify discrepancies
#' between summary statistics and LD matrix extracted from a reference panel. The model can be executed chromosome-wise to increase power. 
#'@param z.list Input list of summary statistics at the testing loci; each element of the list is the summary statistics at each individual locus.
#'@param ld.list Input list of LD correlation matrix at the testing loci; each element of the list is the LD matrix at each individual locus.
#'@param w.list Input list of the functional annotations at the testing loci; each element of the list is the functional annotation matrix at each individual locus.
#'@param lambda.list Input list of the hyper-parameter \eqn{\eta} at the testing loci; each element of the list is the hyper-parameter of each individual locus.
#'@param label.list Input list of the names at the testing loci. Default is NULL. 
#'@param effect.size.prior The prior of the effect size. The choice are 'Cauchy' and 'Spike-slab' priors, where the 'Spike-slab' prior is the default prior.
#'@param input.alpha The elastic net mixing parameter, where \eqn{0\le}\eqn{\alpha}\eqn{\le 1}.
#'@param y.var The input variance of the summary statistics, the default value is 1 as the summary statistics are standardized. 
#'@param rho.index A number between 0 and 1 specifying \eqn{\rho} for the estimated credible sets.
#'@param BF.index  The threshold of the Bayes factor of the estimated credible models. The default setting is 10.
#'@param outlier.BF.index  The Bayes threshold for the Bayesian hypothesis test for outlier detection. The default setting is 1/3.2. 
#'@param outlier.switch The indicator variable for outlier detection. We suggest that the detection should always be turned on if using external LD matrix. 
#'@param num.causal The maximum number of causal variants assumed per locus, which is 10 causal SNPs per locus by default.
#'@param Max.Model.Dim Maximum number of the top candidate models based on the unnormalized posterior probability. 
#'@param all.inner.iter Maximum iterations for Shotgun algorithm to run per iteration within EM algorithm.
#'@param all.iter Maximum iterations for EM algorithm to run.
#'@param tau The prior precision parameter of the effect size. The default value is computed based on the scenario of n=10000 and 1% heritability.
#'@param output.labels Output directory where output will be written while CARMA is running. Default is the OS root directory ".".
#'@param epsilon.threshold Convergence threshold measured by average of Bayes factors.
#'@param printing.log Whether print the running log while running CARMA.
#'@param EM.dist The distribution used to model the prior probability of being causal as a function of functional annotations. The default distribution is logistic distribution. 
#'@return The return is a list, for each list:
#'\itemize{
#'\item pip - The posterior inclusion probability of each individual locus.
#'\item Credibleset - The information on the credible set given a threshold \eqn{\rho}.
#'\item Credible model - The information on the credible model given a threshold  of the Bayes factor.
#'\item Outliers - The information regarding the detected outliers.
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
#'lambda.list<-list()
#'lambda.list[[1]]<-1/sqrt(p)
#'CARMA.result<-CARMA(z.list,ld.list=ld.list,
#'lambda.list = lambda.list,effect.size.prior='Hyper-g')
CARMA<-function(z.list,ld.list,w.list=NULL,lambda.list=NULL,output.labels='.',label.list=NULL,
                                 effect.size.prior='Spike-slab',rho.index=0.99,BF.index=10,EM.dist='Logistic',
                                Max.Model.Dim=2e+5,all.iter=3,all.inner.iter=10,input.alpha=0,epsilon.threshold=1e-5,printing.log=F,
                                 num.causal=10,y.var=1,tau=0.04,outlier.switch=T,outlier.BF.index=1/3.2,prior.prob.computation='Logistic'){
                                     

  Sys.setenv("PKG_CXXFLAGS"="-std=c++11")    #compile functions that use C++11 in R
  ##########  Feature learning step for the CARMA algorithm, such as learning the total number of input loci, the total number of variants at each locus, etc.#########
  ##########  Additionally, CARMA defines the lists of the model spaces and the likelihood of all input loci###########
  {
    log.2pi<-log(2*pi)
    L<-length(z.list)
    p.list<-list()
    for(i in 1:L){
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
        invariant.var.index<-which((apply(w.list[[i]][,-1],2,sd))==0)
        if(length(invariant.var.index)!=0){
          invariant.var<-w.list[[i]][,invariant.var.index+1]
          w.list[[i]]<-as.matrix(cbind(1,scale(w.list[[i]][,-1])))
          w.list[[i]][,invariant.var.index+1]<-invariant.var
        }else{
          w.list[[i]]<-as.matrix(cbind(1,scale(w.list[[i]][,-1])))
        }
        
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
    

  }
  ###The M-step of the EM algorithm for incorporating functional annotations
  EM.M.step.func<-function(input.response=NULL,w=w,input.alpha=0.5,EM.dist='Logistic'){
    if(EM.dist=='Poisson'){
      count.index<-input.response
      cv.poisson<-cv.glmnet(w,count.index,family = 'poisson',alpha=input.alpha,type.measure='deviance' )
      cv.index<-which(cv.poisson$lambda==cv.poisson$lambda.min)
      glm.beta<-as.matrix(c(cv.poisson$glmnet.fit$a0[cv.index],cv.poisson$glmnet.fit$beta[-1,cv.index]))
    }
    if(EM.dist=='Logistic'){
      response.matrix<-matrix(c(1-input.response, input.response),length(input.response),2)
      cv.logistic<-cv.glmnet(x=w,y=response.matrix, family='binomial',alpha=input.alpha,type.measure = 'deviance')
      cv.index<-which(cv.logistic$lambda==cv.logistic$lambda.min)
      glm.beta<-as.matrix(c(cv.logistic$glmnet.fit$a0[cv.index],cv.logistic$glmnet.fit$beta[-1,cv.index]))
      
    }
    
    return(glm.beta=glm.beta)
  }
  ###The computation of the credible set based on the final results of the fine-mapping step.
  credible.set.fun.improved<-function(pip,ld,true.beta=NULL,rho=0.99){
    candidate.r<-c((seq(from=.5,0.95,by=0.05)),seq(0.96,0.99,0.01))
    
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
        
        cor.group<-as.numeric(names(which(abs(working.ld[which(pip.order[s]==as.numeric(colnames(working.ld))),])>cor.threshold)))
        if(sum(pip[cor.group])>rho){
          group.pip<- pip[cor.group]
          snp.index<- cor.group[order(group.pip,decreasing = T)][1: min(which(cumsum( sort(group.pip,decreasing = T))>rho))]
          snp.list[[r]][[group.index]]<-snp.index
          group.index<-group.index+1
          pip.order<-pip.order[-match(snp.index,pip.order)]
          working.ld<-working.ld[-match(snp.index,as.numeric(colnames(working.ld))),-match(snp.index,as.numeric(colnames(working.ld)))]
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
  ###The computation of the credible models based on the final results of the fine-mapping step.
  credible.model.fun<-function(likelihood,model.space,bayes.threshold=10){
    na.index<-which(is.na(likelihood))
    if(length(na.index)!=0){
      likelihood<-likelihood[-na.index]
      model.space<-model.space[-na.index,]
    }
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

  
  
#########The module function of the CARMA fine-mapping step for each locus included in the analysis##########
  
  Module.Cauchy.Shotgun<-function(z,ld.matrix,Max.Model.Dim=1e+4,input.S=NULL,lambda,label,printing.log=F,
                                              num.causal=10,output.labels,y.var=1,effect.size.prior=effect.size.prior,model.prior=model.prior,
                                              outlier.switch,input.conditional.S.list=NULL,tau=1/0.05^2,
                                              C.list=NULL,prior.prob=NULL,epsilon=1e-3,inner.all.iter=10){
       {
     #######The prior distributions on the model space#########
         prob.list<-list()
         p<-nrow(z);
         if(model.prior=='input.prob'){
            posi.log.pro<-log(prior.prob)
            nega.log.pro<-log(1-prior.prob) 
            input.prior.dist<-function(x){
              variable.index<-which(x==1)
              if(any(variable.index)){
                return( sum(posi.log.pro[variable.index])+sum(nega.log.pro[(1:p)[-variable.index]])-sum(nega.log.pro))
              }else{
                return(sum(nega.log.pro))
              }
            }
            prior.dist<-input.prior.dist
          }
       
      
       if(model.prior=='Poisson'){
        Poisson.prior.dist<-function(t){
          dim.model<-sum(t)
          result<-dim.model*log(lambda)+lfactorial(p-dim.model)-lfactorial(p)
          return(result)
        }
        prior.dist<-Poisson.prior.dist
       }
      if(model.prior=='beta-binomial'){
        beta.binomial.dist<-function(t){
          dim.model<-sum(t)
          result<- lbeta(dim.model+1,p-dim.model+9)-lbeta(1,p+9)
          return(result)
        }
        prior.dist<-beta.binomial.dist
      }
      
      ###### The marginal likelihood defined by the prior distribution of the effect size#########
      
         if(effect.size.prior=='Cauchy'){
           marginal_likelihood=Cauchy_fixed_sigma_marginal
           tau.sample<-rgamma(1e+5,0.5,rate=0.5)
           if(outlier.switch){
             outlier_likelihood=outlier_Cauchy_fixed_sigma_marginal
             outlier.tau=tau.sample
           }
         }
         if(effect.size.prior=='Hyper-g'){
           marginal_likelihood=hyper_g_fixed_sigma_marginal
           tau.sample<-rgamma(1e+5,0.5,rate=0.5)
         }
         if(effect.size.prior=='Normal'){
           marginal_likelihood=Normal_fixed_sigma_marginal
           tau.sample<-tau
           if(outlier.switch){
             outlier_likelihood=outlier_Normal_fixed_sigma_marginal
             outlier.tau=tau.sample
           }
         }
         if(effect.size.prior=='Spike-slab'){
           marginal_likelihood=ind_Normal_fixed_sigma_marginal
           tau.sample<-tau
           if(outlier.switch){
         
               outlier_likelihood=outlier_ind_Normal_marginal
               outlier.tau=tau.sample
             
           }
         }
   
    ########Feature learning for the fine-mapping step, such as learning the visited model space from the previous iterations#######

         p<-nrow(z);
      log.2pi<-log(2*pi)

      B<-Max.Model.Dim
      stored.result.prob<-rep(0,p)
      stored.bf<-0
      Sigma<-as.matrix(ld.matrix)
      
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
           conditional.S<-input.conditional.S.list$Index
           conditional.S<-unique(conditional.S)
           S<-conditional.S
         }
    }
    
   
    ###The function that defines neighborhood model space
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
    ####Function that computes posterior inclusion probability based on the marginal likelihood and model space
    PIP.func<-function(likeli,model.space){
      infi.index<-which(is.infinite(likeli))
      if(length(infi.index)!=0){
        likeli<-likeli[-infi.index]
        model.space<-model.space[-infi.index,]
      }
      na.index<-which(is.na(likeli))
      if(length(na.index)!=0){
        likeli<-likeli[-na.index]
        model.space<-model.space[-na.index,]
      }
      aa<-likeli-max(likeli,na.rm=T)
      prob.sum<-sum(exp(aa))
      result.prob<-rep(NA,p)
      for(i in 1:p){
        result.prob[i]<-sum(exp(aa[ which(model.space[,i]==1)]))/prob.sum
      }
      return(result.prob)
    }
    index.fun.inner<-function(x){
      m<-as(as(as(matrix(0,nrow=nrow(x),ncol=p), "dMatrix"), "generalMatrix"), "TsparseMatrix")
      m@i<-as.integer(rep(1:nrow(x)-1,each=ncol(x)))
      m@j<-as.integer(c(t(x))-1)
      m@x=rep(1,nrow(x)*ncol(x))
      m<-as(m,"CsparseMatrix")
      return(m)
    }
    index.fun<-function(outer.x,Max.Model.Dimins=10){
     if(nrow(outer.x)>1000){
    index.bins<-which((1:nrow(outer.x))%%floor(nrow(outer.x)/Max.Model.Dimins)==0)
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
    ridge.fun<-function(x){
      temp.ld.S<-x*modi.ld.S+(1-x)*diag(nrow(modi.ld.S))
      temp.Sigma[test.S,test.S]<-temp.ld.S
      return( outlier_likelihood(test.S,temp.Sigma,z,outlier.tau,length(test.S),1) )
    }
    
    ######################################################
    for(l in 1:inner.all.iter){
      for(h in 1:10){
        ##############Shotgun COMPUTATION ############
        { 
          set.gamma<-set.gamma.func(S,conditional.S)  
          if(is.null(conditional.S)){
            working.S=S
            base.model<-null.model
            base.model.margin<-null.margin
          }else{
            working.S=S[-match(conditional.S,S)]
            if(length(working.S)!=0){
            base.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'CsparseMatrix')
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
            S.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'CsparseMatrix')
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
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
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
            }
            
            add.B<-list()
            add.B[[1]]<-c(set.gamma.margin[[1]],
                          set.gamma.margin[[2]],
                          set.gamma.margin[[3]])
            add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
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
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
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
            add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
            for(i in 1:3){
              add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[i]])
            }
           }
          if(length(working.S)==0){
            
            for(i in 2){
              matrix.gamma[[i]]<-index.fun(set.gamma[[i]])
              
              if(length(C.list[[2]])<ncol(set.gamma[[i]])){
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
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
                set.gamma.margin[[i]]<- set.gamma.margin[[i]]+ apply(matrix.gamma[[i]],1,prior.dist)
              }
            }
            
            add.B<-list()
            add.B[[1]]<-c(set.gamma.margin[[2]])
            add.B[[2]]<-matrix.gamma[[2]]
          }
          ########## add visited models into the storage space of models###############
          
          
          add.index<-match.dgCMatrix(B.list[[2]],add.B[[2]])
          if(length(which(!is.na(add.index)))>10){
            check.index<-sample(which(!is.na(add.index)),10)
 
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
          ###################Select next visiting model###############
          
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
            
        #######The Bayesian hypothesis testing for Z-scores/LD discrepancies########
            if(outlier.switch){
            for(i in 2:length(set.gamma)){
              repeat{
                
              aa<-set.gamma.margin[[i]]-current.log.margin
              aa<-aa-aa[which.max(aa)]
              if(length(which(is.nan(aa)))!=0){
                aa[which(is.nan(aa))]<-min(aa[!is.na(aa)])
              }
            
              set.star$gamma.set.index[i]<-c(sample((1:length(set.gamma.margin[[i]])),
                                                    1,prob=exp(aa)))
              set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
     
              test.S<-set.gamma[[i]][set.star$gamma.set.index[i],]
              
              modi.Sigma<-Sigma
              temp.Sigma<-Sigma
              if(length(test.S)>1){
              
              modi.ld.S<- modi.Sigma[test.S,test.S]
              
                opizer<-optimize(ridge.fun,interval=c(0,1),maximum = T)
                modi.ld.S<-opizer$maximum*modi.ld.S+(1-opizer$maximum)*diag(nrow(modi.ld.S)) 
             
              
              modi.Sigma[test.S,test.S]<-modi.ld.S
        
              test.log.BF<-outlier_likelihood(test.S,Sigma,z,outlier.tau,length(test.S),1)-outlier_likelihood(test.S,modi.Sigma,z,outlier.tau,length(test.S),1)
              test.log.BF<--abs(test.log.BF)
              if(printing.log==T){
              print(paste0('Outlier BF: ', test.log.BF))
              print(test.S)
              print(paste0('This is xi hat: ', opizer))
              }
              }
                
              if(exp(test.log.BF)<outlier.BF.index){
                set.gamma[[i]]<-set.gamma[[i]][-set.star$gamma.set.index[i],]
                set.gamma.margin[[i]]<-set.gamma.margin[[i]][-set.star$gamma.set.index[i]]
                conditional.S<-c(conditional.S,test.S[is.na(match(test.S,working.S))])
                conditional.S<-unique(conditional.S)
              }else{
                break
              }
              }
              rm(aa)
            }
            }else{
              for(i in 2:length(set.gamma)){
                 
                  aa<-set.gamma.margin[[i]]-current.log.margin
                  aa<-aa-aa[which.max(aa)]
                  if(length(which(is.nan(aa)))!=0){
                    aa[which(is.nan(aa))]<-min(aa[!is.na(aa)])
                  }
                  
                  set.star$gamma.set.index[i]<-c(sample((1:length(set.gamma.margin[[i]])),
                                                        1,prob=exp(aa)))
                  set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
                rm(aa)
              }
            }
            if(printing.log==T){
           print(set.star)
          }
            if(length(working.S)==num.causal){
              set.star<-set.star[-2,]
              aa<-set.star$margin-current.log.margin-max(set.star$margin-current.log.margin)
              sec.sample<-sample(c(1,3),1,prob=exp(aa))
              S<-set.gamma[[sec.sample]][set.star$gamma.set.index[[which(sec.sample==set.star$set.index)]] ,]
            }else{
              aa<-set.star$margin-current.log.margin-max(set.star$margin-current.log.margin)
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
                         
                         set.star$gamma.set.index[2] <-c(sample((1:length(set.gamma.margin[[2]]))[order(exp(aa),decreasing = T)[1:(min(length(aa),floor(p/2)))]],
                                                             1,prob=exp(aa)[order(exp(aa),decreasing = T)[1:(min(length(aa),floor(p/2)))]]))
                         set.star$margin[2]<-set.gamma.margin[[2]][  set.star$gamma.set.index[2]]
                         
                         S<-set.gamma[[2]][set.star$gamma.set.index[2],]
                         if(printing.log==T){
                           print(set.star)
                         }
            }
          if(printing.log==T){
          print(paste0('this is running S: ',paste0(S,collapse = ',')))
          }
          S<-unique(c(S,conditional.S))
        }
        
     
        }
      ######Output of the results of the module function######
    
      
      result.B.list<-list()
      if(!is.null(conditional.S)){
      all.c.index<-c()
 
     
        for(tt in conditional.S){
          c.index<-(B.list[[2]]@i[min(length(B.list[[2]]@i),(B.list[[2]]@p[tt]+1)):B.list[[2]]@p[tt+1]])+1
          all.c.index<-c(all.c.index,c.index)
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
      
      if(num.causal==1){
        single.set<-matrix(1:p,p,1)
        single.marginal<-apply(single.set,1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
        aa<-single.marginal-max(single.marginal,na.rm=T)
        prob.sum<-sum(exp(aa))
        result.prob<-(exp(aa))/prob.sum
      }else{
        result.prob<-PIP.func(result.B.list[[1]],result.B.list[[2]])
      }
      conditional.S.list<-data.frame(Index=conditional.S,Z=z[conditional.S,])
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
     #  print(difference)
       if(difference<epsilon){
         break
       }else{
         stored.bf<-mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])
       }
    }
    return(list(result.B.list,C.list,result.prob,conditional.S.list,prob.list))
  }
  ######## Burning step###########
  previous.result<-list()
  ########Run fine-mapping step (module function) for each locus included in the analysis
  for(i in 1:L){
    t0=Sys.time()
   all.C.list[[i]]<-Module.Cauchy.Shotgun(z.list[[i]],ld.list[[i]],epsilon=epsilon.list[[i]],
                                           Max.Model.Dim=Max.Model.Dim,lambda = lambda.list[[i]],
                                          outlier.switch=outlier.switch,tau=tau,
                                           num.causal = num.causal,y.var=y.var,printing.log=printing.log,
                                           label = label.list[[i]],output.labels = output.labels,
                                           effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
   t1=Sys.time()-t0
   print(paste0('This is locus ',i,' burning time'))
   print((t1))
  
  }
  ########Running CARMA######## 
  for(g in 1:all.iter){ 
    if(outlier.switch){
      delete.list<-list()
      for(i in 1:L){
        delete.list[[i]]<-integer(0)
        if(nrow(all.C.list[[i]][[4]])!=0){
        temp.delete.list<-c(all.C.list[[i]][[4]]$Index)
        delete.list[[i]]<-temp.delete.list
        }
      }
    }else{
      delete.list<-list()
      for(i in 1:L){
      delete.list[[i]]<-integer(0)
      }
      }
 
    #########If the list of annotations is non-empty, then the PIPs and functional annotations at all loci are aggregated for the M-step of the EM algorithm
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
    if(!is.null(w.list)){
      
    if(EM.dist=='Poisson'){
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
    M.step.response=model.space.count
   
    }
   ######The M step of the EM algorithm
    if(EM.dist=='Logistic'){
      M.step.response<-c()
      for(i in 1:L){
        if(length(delete.list[[i]])!=0){
          indi.pip<-(all.C.list[[i]][[3]][-delete.list[[i]]])
         }else{
          indi.pip<-(all.C.list[[i]][[3]])
        }
        M.step.response<-c(M.step.response,indi.pip)
      }
    }

    try.index<-try(glm.beta<-EM.M.step.func(input.response  = M.step.response ,w=w,input.alpha=input.alpha,EM.dist=EM.dist))
    prior.prob.list<-list()
    
    if(class(try.index)[1]!='try-error'){
      for(i in 1:L){
        if(prior.prob.computation=='Intercept.approx'){
          glm.beta[1]=log((min(Max.Model.Dim,nrow(all.C.list[[i]][[1]][[2]])))*lambda.list[[i]]/(lambda.list[[i]]+p.list[[i]]))
          prior.prob.list[[i]]<-(exp(w.list[[i]]%*%glm.beta)/(max(1+max(exp(w.list[[i]]%*%glm.beta)),min(Max.Model.Dim,nrow(all.C.list[[i]][[1]][[2]])))))
      #    print(sort(prior.prob.list[[i]],decreasing = T)[1:10])
        }
        
        if(prior.prob.computation=='Logistic'){
          prior.prob.list[[i]]<-plogis(w.list[[i]]%*%glm.beta)
      #    print(sort(prior.prob.list[[i]],decreasing = T)[1:10])
        }
        if(!is.null(output.labels)){
        write.table((glm.beta),file=paste0(output.labels,'/post_', label.list[[i]],'_theta.txt'),row.names = F,append = F,col.names = F)
        }
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
   #######Fine-mapping step for each locus, i.e., the E-step in the EM algorithm
    for(i in 1:L){
      t0=Sys.time()
      all.C.list[[i]]<-Module.Cauchy.Shotgun(z=z.list[[i]],ld.list[[i]],input.conditional.S.list = all.C.list[[i]][[4]],
                                             Max.Model.Dim=Max.Model.Dim,y.var=y.var,num.causal = num.causal,epsilon=epsilon.list[[i]],
                                             C.list = all.C.list[[i]][[2]],prior.prob = prior.prob.list[[i]],
                                             outlier.switch=outlier.switch,tau=tau,printing.log=printing.log,
                                             lambda = lambda.list[[i]], label = label.list[[i]],output.labels = output.labels,
                                             effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
      t1=Sys.time()-t0
      print(paste0('This is locus ',i,' computing time'))
      print((t1))
    }

    difference<-0
     for(i in 1:L){
       difference<-difference+abs(previous.result[[i]]-mean(all.C.list[[i]][[1]][[1]][1:round(quantile(1:length(all.C.list[[i]][[1]][[1]]),probs = 0.25))]))
     }
    if(printing.log==T){
      print(paste0('This is difference; ',difference))
    }
   
     if(difference<all.epsilon.threshold){
       break
     }
  }
  ### Output of the results of CARMA
  results.list<-list()
  for(i in 1:L){
    results.list[[i]]<-list()
    pip=all.C.list[[i]][[3]]
    credible.set<-credible.set.fun.improved(pip,ld.list[[i]],rho=rho.index)
    credible.model<-credible.model.fun(all.C.list[[i]][[1]][[1]],all.C.list[[i]][[1]][[2]],bayes.threshold = BF.index)
    results.list[[i]][[1]]<-pip
    results.list[[i]][[2]]<-credible.set
    results.list[[i]][[3]]<-credible.model
    results.list[[i]][[4]]<-all.C.list[[i]][[4]]
    names(results.list[[i]])<-c('PIPs','Credible set','Credible model','Outliers')
    
    }
  return(results.list)
}

