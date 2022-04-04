#' CAusal Robust Mapping Method with Annotations(CARMA) 
#' 
#' Performs a Bayesian fine-mapping model in order to identify putative causal variants at GWAS loci. This function requires the summary statistics
#' of the SNPs in the testing loci, and the corresponding simple-LD matrices for fine-mapping. Functional annotations can be included as the prior 
#' information of the causality of the testing SNPs. The model can be executed chromosome-wise to increase power. 
#'@param z.list Input list of the summary statistics of the testing loci, and each element of the list is the summary statistics of each individual locus.
#'@param ld.list Input list of the LD correlation matrix of the testing loci, and each element of the list is the LD matrix of each individual locus.
#'@param w.list Input list of the functional annotations of the testing loci, and each element of the list is the functional annotation matrix of each individual locus. 
#'@param lambda.list Input list of the hyper-parameter \eqn{\eta} of the testing loci, and each element of the list is the hyper-parameter of each individual locus.
#'@param label.list Input list of the names of the testing loci. Default is NULL. 
#'@param effect.size.prior The prior of the effect size. The choice are 'Cauchy' and 'Hyper-g' priors, where the Cauchy prior is the default prior.
#'@param input.alpha The elastic net mixing parameter, where \eqn{0\le}\eqn{\alpha}\eqn{\le 1}.
#'@param rho.index A number between 0 and 1 specifying \eqn{\rho} of the estimated credible sets.
#'@param BF.index  A number greater than 1 to specifying the threshold of the Bayes factor of the estimated credible models.
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
#'CARMA.result<-CARMA(z.list,ld.list=ld.list,lambda.list = lambda.list,effect.size.prior='Hyper-g')


CARMA<-function(z.list,ld.list,w.list=NULL,lambda.list,
                                 effect.size.prior='Cauchy',rho.index=0.99,BF.index=10,
                                 Max.Model.Dim=1e+4,all.iter=10,all.inner.iter=10,label.list=NULL,
                                output.labels=NULL,input.alpha=0.5,epsilon.threshold=1e-3){
                                      
  Sys.setenv("PKG_CXXFLAGS"="-std=c++11")    
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
  
  
  ########## Input data#########
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
      w.list[[i]]<-as.matrix(cbind(1,scale(w.list[[i]][,-1])))
      }
    }
    Sigma.list<-list()
    Sigma.inv.list<-list()
    zSigmaz.list<-list()
    for(i in 1:L){
    Sigma.list[[i]]<-as.matrix(ld.list[[i]])
    Sigma.inv.list[[i]]<-ginv(Sigma.list[[i]])
    zSigmaz.list[[i]]<-as.numeric(t(z.list[[i]])%*%Sigma.inv.list[[i]]%*%z.list[[i]])
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
    credible.set.list<-list()
    credible.model.list<-list()
    if(is.null(label.list)){
      for(i in 1:L){
        label.list[[i]]=paste0('locus_',i)
      }
    }
        all.epsilon.threshold<-0
    epsilon.list<-list()
    for(i in 1:L){
      epsilon.list[[i]]<-epsilon.threshold*p.list[[i]]
      all.epsilon.threshold<-all.epsilon.threshold+  epsilon.threshold*p.list[[i]]
    }
    model.prior='Poisson'
    standardize.model.space=T
    
   ####### each module CS compute their own null margin.
  }
  
#########Module model##########
  
  Module.Cauchy.Shotgun<-function(z,ld.matrix,Sigma.inv.input,zSigmaz.input,Max.Model.Dim=1e+4,input.S=NULL,lambda,label,
                                              L=NULL,output.labels,effect.size.prior=effect.size.prior,model.prior=model.prior,
                                              C.list=NULL,prior.prob=NULL,epsilon=1e-3,inner.all.iter=10){ 
    {
      ######Define prior distribution ################
  
       if(model.prior=='input.prob'){
         posi.log.pro<-log(prior.prob)
         nega.log.pro<-log(1-prior.prob) 
         input.prior.dist<-function(x){
           variable.index<-which(x==1)
           if(any(variable.index)){
             return( sum(posi.log.pro[variable.index])+sum(nega.log.pro[(1:p)[-variable.index]]))
           }else{
             return(sum(nega.log.pro))
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
        marginal_likelihood=Cauchy_marginal
      }
      if(effect.size.prior=='Hyper-g'){
        marginal_likelihood=hyper_g_marginal
      }
      
      p<-nrow(z);
      log.2pi<-log(2*pi)
      tau.sample<-rgamma(1e+5,0.5,rate=0.5)
      B<-Max.Model.Dim
      stored.bf<-0
      Sigma<-as.matrix(ld.matrix)
      Sigma.inv<-Sigma.inv.input
      zSigmaz<-zSigmaz.input
      #########
      if(!is.null(input.S)){
        S<-input.S
      }else{
        S<-integer(0)
      }
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
    }
    
    
    set.gamma.func<-function(S){
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
        #set of gamma+
        set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
        set.gamma[[2]]<-t(set.gamma[[2]])
        #set of gamma=
        set.gamma[[3]]<-t(add.function(set.gamma[[1]][1,]))
      }
      return(set.gamma)
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
    
    index.fun<-function(x){
      m<-as(matrix(0,nrow=nrow(x),ncol=p),'dgTMatrix')
      m@i<-as.integer(rep(1:nrow(x)-1,each=ncol(x)))
      m@j<-as.integer(c(t(x))-1)
      m@x=rep(1,nrow(x)*ncol(x))
      m<-as(m,"dgCMatrix") 
      return(m)
    }
    ######################################################
    for(l in 1:inner.all.iter){
      for(h in 1:10){
        ##############COMPUTATION ############
        {
          set.gamma<-set.gamma.func(S)  
          set.gamma.margin<-list()
          set.gamma.prior<-list()
          matrix.gamma<-list()
          if(length(S)!=0){
            S.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'dgCMatrix')
            S.model[,S]<-1
            current.log.margin<-marginal_likelihood(S,Sigma,z,zSigmaz,tau=tau.sample,p=p,p_S=p_S)+prior.dist(S.model)
          }else{
            current.log.margin<-prior.dist(null.model)
          }
          
          if(length(S)>1){
            
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
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,zSigmaz=zSigmaz,tau=tau.sample,p=p,p_S=p_S)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,zSigmaz=zSigmaz,tau=tau.sample,p=p,p_S=p_S)
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
          
          if(length(S)==1){
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
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,zSigmaz=zSigmaz,tau=tau.sample,p=p,p_S=p_S)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,zSigmaz=zSigmaz,tau=tau.sample,p=p,p_S=p_S)
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
          if(length(S)==0){
            
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
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,zSigmaz=zSigmaz,tau=tau.sample,p=p,p_S=p_S)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,zSigmaz=zSigmaz,tau=tau.sample,p=p,p_S=p_S)
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
          if(length(S)!=0){
            set.star<-data.frame(set.index=1:3,gamma.set.index=rep(NA,3),margin=rep(NA,3))
            for(i in 1){
              aa<-set.gamma.margin[[i]]-current.log.margin
              aa<-aa-aa[which.max(aa)]
              aa[which(is.nan(aa))]<-aa[which.min(aa)]
              set.star$gamma.set.index[i] <-c(sample(1:length(set.gamma.margin[[i]]),1,prob=exp(aa)))
              set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
              rm(aa)
            }
            for(i in 2:length(set.gamma)){
              aa<-set.gamma.margin[[i]]-current.log.margin
              aa<-aa-aa[which.max(aa)]
              aa[which(is.nan(aa))]<-aa[which.min(aa)]
              set.star$gamma.set.index[i]<-c(sample((1:length(set.gamma.margin[[i]]))[order(exp(aa),decreasing = T)[1:(floor(length(set.gamma.margin[[i]])/2))]],
                                                    1,prob=exp(aa)[order(exp(aa),decreasing = T)[1:(floor(length(set.gamma.margin[[i]])/2))]]))
              set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
              rm(aa)
            }
            
            
           print(set.star)
            ##second sample
            aa<-set.star$margin-current.log.margin-max(set.star$margin-current.log.margin)
            sec.sample<-sample(1:3,1,prob=exp(aa) )
            S<-set.gamma[[sec.sample]][set.star$gamma.set.index[[sec.sample]] ,]
          }else{
            set.star<-data.frame(set.index=1:1,gamma.set.index=rep(NA,1),margin=rep(NA,1))
            
            set.star$gamma.set.index <-c(sample((1:length(set.gamma.margin[[2]]))[order(exp(set.gamma.margin[[2]]-current.log.margin),decreasing = T)[1:(floor(p/2))]],
                                                1,prob=exp(set.gamma.margin[[2]]-current.log.margin)[order(exp(set.gamma.margin[[2]]-current.log.margin),decreasing = T)[1:(floor(p/2))]]))
            set.star$margin<-set.gamma.margin[[2]][  set.star$gamma.set.index[2]]
            S<-set.gamma[[2]][set.star$gamma.set.index]
            print(set.star)
          }
          
          print(S)
          #print(h)
          
        }
      }
 
      ######Result
      
      result.B.list<-list()
      result.B.list[[1]]<-B.list[[1]][1:min(length(B.list[[1]]),B)]
      result.B.list[[2]]<-B.list[[2]][1:min(length(B.list[[1]]),B),]
      result.prob<-PIP.func(result.B.list[[1]],result.B.list[[2]])
      
      if(!is.null(output.labels)){
        if(dir.exists(output.labels)==F ){
          dir.create(output.labels,recursive = T)
        }
      write.table(result.B.list[[1]],file=paste0(output.labels,'/post_',label,'_poi_likeli_','.txt'),row.names = F,col.names = F)
      writeMM(result.B.list[[2]],file=paste0(output.labels,'/post_',label,'_poi_gamma_','.mtx'))
      write.table((result.prob),file=paste0(output.labels,'/post_', label,'_','.txt'),row.names = F,append = F,col.names = F)
      }
      
      difference<-abs(mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])-stored.bf)
      print(difference)
      if(difference<epsilon){
        break
      }else{
        stored.bf<-mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])
      }
    }
    
    return(list(result.B.list,C.list,result.prob))
  }
  
  ######## All burning###########
  previous.result<-list()
  for(i in 1:L){
    t0=Sys.time()
   all.C.list[[i]]<-Module.Cauchy.Shotgun(z.list[[i]],ld.list[[i]],Sigma.inv.input=Sigma.inv.list[[i]],zSigmaz.input = zSigmaz.list[[i]],
                                           Max.Model.Dim=Max.Model.Dim, input.S = S.list[[i]],epsilon=epsilon.list[[i]],
                                           lambda = lambda.list[[i]],
                                           label = label.list[[i]],output.labels = output.labels,
                                           effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
   t1=Sys.time()-t0
   print(paste0('This is locus ',i,' burning time'))
   print((t1))
  }

  if(!is.null(w.list)){
  w<-matrix(NA,nrow=0,ncol=ncol(w.list[[1]]))
  colnames(w)<-colnames(w.list[[1]])
  for(i in 1:L){
    w<-rbind(w,w.list[[i]])
  }
  }
  for(g in 1:all.iter){ 

    ################################
    if(!is.null(w.list)){
    if(!standardize.model.space){
      model.space.count<-c()
      for(i in 1:L){
        model.space.count<-c(model.space.count,colSums(all.C.list[[i]][[1]][[2]]))
      }
    }else{
      model.space.count<-c()
      for(i in 1:L){
        indi.count<-colSums(all.C.list[[i]][[1]][[2]])
        indi.count<-floor(colSums(all.C.list[[i]][[1]][[2]])/nrow(all.C.list[[i]][[1]][[2]])*Max.Model.Dim)
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
        write.table((glm.beta),file=paste0(output.labels,'/post_', label.list[[i]],'_','_theta.txt'),row.names = F,append = F,col.names = F)
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
    ########    

    for(i in 1:L){
      previous.result[[i]]<-mean(all.C.list[[i]][[1]][[1]][1:round(quantile(1:length(all.C.list[[i]][[1]][[1]]),probs = 0.25))])
      t0=Sys.time()
      all.C.list[[i]]<-Module.Cauchy.Shotgun(z.list[[i]],ld.list[[i]],Sigma.inv.input=Sigma.inv.list[[i]],zSigmaz.input = zSigmaz.list[[i]],
                                             Max.Model.Dim=Max.Model.Dim,
                                             C.list = all.C.list[[i]][[2]],prior.prob = prior.prob.list[[i]],
                                             input.S = S.list[[i]],lambda = lambda.list[[i]],epsilon=epsilon.list[[i]],
                                             label = label.list[[i]],output.labels = output.labels,
                                             effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
      t1=Sys.time()-t0
      print(paste0('This is locus ',i,' computing time'))
      print((t1))
      
    
    }
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
    names(results.list[[i]])<-c('PIPs','Credible set','Credible model')
    
  }
  return(results.list)
  
}
