% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CARMA.R
\name{CARMA}
\alias{CARMA}
\title{CARMA}
\usage{
CARMA(
  z.list,
  ld.list,
  w.list = NULL,
  lambda.list = NULL,
  output.labels = ".",
  label.list = NULL,
  effect.size.prior = "Spike-slab",
  rho.index = 0.99,
  BF.index = 10,
  EM.dist = "Logistic",
  Max.Model.Dim = 2e+05,
  all.iter = 3,
  all.inner.iter = 10,
  input.alpha = 0,
  epsilon.threshold = 1e-05,
  printing.log = F,
  num.causal = 10,
  y.var = 1,
  tau = 0.04,
  outlier.switch = T,
  outlier.BF.index = 1/3.2,
  prior.prob.computation = "Logistic"
)
}
\arguments{
\item{z.list}{Input list of summary statistics at the testing loci; each element of the list is the summary statistics at each individual locus.}

\item{ld.list}{Input list of LD correlation matrix at the testing loci; each element of the list is the LD matrix at each individual locus.}

\item{w.list}{Input list of the functional annotations at the testing loci; each element of the list is the functional annotation matrix at each individual locus.}

\item{lambda.list}{Input list of the hyper-parameter \eqn{\eta} at the testing loci; each element of the list is the hyper-parameter of each individual locus.}

\item{output.labels}{Output directory where output will be written while CARMA is running. Default is the OS root directory ".".}

\item{label.list}{Input list of the names at the testing loci. Default is NULL.}

\item{effect.size.prior}{The prior of the effect size. The choice are 'Cauchy' and 'Spike-slab' priors, where the 'Spike-slab' prior is the default prior.}

\item{rho.index}{A number between 0 and 1 specifying \eqn{\rho} for the estimated credible sets.}

\item{BF.index}{The threshold of the Bayes factor of the estimated credible models. The default setting is 10.}

\item{EM.dist}{The distribution used to model the prior probability of being causal as a function of functional annotations. The default distribution is logistic distribution.}

\item{Max.Model.Dim}{Maximum number of the top candidate models based on the unnormalized posterior probability.}

\item{all.iter}{Maximum iterations for EM algorithm to run.}

\item{all.inner.iter}{Maximum iterations for Shotgun algorithm to run per iteration within EM algorithm.}

\item{input.alpha}{The elastic net mixing parameter, where \eqn{0\le}\eqn{\alpha}\eqn{\le 1}.}

\item{epsilon.threshold}{Convergence threshold measured by average of Bayes factors.}

\item{printing.log}{Whether print the running log while running CARMA.}

\item{num.causal}{The maximum number of causal variants assumed per locus, which is 10 causal SNPs per locus by default.}

\item{y.var}{The input variance of the summary statistics, the default value is 1 as the summary statistics are standardized.}

\item{tau}{The prior precision parameter of the effect size. The default value is computed based on the scenario of n=10000 and 1% heritability.}

\item{outlier.switch}{The indicator variable for outlier detection. We suggest that the detection should always be turned on if using external LD matrix.}

\item{outlier.BF.index}{The Bayes threshold for the Bayesian hypothesis test for outlier detection. The default setting is 1/3.2.}
}
\value{
The return is a list, for each list:
\itemize{
\item pip - The posterior inclusion probability of each individual locus.
\item Credibleset - The information on the credible set given a threshold \eqn{\rho}.
\item Credible model - The information on the credible model given a threshold  of the Bayes factor.
\item Outliers - The information regarding the detected outliers.
}
}
\description{
Performs a Bayesian fine-mapping model in order to identify putative causal variants at GWAS loci. The model requires the summary statistics
for the SNPs at the testing loci, the corresponding LD matrices for fine-mapping, and an estimated trait variance. Functional annotations can be included as prior 
information on the causality of the testing SNPs. The model also provides a procedure of outlier detection, which aims to identify discrepancies
between summary statistics and LD matrix extracted from a reference panel. The model can be executed chromosome-wise to increase power.
}
\details{
The function performs a Bayesian fine-mapping method.
}
\examples{
# Example 
set.seed(1)
n = 400
p = 500
beta = rep(0,p)
beta[1] = 1
X = matrix(rnorm(n*p),nrow = n,ncol = p)
X = scale(X,center = TRUE,scale = TRUE)
y = drop(X \%*\% beta + rnorm(n))
SS=compute_summary_statistics(X,y)
z.list<-list()
z.list[[1]]<-(SS$betahat/SS$sebetahat)
ld.list<-list()
ld.list[[1]]<-cov(X)
lambda.list<-list()
lambda.list[[1]]<-1/sqrt(p)
CARMA.result<-CARMA(z.list,ld.list=ld.list,
lambda.list = lambda.list,effect.size.prior='Hyper-g')
}
