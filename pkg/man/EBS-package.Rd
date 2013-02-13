\name{EBS-package}
\alias{EBS-package}
\alias{EBS}
\docType{package}
\title{
Contains funcions that run exact bayesian changepoint methods and return changepoint probabilities and ICL criteria for model selection
}
\description{
Implements changepoint method in an exact baysian framework for finding single and multiple changepoints within data.  Retrieves each changepoint probabilities for segmentations in 1 to Kmax segments. Chooses the optimal number of segments according to the ICL criterion. Compares change-point location between profiles using credibility intervals or likelihood ratios.
}
\details{
\tabular{ll}{
Package: \tab EBS\cr
Type: \tab Package\cr
Version: \tab 2.0 \cr
Date: \tab 2012-11-26\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Alice Cleynen

Maintainer: Alice Cleynen <alice.cleynen@agroparistech.fr>
}

\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems

Johnson, Kotz & Kemp: Univariate Discrete Distributions

Hall, Kay & Titterington: Asymptotically optimal difference-based estimation of variance in non-parametric regression
}

\keyword{changepoint}
\keyword{segmentation}
\keyword{exact}
\keyword{Bayesian}

\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
out <- EBSegmentation(x,model=1,Kmax=20)
bic <- EBSBIC(out)
print(bic$NbBIC)
icl <- EBSICL(out)
print(icl$NbICL)
plot(bic$BIC,type='b',pch=1,col='blue',ylim=c(0,1000))
lines(icl$ICL,type='b',pch=2,col='red')
EBSPlotProba(out, 5, data=TRUE)

# changes for Negative Binomial model, comparison of two profiles
set.seed(1)
x1<-c(rnbinom(125,size=0.2,prob=0.8),rnbinom(100,size=0.2, prob=0.1),
rnbinom(50,size=0.2,prob=0.6),rnbinom(75,size=0.2, prob=0.95),
rnbinom(50,size=0.2,prob=0.25))
x2<-c(rnbinom(125,size=0.15,prob=0.75),rnbinom(75,size=0.15,prob=0.2),
rnbinom(75,size=0.15,prob=0.9),rnbinom(125,size=0.15,prob=0.1))
M<-rbind(x1,x2)
E <- EBSProfiles(M,model=3,K=10,homoscedastic=TRUE) 

# Computes probabilities for both profile assuming independance but common 
#overdispersion
EBSPlotProbaProfiles(E,K=c(5,4)) 

# Plots posterior distribution of each change points of the two profiles, 
#the first into 5 segments, the second into 4.
mass<-CompCredibility(E,Conditions=c(1,2),Tau=c(1,1),K=c(5,4)) 

# Computes the distribution and credibility interval of the difference of 
#location of the first change point of the two profiles, 
#the first being devided into 5 segments, the second into 4
mass$massto0
DecisionStatistic<-EBSStatistic(E,Conditions=c(1,2),Tau=c(1,1)) 

# Computes the likelihood ratio of the profiles having same first 
#change-point versus complementary.

}
