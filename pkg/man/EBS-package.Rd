\name{EBS-package}
\alias{EBS-package}
\alias{EBS}
\docType{package}
\title{
Contains funcions that run exact bayesian changepoint methods and return changepoint probabilities and ICL criteria for model selection
}
\description{
Implements changepoint method in an exact baysian framework for finding single and multiple changepoints within data.  Retrieves each changepoint probabilities for segmentations in 1 to Kmax segments. Chooses the optimal number of segments according to the ICL criterion.
}
\details{
\tabular{ll}{
Package: \tab EBS\cr
Type: \tab Package\cr
Version: \tab 1.0 \cr
Date: \tab 2012-05-15\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Alice Cleynen, Guillem Rigaill, Stephane Robin
Maintainer: Alice Cleynen <alice.cleynen@agroparistech.fr>
}

\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}

\keyword{changepoint}
\keyword{segmentation}
\keyword{exact}
\keyword{Bayesian}

\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
out <- EBSegmentation(x,Kmax=20)
bic <- EBSBIC(out)
print(bic$NbBIC)
icl <- EBSICL(out)
print(icl$NbICL)
plot(bic$BIC,type='b')
lines(icl$ICL,type='b')
EBSPlotProba(out, icl$NbICL, data=TRUE, file="my-segmentation.pdf")

}
