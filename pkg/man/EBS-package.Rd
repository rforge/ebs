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
Alice Cleynen
Maintainer: Alice Cleynen <alice.cleynen@agroparistech.fr>
}

\references{
Guillem Rigaill and Emilie Lebarbier and Stephane Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}

}
\keyword{changepoint}
\keyword{segmentation}
\keyword{exact}
\keyword{Bayesian}

\examples{
# changes for Poisson model
x=c(rpois(100,1),rpois(50,10),rpois(150,3))
ans=EBSegmentation(x) 
prob=EBSDistrib(ans,1,3) #finds the distribution of the first breakpoint of a segmentation in 2 segments
plot(prob, type='l')
bestk=EBSICL(ans)
print(bestk) #finds optimal number of segments according to ICL criterion

}
