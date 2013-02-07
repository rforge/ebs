\name{EBSBIC}
\alias{EBSBIC}
\alias{EBSBIC.default}
\title{
Model Selection by BIC Criterion
}
\description{
 Computes the exact BIC criterion: -Loglikelihood (data,K) and chooses the optimal number of segments as k= argmin(BIC)
}
\usage{
EBSBIC(x, prior=numeric())
}
\arguments{
  \item{x}{An object of class EBS returned by function EBSegmentation applied to data of interest.}
  \item{prior}{A vector of size Kmax giving prior probabilities for segment numbers.}
}

\details{
	This function is used to choose the optimal K according to the BIC criteria. 
}
\value{
	\item{NbBIC}{An integer containing the choice of the optimal number of segments.} 
	\item{BIC}{A vector of length Kmax returning -Loglikelihood (data,K).}  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSegmentation}}, \code{\link{EBSICL}}, \code{\link{EBSPostK}}
}

\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
out <- EBSegmentation(x,model=1,Kmax=20)
bestKBIC=EBSBIC(out)$NbBIC
print(bestKBIC)
}

