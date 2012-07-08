\name{EBSBIC}
\alias{EBSBIC}
\alias{EBSBIC.default}
\title{
Model Selection by BIC Criterion
}
\description{
 Computes the exact BIC criterion: -Loglikelihood (data,K)) and chooses the optimal number of segments as k= argmin(BIC)
}
\usage{
EBSBIC(x)
}
\arguments{
  \item{x}{
	An object of class EBS returned by function EBSegmentation applied to data of interest.
}
}

\details{
	This function is used to choose the optimal K according to the BIC criteria. 
}
\value{
	\item{NbBIC}{An integer containing the choice of the optimal number of segments.} 
	\item{BIC}{Vector of -Loglikelihood (data,K).}  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{
Alice Cleynen, Guillem Rigaill, Stephane Robin
}


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

