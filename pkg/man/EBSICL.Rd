\name{EBSICL}
\alias{EBSICL}
\alias{EBSICL.default}
\title{
Model Selection by Integrated Completed Likelihood criterion
}
\description{
 Computes the exact ICL criterion: -Loglikelihood (data,K) + H(m|K) where H is the entropy of the segmentation, and chooses the optimal number of segments as k= argmin(ICL)
}
\usage{
EBSICL(x)
}
\arguments{
  \item{x}{
	An object of class EBS returned by function EBSegmentation applied to data of interest.
}
}

\details{
	This function is used to compute the entropy of the segmentation in k segments (for k in 1 to Kmax) and choose the optimal K according to the ICL criteria. 
}
\value{
	\item{NbICL}{An integer containing the choice of the optimal number of segments. }
	\item{ICL}{Vector of length x$Kmax containing the ICL values.} 
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{
Alice Cleynen, Guillem Rigaill, Stephane Robin
}


\seealso{
\code{\link{EBSegmentation}}, \code{\link{EBSBIC}}
}

\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
out <- EBSegmentation(x,model=1,Kmax=20)
bestK=EBSICL(out)$NbICL
print(bestK)
}

