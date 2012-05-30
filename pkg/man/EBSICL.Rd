\name{EBSICL}
\alias{EBSICL}
\alias{EBSICL.default}
\title{
Model Selection by Integrated Completed Likelihood criterion
}
\description{
Chooses the optimal number of segments according to the exact ICL criterion: k= argmin(-Loglikelihood (data,K) + H(m|K) ) where H is the entropy of the segmentation.
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
	An integer containing the choice of the optimal number of segments.  
}
\references{
Guillem Rigaill and Emilie Lebarbier and Stephane Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{
Alice Cleynen
}


\seealso{
\code{\link{EBSegmentation}}, \code{\link{EBSDistrib}}
}
\examples{
# changes for Poisson model
x=c(rpois(15,1),rpois(60,4),rpois(200,0.8), rpois(85,3.1))
ans=EBSegmentation(x,Kmax=20) 
bestK=EBSICL(ans)
print(bestK)
}

