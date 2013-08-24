
\name{EBSPriorProfiles}
\alias{EBSPriorProfiles}
\alias{EBSPriorProfiles.default}
\title{
Matrix of prior values for each interval
}
\description{
Computes the matrix of indexes nr (values associated to prior on segmentation) for each segment r. 
}
\usage{
EBSPriorProfiles(n=numeric(), K = 3, unif=TRUE)
}
\arguments{
  \item{n}{
	The size of the series.
}
  \item{K}{
	A vector of size the number of series which elements are the maximum number of segments that will be considered for each profile.
}
  \item{unif}{
	A boolean stating whether the uniform prior will be used for the segmentation. 
}


}
\details{
	This function is used to compute the values associated with the prior on the segmentation to use in computations such as ICL.  
}
\value{
	An object of class EBSProfiles with values associated to prior on segmentation.  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSProfiles}}, \code{\link{EBSPrior}}
}
\examples{
outPrior <- EBSPriorProfiles(1000,K=c(4,5,5),unif=FALSE)
}

