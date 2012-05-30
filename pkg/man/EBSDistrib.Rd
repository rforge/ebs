\name{EBSDistrib}
\alias{EBSDistrib}
\alias{EBSDistrib.default}
\title{
Probability of kth changepoint in a segmentation in K segments
}
\description{
Calculates the bayesian posterior probability of a changepoint of a segmentation in K segments. 
}
\usage{
EBSDistrib(x, k, Kk,...)
}
\arguments{
  \item{x}{
	An object of class EBS returned by function EBSegmentation applied to data of interest.
}
  \item{k}{
	The rank of the breakpoint for which the posterior distribution is wanted. Must have 0<k<K
}
  \item{K}{
	The number of segments for the segmentation of interest. Must have 1<K<=x$Kmax
}


}
\details{
	This function is used to compute the posterior distribution of kth changepoint for a segmentation in K segments. 
}
\value{
	An vector containing distibution of kth changepoint in a segmentation in K segments.  
}
\references{
Guillem Rigaill and Emilie Lebarbier and Stephane Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{
Alice Cleynen
}


\seealso{
\code{\link{EBSegmentation}}, \code{\link{EBSICL}}
}
\examples{
# changes for Poisson model
x=c(rpois(150,1),rpois(50,3),rpois(20,0.1), rpois(80,2.5))
ans=EBSegmentation(x,Kmax=10) 
y1=EBSDistrib(ans,1,4)
plot(y1,type='l')
}

