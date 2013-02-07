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
EBSDistrib(x, k, Kk)
}
\arguments{
  \item{x}{
	An object of class EBS returned by function EBSegmentation applied to data of interest.
}
  \item{k}{
	The rank of the breakpoint for which the posterior distribution is wanted. Must have 0<k<Kk
}
  \item{Kk}{
	The number of segments for the segmentation of interest. 

Must have 2<Kk<=x$Kmax
}


}
\details{
	This function is used to compute the posterior distribution of kth changepoint for a segmentation in K segments. 
}
\value{
	A vector containing distibution of kth changepoint in a segmentation in Kk segments.  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSegmentation}}, \code{\link{EBSPlotProba}}
}
\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
out <- EBSegmentation(x,model=1,Kmax=20)
y1=EBSDistrib(out,1,5)
plot(y1,type='l')
}

