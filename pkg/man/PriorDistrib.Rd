\name{PriorDistrib}
\alias{PriorDistrib}
\title{
Prior distribution of change-point when uniform prior on segmentation
}
\description{
 Computes the prior distribution of a given change-point when using a uniform prior on segmentation with known K
}
\usage{
PriorDistrib(n,k,K)
}
\arguments{
  \item{n}{Integer giving length of dataset.}
  \item{k}{Integer of index of given change-point.}
  \item{K}{Integer giving number of segments in segmentation.}
}

\details{
	This function is used to compare prior and posterior change-point distributions.
}
\value{
	A vector of length n with the change-point distribution.  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSegmentation}}}


\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
y=PriorDistrib(length(x),1,5)
plot(y,type='l')
}

