\name{TruncPois}
\alias{TruncPois}
\title{
Truncated Poisson distribution
}
\description{
 Given a parameter lambda and a number Kmax, computes the distribution of the truncated Poisson(lambda) up to Kmax.
}
\usage{
TruncPois(lambda,Kmax)
}
\arguments{
  \item{lambda}{Value of the Poisson parameter wanted}
  \item{Kmax}{Maximum number of segments considered}
}

\details{
	Given a parameter lambda and a number Kmax, computes the distribution of the truncated Poisson(lambda) up to Kmax.
}
\value{
	A vector of probabilities of size Kmax with truncated poisson probability.  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSPostK}} }


\examples{
# truncated Poisson with parameter 3.5 and Kmax=20
priorK<-TruncPois(3.5,20)
}

