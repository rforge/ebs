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
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSPostK}} }


\examples{
# truncated Poisson with parameter 3.5 and Kmax=20
priorK<-TruncPois(3.5,20)
}

