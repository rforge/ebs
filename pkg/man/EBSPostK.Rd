\name{EBSPostK}
\alias{EBSPostK}
\alias{EBSPostK.default}
\title{
Posterior probabilities of the number of segments
}
\description{
For a single profile, computes the posterior probabilities of the number of segments given its prior.
}
\usage{
EBSPostK(x, prior=numeric())
}
\arguments{
  \item{x}{An object of class EBS returned by function EBSegmentation applied to data of interest.}
  \item{prior}{A vector of length Kmax giving prior probabilities on the value of K. Default value is uniform on 1:Kmax.}
}

\details{This function computes the posterior probabilities of the number of segments, so that the user can use Bayesian Model Averaging.}

\value{
	A vector of size Kmax containing the probability of each of the number of segments.  
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}

\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
out <- EBSegmentation(x,model=1,Kmax=20)
Posterior<-EBSPostK(out)
plot(Posterior, type='b')
}

\seealso{\code{\link{EBSegmentation}}, \code{\link{EBSBIC}}, \code{\link{EBSICL}} }
