\name{EBSPostK}
\alias{EBSPostK}
\alias{EBSPostK.default}
\title{
Posterior probabilities of the number of segments
}
\description{
Computes the posterior probabilities of the number of segments assuming the prior is uniform on 1:Kmax
}
\usage{
EBSPostK(x)
}
\arguments{
  \item{x}{
	An object of class EBS returned by function EBSegmentation applied to data of interest.
}
}

\details{
	This function computes the posterior probabilities of the number of segments assuming the prior is uniform on 1:Kmax, so that the user can use Bayesian Model Averaging.
}
\value{
	A vector of size Kmax containing the probability of each of the number of segments.  
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
Posterior=EBSPostK(out)
plot(Posterior, type='b')
}

