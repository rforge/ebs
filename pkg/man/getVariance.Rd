\name{getVariance}
\alias{getVariance}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getVariance(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot Variance of an object of class EBS
}
\value{
	If model is Gaussian homoscedastic, the value of the variance used for the analysis.
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Variance}} }


\examples{
x=new("EBS") # new EBS object
getVariance(x) # retrieves the Variance slot from x}

\keyword{methods}
\keyword{Variance}
