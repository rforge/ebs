\name{Variance}
\alias{Variance}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Variance(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot Variance of an object of class EBSProfiles
}
\value{
	If model is Gaussian homoscedastic, the value of the variance used for each profile in the analysis.
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Variance}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Variance(x) # retrieves the Variance slot from x
}

\keyword{methods}
\keyword{Variance}
