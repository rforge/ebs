\name{HyperParameters}
\alias{HyperParameters}

\title{ Generic function}
\description{
 Generic function
}
\usage{
HyperParameters(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot HyperParameters of an object of class EBSProfiles
}
\value{
	A vector of HyperParameters used for each profile 
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getHyperParameters}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
HyperParameters(x) # retrieves the HyperParameters slot from x
}

\keyword{methods}
\keyword{HyperParameters}
