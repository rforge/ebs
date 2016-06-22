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
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
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
