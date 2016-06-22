\name{getHyperParameters}
\alias{getHyperParameters}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getHyperParameters(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot HyperParameters of an object of class EBSProfiles
}
\value{
	A vector of HyperParameters used for the analysis of the dataset.
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{HyperParameters}} }


\examples{
x=new("EBS") # new EBS object
getHyperParameters(x) # retrieves the HyperParameters slot from x
}

\keyword{methods}
\keyword{HyperParameters}
