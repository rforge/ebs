\name{NbConditions}
\alias{NbConditions}

\title{ Generic function}
\description{
 Generic function
}
\usage{
NbConditions(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the number of profiles analyzed from an object of class EBSProfiles
}
\value{
	An integer with the number of profiles.  
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getLength}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
NbConditions(x) # retrieves the number of profiles from x
}

\keyword{methods}
\keyword{length}
