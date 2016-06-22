\name{Kmax}
\alias{Kmax}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Kmax(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot K of an object of class EBSProfiles
}
\value{
	An integer of the maximum number of segments considered for the analysis of the dataset.
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getKmax}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Kmax(x) # retrieves the Kmax slot from x
}

\keyword{methods}
\keyword{Kmax}
