\name{Length}
\alias{Length}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Length(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the length of the dataset used from an object of class EBSProfiles
}
\value{
	An integer with size of data used.  
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getData}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Length(x) # retrieves the length of data from x
}

\keyword{methods}
\keyword{length}
