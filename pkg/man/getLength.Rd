\name{getLength}
\alias{getLength}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getLength(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the length of the dataset used from an object of class EBS
}
\value{
	An integer with size of a profile.  
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Length}} }


\examples{
x=new("EBS") # new EBS object
getLength(x) # retrieves the length of data from x
}

\keyword{methods}
\keyword{length}
