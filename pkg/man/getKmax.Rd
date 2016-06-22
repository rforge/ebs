\name{getKmax}
\alias{getKmax}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getKmax(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot Kmax of an object of class EBS
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
\code{\link{Kmax}} }


\examples{
x=new("EBS") # new EBS object
getKmax(x) # retrieves the Kmax slot from x
}

\keyword{methods}
\keyword{Kmax}
