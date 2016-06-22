\name{getPriorm}
\alias{getPriorm}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getPriorm(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot unif of an object of class EBS
	}
\value{
	A boolean stating whether prior on segmentation is uniform. 
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Col}},\code{\link{getLi}} }


\examples{
x=new("EBS") # new EBS object
getPriorm(x) # retrieves the unif slot from x
}

\keyword{methods}
\keyword{unif}
