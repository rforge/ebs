\name{getCol}
\alias{getCol}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getCol(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot Col of an object of class EBS
}
\value{
	A list where each element is a the Col matrix of size (n+1)*Kmax of a profile. 
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
getCol(x) # retrieves the Col slot from x
}

\keyword{methods}
\keyword{Col}
