\name{getP}
\alias{getP}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getP(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot P of an object of class EBS
}
\value{
	A list of the matrix P of each profiles.
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{matProba}} }


\examples{
x=new("EBS") # new EBS object
getP(x) # retrieves the P slot from x
}

\keyword{methods}
\keyword{P}
