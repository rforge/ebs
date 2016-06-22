\name{Li}
\alias{Li}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Li(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot Li of an object of class EBSProfiles
}
\value{
	A matrix of size Kmax*(n+1)  
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Col}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Li(x) # retrieves the Li slot from x
}

\keyword{methods}
\keyword{Li}
