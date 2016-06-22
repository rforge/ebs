\name{matProba}
\alias{matProba}

\title{ Generic function}
\description{
 Generic function
}
\usage{
matProba(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot P of an object of class EBSProfiles
}
\value{
	A matrix of size (n+1)*(n+1) containing segment probabilities.
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getP}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
matProba(x) # retrieves the P slot from x
}

\keyword{methods}
\keyword{matProba}
