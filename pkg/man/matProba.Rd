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
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
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
