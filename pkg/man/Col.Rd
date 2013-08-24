\name{Col}
\alias{Col}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Col(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot Col of an object of class EBSProfiles
}
\value{
	A matrix of size (n+1)*Kmax  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Li}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Col(x) # retrieves the Col slot from x
}

\keyword{methods}
\keyword{Col}
