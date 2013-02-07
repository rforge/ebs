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
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot Li of an object of class EBS
}
\value{
	A matrix of size Kmax*(n+1)  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Col}} }


\examples{
x=new("EBS") # new EBS object
Li(x) # retrieves the Li slot from x
}

\keyword{methods}
\keyword{Li}
