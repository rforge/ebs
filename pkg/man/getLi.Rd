\name{getLi}
\alias{getLi}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getLi(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot Li of an object of class EBS
}
\value{
	A list where each element is the Li matrix of size Kmax*(n+1) of a profile 
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getCol}}, \code{\link{Li}} }


\examples{
x=new("EBS") # new EBS object
getLi(x) # retrieves the Li slot from x
}

\keyword{methods}
\keyword{Li}
