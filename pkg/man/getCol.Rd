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
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
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
