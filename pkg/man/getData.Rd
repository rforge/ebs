\name{getData}
\alias{getData}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getData(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the dataset used from an object of class EBS.
}
\value{
	A vector with the dataset of a profile analyzed.  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{dataLength}}}


\examples{
x=new("EBS") # new EBS object
getData(x) # retrieves the Data from x
}

\keyword{methods}
\keyword{Data}
