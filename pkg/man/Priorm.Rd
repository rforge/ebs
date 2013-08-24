\name{Priorm}
\alias{Priorm}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Priorm(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot unif of an object of class EBSProfiles
	}
\value{
	A boolean stating whether prior on segmentation is uniform. 
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Col}},\code{\link{getLi}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Priorm(x) # retrieves the unif slot from x
}

\keyword{methods}
\keyword{unif}
