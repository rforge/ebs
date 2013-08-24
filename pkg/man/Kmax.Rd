\name{Kmax}
\alias{Kmax}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Kmax(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot K of an object of class EBSProfiles
}
\value{
	An integer of the maximum number of segments considered for the analysis of the dataset.
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getKmax}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Kmax(x) # retrieves the Kmax slot from x
}

\keyword{methods}
\keyword{Kmax}
