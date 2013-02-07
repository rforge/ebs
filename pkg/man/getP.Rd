\name{getP}
\alias{getP}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getP(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot P of an object of class EBSProfiles
}
\value{
	A list of the matrix P of each profiles.
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{MatProba}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
getP(x) # retrieves the P slot from x
}

\keyword{methods}
\keyword{P}
