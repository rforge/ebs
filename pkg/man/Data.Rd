\name{Data}
\alias{Data}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Data(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the datasets used from an object of class EBSProfiles.
}
\value{
	A matrix where each row is the dataset of a profile analyzed.  
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Length}}, \code{\link{NbConditions}}}


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Data(x) # retrieves the Data from x
}

\keyword{methods}
\keyword{Data}
