\name{Model}
\alias{Model}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Model(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot Model of an object of class EBSProfiles
}
\value{
	An object of class string returning the model used for the analysis of the datasets.
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getModel}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Model(x) # retrieves the Model slot from x
}

\keyword{methods}
\keyword{Model}
