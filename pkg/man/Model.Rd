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
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
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
