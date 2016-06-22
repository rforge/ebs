\name{getModel}
\alias{getModel}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getModel(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot Model of an object of class EBS
}
\value{
	An object of class string returning the model used for the analysis of the dataset.
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Model}} }


\examples{
x=new("EBS") # new EBS object
getModel(x) # retrieves the Model slot from x
}

\keyword{methods}
\keyword{Model}
