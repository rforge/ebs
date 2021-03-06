\name{getOverdispersion}
\alias{getOverdispersion}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getOverdispersion(object)
}
\arguments{
  \item{object}{An object of class EBS}
}

\details{
	Returns the slot Overdispersion of an object of class EBS
}
\value{
	If model is Negative Binomial, the value of the overdispersion used for the analysis.
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{Overdispersion}} }


\examples{
x=new("EBS") # new EBS object
getOverdispersion(x) # retrieves the P slot from x
}

\keyword{methods}
\keyword{Overdispersion}
