\name{Overdispersion}
\alias{Overdispersion}

\title{ Generic function}
\description{
 Generic function
}
\usage{
Overdispersion(object)
}
\arguments{
  \item{object}{An object of class EBSProfiles}
}

\details{
	Returns the slot Overdispersion of an object of class EBSProfiles
}
\value{
	If model is Negative Binomial, the value of the overdispersion used for each profile in the analysis.
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getOverdispersion}} }


\examples{
x=new("EBSProfiles") # new EBSProfiles object
Overdispersion(x) # retrieves the P slot from x
}

\keyword{methods}
\keyword{Overdispersion}
