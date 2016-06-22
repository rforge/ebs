\name{EBSPrior}
\alias{EBSPrior}
\alias{EBSPrior.default}
\title{
Matrix of prior values for each interval
}
\description{
Computes the matrix of indexes nr (values associated to prior on segmentation) for each segment r. 
}
\usage{
EBSPrior(n=numeric(), Kmax = 15, unif = TRUE)
}
\arguments{
  \item{n}{
	The size of the series.
}
  \item{Kmax}{
	The maximum number of segments that will be considered.
}
  \item{unif}{
	A boolean stating whether the uniform prior will be used for the segmentation. 
}


}
\details{
	This function is used to compute the values associated with the prior on the segmentation to use in computations such as ICL.  
}
\value{
	An object of class EBS with values associated to prior on segmentation.  
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSegmentation}}
}
\examples{
outPrior <- EBSPrior(1000,Kmax=20,unif=FALSE)
}

