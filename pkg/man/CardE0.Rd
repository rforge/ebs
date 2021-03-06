\name{CardE0}
\alias{CardE0}
\title{
Prior probability of E0
}
\description{
 Probability under uniform that profiles share same change-point location
}
\usage{
CardE0(n,k,K,unif=TRUE)
}
\arguments{
  \item{n}{Integer giving length of datasets.}
  \item{k}{Vector of integers giving change-point number of datasets.}
  \item{K}{Vector of integers giving number of segments in segmentation of each profile.}
  \item{unif}{Boolean stating whether prior on segmentation is uniform given K. If false then prior favors segments of equal length.}
}

\details{
	Returns the probability under the uniform that profiles of length n have their k[i]th change-point at same location when segmented into in K[i] segments.
}
\value{
	A numeric object giving the prior probability of E0.  
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{PriorDistrib}} }


\examples{
# probability that two profiles of size 100 have their 1st and second change-point at same location
# when segmented into 5 segments: 
n<-100
k<-c(1,3)
K<-c(5,5)
CardE0(n,k,K,TRUE)
}

