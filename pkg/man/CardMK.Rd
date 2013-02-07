\name{CardMK}
\alias{CardMK}
\title{
Number of segmentations in K segments
}
\description{
 Number of all possible segmentations in K segments
}
\usage{
CardMK(n,K)
}
\arguments{
  \item{n}{Integer giving length of dataset.}
  \item{K}{Integer giving number of segments in segmentation.}
}

\details{
	Returns the number of possible segmentations of [1,n] in K segments
}
\value{
	A numeric object giving the number of possible segmentations of [1,n] in K segments.  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{PriorDistrib}} }


\examples{
# number of partitions of {1, 100} in 5 segments: 
n<-100
K<-5
CardMK(n,K)
}

