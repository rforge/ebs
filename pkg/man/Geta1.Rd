\name{Geta1}
\alias{Geta1}
\title{
Probability of segmentation from E_1 for new measure P
}
\description{
 Given a probability for segmentation from E_0, computes the associated probability for segmentation from E_1
}
\usage{
Geta1(n,k,K,a0)
}
\arguments{
  \item{n}{Integer giving length of dataset.}
  \item{k}{Integer of index of given change-point.}
  \item{K}{Integer giving number of segments in segmentation.}
  \item{a0}{Double giving ratio of probability of segmentation from E_0 under uniform and re-scaled measures.}
}

\details{
	Given a probability for segmentation from E_0, this function computes the associated probability for segmentation from E_1
}
\value{
	A numeric object giving a_1.  
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSProfiles}},\code{\link{EBSStatistic}} }


\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
a0=2
a1=Geta1(length(x),1,5,a0)
}

