\name{EBSICLProfiles}
\alias{EBSICLProfiles}
\alias{EBSICLProfiles.default}
\title{
Model Selection by Integrated Completed Likelihood criterion
}
\description{
 For each profile, computes the exact ICL criterion: -Loglikelihood (data,K) + H(m|K) where H is the entropy of the segmentation, and chooses the optimal number of segments as k= argmin(ICL)
}
\usage{
EBSICLProfiles(x, prior=numeric())
}
\arguments{
  \item{x}{An object of class EBSPofiles returned by function EBSProfiles applied to matrix of profiles of interest.}
  \item{prior}{A vector of length Kmax giving prior probabilities on the value of K. Default value is uniform on 1:Kmax.}
}

\details{
	For each condition, this function is used to compute the entropy of the segmentation in k segments (for k in 1 to Kmax) and choose the optimal K according to the ICL criteria. 
}
\value{
	\item{NbICL}{A vector containing the choice of the optimal number of segments for each profile. }
	\item{ICL}{A list of vector (one for each condition) of length getK(x)[l] containing the ICL values.} 
}
\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}
\author{Alice Cleynen}


\seealso{
\code{\link{EBSProfiles}}, \code{\link{EBSICL}}
}

\examples{
# changes for Poisson model
set.seed(1)
x1<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
x2<-c(rpois(100,1),rpois(100,3),rpois(75,2),rpois(125,0.5))
M<-rbind(x1,x2)
E<-EBSProfiles(M, K=8)
out <- EBSICLProfiles(E)
print(out$NbICL)
}

