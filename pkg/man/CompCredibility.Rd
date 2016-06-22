\name{CompCredibility}
\alias{CompCredibility}
\alias{CompCredibility.default}
\alias{print.Credibility}
\alias{plot.Credibility}
\title{Comparison of two profiles with credibility intervals}

\description{
Given two conditions, computes the posterior distribution of the difference of change-point locations, and its credibility interval.
}
\usage{
CompCredibility(x, Conditions, Tau = numeric(), K = numeric()) 
}

\arguments{
  \item{x}{An object of class EBSProfiles.}
  \item{Conditions}{A vector of length 2 containing the index of the two conditions to compare.}
  \item{Tau}{A vector of length 2 containing the index of the change-point of interest of the two conditions to compare.}
  \item{K}{A vector of length 2 containing the maximum number of segments for the segmentation of each of the two profiles to compare.}
}

\details{This function is used to compute the posterior credibility interval of the difference of change-point locations between two profiles assumed to be independant.}

\value{  
	\item{Distribution}{Posterior distribution of the difference between the location of change-points of interest for the two conditions}
	\item{masswith0}{Mass of the smallest credibility interval up to and including data-point 0}
	\item{massto0}{Mass of the smallest credibility interval up to but excluding data-point 0}
}

\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Johnson, Kotz & Kemp: Univariate Discrete Distributions

Hall, Kay & Titterington: Asymptotically optimal difference-based estimation of variance in non-parametric regression
}

\author{Alice Cleynen}


\seealso{\code{\link{EBSProfiles}}, \code{\link{EBSStatistic}}}
\examples{
# changes for Poisson model
set.seed(1)
x1<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
x2<-c(rpois(125,3),rpois(75,4),rpois(75,1),rpois(125,8))
M<-rbind(x1,x2)
E <- EBSProfiles(M,model=1,K=10)
out<-CompCredibility(E,Conditions=c(1,2),Tau=c(1,1),K=c(5,4))
}

