\name{EBSPlotProbaProfiles}
\alias{EBSPlotProbaProfiles}
\alias{EBSPlotProbaProfiles.default}

\title{Plot distribution of changepoints of each profile}

\description{For the set of profiles and their number of segments, plots the posterior distribution of each of the change-points}

\usage{
EBSPlotProbaProfiles(x,K=numeric(),data=FALSE)
}

\arguments{
  \item{x}{An object of class EBSProfiles returned by function EBSProfiles applied to the matrix of profiles of interest.}
  \item{K}{The vector of number of segments of the segmentation for which the posterior distributions are wanted (one value for each profile). Must have 2<K[l]<getK(x)[l] for all profile l}
  \item{data}{A logical vector. If TRUE, the data is plotted as well as the posterior distributions. By default, data=FALSE.}
}

\details{
	This function is used to plot the posterior distribution of all changepoints of a segmentation in K[l] segments for all profiles l. Graph is subdivided into NbConditions plots, one for each profile. 
}

\value{
	A plot of the posterior distributions, for each profile.   
}

\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}
}

\author{Alice Cleynen}


\seealso{
\code{\link{EBSegmentation}}, \code{\link{EBSPlotProba}}
}
\examples{
# changes for Poisson model
set.seed(1)
x1<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
x2<-c(rpois(125,3),rpois(75,4),rpois(75,1),rpois(125,8))
M<-rbind(x1,x2)
E <- EBSProfiles(M,model=1,K=10)
EBSPlotProbaProfiles(E, K=c(5,4), data=TRUE)
}

