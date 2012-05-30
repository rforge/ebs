\name{EBSPlotProba}
\alias{EBSPlotProba}
\alias{EBSPlotProba.default}

\title{
Plot distribution of changepoints
}

\description{
Given a number of segments, plots the posterior distribution of each of the changepoints
}

\usage{
EBSPlotProba(x, k, data=FALSE, file=character(), type='pdf',...)
}

\arguments{
  \item{x}{
	An object of class EBS returned by function EBSegmentation applied to data of interest.
}
  \item{k}{
	The number of segments of the segmentation for which the posterior distribution is wanted. Must have 2<k<x$Kmax
}

  \item{data}{
	A logical vector. If TRUE, the data is plotted as well as the posterior distributions. By default, data=FALSE.
}
  \item{file}{
	An object of type string. If filled, the plot is saved in a file which name is given by the file argument, and which type is given by argument type.
}
  \item{type}{
	An object of type string. If file is filled, argument type determines the type of the file saved. Possible values are 'pdf', 'png' and 'ps'. By default, 'pdf' is used.
}
}

\details{
	This function is used to plot the posterior distribution of all changepoints of a segmentation in k segments. 
}

\value{
	A plot of the posterior distributions.   
}

\references{
Guillem Rigaill and Emilie Lebarbier and Stephane Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}

\author{
Alice Cleynen
}


\seealso{
\code{\link{EBSegmentation}}, \code{\link{EBSDistribution}}
}
\examples{
# changes for Poisson model
x=c(rpois(150,1),rpois(50,3),rpois(20,0.1), rpois(80,2.5))
ans=EBSegmentation(x,Kmax=10) 
EBSPlotProba(ans,4)
EBSPlotProba(ans,4,data=TRUE,file="mysegmentation.png",type='png')
}

