\name{EBSPlotProba}
\alias{EBSPlotProba}
\alias{EBSPlotProba.default}

\title{Plot distribution of changepoints of one profile}

\description{Given a profile and its number of segments, plots the posterior distribution of each of the change-points}

\usage{
EBSPlotProba(x,K,data=FALSE, file=character(), type='pdf')
}

\arguments{
  \item{x}{An object of class EBS returned by function EBSegmentation applied to data of interest.}
  \item{K}{The number of segments of the segmentation for which the posterior distributions are wanted. Must have 2<K<x$Kmax}
  \item{data}{A logical vector. If TRUE, the data is plotted as well as the posterior distributions. By default, data=FALSE.}
  \item{file}{An object of type string. If filled, the plot is saved in a file which name is given by the file argument, and which type is given by argument type.}
  \item{type}{An object of type string. If file is filled, argument type determines the type of the file saved. Possible values are 'pdf', 'png' and 'ps'. By default, 'pdf' is used.}
}

\details{
	For a single profile, this function is used to plot the posterior distribution of all changepoints of a segmentation in K segments. 
}

\value{
	A plot of the posterior distributions.   
}

\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems 
}

\author{Alice Cleynen}


\seealso{
\code{\link{EBSegmentation}}, \code{\link{EBSDistrib}}
}
\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
out <- EBSegmentation(x,model=1,Kmax=20)
EBSPlotProba(out,4)
EBSPlotProba(out,4,data=TRUE,file="mysegmentation.png",type='png')
}

