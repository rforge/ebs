\name{GetCondition}
\alias{GetCondition}
\alias{GetCondition.default}
\title{Segmentation information for one profile out of the set}
\description{Retrieves the segmentation information of one particular profile out of the set given}
\usage{GetCondition(x, Condition = numeric())}
\arguments{
  \item{x}{An object of class EBSProfiles}
  \item{Condition}{the index of the profile of interest}
}

\details{This function is used to retrieve the segmentation information of a profile when a set of conditions have been analysed together. }
\value{ An object of class "EBS".  }

\references{
Rigaill, Lebarbier & Robin (2012): Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \emph{Statistics and Computing}

Cleynen & Robin (2014): Comparing change-point location in independent series \emph{Statistics and Computing}

Johnson, Kotz & Kemp: Univariate Discrete Distributions

Hall, Kay & Titterington: Asymptotically optimal difference-based estimation of variance in non-parametric regression
}

\author{Alice Cleynen}


\seealso{\code{\link{EBS-class}},\code{\link{EBSProfiles-class}},\code{\link{EBSProfiles}},\code{\link{EBSegmentation}}}

\examples{
# changes for Poisson model
set.seed(1)
x1<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
x2<-c(rpois(100,1),rpois(100,3),rpois(75,2),rpois(50,0.5),rpois(75,3))
M<-rbind(x1,x2)
E<-EBSProfiles(M)
C1 <- GetCondition(E,1)
}

