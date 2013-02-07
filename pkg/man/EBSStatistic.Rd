\name{EBSStatistic}
\alias{EBSStatistic}
\alias{EBSStatistic.default}
\title{Statistic for Profile Comparison} 

\description{Posterior probability of profiles having same change-point location}

\usage{EBSStatistic(x, Conditions = numeric(), Tau = numeric(), K = numeric(), a0=1)}

\arguments{
  \item{x}{An object of class EBSProfiles, outcome of function EBSProfiles applied to matrix of profiles of interest}
  \item{Conditions}{A vector containing the set of index of profiles to compare using posterior probabilities.}
  \item{Tau}{The vector of index of the change-point of interest for each profile. If field is left empty, default value is $1$ for all profiles. If only one value is entered, this value is used for all profiles.}
  \item{K}{The vector of number of segments of the segmentation for each profile. If field is left empty, function calls EBSICLProfiles to choose value of K. If only one value is entered, this value is used for all profiles.}
  \item{a0}{The ratio of prior probability of having same change-point under uniform segmentation and alternative neasure (see details). If field is left empty, default value is 1.}
}

\value{ The value of the likelihood ratio (see details).}

\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems 

Johnson, Kotz & Kemp: Univariate Discrete Distributions

Hall, Kay & Titterington: Asymptotically optimal difference-based estimation of variance in non-parametric regression
}

\author{Alice Cleynen}


\seealso{\code{\link{EBSProfiles-class}}, \code{\link{EBSProfiles}}, \code{\link{EBSegmentation}}, \code{\link{EBSICLProfiles}}}

\examples{
# Comparison of 3rd change-point of 2 profiles segmented in 5 and 4 segments with Poisson model
set.seed(1)
x1<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
x2<-c(rpois(100,1),rpois(100,3),rpois(75,2),rpois(125,0.5))
M<-rbind(x1,x2)
E<-EBSProfiles(M, K=8)
out <- EBSStatistic(E,Conditions=c(1,2),Tau=c(3,3),K=c(5,4))
}

\details{This function returns \eqn{p(\mathcal{E}_0|Y,K)}{p(E_0|Y,K)} where:

- \eqn{Y}{Y} is the matrix of data: \eqn{Y=(Y^1, \dots, Y^I)}{Y=(Y^1, \dots, Y^I)},

- \eqn{m_l}{m_l} is the segmentation of profile \eqn{l}{l},

- \eqn{k_l}{k_l} is the index of the change-point of interest in profile \eqn{l}{l}, and \eqn{\tau_{k_l}}{tau_(k_l)} is the corresponding change-point,

- \eqn{\mathcal{E}_0}{E_0} denotes the event \eqn{\tau_{k_1} = \dots = \tau_{k_I}}{tau_(k_1)=...=tau_(k_I)},

- \eqn{q}{q} denotes the uniform prior on segmentation \eqn{m}{m},

- \eqn{p}{p} denotes a probability measure chosen by user for which:

- \eqn{a_0=p(\mathcal{E}_0|K)/q(\mathcal{E}_0|K)}{a_0=p(E_0|K)/q(E_0|K)} 

- \eqn{a_1=p(\mathcal{E}_1|K)/q(\mathcal{E}_1|K)}{a_1=p(E_1|K)/q(E_1|K)}


Then the function returns \deqn{S(Y)= \frac{a_0 q(\mathcal{E}_0|Y,K) }{ (a_0-a_1)q(\mathcal{E}_0|Y,K) +a_1}}{S(Y)=a_0 q(E_0|Y,K) /[(a_0-a_1)q(E_0|Y,K) +a_1]}
}

