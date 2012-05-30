\name{EBSegmentation}
\alias{EBSegmentation}
\alias{EBSegmentation.default}
\alias{print.EBS}
\title{
Exact Bayesian Segmentation for Poisson, Negative Binomial and Normal models
}
\description{
Calculates the bayesian probability of each segmentation in 1 to Kmax segments (assuming the data is poisson, normal or negative binomial distributed) and returns object of class EBS.
}
\usage{
EBSegmentation(data=numeric(), model=1, Kmax = 15, hyper = numeric(), theta=0.15, uniform=TRUE)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{model}{
	Model under which data is assumed to be distributed. Possible values are 1 for Poisson, 2 for Normal and 3 for Negative Binomial.
}
  \item{Kmax}{
	The maximum number of segments for the segmentation. Function will find explore the set of all possible segmentations in k segments for k in 1 to Kmax.
}
  \item{hyper}{
	The set of hyper-parameters for the prior on the data-distribution. If model is Poisson the conjugate law is Gamma and 2 parameters are needed. If model is Negative Binomial the conjugate is Beta and 2 parameters are needed. If model is Normal the prior on the mean is normal and the prior on the inverse variance is Gamma. Four parameters are needed. By default, if uniform is true and the parameters are set to use non informative priors. If parameters are user-defined, uniform has to be set to FALSE.
}

  \item{theta}{
	If model=3 (Negative binomial), the value of the inverse of the overdispersion parameter. By default set to 0.15 as a classical value for NGS data.
}
  \item{uniform}{
	Logical, default is TRUE.  If true then default parameters are used and set to non-informative priors. If False user HAS to enter parameters.
}

}
\details{
	This function is used to compute the matrix of segment probabilities assuming data is poisson, normal or negative binomial distributed. The probability of each interval being devided in k segments (k in 1 to Kmax) is computed. 
}
\value{
	An object of class "EBS".  
	\item{model}{Emission distribution (Poisson, Normal, or Negative Binomial)}
	\item{length}{the length of the data-set}
	\item{Kmax}{the maximum number of segments for the segmentation}
	\item{HyperParameters}{The hyperparameters used for the prior on the data distribution}
	\item{Li}{a matrix of size Kmax*(length+1). Element [i,j] is the log-probability of interval [1,j[ being segmented in j segments}
	\item{Col}{a matrix of size (length+1)*Kmax. Element [i,j] is the log-probability of interval [i,n] being segmented in i segments}
	\item{matProba}{a matrix of size (length+1)*(length+1). Element [i,j] is the log-probability of interval [i,j[}
}
\references{
Guillem Rigaill and Emilie Lebarbier and Stephane Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}
}
\author{
Alice Cleynen
}


\seealso{
\code{\link{EBSICL}},\code{\link{EBSDistrib}}
}
\examples{
# changes for Normal model
x=c(rnorm(50,1,1),rnorm(50,3,0.5),rnorm(100,3,2))
ans=EBSegmentation(x,model=2) 
}

