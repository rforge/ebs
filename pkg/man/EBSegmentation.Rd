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
EBSegmentation(data=numeric(), model=1, Kmax = 15, hyper = numeric(), theta = numeric(), var = numeric())
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{model}{
	Model under which data is assumed to be distributed. Possible values are 1 for Poisson, 2 for Normal Homoscedastic, 3 for Negative Binomial and 4 for Normal Heteroscedastic.
}
  \item{Kmax}{
	The maximum number of segments for the segmentation. Function will find explore the set of all possible segmentations in k segments for k in 1 to Kmax.
}
  \item{hyper}{
	The set of hyper-parameters for the prior on the data-distribution. If model is Poisson the conjugate law is Gamma and 2 parameters are needed. If model is Negative Binomial the conjugate is Beta and 2 parameters are needed. If model is Normal the prior on the mean is normal, and if it is heteroscedastic the prior on the inverse variance is Gamma, so that 4 parameters are needed. If the user does not give his own hyperparameters, the package uses the following default values:

For the Poisson model, Gamma(1,1) is used. For Negative Binomial model, Jeffreys' prior, Beta(1/2,1/2) is used. For the Normal Homoscedastic, N(0,1) is used for a prior on the mean. Finally, for the Normal Heteroscedastic, the package computes the MAD on the data and fits an inverse-gamma distribution on the result. The parameters are used for the prior on the variance: IG(alpha,beta), and the prior on the mean is N(0,2*beta).
}

  \item{theta}{
	If model=3 (Negative binomial), the value of the inverse of the overdispersion parameter. If the user does not give his own hyperparameters, the package uses a modified version of Johnson and Kotz's estimator where the mean is replaced by the median. 
}
  \item{var}{
	If model=2 (Normal Homoscedastic), the value of the variance. If the user does not give his own hyperparameters, the package uses Hall's estimator whith d=4.
}

}
\details{
	This function is used to compute the matrix of segment probabilities assuming data is poisson, normal or negative binomial distributed. The probability of each interval being devided in k segments (k in 1 to Kmax) is computed. 
}
\value{
	An object of class "EBS".  
	\item{model}{Emission distribution (Poisson, Normal Homoscedastic, Negative Binomial or Normal Heteroscedastic)}
	\item{length}{the length of the data-set}
	\item{Kmax}{the maximum number of segments for the segmentation}
	\item{HyperParameters}{The hyperparameters used for the prior on the data distribution}
	\item{Li}{a matrix of size Kmax*(length+1). Element [i,j] is the log-probability of interval [1,j[ being segmented in j segments}
	\item{Col}{a matrix of size (length+1)*Kmax. Element [i,j] is the log-probability of interval [i,n] being segmented in i segments}
	\item{matProba}{a matrix of size (length+1)*(length+1). Element [i,j] is the log-probability of interval [i,j[}
}
\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}

Johnson, Kotz & Kemp: Univariate Discrete Distributions

Hall, Kay & Titterington: Asymptotically optimal difference-based estimation of variance in non-parametric regression

}
\author{
Alice Cleynen, Guillem Rigaill, Stephane Robin
}


\seealso{
\code{\link{EBSICL}},\code{\link{EBSDistrib}}
}
\examples{
# changes for Poisson model
set.seed(1)
x<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
out <- EBSegmentation(x,model=1,Kmax=20)
}

