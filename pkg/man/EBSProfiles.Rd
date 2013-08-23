\name{EBSProfiles}
\alias{EBSProfiles}
\alias{EBSProfiles.default}
\title{Exact Bayesian Segmentation for multiple profiles}

\description{
For each profile i, calculates the bayesian probability of each segmentation in 1 to K[i] segments (assuming the data is poisson, normal or negative binomial distributed) and returns object of class EBSProfiles.
}
\usage{
EBSProfiles(data=numeric(), model=1, K = 3, hyper = numeric(), 

theta = numeric(), var = numeric(), homoscedastic = FALSE, unif= TRUE)
}

\arguments{
  \item{data}{A matrix where each line contains the data of one profile within which you wish to find changepoints.}
  \item{model}{Model under which each profile is assumed to be distributed. Possible values are 1 for Poisson, 2 for Normal Homoscedastic, 3 for Negative Binomial and 4 for Normal Heteroscedastic.}
  \item{K}{A vector containing the maximum number of segments for the segmentation of each profile. Function will explore the set of all possible segmentations in k segments for k in 1 to K[i]. If length(K)=1, the same value of K will be used for each profile.}
  \item{hyper}{The set of hyper-parameters for the prior on the data-distribution. If model is Poisson the conjugate law is Gamma and 2 parameters are needed for each profile (ie vector of length 2*(number of profiles)). If model is Negative Binomial the conjugate is Beta and 2 parameters are needed for each profile (ie vector of length 2*(number of profiles)). If model is Normal the prior on the mean is normal, and if it is heteroscedastic the prior on the inverse variance is Gamma, so that 4 parameters are needed for each profile (ie vector of length 4*(number of profiles)). The first two are the mean hyperparameters, the last two are the variance's. If the user does not give his own hyperparameters, the package uses the following default values:

For the Poisson model, Gamma(1,1) is used. For Negative Binomial model, Jeffreys' prior, Beta(1/2,1/2) is used. For the Normal Homoscedastic, N(0,1) is used for a prior on the mean. Finally, for the Normal Heteroscedastic, the package computes the MAD on the data and fits an inverse-gamma distribution on the result. The parameters are used for the prior on the variance: IG(alpha,beta), and the prior on the mean is N(0,2*beta).}

  \item{theta}{If model=3 (Negative binomial), the vector of values of the inverse of the overdispersion parameter for each profile. If the user does not give his own hyperparameters, the package uses a modified version of Johnson and Kotz's estimator where the mean is replaced by the median. If homoscedastic is TRUE, the median is taken over all profiles, else one value per profile is computed. }
  \item{var}{If model=2 (Normal Homoscedastic), the vector of values of the variance. If the user does not give his own hyperparameters, the package uses Hall's estimator whith d=4. If homoscedastic is TRUE, the mean of the estimate over all profile is used, else one value per profile is computed.}
	
	 \item{homoscedastic}{If model=2 (Normal Homoscedastic) or model=3, indicates whether the fixed parameter (variance or overidspersion) is common for all profiles or is profile-specific.}
	 
	 \item{unif}{A boolean stating whether prior on segmentation is uniform given number of segments. If false, then the prior favors segmentation with segments of equal length, i.e. n_r is proportional to the inverse of segment length.}
}

\details{This function is used to compute the matrix of segment probabilities assuming data is poisson, normal or negative binomial distributed. The probability of each interval being divided in k segments (k in 1 to Kmax) is computed. }

\value{
	An object of class "EBSProfiles".  
	\item{model}{Emission distribution (Poisson, Normal Homoscedastic, Negative Binomial or Normal Heteroscedastic)}
	\item{length}{the length of each profile}
	\item{NbConditions}{the number of profiles}
	\item{K}{the maximum number of segments for the segmentation for each profile}
	\item{HyperParameters}{The hyperparameters used for the prior on the data distribution for each profile}
	\item{Variance}{the vector of variances if model is Normal Homoscedastic}
	\item{overdispersion}{the vector of overdispersions if model is negative Binomial}
	\item{Li}{a list (one element per profile) of matrix of size Kmax*(length+1). Element [i,j] is the log-probability of interval [1,j[ being segmented in j segments}
	\item{Col}{a list (one element per profile) of matrix of size (length+1)*Kmax. Element [i,j] is the log-probability of interval [i,n] being segmented in i segments}
	\item{P}{a list (one element per profile) of matrix of size (length+1)*(length+1). Element [i,j] is the log-probability of interval [i,j[}
}

\references{
Rigaill, Lebarbier & Robin: Exact posterior distributions over the segmentation space and model selection for multiple change-point detection problems \url{Arxiv:1004.4347}

Johnson, Kotz & Kemp: Univariate Discrete Distributions

Hall, Kay & Titterington: Asymptotically optimal difference-based estimation of variance in non-parametric regression
}

\author{Alice Cleynen}


\seealso{
\code{\link{EBSegmentation}}}
\examples{
# changes for Poisson model
set.seed(1)
x1<-c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
x2<-c(rpois(125,3),rpois(75,4),rpois(75,1),rpois(125,8))
M<-rbind(x1,x2)
out <- EBSProfiles(M,model=1,K=10)
}

