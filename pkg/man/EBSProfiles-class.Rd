\name{EBSProfiles-class}
\Rdversion{1.1}
\docType{class}
\alias{EBSProfiles-class}
\alias{EBSProfiles-method, show}
\alias{EBSProfiles-method,Model}
\alias{EBSProfiles-method,Li}
\alias{EBSProfiles-method,Col}
\alias{EBSProfiles-method,Length}
\alias{EBSProfiles-method,matProba}
\alias{EBSProfiles-method,Overdispersion}
\alias{EBSProfiles-method,Data}
\alias{EBSProfiles-method,Kmax}
\alias{EBSProfiles-method,HyperParameters}
\alias{EBSProfiles-method,Variance}
\alias{EBSProfiles-method,Priorm}


\title{Class "EBSProfiles"}
\description{
A class for Bayesian Segmentation and comparison of multiple profiles.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EBSProfiles", ...)}.
  \describe{
    \item{\code{new("EBSProfiles", ...)}:}{creates a new object with class EBSProfiles }
  }
}


         
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"character"}, the assumed (identic) distribution class of each profile }
    \item{\code{data}:}{Object of class \code{"numeric"}, the matrix of profiles (lines) to be segmented }
    \item{\code{length}:}{Object of class \code{"numeric"}, the length of each profile }
    \item{\code{NbConditions}:}{Object of class \code{"numeric"}, the number of profiles }
    \item{\code{K}:}{Object of class \code{"numeric"}, the vector of maximum number of segments considered for each profile}
    \item{\code{HyperParameters}:}{Object of class \code{"numeric"}, the vector of hyperparameters used for the prior distribution on the model parameters for each profile}
    \item{\code{Variance}:}{Object of class \code{"numeric"}, if model = Normal Homoscedastic, the vector of variances used in the analysis}
    \item{\code{overdispersion}:}{Object of class \code{"numeric"}, if model = Negative Binomial, the vector of overdispersions used in the analysis}
    \item{\code{Li}:}{Object of class \code{"numeric"}, a list of the matrix of each profile, of size Kmax*(length+1) where element [i,j] is the log-probability of interval [1,j[ being segmented in j segments}
    \item{\code{Col}:}{Object of class \code{"numeric"}, a list of the matrix of each profile, of size (length+1)*Kmax where element [i,j] is the log-probability of interval [i,n] being segmented in i segments}
    \item{\code{P}:}{a list of the matrix of each profile, of size (length+1)*(length+1) where element [i,j] is the log-probability of interval [i,j[}
    \item{\code{unif}:}{a boolean stating whether the prior on the segmentation is uniform given the number of segments.}
  }
}
\section{Methods}{
  \describe{
    \item{Model}{\code{signature(object = "EBSProfiles")}: retrieves model slot }
    \item{Data}{\code{signature(object = "EBSProfiles")}: retrieves data slot }
    \item{Length}{\code{signature(object = "EBSProfiles")}: retrieves length slot }
    \item{NbConditions}{\code{signature(object = "EBSProfiles")}: retrieves NbConditions slot }
    \item{Kmax}{\code{signature(object = "EBSProfiles")}: retrieves K slot }
    \item{HyperParameters}{\code{signature(object = "EBSProfiles")}: retrieves HyperParameters slot }
    \item{Variance}{\code{signature(object = "EBSProfiles")}: retrieves Variance slot }
    \item{Overdispersion}{\code{signature(object = "EBSProfiles")}: retrieves overdispersion slot }
    \item{Li}{\code{signature(object = "EBSProfiles")}: retrieves Li slot }
    \item{Col}{\code{signature(object = "EBSProfiles")}: retrieves Col slot }
    \item{matProba}{\code{signature(object = "EBSProfiles")}: retrieves P slot }
    \item{Priorm}{\code{signature(object = "EBSProfiles")}: retrieves unif slot }
  }
}

\author{
Alice Cleynen
}

\seealso{
\code{\link{EBSegmentation}},\code{\link{EBSProfiles}}, \code{\link{Classes}}
}
\examples{
showClass("EBS") # shows the structure of the cpt class

x1<-c(rpois(100,2),rpois(200,5))
x2<-c(rpois(100,3),rpois(150,8),rpois(50,2))
data<-rbind(x1,x2)
# creates a new EBSProfiles object containing the segmentation of 
#profiles x1 and x2
E<-EBSProfiles(data,K=c(2,3))
class(E) # verifies the class of E
Model(E) # retrieves model of the segmentation
# retrieves the maximal number of segments considered for profile x1
Kmax(E)[1]
}
\keyword{classes}
\keyword{EBSProfiles}
