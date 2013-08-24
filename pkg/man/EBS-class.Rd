\name{EBS-class}
\Rdversion{1.1}
\docType{class}
\alias{EBS-class}
\alias{EBS-method, show}
\alias{EBS-method,getModel}
\alias{EBS-method,getLi}
\alias{EBS-method,getCol}
\alias{EBS-method,getLength}
\alias{EBS-method,getP}
\alias{EBS-method,getOverdispersion}
\alias{EBS-method,getData}
\alias{EBS-method,getKmax}
\alias{EBS-method,getHyperParameters}
\alias{EBS-method,getVariance}
\alias{EBS-method,getPriorm}

\title{Class "EBS"}
\description{
A class for Bayesian Segmentation objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EBS", ...)}.
  \describe{
    \item{\code{new("EBS", ...)}:}{creates a new object with class EBS }
  }
}

         
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"character"}, the assumed distribution of the data }
    \item{\code{data}:}{Object of class \code{"numeric"}, the data to be segmented }
    \item{\code{length}:}{Object of class \code{"numeric"}, the length of the profile }
    \item{\code{Kmax}:}{Object of class \code{"numeric"}, the maximum number of segments considered for the segmentation}
    \item{\code{HyperParameters}:}{Object of class \code{"numeric"}, the hyperparameters used for the prior distribution on the model parameters}
    \item{\code{Variance}:}{Object of class \code{"numeric"}, if model = Normal Homoscedastic, the variance used in the analysis}
    \item{\code{overdispersion}:}{Object of class \code{"numeric"}, if model = Negative Binomial, the overdispersion used in the analysis}
    \item{\code{Li}:}{Object of class \code{"numeric"}, matrix of size Kmax*(length+1) where element [i,j] is the log-probability of interval [1,j[ being segmented in i segments}
    \item{\code{Col}:}{Object of class \code{"numeric"}, a matrix of size (length+1)*Kmax where element [i,j] is the log-probability of interval [i,n] being segmented in j segments}
    \item{\code{matProba}:}{a matrix of size (length+1)*(length+1)where element [i,j] is the log-probability of interval [i,j[}
    \item{\code{unif}:}{a boolean stating whether prior on segmentation is uniform given number of segments}
  }
}
\section{Methods}{
  \describe{
    \item{getModel}{\code{signature(object = "EBS")}: retrieves model slot }
    \item{getData}{\code{signature(object = "EBS")}: retrieves data slot }
    \item{getLength}{\code{signature(object = "EBS")}: retrieves length slot }
    \item{getKmax}{\code{signature(object = "EBS")}: retrieves Kmax slot }
    \item{getHyperParameters}{\code{signature(object = "EBS")}: retrieves HyperParameters slot }
    \item{getVariance}{\code{signature(object = "EBS")}: retrieves Variance slot }
    \item{getOverdispersion}{\code{signature(object = "EBS")}: retrieves overdispersion slot }
    \item{getLi}{\code{signature(object = "EBS")}: retrieves Li slot }
    \item{getCol}{\code{signature(object = "EBS")}: retrieves Col slot }
    \item{getP}{\code{signature(object = "EBS")}: retrieves matProba slot }
    \item{getPriorm}{\code{signature(object = "EBS")}: retrieves unif slot }
  }
}

\author{
Alice Cleynen
}

\seealso{
\code{\link{EBSegmentation}}
}
\examples{
showClass("EBS") # shows the structure of the cpt class

data<-c(rpois(100,2),rpois(100,5))
# creates a new EBS object containing the segmentation of x
E<-EBSegmentation(data) 
class(E) # verifies the class of E
getModel(E) # retrieves model of the segmentation
getLength(E) # retrieves the length of the signal
}
\keyword{classes}
\keyword{EBS}
