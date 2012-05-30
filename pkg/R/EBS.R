EBSegmentation <- function(data=numeric(), model=1, Kmax = 15, hyper = numeric(), theta=0.15, uniform=TRUE) UseMethod("EBSegmentation")

EBSegmentation.default <-function(data=numeric(), model=1, Kmax = 15, hyper = numeric(), theta=0.15, uniform=TRUE)
{
  if ((model!=1)&(model!=2)&(model!=3))
    stop("Choose model=1 (Poisson), 2 (normal) or 3 (Negative Binomial)")
  if (length(data)==0)
    stop("Give me a vector of data to segment")

  if(uniform)
    if(model==1)
	hyper=c(2,2) else if(model==3)
	hyper=c(1,1) else
	hyper=c(2,0,2,1)
  
  hyper=as.vector(hyper)

  if (((model==1)|(model==3)) & (length(hyper)!=2))
    stop("for Poisson or Negative Binomial models two hyper-parameters are needed")
  if ((model==2) & (length(hyper)!=4))
    stop("for Normal model four hyper-parameters are needed")

  n=length(data)
  Li=matrix(0,nrow=Kmax,ncol=(n+1))
  Li = as.vector(Li)
  Col=matrix(0,ncol=Kmax,nrow=(n+1))
  Col = as.vector(Col)
  P=matrix(0,nrow=(n+1),ncol=(n+1))
  P=as.vector(P)

  if (model==1)
    Rep<-.C("SegmentPoisson", Size = as.integer(n),KMax = as.integer(Kmax), hyper = as.double(hyper), Data = as.integer(data), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS") else if (model==3)
    Rep<-.C("SegmentBinNeg", Size = as.integer(n),KMax = as.integer(Kmax), hyper = as.double(hyper), theta = as.double(theta), Data = as.integer(data), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS") else if (model==2)
    Rep<-.C("SegmentGaussienne", Size = as.integer(n),KMax = as.integer(Kmax), hyper = as.double(hyper), Data = as.integer(data), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS")

  resLi=matrix(Rep$Li,ncol=Kmax)
  resLi<-t(resLi)
  resCol=matrix(Rep$Col,ncol=Kmax)
  resP=matrix(Rep$P,ncol=(n+1))
  resP = t(resP)
  if (model==1) 
  {
      model.dist="Poisson"
      EBSegmentation.res=list(model=model.dist, data=data, length=n, Kmax=Kmax, HyperParameters = hyper, Li=resLi, Col=resCol, matProba = resP)
  }

  if (model==3) 
  {
      model.dist="Negative Binomial"
      EBSegmentation.res=list(model=model.dist, data=data, length=n, Kmax=Kmax, HyperParameters = hyper, overdispersion = theta, Li=resLi, Col=resCol, matProba = resP)
  }

  if (model==2) 
  {
      model.dist="Normal"
      EBSegmentation.res=list(model=model.dist, data=data, length=n, Kmax=Kmax, HyperParameters = hyper, Li=resLi, Col=resCol, matProba = resP)
  }

  class(EBSegmentation.res)= "EBS"
  EBSegmentation.res
}

EBSDistrib<-function(x, k, Kk,...) UseMethod("EBSDistrib")

EBSDistrib.default<-function(x, k, Kk,...)
{
  if(k<1)
    stop("k has to be >0")
  if(Kk>x$Kmax)
    stop("I only know the segmentation up to Kmax, chose K<=Kmax")

  DataLi=as.vector(t(x$Li))
  DataCol=as.vector(x$Col)
  n=x$length
  Dis=rep(0,n+1)
  Dis=as.vector(Dis)

  Rep<-.C("Distribution",Siz = as.integer(n+1), k=as.integer(k), Kk=as.integer(Kk), Col=as.double(DataCol), Li=as.double(DataLi), Dist=as.double(Dis))

  EBSDistrib.res = Rep$Dist
  EBSDistrib.res
}

EBSICL<-function(x,...) UseMethod("EBSICL")

EBSICL.default<-function(x,...)
{

  DataLi=as.vector(t(x$Li))
  DataCol=as.vector(x$Col)
  DataP=as.vector(t(x$matProba))
  n=x$length
  Kmax=x$Kmax
  k=0

  Rep<-.C("ChooseICL",Siz = as.integer(n+1), Kmax=as.integer(Kmax), Col=as.double(DataCol), Li=as.double(DataLi), P = as.double(DataP), kICL = as.integer(k))

  EBSICL.res = Rep$kICL
  EBSICL.res
}

EBSBIC<-function(x,...) UseMethod("EBSBIC")

EBSBIC.default<-function(x,...)
{

  DataCol=as.vector(x$Col)
  n=x$length
  Kmax=x$Kmax
  k=0

  Rep<-.C("ChooseBIC",Siz = as.integer(n+1), Kmax=as.integer(Kmax), Col=as.double(DataCol), kBIC = as.integer(k))

  EBSBIC.res = Rep$kBIC
  EBSBIC.res
}

EBSPostK<-function(x,...) UseMethod("EBSPostK")

EBSPostK.default<-function(x,...)
{

  DataCol=as.vector(x$Col)
  n=x$length
  Kmax=x$Kmax
  post=as.vector(rep(0,Kmax))

  Rep<-.C("PostK",Siz = as.integer(n+1), Kmax=as.integer(Kmax), Col=as.double(DataCol), Post = as.double(post))

  EBSPostK.res = Rep$Post
  EBSPostK.res
}

EBSPlotProba<-function(x,k,data=FALSE, file=character(), type='pdf',...) UseMethod("EBSPlotProba")

EBSPlotProba.default<-function(x,k,data=FALSE, file=character(), type='pdf',...)
{

  if(k<2)
    stop("k has to be >1")
  if(k>x$Kmax)
    stop("I only know the segmentation up to Kmax, chose k<=Kmax")

  y<-list()
  a<-rep(0,(k-1))
  for (i in 1:(k-1))
  {
    y[[i]]<-EBSDistrib(x,i,k)
    a[i] = max(y[[i]])
  }
  b= max(a)

  if (data)
  {
    if(length(file)==0)
    {
      par(ann=FALSE)
      plot(x$data, pch=1)
      par(ann=FALSE,new=TRUE)
      plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
      axis(4,col='blue')
      if (k>2)
        for (i in 2:(k-1))
          lines(y[[i]],col='blue')
    } else
    {
      if (type=='pdf')
      {
	pdf(file)
        par(ann=FALSE)
        plot(x$data, pch=1)
        par(ann=FALSE,new=TRUE)
        plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
        axis(4,col='blue')
        if (k>2)
          for (i in 2:(k-1))
            lines(y[[i]],col='blue')
	dev.off()
      }

      if (type=='png')
      {
	png(file)
        par(ann=FALSE)
        plot(x$data, pch=1)
        par(ann=FALSE,new=TRUE)
        plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
        axis(4,col='blue')
        if (k>2)
          for (i in 2:(k-1))
            lines(y[[i]],col='blue')
	dev.off()
      }
      if (type=='ps')
      {
	postscript(file)
        par(ann=FALSE)
        plot(x$data, pch=1)
        par(ann=FALSE,new=TRUE)
        plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
        axis(4,col='blue')
        if (k>2)
          for (i in 2:(k-1))
            lines(y[[i]],col='blue')
	dev.off()
      }
    }
  } else
  {
    if (length(file)==0)
    {
      plot(y[[1]], type='l', ylim=c(0,b), col='blue')
      if (k>2)
        for (i in 2:(k-1))
          lines(y[[i]],col='blue')
    } else
    {
      if (type=='pdf')
      {
	pdf(file)
        plot(y[[1]], type='l', ylim=c(0,b), col='blue')
        if (k>2)
          for (i in 2:(k-1))
            lines(y[[i]],col='blue')
	dev.off()
      }
      if (type=='png')
      {
	png(file)
        plot(y[[1]], type='l', ylim=c(0,b), col='blue')
        if (k>2)
          for (i in 2:(k-1))
            lines(y[[i]],col='blue')
	dev.off()
      }
      if (type=='ps')
      {
	postscript(file)
        plot(y[[1]], type='l', ylim=c(0,b), col='blue')
        if (k>2)
          for (i in 2:(k-1))
            lines(y[[i]],col='blue')
	dev.off()
      }
    }
  }
}

EBSPostMean<-function(x,k,...) UseMethod("EBSPostMean")

EBSPostMean.default<-function(x,k,...)
{

  if(k<2)
    stop("k has to be >1")
  if(k>x$Kmax)
    stop("I only know the segmentation up to Kmax, chose k<=Kmax")

  DataLi=as.vector(t(x$Li))
  DataCol=as.vector(x$Col)
  DataP=as.vector(t(x$matProba))
  Dat=as.vector(x$data)
  n=x$length
  Pos=as.vector(rep(0,n))

  Rep<-.C("GetPostMean",Siz = as.integer(n+1), Kk=as.integer(k), Data=as.double(Dat), Col=as.double(DataCol), Li=as.double(DataLi), P = as.double(DataP), Post = as.double(Pos))
  EBSPostMean.res = Rep$Post
  EBSPostMean.res
}


print.EBS <- function(x,...)
{
  cat("\n Model used for the segmentation: \n")
  print(x$model)
  cat("\n Length of data: \n")
  print(x$length)
  cat("\n data: \n")
  str(x$data)
  cat("\n Maximum number of segments considered for the segmentation \n")
  print(x$Kmax)
  cat("\n Hyper-parameters used for prior on data distribution: \n")
  if(x$model=="Poisson")
  {
    alpha = x$HyperParameters[1]
    beta = x$HyperParameters[2]
    cat("Gamma(alpha,beta): \n alpha=") 
    print(alpha[1])
    cat("\n beta= ") 
    print(beta)
    cat(" \n") 
  } else if (x$model=="Negative Binomial")
  {
    alpha = x$HyperParameters[1]
    beta = x$HyperParameters[2]
    cat("Beta(alpha, beta): \n alpha= ") 
    print(alpha[1])
    cat("\n beta= ") 
    print(beta)
    cat(" \n") 
  } else 
  {
    cat("for inverse of variance: Gamma(alpha, beta): \n 2 *alpha = ") 
    print(x$HyperParameters[1])
    cat("\n 2* beta= ") 
    print(x$HyperParameters[3])
    cat("\n for mean: Normal(mu,sigma): \n mu=") 
    print(x$HyperParameters[2])
    cat("\n variance/sigma = ") 
    print(x$HyperParameters[4])
    cat(" \n") 
  }
  if(x$model=="Negative Binomial")
  {
    cat("\n used value for inverse of overdispersion \n")
    print(x$overdispersion)
  }  
  cat("\n Log-proba [1,i[ in j segments: ($Li)")
  str(x$Li)
  cat("\n Log-proba [i,n+1[ in j segments: ($Col)")
  str(x$Col)
  cat("\n Log-proba [i,j[: ($matProba)")
  str(x$matProba)

}
