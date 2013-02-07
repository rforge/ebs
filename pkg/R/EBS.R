EBSegmentation <- function(data=numeric(), model=1, Kmax = 15, hyper = numeric(), theta = numeric(), var = numeric()) UseMethod("EBSegmentation")

EBSegmentation.default <-function(data=numeric(), model=1, Kmax = 15, hyper = numeric(), theta = numeric(), var = numeric())
{
  if ((model!=1)&(model!=2)&(model!=3)&(model!=4))
    stop("Choose model=1 (Poisson), 2 (Normal Homoscedastic), 3 (Negative Binomial) or 4 (Normal Heteroscedastic)")
  if (length(data)==0)
    stop("Give me a vector of data to segment")
  n=length(data)

  if ((model==2)&(length(var)==0))
  {
		data1<-data[-(1:3)]
		data2<-data[-c(1,2,n)]
		data3<-data[-c(1,n-1,n)]
		data4<-data[-c(n-2,n-1,n)]
		d<-c(0.1942, 0.2809, 0.3832, -0.8582)
		v2<-d[1]*data1+d[2]*data2+d[3]*data3+d[4]*data4
		v2<-v2*v2
		var<-sum(v2)/(n-3)
  }

  if ((model==3)&(length(theta)==0))
  {
    i<-1
    pold<-0
    while (i<n)
    {
      pnew=0
      while((i<n) & (data[i]!=0)) i<-i+1
      while((i<n) & (data[i]==0) )
      { 
        pnew=pnew+1
        i<-i+1
      }
      if(pnew>pold)
        pold=pnew
    }
    h<-max(2*pold,15)
		Xcum = cumsum(data)
		X2cum = cumsum(data^2)
		M = (Xcum[h:n] - c(0, Xcum[1:(n-h)])) / h
		S2 = (X2cum[h:n] - c(0, X2cum[1:(n-h)])) / (h-1) - h/(h-1)*M^2
		K = M^2 / (S2-M)
    theta = median(K)
    while ((theta<0)&(h<(n/2)))
    {
    	h<-2*h
			M = (Xcum[h:n] - c(0, Xcum[1:(n-h)])) / h
			S2 = (X2cum[h:n] - c(0, X2cum[1:(n-h)])) / (h-1) - h/(h-1)*M^2
			K = M^2 / (S2-M)
		  theta = median(K)   	
    }
  }

  if ((model==4) & (length(hyper)==0))
  {
		me<-median(data)
		int<-abs(data-me)
		OK<-which(int!=0)
		inverse<-1/int[OK]
		y<-fitdistr(inverse,"gamma")
  }

  if(length(hyper)==0)
    if(model==1)
			hyper=c(1,1) else if(model==3)
			hyper=c(1/2,1/2) else if(model==2)
			hyper=c(0,1) else
			hyper=c(0,1,y$estimate[1],y$estimate[2])
  
  hyper=as.vector(hyper)

  if (((model==1)|(model==2)|(model==3)) & (length(hyper)!=2))
    stop("for Poisson, Normal Homoscedastic and Negative Binomial models, two hyper-parameters are needed")
  if ((model==4) & (length(hyper)!=4))
    stop("for Normal Heteroscedastic model four hyper-parameters are needed")


  Li=matrix(0,nrow=Kmax,ncol=(n+1))
  Li = as.vector(Li)
  Col=matrix(0,ncol=Kmax,nrow=(n+1))
  Col = as.vector(Col)
  P=matrix(0,nrow=(n+1),ncol=(n+1))
  P=as.vector(P)

  if (model==1)
    Rep<-.C("SegmentPoisson", Size = as.integer(n),KMax = as.integer(Kmax), hyper = as.double(hyper), Data = as.integer(data), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS") else if (model==3)
    Rep<-.C("SegmentBinNeg", Size = as.integer(n),KMax = as.integer(Kmax), hyper = as.double(hyper), theta = as.double(theta), Data = as.integer(data), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS") else if (model==2)
    Rep<-.C("SegmentGaussienneHomo", Size = as.integer(n),KMax = as.integer(Kmax), hyper = as.double(hyper), Var = as.double(var), Data = as.double(data), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS") else if (model==4)
    Rep<-.C("SegmentGaussienne", Size = as.integer(n),KMax = as.integer(Kmax), hyper = as.double(hyper), Data = as.double(data), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS")

  resLi=matrix(Rep$Li,ncol=Kmax)
  resLi<-t(resLi)
  resCol=matrix(Rep$Col,ncol=Kmax)
  resP=matrix(Rep$P,ncol=(n+1))
  resP = t(resP)
  if (model==1) 
  {
      model.dist="Poisson"
      EBSegmentation.res=new("EBS", model=model.dist, data=data, length=n, Kmax=Kmax, HyperParameters = hyper, Li=resLi, Col=resCol, matProba = resP)
  }

  if (model==3) 
  {
      model.dist="Negative Binomial"
      EBSegmentation.res=new("EBS",model=model.dist, data=data, length=n, Kmax=Kmax, HyperParameters = hyper, overdispersion = theta, Li=resLi, Col=resCol, matProba = resP)
  }

  if (model==2) 
  {
      model.dist="Normal Homoscedastic"
      EBSegmentation.res=new("EBS",model=model.dist, data=data, length=n, Kmax=Kmax, HyperParameters = hyper, Variance = var, Li=resLi, Col=resCol, matProba = resP)
  }

  if (model==4) 
  {
      model.dist="Normal Heteroscedastic"
      EBSegmentation.res=new("EBS",model=model.dist, data=data, length=n, Kmax=Kmax, HyperParameters = hyper, Li=resLi, Col=resCol, matProba = resP)
  }

  EBSegmentation.res
}


EBSDistrib<-function(x, k, Kk) UseMethod("EBSDistrib")
EBSDistrib.default<-function(x, k, Kk)
{
	if(class(x)!="EBS")
		stop('object x must be of class EBS')
  if(k<1)
    stop("k has to be >0")
  if(Kk>getKmax(x))
    stop("I only know the segmentation up to Kmax, chose K<=Kmax")

	n=dataLength(x)
	t<-(k+1):(n-(Kk-k-1))
	Aux<-Li(x)[k,t]+Col(x)[t,Kk-k]
	ma<-max(Aux)
	Aux2<-exp(Aux-ma)
	EBSDistrib.res<-c(rep(0,k),Aux2/sum(Aux2),rep(0,Kk-k-1))
  EBSDistrib.res
}

EBSICL<-function(x, prior=numeric()) UseMethod("EBSICL")
EBSICL.default<-function(x, prior=numeric())
{
	if(class(x)!="EBS")
		stop('object x must be of class EBS')
  DataLi=as.vector(t(Li(x)))
  DataCol=as.vector(Col(x))
  DataP=as.vector(t(MatProba(x)))
  n=dataLength(x)
  Kmax=getKmax(x)
  k=0
  icl=as.vector(rep(0,Kmax))
  if (length(prior)==0)
    prior=rep(1/Kmax, Kmax)
  if (length(prior)!=Kmax)
    stop("Prior size must be equal to Kmax")
  if (sum(prior)!=1)
    stop("Sum of elements of prior on K must be equal to one")

  Rep<-.C("ChooseICL",Siz = as.integer(n+1), Kmax=as.integer(Kmax), PriorK=as.double(prior), Col=as.double(DataCol), Li=as.double(DataLi), P = as.double(DataP), ICL=as.double(icl), kICL = as.integer(k))

  EBSICL.res = list(ICL=Rep$ICL, NbICL=Rep$kICL)
  EBSICL.res
}

EBSBIC<-function(x, prior=numeric()) UseMethod("EBSBIC")
EBSBIC.default<-function(x, prior=numeric())
{

  DataCol=as.vector(Col(x))
  n=dataLength(x)
  Kmax=getKmax(x)
  k=0
  bic=as.vector(rep(0,Kmax))
  if (length(prior)==0)
    prior=rep(1/Kmax, Kmax)
  if (length(prior)!=Kmax)
    stop("Prior size must be equal to Kmax")
  if (sum(prior)!=1)
    stop("Sum of elements of prior on K must be equal to one")

  Rep<-.C("ChooseBIC",Siz = as.integer(n+1), Kmax=as.integer(Kmax), PriorK=as.double(prior), Col=as.double(DataCol), BIC = as.double(bic), kBIC = as.integer(k))

  EBSBIC.res = list(BIC=Rep$BIC, NbBIC=Rep$kBIC)
  EBSBIC.res
}

EBSPostK<-function(x, prior=numeric()) UseMethod("EBSPostK")
EBSPostK.default<-function(x, prior=numeric())
{
	if(class(x)!="EBS")
		stop('object x must be of class EBS')
  DataCol=as.vector(Col(x))
  n=dataLength(x)
  Kmax=getKmax(x)
  post=as.vector(rep(0,Kmax))
  if (length(prior)==0)
    prior=rep(1/Kmax, Kmax)
  if (length(prior)!=Kmax)
    stop("Prior size must be equal to Kmax")
  if (sum(prior)!=1)
    stop("Sum of elements of prior on K must be equal to one")

  EBSPostK.res = exp(-EBSBIC(x,prior)$BIC+min(EBSBIC(x,prior)$BIC))/sum(exp(-EBSBIC(x,prior)$BIC+min(EBSBIC(x,prior)$BIC)))
  EBSPostK.res
}

EBSPlotProba<-function(x,K,data=FALSE, file=character(), type='pdf') UseMethod("EBSPlotProba")
EBSPlotProba.default<-function(x,K,data=FALSE, file=character(), type='pdf')
{
	if(class(x)!="EBS")
		stop('object x must be of class EBS')
  if(K<2)
    stop("K has to be >1")
  if(K>getKmax(x))
    stop("I only know the segmentation up to Kmax, chose K<=Kmax")

  y<-list()
  a<-rep(0,(K-1))
  for (i in 1:(K-1))
  {
    y[[i]]<-EBSDistrib(x,i,K)
    a[i] = max(y[[i]])
  }
  b= max(a)

  if (data)
  {
    if(length(file)==0)
    {
      par(ann=FALSE)
      plot(getData(x), pch=1)
      par(ann=FALSE,new=TRUE)
      plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
      axis(4,col='blue')
      if (K>2)
        for (i in 2:(K-1))
          lines(y[[i]],col='blue')
    } else
    {
      if (type=='pdf')
      {
				pdf(file)
        par(ann=FALSE)
        plot(getData(x), pch=1)
        par(ann=FALSE,new=TRUE)
        plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
        axis(4,col='blue')
        if (K>2)
          for (i in 2:(K-1))
            lines(y[[i]],col='blue')
				dev.off()
      }

      if (type=='png')
      {
				png(file)
        par(ann=FALSE)
        plot(getData(x), pch=1)
        par(ann=FALSE,new=TRUE)
        plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
        axis(4,col='blue')
        if (K>2)
          for (i in 2:(K-1))
            lines(y[[i]],col='blue')
				dev.off()
      }
      if (type=='ps')
      {
				postscript(file)
        par(ann=FALSE)
        plot(getData(x), pch=1)
        par(ann=FALSE,new=TRUE)
        plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
        axis(4,col='blue')
        if (K>2)
          for (i in 2:(K-1))
            lines(y[[i]],col='blue')
				dev.off()
      }
    }
  } else
  {
    if (length(file)==0)
    {
      plot(y[[1]], type='l', ylim=c(0,b), col='blue')
      if (K>2)
        for (i in 2:(K-1))
          lines(y[[i]],col='blue')
    } else
    {
      if (type=='pdf')
      {
				pdf(file)
        plot(y[[1]], type='l', ylim=c(0,b), col='blue')
        if (K>2)
          for (i in 2:(K-1))
            lines(y[[i]],col='blue')
				dev.off()
      }
      if (type=='png')
      {
				png(file)
        plot(y[[1]], type='l', ylim=c(0,b), col='blue')
        if (K>2)
          for (i in 2:(K-1))
            lines(y[[i]],col='blue')
				dev.off()
      }
      if (type=='ps')
      {
				postscript(file)
        plot(y[[1]], type='l', ylim=c(0,b), col='blue')
        if (K>2)
          for (i in 2:(K-1))
            lines(y[[i]],col='blue')
				dev.off()
      }
    }
  }
}




TruncPois <-function(lambda,Kmax)
{
  Vector<-rep(0,Kmax)
  sum<-0
  for (i in 1:Kmax)
  {
    Vector[i]=exp(-lambda+i*log(lambda)-lgamma(i+1))
    sum = sum + Vector[i]
  }
  for (i in 1:Kmax)
    Vector[i] = Vector[i]/sum
  Vector    
}
