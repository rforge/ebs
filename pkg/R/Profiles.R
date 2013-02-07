EBSProfiles <- function(data=numeric(), model=1, K = 3, hyper = numeric(), theta = numeric(), var = numeric(), homoscedastic = FALSE) UseMethod("EBSProfiles")

EBSProfiles.default <-function(data = numeric(), model=1, K=3, hyper = numeric(), theta = numeric(), var = numeric(), homoscedastic = FALSE)
{
  if ((model!=1)&(model!=2)&(model!=3)&(model!=4))
    stop("Choose model=1 (Poisson), 2 (Normal Homoscedastic), 3 (Negative Binomial) or 4 (Normal Heteroscedastic)")
  if (length(data)==0)
    stop("Give me a matrix of profiles to segment")
  n=ncol(data)
 	NbConditions<-length(data)/n
  hyp4<-rep(0,(4*NbConditions))
  
  if (length(K)==1)
  	K<-rep(K,NbConditions)
  if (length(K)==0)
  	K<-rep(3,NbConditions)
  if (length(K)!=NbConditions)
  	stop("Give me either one K per profile, or one common for all")
  	
  if (length(theta)==1)
  	theta<-rep(theta,NbConditions)
  if (length(var)==1)
  	var<-rep(var,NbConditions)

  if ((model==2)&(length(var)==0))
  {
  	if (homoscedastic)
  		var<-0 else
			var<-rep(0,NbConditions)
		for (i in 1:NbConditions)
		{
			data1<-data[i,-(1:3)]
			data2<-data[i,-c(1,2,n)]
			data3<-data[i,-c(1,n-1,n)]
			data4<-data[i,-c(n-2,n-1,n)]
			d<-c(0.1942, 0.2809, 0.3832, -0.8582)
			v2<-d[1]*data1+d[2]*data2+d[3]*data3+d[4]*data4
			v2<-v2*v2
			if (homoscedastic)
  			var<-var+sum(v2)/(4*n-12) else
				var[i]<-sum(v2)/(n-3)
	  }
	  if (homoscedastic)
	  	var<-rep(var, NbConditions)
  }

  if ((model==3)&(length(theta)==0))
  {
  	if (homoscedastic)
  		tall<-NULL else
  		theta<-rep(0,NbConditions)
  	for(l in 1:NbConditions)
  	{
		  i<-1
		  pold<-0
		  while (i<n)
		  {
		    pnew=0
		    while((i<n) & (data[l,i]!=0)) i<-i+1
		    while((i<n) & (data[l,i]==0) )
		    { 
		      pnew=pnew+1
		      i<-i+1
		    }
		    if(pnew>pold)
		      pold=pnew
		  }
		  h2<-max(2*pold,15)
			Xcum = cumsum(data[l,])
			X2cum = cumsum(data[l,]^2)
			MA = (Xcum[h2:n] - c(0, Xcum[1:(n-h2)])) / h2
			S2 = (X2cum[h2:n] - c(0, X2cum[1:(n-h2)])) / (h2-1) - h2/(h2-1)*MA^2
			K2 = MA^2 / (S2-MA)
			if (homoscedastic)
			{
				tall<-c(tall,K2)
			} else
			{
		  	theta[l] = median(K2)
		  	while ((theta[l]<0)&(h2<(n/2)))
				{
					h2<-2*h2
					MA = (Xcum[h2:n] - c(0, Xcum[1:(n-h2)])) / h2
					S2 = (X2cum[h2:n] - c(0, X2cum[1:(n-h2)])) / (h2-1) - h2/(h2-1)*MA^2
					K2 = MA^2 / (S2-MA)
					theta[l] = median(K2)   	
				}
		 	}
		  
  	}
  	if (homoscedastic)
  		theta<-rep(median(tall),NbConditions)
  }

  if ((model==4) & (length(hyper)==0))
  	for(i in 1:NbConditions)
  	{
    	me<-median(data[1,])
    	int<-abs(data[1,]-me)
    	OK<-which(int!=0)
    	inverse<-1/int[OK]
    	y<-fitdistr(inverse,"gamma")
    	hyp4[4*(i-1)+2]<-1
    	hyp4[4*(i-1)+3]<-y$estimate[1]
    	hyp4[4*(i-1)+4]<-y$estimate[2]
    }

  if(length(hyper)==0)
    if(model==1)
			hyper=rep(c(1,1),NbConditions) else if(model==3)
			hyper=rep(c(1/2,1/2),NbConditions) else if(model==2)
			hyper=rep(c(0,1),NbConditions) else
			hyper=hyp4
  hyper=as.vector(hyper)
  
  if (((model==1)|(model==2)|(model==3)) & (length(hyper)!=(2*NbConditions)))
    stop("for Poisson, Normal Homoscedastic and Negative Binomial models, two hyper-parameters are needed for each condition")
  if ((model==4) & (length(hyper)!=(4*NbConditions)))
    stop("for Normal Heteroscedastic model four hyper-parameters are needed for each condition")
    
  Lin<-list()
  Coln<-list()
  Pn<-list()
  if ((model==3) &(length(theta)!=NbConditions))
  	stop("Give me either one theta per profile, or one common for all")
  if ((model==2) &(length(var)!=NbConditions))
  	stop("Give me either one variance per profile, or one common for all")
	
	for (i in 1:NbConditions)
	{
		Km<-K[i]
		Li=matrix(0,nrow=K[i],ncol=(n+1))
		Li = as.vector(Li)
		Col=matrix(0,ncol=K[i],nrow=(n+1))
		Col = as.vector(Col)
		P=matrix(0,nrow=(n+1),ncol=(n+1))
		P=as.vector(P)
		M=as.vector(data[i,])		
		if (model==4)
			hyp<-c(hyper[4*(i-1)+1],hyper[4*(i-1)+2],hyper[4*(i-1)+3],hyper[4*(i-1)+4])  else hyp<-c(hyper[2*(i-1)+1],hyper[2*(i-1)+2]) 
		if (model==3) thei<-theta[i]
		if (model==2)	vari<-var[i]

		if (model==1)
    	Rep<-.C("SegmentPoisson", Size = as.integer(n),KMax = as.integer(Km), hyper = as.double(hyp), Data = as.integer(M), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS") else if (model==3)
    	Rep<-.C("SegmentBinNeg", Size = as.integer(n),KMax = as.integer(Km), hyper = as.double(hyp), theta = as.double(thei), Data = as.integer(M), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS") else if (model==2)
    	Rep<-.C("SegmentGaussienneHomo", Size = as.integer(n),KMax = as.integer(Km), hyper = as.double(hyp), Var = as.double(vari), Data = as.double(M), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS") else if (model==4)
    	Rep<-.C("SegmentGaussienne", Size = as.integer(n),KMax = as.integer(Km), hyper = as.double(hyp), Data = as.double(M), Col = as.double(Col), Li = as.double(Li), P = as.double(P),  PACKAGE="EBS")


		Lin[[i]]=t(matrix(Rep$Li,ncol=Km))
		Coln[[i]]=matrix(Rep$Col,ncol=Km)
		Pn[[i]]=t(matrix(Rep$P,ncol=(n+1)))
	}
	
	Datasets<-data
	rnames<-NULL
	for (i in 1:NbConditions)
	{
		rnames<-c(rnames,paste("Profile ",i,sep=""))
	}
	row.names(Datasets)<-rnames
  if (model==1) 
  {
      model.dist="Poisson"
      EBSProfiles.res=new("EBSProfiles", model=model.dist, data=Datasets, length=n, NbConditions = NbConditions, K=K, HyperParameters = hyper, Li=Lin, Col=Coln, P=Pn)
  }

  if (model==3) 
  {
      model.dist="Negative Binomial"
      EBSProfiles.res=new("EBSProfiles", model=model.dist, data=Datasets, length=n, NbConditions = NbConditions, K=K, HyperParameters = hyper, overdispersion = theta, Li=Lin, Col=Coln, P=Pn)
  }

  if (model==2) 
  {
      model.dist="Normal Homoscedastic"
      EBSProfiles.res=new("EBSProfiles", model=model.dist, data=Datasets, length=n, NbConditions = NbConditions, K=K, HyperParameters = hyper, Variance = var, Li=Lin, Col=Coln, P=Pn)
  }

  if (model==4) 
  {
      model.dist="Normal Heteroscedastic"
      EBSProfiles.res=new("EBSProfiles", model=model.dist, data=Datasets, length=n, NbConditions = NbConditions, K=K, HyperParameters = hyper, Li=Lin, Col=Coln, P=Pn)
  }

  EBSProfiles.res
}

##

GetCondition<-function(x, Condition = numeric()) UseMethod("GetCondition")
GetCondition.default<-function(x, Condition = numeric())
{
	if (length(Condition)!=1)
		stop('I need one and only one number of profile to return')
 	Lii<-getLi(x)[[Condition]]
  Coli<-getCol(x)[[Condition]]
  Pi<-getP(x)[[Condition]]
  if (Model(x)=="Poisson")
  {
  	hyper<-c(HyperParameters(x)[2*(Condition-1)+1],HyperParameters(x)[2*(Condition-1)+2])
		GetCondition.res=new("EBS", model=Model(x), data=Data(x)[Condition,], length=getLength(x), Kmax=getK(x)[Condition], HyperParameters=hyper, Li=Lii, Col= Coli, matProba=Pi)
	} else if (Model(x)=="Negative Binomial")
	{
  	hyper<-c(HyperParameters(x)[2*(Condition-1)+1],HyperParameters(x)[2*(Condition-1)+2])
  	theta<-Overdispersion(x)[Condition]
		GetCondition.res=new("EBS", model = Model(x), data = Data(x)[Condition,], length = getLength(x), Kmax = getK(x)[Condition], HyperParameters = hyper, overdispersion = theta, Li = Lii, Col = Coli, matProba = Pi)
	} else if (Model(x)=="Normal Homoscedastic")
	{
  	hyper<-c(HyperParameters(x)[2*(Condition-1)+1],HyperParameters(x)[2*(Condition-1)+2])
  	var<-Variance(x)[Condition]
		GetCondition.res=new("EBS", model = Model(x), data = Data(x)[Condition,], length = getLength(x), Kmax = getK(x)[Condition], HyperParameters = hyper, Variance = var, Li = Lii, Col = Coli, matProba = Pi)
	} else
	{
		hyper<- c(HyperParameters(x)[4*(Condition-1)+1], HyperParameters(x)[4*(Condition-1)+2], HyperParameters(x)[4*(Condition-1)+3], HyperParameters(x)[4*(Condition-1)+4])
		GetCondition.res=new("EBS", model = Model(x), data = Data(x)[Condition,], length = getLength(x), Kmax = getK(x)[Condition], HyperParameters = hyper, Li = Lii, Col = Coli, matProba = Pi)  	
  }
	GetCondition.res
}

##

EBSStatistic <- function(x, Conditions = numeric(), Tau = numeric(), K = numeric(),a0=1) UseMethod("EBSStatistic")
EBSStatistic.default <- function(x, Conditions = numeric(), Tau = numeric(), K = numeric(),a0=1)
{
	if (class(x)!="EBSProfiles")
		stop('x must be an object of class EBSProfiles')
	I<-length(Conditions)
	if (I==0)
		stop('please specify which profiles to compare')
	if (I>NbConditions(x))
		stop('Must compare a maximum of NbConditions(x) profiles')
	if (length(K)==1)
		K<-rep(K,I)
	if (length(K)==0)
		K<-getK(x)[Conditions]
	if (length(Tau)==1)
		Tau<-rep(Tau,I)
	if (length(Tau)==0)
		Tau<-rep(1,I)
	if (length(K)!=I)
		stop('need a value of K per profile')
	if (length(Tau)!=I)
		stop('need a change-point number per profile')
	for (i in 1:I)
		if(K[i]>getK(x)[Conditions[i]])
			stop('number of segments can not be larger than that given to function EBSProfiles')
	for (i in 1:I)
		if(Tau[i]>=K[i])
			stop('number of change-point has to be strictly inferior to number of segment in profile')
	
	n<-getLength(x)		
	y<-EBSDistrib(GetCondition(x,Conditions[1]),Tau[1],K[1])
	for (i in 2:I)	
		y<-y*EBSDistrib(GetCondition(x,Conditions[i]),Tau[i],K[i])
	y0<-sum(y)
	a1<-Geta1(n,Tau,K,a0)
	EBSStatistic.res<-a0*y0/((a0-a1)*y0+a1)
	EBSStatistic.res
}


###

CompCredibility<-function(x, Conditions, Tau = numeric(), K = numeric()) UseMethod("CompCredibility")
CompCredibility.default<-function(x, Conditions, Tau = numeric(), K = numeric()) 
{
	if (class(x)!="EBSProfiles")
		stop('object x must be of class EBSProfiles')
	if (length(Conditions)!=2)
		stop('Only two conditions can be compared with credibility intervals')
	if (max(Conditions)>NbConditions(x))
		stop('at least one condition you which to compare is not contained in x')
	if (length(K)==1)
		K<-rep(K,2)
	if (length(K)==0)
		K<-rep(3,2)
	if (length(Tau)==1)
		Tau<-rep(Tau,2)
	if (length(Tau)==0)
		Tau<-rep(1,2)
	if (length(K)!=2)
		stop('need a value of K per profile')
	if (length(Tau)!=2)
		stop('need a change-point number per profile')
	for (i in 1:2)
		if(K[i]>getK(x)[Conditions[i]])
			stop('number of segments can not be larger than that given to function EBSProfiles')
	for (i in 1:2)
		if(Tau[i]>=K[i])
			stop('number of change-point has to be strictly inferior to number of segment in profile')
	y1<-EBSDistrib(GetCondition(x,Conditions[1]),Tau[1],K[1])
	y2<-EBSDistrib(GetCondition(x,Conditions[2]),Tau[2],K[2])
	y<-NULL
	for (d in (2-getLength(x)):0)
	{
	  rab<-rep(0,abs(d))
	  a1<-c(rab,y1)
	  b1<-c(y2,rab)
	  temp1<-sum(a1*b1)
	  y<-rbind(y,c(d,temp1))
	}
	for (d in 1:(getLength(x)-2))
	{
	  rab<-rep(0,abs(d))
	  a1<-c(rab,y2)
	  b1<-c(y1,rab)
	  temp1<-sum(a1*b1)
	  y<-rbind(y,c(d,temp1))
	}
	aux<-sort(y[,2],decreasing=TRUE,index.return=TRUE)
	a<-sum(aux$x[1:which(aux$ix==(getLength(x)-1))])
	b<-a-aux$x[which(aux$ix==(getLength(x)-1))]
	cred<-list(Distribution=y,masswith0=a,massto0=b)
	class(cred)<-"Credibility"	
	cred			
}

print.Credibility<-function(x, ...)
{
	cat('Distribution of the difference Tau_i - Tau_j : ($Distribution) \n')
	str(x$Distribution)
	cat('Mass of credibility interval reaching and including 0 : ($masswith0) \n')
	print(x$masswith0)
	cat('Mass of credibility interval before reaching 0 : ($massto0) \n')
	print(x$massto0)
}

###

EBSICLProfiles<-function(x, prior=numeric()) UseMethod("EBSICLProfiles")
EBSICLProfiles.default<-function(x, prior=numeric())
{
	if (class(x)!="EBSProfiles")
		stop('object x must be of class EBSProfiles')
	I<-NbConditions(x)
	KMAX=max(getK(x))
	if (length(prior)==0)
	{
		prior<-matrix(0,nrow=I,ncol=KMAX)
		for (i in 1:I)
			for (k in 1:getK(x)[i])
				prior[i,k]<-1/getK(x)[i]
	}
	if (nrow(prior)!=I)
		stop('Need a prior on the number of segments for each profile')
	icl<-list()
	K<-rep(0,I)
	for (i in 1:I)
	{
		aux<-EBSICL(GetCondition(x,i),prior=prior[i,])
		icl[[i]]<-aux$ICL
		K[i]<-aux$NbICL
	}
	res<-list(ICL=icl, NbICL=K)
	res
}



###

EBSPlotProbaProfiles<-function(x,K=numeric(),data=FALSE) UseMethod("EBSPlotProbaProfiles")
EBSPlotProbaProfiles.default<-function(x,K=numeric(),data=FALSE)
{
	if(class(x)!="EBSProfiles")
		stop('object x must be of class EBSProfiles')
	I<-NbConditions(x)
	if (length(K)==0)
		K<-EBSICLProfiles(x)$NbICL
	if (length(K)==1)
		K<-rep(K,I)
	for (i in 1:I)
	{
  	if(K[i]<2)
    	stop("K has to be >1")
  	if(K[i]>getK(x)[i])
    	stop("I only know the segmentation up to Kmax, chose K<=Kmax")
  }
	par(mfrow=c(I,1))
	for (j in 1:I)
	{
	  y<-list()
  	a<-rep(0,(K[j]-1))
	  for (i in 1:(K[j]-1))
  	{
  	  y[[i]]<-EBSDistrib(GetCondition(x,j),i,K[j])
  	  a[i] = max(y[[i]])
  	}
  	b= max(a)
	  if (data)
  	{
      par(ann=FALSE)
      plot(Data(x)[j,], pch=1)
      par(ann=FALSE,new=TRUE)
      plot.default(y[[1]], type='l', col='blue',ylim=c(0,b),axes=FALSE)
      axis(4,col='blue')
      if (K[j]>2)
        for (i in 2:(K[j]-1))
          lines(y[[i]],col='blue')
    } else
  	{
    	plot(y[[1]], type='l', ylim=c(0,b), col='blue')
      if (K[j]>2)
        for (i in 2:(K[j]-1))
          lines(y[[i]],col='blue')
    }
  } 
}
######



CardMK<-function(n,K)
{
	res<-1
	I<-length(K)
	for (i in 1:I)
		res=res*choose(n-1,K[i]-1)
	return(res)
}

PriorDistrib<-function(n,k,K)
{
	if (k>=K)
		stop('k must be smaller than K')
		I<-length(K)
		MK<-CardMK(n,K)
		t<-(k+1):(n-(K-1-k))
	dis<-c(rep(0,k),choose(t-2,k-1)*choose(n-t-1,K-k-1),rep(0,K-1-k))/MK
	dis
}

Geta1<-function(n,k,K,a0)
{
	if (length(k)!=length(K))
		stop('k and K must be vectors of same length')
	I<-length(K)
	y<-rep(1,n)
	for (i in 1:I)
		y<-y*PriorDistrib(n,k[i],K[i])
	card<-sum(y)
	a1<-(1-a0*card)/(1-card)
	a1
}

