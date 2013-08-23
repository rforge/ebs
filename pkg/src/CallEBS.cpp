#include "CallEBS.h"
#include "Constants.h"
#include "LogPoisson.h"
#include "LogBinNeg.h"
#include "LogGaussienne.h"
#include "LogGaussienneHomo.h"
#include "PriorUnif.h"
#include "PriorSize.h"



void CallEBSPoisson(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P, bool* u)
{
    int K= *KMax;
    int n = *Size;
    double Hyper1 = hyper[0];
    double Hyper2 = hyper[1];
    MyVector<int> MyData(n,0);
    for (int i=0; i<n; i++)
	  MyData[i] = Data[i];
    Observations<int> LesObservations(MyData);
    LogPoisson BN(LesObservations, Hyper1, Hyper2);
    Distributions<LogPoisson> Dis(n+1, K, &BN, u);
    for (int k=0; k<K; k++)
	  for (int i=0; i<(n+1); i++)
	  {
	    Li[k*(n+1)+i]=Dis.GetDataLi()[k][i];
	    Col[k*(n+1)+i]=Dis.GetDataCol()[i][k];
	  }

    for (int i=0; i<(n+1); i++)
	  for (int j=0; j<(n+1); j++)
	    P[i*(n+1)+j]=Dis.P.Data[i][j];

    return;
}

void CallPriorSize(int *Size, int *KMax, double* Col, double* Li, double* P, bool* u)
{
    int K= *KMax;
    int n = *Size;
    double r = 0;
   	PriorSize PS(r);
    Distributions<PriorSize> Dis(n+1, K, &PS, u);
    for (int k=0; k<K; k++)
	  for (int i=0; i<(n+1); i++)
	  {
	    Li[k*(n+1)+i]=Dis.GetDataLi()[k][i];
	    Col[k*(n+1)+i]=Dis.GetDataCol()[i][k];
	  }

    for (int i=0; i<(n+1); i++)
	  for (int j=0; j<(n+1); j++)
	    P[i*(n+1)+j]=Dis.P.Data[i][j];

    return;
}

void CallPriorUnif(int *Size, int *KMax, double* Col, double* Li, double* P, bool* u)
{
    int K= *KMax;
    int n = *Size;
    double r=0;
   	PriorUnif PS(r);
    Distributions<PriorUnif> Dis(n+1, K, &PS, u);
    for (int k=0; k<K; k++)
	  for (int i=0; i<(n+1); i++)
	  {
	    Li[k*(n+1)+i]=Dis.GetDataLi()[k][i];
	    Col[k*(n+1)+i]=Dis.GetDataCol()[i][k];
	  }

    for (int i=0; i<(n+1); i++)
	  for (int j=0; j<(n+1); j++)
	    P[i*(n+1)+j]=Dis.P.Data[i][j];

    return;
}

void CallEBSBinNeg(int *Size, int *KMax, double* hyper, double* theta, int* Data, double* Col, double* Li, double* P, bool* u)
{
    int K= *KMax;
    int n = *Size;
    double Hyper1 = hyper[0];
    double Hyper2 = hyper[1];
    double Theta = *theta;
    MyVector<int> MyData(n,0);
    for (int i=0; i<n; i++)
	  MyData[i] = Data[i];
    Observations<int> LesObservations(MyData);
    LogBinNeg BN(LesObservations, Theta, Hyper1, Hyper2);
    Distributions<LogBinNeg> Dis(n+1, K, &BN, u);
    for (int k=0; k<K; k++)
	  for (int i=0; i<(n+1); i++)
	  {
	    Li[k*(n+1)+i]=Dis.GetDataLi()[k][i];
	    Col[k*(n+1)+i]=Dis.GetDataCol()[i][k];
	  }

    for (int i=0; i<(n+1); i++)
	  for (int j=0; j<(n+1); j++)
	    P[i*(n+1)+j]=Dis.P.Data[i][j];
    return;
}

void CallEBSGaussienne(int *Size, int *KMax, double* hyper, double* Data, double* Col, double* Li, double* P, bool* u)
{
    int K= *KMax;
    int n = *Size;
    double Hyper1 = hyper[0];
    double Hyper2 = hyper[1];
    double Hyper3 = hyper[2];
    double Hyper4 = hyper[3];
    MyVector<double> MyData(n,0);
    for (int i=0; i<n; i++)
	MyData[i] = Data[i];
    Observations<double> LesObservations(MyData);
    LogGaussienne G(LesObservations, Hyper1, Hyper2, Hyper3, Hyper4);
    Distributions<LogGaussienne> Dis(n+1, K, &G, u);
    for (int k=0; k<K; k++)
	for (int i=0; i<(n+1); i++)
	{
	  Li[k*(n+1)+i]=Dis.GetDataLi()[k][i];
	  Col[k*(n+1)+i]=Dis.GetDataCol()[i][k];
	}

    for (int i=0; i<(n+1); i++)
	for (int j=0; j<(n+1); j++)
	  P[i*(n+1)+j]=Dis.P.Data[i][j];

    return;
}

void CallEBSGaussienneHomo(int *Size, int *KMax, double* hyper, double *Var, double* Data, double* Col, double* Li, double* P, bool* u)
{
    int K= *KMax;
    int n = *Size;
    double v = Var[0];
    double Hyper1 = hyper[0];
    double Hyper2 = hyper[1];
    MyVector<double> MyData(n,0);
    for (int i=0; i<n; i++)
	MyData[i] = Data[i];
    Observations<double> LesObservations(MyData);
    LogGaussienneHomo G(LesObservations, v, Hyper1, Hyper2);
    Distributions<LogGaussienneHomo> Dis(n+1, K, &G, u);
    for (int k=0; k<K; k++)
	for (int i=0; i<(n+1); i++)
	{
	  Li[k*(n+1)+i]=Dis.GetDataLi()[k][i];
	  Col[k*(n+1)+i]=Dis.GetDataCol()[i][k];
	}

    for (int i=0; i<(n+1); i++)
	for (int j=0; j<(n+1); j++)
	  P[i*(n+1)+j]=Dis.P.Data[i][j];

    return;
}

void BreakDistrib(int *Siz, int *kk, int*KK, double* Col, double* Li, double* Dist)
{
    int k = *kk;
    int Kk = *KK;
    int Size = *Siz;
    double* Aux = new double[Size];
    double Norm = 0;
    for (int i=0; i<Size; i++)
	Aux[i] = Li[(k-1)*Size+i] + Col[i+(Kk-k-1)*Size];
    double Max = Aux[0];
    for (int i = 0; i < Size; i++)
      if (Aux[i] > Max)
     	Max = Aux[i];
    for (int i=0; i<Size; i++)
	Aux[i] -= Max;
    for (int i=0; i<Size; i++)
    {
	Dist[i] = exp(Aux[i]);
	Norm += Dist[i];
    }
    for (int i=0; i<Size; i++)
		Dist[i] /= Norm;
    delete[] Aux;
    return;
}

void ProbaSegment(int *Siz, int*Kk, double* Col, double* Li, double *P, double* Pseg)
{

  int Size = *Siz;
  int k = *Kk;
  for (int j=1; j<(Size-k); j++)
	Pseg[j] = exp(P[j]-Col[(k-1)*Size]+Col[j+(k-2)*Size]);
  for (int i=(k-1); i<(Size-1); i++)
	Pseg[i*Size+Size-1] = exp(P[i*Size+Size-1]-Col[(k-1)*Size]+Li[(k-2)*Size+i]);
  if(k>2)
  {
    for (int i=1; i<(Size-1); i++)
      for (int j=i+1; j<(Size-1); j++)
      {
	double Aux = MINUS_INFINITY;
	for (int p=2; p<k; p++)
	  if ((i>p-2) & (i<Size-2-k+p) & (j<Size-k+p))
	    Aux = lsum(Aux,sumoflogs(Li[(p-2)*Size+i],Col[j+(k-p-1)*Size]));
	if ((i>0) & (j<(Size-1)))
	    Pseg[i*Size+j] = exp(P[i*Size+j]-Col[(k-1)*Size] + Aux);
      }
  }
  return;
}

void GroupSegment(int *Siz, int*Kk, double* Col, double* Li, double *P, double* Pseg) // pour calculer la matrice des S_K(i,j)*log(f(i,j))
{

  int Size = *Siz; // la taille de la matrice
  int k = *Kk; // la valeur de K-1 (attention decalage d'indice)
  for (int i=0; i<(Size*Size-1); i++)
  	Pseg[i] = 0;
  for (int j=1; j<(Size-k); j++) // les termes quand i=1 (zero en valeur c++) (necessairement le premier segment)
		Pseg[j] = P[j] * exp(P[j]-Col[(k-1)*Size]+Col[j+(k-2)*Size]);
  for (int i=(k-1); i<(Size-1); i++) // les termes quand j=n+1 (n en valeur c++) (necessairement le dernier segment)
		Pseg[i*Size+Size-1] = P[i*Size+Size-1] * exp(P[i*Size+Size-1]-Col[(k-1)*Size]+Li[(k-2)*Size+i]);
  if(k>2)
    for (int i=1; i<(Size-1); i++)
      for (int j=i+1; j<(Size-1); j++)
      {
				double Aux = MINUS_INFINITY;
				for (int p=2; p<k; p++) // les valeurs possibles du numero de segment
					if ((i>p-2) & (i<Size-2-k+p) & (j<Size-k+p))
						Aux = lsum(Aux,sumoflogs(Li[(p-2)*Size+i],Col[j+(k-p-1)*Size]));
				Pseg[i*Size+j] = P[i*Size+j] * exp(P[i*Size+j]-Col[(k-1)*Size] + Aux);
      }
  return;
}



void GetICL(int *Siz, int *Kmax, double* PriorK, double* Col, double* Li, double *P, double* ICL, int* kICL, bool* u)
{
  int K = *Kmax;
  int Size = *Siz;
  int Size2 = Size+1;
  ICL[0]=-Col[0]-PriorK[0];
  double* cons = new double[K];
  	for (int l=0; l<K; l++)
  		cons[l] = 0;
  if (!(*u))
  {
  	double* pr= new double[Size2*Size2]; 
			for (int l=0; l<(Size2*Size2); l++)
				pr[l] = 0;
		double* col = new double[Size2*K];
		double* li = new double[Size2*K];
			for (int l=0; l<(Size2*K); l++)
			{
				col[l] = 0;
				li[l] = 0;
			}
		CallPriorSize(Siz, Kmax, col, li, pr,u);
		for (int l=1; l<K; l++)
  		cons[l] = col[(l-1)*Size];
  	delete[] pr;
		delete[] col;
		delete[] li;
	}
  for (int k=1; k<K; k++)
  {
    int* kneeded= new int[1];
    kneeded[0]=k+1;
    double* PS= new double[Size*Size]; // pour chaque valeur de k, construction de la matrice S_k(i,j)*log(f(i,j))
    for (int l=0; l<(Size*Size); l++)
			PS[l] = 0;
    GroupSegment(Siz, kneeded, Col, Li, P, PS);
    
    if (k==1)
      for (int i=1; i<(Size-1); i++)
				ICL[k] -= PS[i] + PS[i*Size+Size-1];
    else
      for (int i=0; i<(Size-1); i++)
					for (int j=i+1; j<std::min(Size-1,Size-k+i); j++)
							ICL[k] -= PS[i*Size+j];
    delete[] PS;
    delete[] kneeded;
    if (!(*u))
    {
		  if (PriorK[k]!=0)
		    ICL[k] +=cons[k] -log(PriorK[k]);
		  else
		  	ICL[k] = PLUS_INFINITY;
    } else
    {
		  if (PriorK[k]!=0)
		    ICL[k] += Norma(0,Size-1,k+1)-log(PriorK[k]);
		  else
		  	ICL[k] = PLUS_INFINITY;    
    }
  }
  delete[] cons;
  int kchoose = 0;
  double mi = PLUS_INFINITY;
  for (int k=1; k<K; k++)
    if (ICL[k]<mi)
    {
			mi = ICL[k];
			kchoose = k+1;
    }
  kICL[0] = kchoose;
  return;
}

void Moyenne(int *Siz, double* Data, double* P)
{
  int Size = *Siz;
  for (int i=0; i<Size; i++)
    for (int j=0; j<Size; j++)
	P[i*Size+j]=0;
  for (int i=0; i<Size; i++)
    for (int j=i; j<Size; j++)
    {
	for (int l=0; l<(j-i); l++)
	  P[i*Size+j] += Data[i+l];
	P[i*Size+j]/=(j-i);
    }
  return;
}

void PosteriorK(int *Siz, int *Kmax, double* PriorK, double* Col, double* PostK, bool* u)
{
  int K = *Kmax;
  int Size = *Siz;
  int Size2 = Size +1;
  double* cons = new double[K];
  	for (int l=0; l<K; l++)
  		cons[l] = 0;
  if (!(*u))
  {
  	double* pr= new double[Size2*Size2]; 
			for (int l=0; l<(Size2*Size2); l++)
				pr[l] = 0;
		double* col = new double[Size2*K];
		double* li = new double[Size2*K];
			for (int l=0; l<(Size2*K); l++)
			{
				col[l] = 0;
				li[l] = 0;
			}
		CallPriorSize(Siz, Kmax, col, li, pr,u);
		for (int l=1; l<K; l++)
  		cons[l] = col[(l-1)*Size];
  	delete[] pr;
		delete[] col;
		delete[] li;
	}
  for (int k=0; k<K; k++)
    if (PriorK[k] != 0)
    {
    	if (!(*u))
      	PostK[k] = -Col[k*Size] + cons[k] -log(PriorK[k]);
    	else
      	PostK[k] = -Col[k*Size] + Norma(0,Size-1,k+1)-log(PriorK[k]);
    }
    else    
      PostK[k] = PLUS_INFINITY;
  double Total = 0;
  double Max = MINUS_INFINITY;
  for (int i = 0; i < K; i++)
    if ((PostK[i] < PLUS_INFINITY) & (PostK[i] > Max))
     	Max = PostK[i];
  for (int i=0; i<K; i++)
    if (PostK[i] < PLUS_INFINITY)
	PostK[i] -= Max;
  for (int k=0; k<K; k++)
    Total += exp(-PostK[k]);
  for (int k=0; k<K; k++)
    PostK[k] = exp(-PostK[k])/Total;
  delete[] cons;
  return;
}

void GetBIC(int *Siz, int *Kmax, double* PriorK, double* Col, double* BIC, int* kBIC, bool* u)
{
  int K = *Kmax;
  int Size = *Siz;
  int Size2 = Size +1;
  double* cons = new double[K];
  	for (int l=0; l<K; l++)
  		cons[l] = 0;
  if (!(*u))
  {
  	double* pr= new double[Size2*Size2]; 
			for (int l=0; l<(Size2*Size2); l++)
				pr[l] = 0;
		double* col = new double[Size2*K];
		double* li = new double[Size2*K];
			for (int l=0; l<(Size2*K); l++)
			{
				col[l] = 0;
				li[l] = 0;
			}
		CallPriorSize(Siz, Kmax, col, li, pr,u);
		for (int l=1; l<K; l++)
  		cons[l] = col[(l-1)*Size];
  	delete[] pr;
		delete[] col;
		delete[] li;
	}
  for (int k=0; k<K; k++)
    if (PriorK[k] != 0)
    {
    	if (!(*u))
    		BIC[k] = -Col[k*Size] + cons[k] -log(PriorK[k]);
    	else
    		BIC[k] = -Col[k*Size] + Norma(0,Size-1,k+1)-log(PriorK[k]);    	
    }
    else    
      BIC[k] = PLUS_INFINITY;
  int kchoose = 0;
  double mi = PLUS_INFINITY;
  for (int k=0; k<K; k++)
    if (BIC[k]<mi)
    {
	  mi = BIC[k];
 	  kchoose = k+1;
    }
  kBIC[0] = kchoose;   
  delete[] cons; 
  return;
}


void PostMean(int *Siz, int *Kk, double* Data, double* Col, double* Li, double *P, double *Post)
{
  int Size = *Siz;
  double* PS= new double[Size*Size];
  for (int l=0; l<(Size*Size); l++)
    PS[l] = 0;
  double* Me= new double[Size*Size];
  for (int l=0; l<(Size*Size); l++)
    Me[l] = 0;
  ProbaSegment(Siz, Kk, Col, Li, P, PS);
  Moyenne(Siz, Data, Me);
  for (int l=0; l<(Size-1); l++)
     for (int i=0; i<(l+1); i++)
	for (int j=(l+1); j<Size; j++)
	  Post[l] += (Me[i*Size+j]*PS[i*Size+j]);
  delete[] PS;
  delete[] Me;
  return;
}


/*
void ProfilesPoisson(int *Siz, int *I, int *Kmax, double* hyper, int* Data, double* LiI, double *ColI, double *RI)
{
	int Ii = *I;
	int n = *Siz;
	int K = *Kmax;
	for (int i=0; i<Ii; i++)
	{
		int* Datai = new int[n];
		for (int l=0; l<n; l++)
			Datai[l] = Data[i*n+l];
		double* hyperi = new double[2];
		hyperi[0] = hyper[2*i];
		hyperi[1] = hyper[2*i+1];
		double* Col = new double[K*(n+1)];
		double* Li = new double[K*(n+1)];
		for (int l=0; l<(K*(n+1)); l++)
		{
			Li[l]=0;
			Col[l]=0;
		}
		double* P = new double[(n+1)*(n+1)];
		for (int l=0; l<(n+1)*(n+1); l++)
			P[l]=0;
		CallEBSPoisson(Siz, Kmax, hyperi, Datai, Col, Li, P);
		for(int j=0; j<(n+1); j++)
		{
			LiI[i*(n+1)+j]=Li[j];
			ColI[i*(n+1)+j]=Col[(K-2)*(n+1)+j];
		}
		RI[i] = Col[(K-1)*(n+1)];
		delete[] Col;
		delete[] Li;
		delete[] Datai;
		delete[] P;
		delete[] hyperi;
	}
}

void ProfilesBinNeg(int *Siz, int *I, int *Kmax, double* hyper, double* theta, int* Data, double* LiI, double *ColI, double *RI)
{
	int Ii = *I;
	int n = *Siz;
	int K = *Kmax;
	for (int i=0; i<Ii; i++)
	{
	  double* thetai = new double[1];
	  thetai[0] = theta[i];
	  double* hyperi = new double[2];
		hyperi[0] = hyper[2*i];
		hyperi[1] = hyper[2*i+1];
		int* Datai = new int[n];
		for (int l=0; l<n; l++)
			Datai[l] = Data[i*n+l];
		double* Col = new double[K*(n+1)];
		double* Li = new double[K*(n+1)];
		for (int l=0; l<(K*(n+1)); l++)
		{
			Li[l]=0;
			Col[l]=0;
		}
		double* P = new double[(n+1)*(n+1)];
		for (int l=0; l<(n+1)*(n+1); l++)
			P[l]=0;
		CallEBSBinNeg(Siz, Kmax, hyperi, thetai, Datai, Col, Li, P);
		for(int j=0; j<(n+1); j++)
		{
			LiI[i*(n+1)+j]=Li[j];
			ColI[i*(n+1)+j]=Col[(K-2)*(n+1)+j];
		}
		RI[i] = Col[(K-1)*(n+1)];
		delete[] Col;
		delete[] Li;
		delete[] Datai;
		delete[] P;
		delete[] thetai;
		delete[] hyperi;
	}
}

void ProfilesGaussienne(int *Siz, int *I, int *Kmax, double* hyper, double* Data, double* LiI, double *ColI, double *RI)
{
	int Ii = *I;
	int n = *Siz;
	int K = *Kmax;
	for (int i=0; i<Ii; i++)
	{
		double* hyperi = new double[4];
		hyperi[0] = hyper[4*i];
		hyperi[1] = hyper[4*i+1];
		hyperi[2] = hyper[4*i+2];
		hyperi[3] = hyper[4*i+3];
		double* Datai = new double[n];
		for (int l=0; l<n; l++)
			Datai[l] = Data[i*n+l];
		double* Col = new double[K*(n+1)];
		double* Li = new double[K*(n+1)];
		for (int l=0; l<(K*(n+1)); l++)
		{
			Li[l]=0;
			Col[l]=0;
		}
		double* P = new double[(n+1)*(n+1)];
		for (int l=0; l<(n+1)*(n+1); l++)
			P[l]=0;
		CallEBSGaussienne(Siz, Kmax, hyperi, Datai, Col, Li, P);
		for(int j=0; j<(n+1); j++)
		{
			LiI[i*(n+1)+j]=Li[j];
			ColI[i*(n+1)+j]=Col[(K-2)*(n+1)+j];
		}
		RI[i] = Col[(K-1)*(n+1)];
		delete[] Col;
		delete[] Li;
		delete[] Datai;
		delete[] P;
		delete[] hyperi;
	}
}

void ProfilesGaussienneHomo(int *Siz, int *I, int *Kmax, double* hyper, double* Var, double* Data, double* LiI, double *ColI, double *RI)
{
	int Ii = *I;
	int n = *Siz;
	int K = *Kmax;
	for (int i=0; i<Ii; i++)
	{
		double* Vari = new double[1];
		Vari[0] = Var[i];
		double* hyperi = new double[2];
		hyperi[0] = hyper[2*i];
		hyperi[1] = hyper[2*i+1];
		double* Datai = new double[n];
		for (int l=0; l<n; l++)
			Datai[l] = Data[i*n+l];
		double* Col = new double[K*(n+1)];
		double* Li = new double[K*(n+1)];
		for (int l=0; l<(K*(n+1)); l++)
		{
			Li[l]=0;
			Col[l]=0;
		}
		double* P = new double[(n+1)*(n+1)];
		for (int l=0; l<(n+1)*(n+1); l++)
			P[l]=0;
		CallEBSGaussienneHomo(Siz, Kmax, hyperi, Vari, Datai, Col, Li, P);
		for(int j=0; j<(n+1); j++)
		{
			LiI[i*(n+1)+j]=Li[j];
			ColI[i*(n+1)+j]=Col[(K-2)*(n+1)+j];
		}
		RI[i] = Col[(K-1)*(n+1)];
		delete[] Col;
		delete[] Li;
		delete[] Datai;
		delete[] P;
		delete[] Vari;
		delete[] hyperi;
	}
}

*/
void LR(int *Siz, int *I, int *ci, int *ni, int* Kmax, double* LiI, double *ColI, double *RI, double* ratio)
{
	int Ii = *I;
	int n = *Siz;
	int K = *Kmax;
	double PH1 = 0;
	for (int i=0; i<(*ni); i++)
		PH1 += RI[ci[i]-1];
	PH1 = exp(PH1);
	double PH0 = 0;
	for (int t=(n+1-K); t>0; t--)
	{
		double aux = 0;
		for (int i=0; i<(*ni); i++)
			aux +=LiI[(ci[i]-1)*(n+1)+t]+ColI[(ci[i]-1)*(n+1)+t];
		PH0 += exp(aux);
	}
	ratio[0] = PH1/PH0;
}


