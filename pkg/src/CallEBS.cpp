#include "CallEBS.h"
#include "LogPoisson.h"
#include "LogBinNeg.h"
#include "LogGaussienne.h"

void CallEBSPoisson(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P)
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
    Distributions<LogPoisson> Dis(n+1, K, &BN);
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

void CallEBSBinNeg(int *Size, int *KMax, double* hyper, double* theta, int* Data, double* Col, double* Li, double* P)
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
    Distributions<LogBinNeg> Dis(n+1, K, &BN);
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

void CallEBSGaussienne(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P)
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
    Distributions<LogGaussienne> Dis(n+1, K, &G);
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

void GroupSegment(int *Siz, int*Kk, double* Col, double* Li, double *P, double* Pseg)
{

  int Size = *Siz;
  int k = *Kk;
  for (int j=1; j<(Size-k); j++)
	Pseg[j] = P[j] * exp(P[j]-Col[(k-1)*Size]+Col[j+(k-2)*Size]);
  for (int i=(k-1); i<(Size-1); i++)
	Pseg[i*Size+Size-1] = P[i*Size+Size-1] * exp(P[i*Size+Size-1]-Col[(k-1)*Size]+Li[(k-2)*Size+i]);
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
	    Pseg[i*Size+j] = P[i*Size+j] * exp(P[i*Size+j]-Col[(k-1)*Size] + Aux);
      }
  }
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

void GetICL(int *Siz, int *Kmax, double* Col, double* Li, double *P, double* ICL, int* kICL)
{
  int K = *Kmax;
  int Size = *Siz;
  ICL[0]=-Col[0];
  for (int k=1; k<K; k++)
  {
    int* kneeded= new int[1];
    kneeded[0]=k+1;
    double* PS= new double[Size*Size];
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
    ICL[k] += Norma(0,Size-1,k+1);
  }
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


void PosteriorK(int *Siz, int *Kmax, double* Col, double* PostK)
{
  int K = *Kmax;
  int Size = *Siz;
  for (int k=0; k<K; k++)
  {
    PostK[k] = -Col[k*Size];
    PostK[k] += Norma(0,Size-1,k+1);
  }
  double Total = 0;
  double Max = PostK[0];
  for (int i = 0; i < K; i++)
    if (PostK[i] > Max)
     	Max = PostK[i];
  for (int i=0; i<K; i++)
	PostK[i] -= Max;
  for (int k=0; k<K; k++)
    Total += exp(-PostK[k]);
  for (int k=0; k<K; k++)
    PostK[k] = exp(-PostK[k])/Total;
  return;
}

void GetBIC(int *Siz, int *Kmax, double* Col, double* BIC, int* kBIC)
{
  int K = *Kmax;
  int Size = *Siz;
  for (int k=0; k<K; k++)
  {
    BIC[k] = -Col[k*Size];
    BIC[k] += Norma(0,Size-1,k+1);
  }
  int kchoose = 0;
  double mi = PLUS_INFINITY;
  for (int k=0; k<K; k++)
    if (BIC[k]<mi)
    {
	  mi = BIC[k];
 	  kchoose = k+1;
    }
  kBIC[0] = kchoose;    
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
