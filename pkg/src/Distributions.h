#ifndef _Distributions_h_
#define _Distributions_h_

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <algorithm>
#include <assert.h>

#include "LogMatrix.h"
#include "LogVector.h"
#include "Observations.h"
#include "Constants.h"
#include "GeneralFunctions.h"

template<typename FunctionTypeName>
class Distributions
{
private:
	int Size;	//Attention, Size est la taille du signal + 1
	int K;
	double **DataCol;
	double **DataLi;

public:
	LogMatrix<FunctionTypeName> P;
	Distributions();
	Distributions(int N, int K);
	Distributions(int N, int K, FunctionTypeName *f);
	void Initialize(int N=0, int K=0, FunctionTypeName *f=NULL);
	int GetSize();
	int GetK();
	double ** GetDataLi();
	double ** GetDataCol();
	double* GetDistBreak(int k, int Kk);	// Rappel: k° rupture=début du segment k+1, Kk=total number of segments
	double** GetProbSeg(int k);	// Proba du segment [i,j[ quand k segments
	//double GetEntropy(int k);	// Entropie de la segmentation en k segments
	int ChooseICL();	// Entropie de la segmentation en k segments

	void operator =(const Distributions<FunctionTypeName> &Original)
	{
		if (this == &Original)
			return;
		if (DataCol != NULL)
		{
			for (int i = 0; i < Size; i++)
				if (DataCol[i] != NULL)
					delete[] DataCol[i];
			delete[] DataCol;
		}
		if (DataLi != NULL)
		{
			for (int i = 0; i < K; i++)
				if (DataLi[i] != NULL)
					delete[] DataLi[i];
			delete[] DataLi;
		}
		P = Original.P;
		Initialize(Original.Size, Original.K);
		for (int i = 0; i < Size; i++)
			for (int j = 0; j < K; j++)
				DataCol[i][j] = Original.DataCol[i][j];
		for (int i = 0; i < K; i++)
			for (int j = 0; j < Size; j++)
				DataLi[i][j] = Original.DataLi[i][j];
		return;
	}

	~Distributions()
	{
		if (DataCol != NULL)
		{
			for (int i = 0; i < Size; i++)
				if (DataCol[i] != NULL)
					delete[] DataCol[i];
			delete[] DataCol;
		}
		if (DataLi != NULL)
		{
			for (int i = 0; i < K; i++)
				if (DataLi[i] != NULL)
					delete[] DataLi[i];
			delete[] DataLi;
		}
	}


};

template<typename FunctionTypeName>
int Distributions<FunctionTypeName>::GetSize()
{
 return Size;
}

template<typename FunctionTypeName>
int Distributions<FunctionTypeName>::GetK()
{
 return K;
}

template<typename FunctionTypeName>
double ** Distributions<FunctionTypeName>::GetDataLi()
{
 return DataLi;
}

template<typename FunctionTypeName>
double ** Distributions<FunctionTypeName>::GetDataCol()
{
 return DataCol;
}

template<typename FunctionTypeName>
Distributions<FunctionTypeName>::Distributions()
{
	Initialize();
}

template<typename FunctionTypeName>
Distributions<FunctionTypeName>::Distributions(int N, int K)
{
	Initialize(N, K);
}

template<typename FunctionTypeName>
Distributions<FunctionTypeName>::Distributions(int N, int K, FunctionTypeName *f)
{
	Initialize(N, K, f);

}

template<typename FunctionTypeName>
void Distributions<FunctionTypeName>::Initialize(int N, int MK, FunctionTypeName *f)
{
	Size = N;
	K = MK;
	DataCol = new double *[Size];
	DataLi = new double *[K];
	for (int i = 0; i < Size; i++)
		DataCol[i] = new double [K];
	for (int i = 0; i < K; i++)
		DataLi[i] = new double [Size];
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < K; j++)
			DataCol[i][j] = MINUS_INFINITY;
	for (int i = 0; i < K; i++)
		for (int j = 0; j < Size; j++)
			DataLi[i][j] = MINUS_INFINITY;
	if (f != NULL)
	{
	  LogMatrix<FunctionTypeName> Aux(N,f);
	  P = Aux;
	  LogVector<FunctionTypeName> Vcol(P, true);
	  LogVector<FunctionTypeName> Vli(P, false);
	  for (int i=0; i<Size; i++)
	  {
		DataCol[i][0] = Vcol.Data[i];
		DataLi[0][i] = Vli.Data[i];
	  }
	  for (int p=1; p<K; p++)
	  {
		Vcol*=P; //On calcule P*Vcol
		Vli *=P; // On calcule Vli*P
	        for (int l=0; l<Size; l++)
		{
			DataCol[l][p] = Vcol.Data[l];
			DataLi[p][l] = Vli.Data[l];
		}
	  }
	}
}


template<typename FunctionTypeName>
double* Distributions<FunctionTypeName>::GetDistBreak(int k,int Kk) //gets proba of breakpoint k for a segmentation in Kk segments
{
	if (k<1)
	{
		std::cerr << "come on, pick at least 2 segments, k has to be >1 " << std::endl;
		exit(222);
	}
	if (K<(k-1))
	{
		std::cerr << "I thought I only had to work up to K = " << K << std::endl;
		exit(223);
	}
	if (K<(Kk))
	{
		std::cerr << "I thought I only had to work up to K = " << K << std::endl;
		exit(223);
	}
	double* Res = new double[Size];
	double* Aux = new double[Size];
	double Norm = 0;
	for (int i=0; i<Size; i++)
		Aux[i] = DataLi[k-1][i] + DataCol[i][Kk-k-1];
	double Max = Aux[0];
  	for (int i = 0; i < Size; i++)
   	  if (Aux[i] > Max)
     	    Max = Aux[i];
	for (int i=0; i<Size; i++)
		Aux[i] -= Max;
	for (int i=0; i<Size; i++)
	{
		Res[i] = exp(Aux[i]);
		Norm += Res[i];
	}
	for (int i=0; i<Size; i++)
		Res[i] /= Norm;
	delete[] Aux;
	return Res;
}


// Cette fonction calcule, pour un segment [i,j[ la sortie S_k([i,j[)*log(f([i,j[) avec les notations de l'article
template<typename FunctionTypeName>
double** Distributions<FunctionTypeName>::GetProbSeg(int k)
{
	if (k<2)
	{
		std::cerr << "come on, pick at least 2 segments, k has to be >=2 " << std::endl;
		exit(224);
	}
	double **Res = new double*[Size];
	for (int i=0; i<Size; i++)
		Res[i] = new double[Size];
	for (int i=0; i<Size; i++)
	  for (int j=0; j<Size; j++)
		Res[i][j] = 0;
	for (int j=1; j<(Size-k); j++)
		Res[0][j] = P.Data[0][j]*exp(P.Data[0][j]-DataCol[0][k-1]+DataCol[j][k-2]);
		//Res[0][j] = P.Data[0][j]*exp(P.Data[0][j]-DataCol[0][k-1]+DataCol[j][k-2]+Norma(0,Size-1,k)-Norma(j,Size-1,k-1));
	for (int i=(k-1); i<(Size-1); i++)
		Res[i][Size-1] = P.Data[i][Size-1]*exp(P.Data[i][Size-1]-DataCol[0][k-1]+DataLi[k-2][i]);
		//Res[i][Size-1] = P.Data[i][Size-1]*exp(P.Data[i][Size-1]-DataCol[0][k-1]+DataLi[k-2][i]+Norma(0,Size-1,k)-Norma(0,i,k-1));
	if(k>2)
	{
	  for (int i=1; i<(Size-1); i++)
	    for (int j=i+1; j<(Size-1); j++)
	    {
	      double Aux = MINUS_INFINITY;
	      for (int p=2; p<k; p++)
		if ((i>p-2) & (i<Size-2-k+p) & (j<Size-k+p))
		  Aux = lsum(Aux,sumoflogs(DataLi[p-2][i],DataCol[j][k-p-1]));
		  //Aux = lsum(Aux,sumoflogs(DataLi[p-2][i]-Norma(0,i,p-1),DataCol[j][k-p-1]-Norma(j,Size-1,k-p)));
		if ((i>0) & (j<(Size-1)))
	      Res[i][j] = P.Data[i][j]*exp(P.Data[i][j]-DataCol[0][k-1] + Aux);
	      //Res[i][j] = P.Data[i][j]*exp(P.Data[i][j]-DataCol[0][k-1] + Aux+Norma(0,Size-1,k));
	    }
	}
	return Res;
}



template<typename FunctionTypeName>
int Distributions<FunctionTypeName>::ChooseICL()
{
	double *ICL=new double[K];
	for (int k=0; k<K; k++)
	  ICL[k]=0;
	for (int k=1; k<K; k++)
	{
	  double** PS=GetProbSeg(k+1);
	  if (k==1)
	    for (int i=1; i<(Size-1); i++)
		ICL[k] -= PS[0][i] + PS[i][Size-1];
	  else
	    for (int i=0; i<(Size-1); i++)
	      for (int j=i+1; j<std::min(Size-1,Size-k+i); j++)
		ICL[k] -= PS[i][j];
	  for (int i=0; i<Size; i++)
		delete[] PS[i];
	  delete[] PS;
	  ICL[k] += Norma(0,Size-1,k+1);
	  std::cout<< "k = " << k+1 << " et ICL(k) = " << ICL[k] << std::endl;
	}
	int kchoose = 0;
	double mi = PLUS_INFINITY;
	for (int k=1; k<K; k++)
	  if (ICL[k]<mi)
	  {
		mi = ICL[k];
		kchoose = k+1;
	  }
	delete[] ICL;
	return kchoose;
}

#endif

