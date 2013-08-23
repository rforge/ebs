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
	Distributions(int N, int K, FunctionTypeName *f, bool *u);
	void Initialize(int N=0, int K=0, FunctionTypeName *f=NULL, bool *u=NULL);
	double ** GetDataLi();
	double ** GetDataCol();
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
Distributions<FunctionTypeName>::Distributions(int N, int K, FunctionTypeName *f, bool *u)
{
	Initialize(N, K, f, u);

}

template<typename FunctionTypeName>
void Distributions<FunctionTypeName>::Initialize(int N, int MK, FunctionTypeName *f, bool *u)
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
	  LogMatrix<FunctionTypeName> Aux(N,f,u);
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


#endif

