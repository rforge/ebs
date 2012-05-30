/*
 *  GeneralMatrix.h
 *  DistributionRuptures
 *
 *  Created by Michel Koskas on 07/07/11.
 *  Copyright 2011 INRA, INA. All rights reserved.
 *
 */


#ifndef _LogMatrix_h_
#define _LogMatrix_h_

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <algorithm>
#include <assert.h>



template<typename FunctionTypeName>
class LogMatrix
{
public:
	int Size;	//Attention, Size est la taille du signal + 1
	double **Data;
	LogMatrix();
	LogMatrix(int N);
	LogMatrix(int N, FunctionTypeName *f);
	LogMatrix( LogMatrix<FunctionTypeName> &M);
	void Initialize(int N = 0, FunctionTypeName *f = NULL);
	void operator =(const LogMatrix<FunctionTypeName> &Original);


	~LogMatrix()
	{
		if (Data != NULL)
		{
			for (int i = 0; i < Size; i++)
				if (Data[i] != NULL)
					delete[] Data[i];
			delete[] Data;
		}
	}


};


template<typename FunctionTypeName>
void LogMatrix<FunctionTypeName>::operator =(const LogMatrix<FunctionTypeName> &Original)
{
	if (this == &Original)
		return;
	if (Data != NULL)
	{
		for (int i = 0; i < Size; i++)
			if (Data[i] != NULL)
				delete[] Data[i];
		delete[] Data;
	}
	Initialize(Original.Size);
	for (int i = 0; i < Size; i++)
		for (int j = i + 1; j < Size; j++)
			Data[i][j] = Original.Data[i][j];
	return;
}

template<typename FunctionTypeName>
LogMatrix<FunctionTypeName>::LogMatrix( LogMatrix<FunctionTypeName> &M)
{
	Initialize(M.Size);
	for (int i = 0; i < Size; i++)
		for (int j = i + 1; j < Size; j++)
			Data[i][j] = M.Data[i][j];
}

template<typename FunctionTypeName>
LogMatrix<FunctionTypeName>::LogMatrix()
{
	Initialize();
}


template<typename FunctionTypeName>
LogMatrix<FunctionTypeName>::LogMatrix(int N)
{
	Initialize(N);
}


template<typename FunctionTypeName>
LogMatrix<FunctionTypeName>::LogMatrix(int N, FunctionTypeName *f)
{
	Initialize(N, f);
}


template<typename FunctionTypeName>
void LogMatrix<FunctionTypeName>::Initialize(int N, FunctionTypeName *f)
{
	Size = N;
	Data = new double *[Size];
	for (int i = 0; i < Size; i++)
		Data[i] = new double[Size];
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size; j++)
			Data[i][j] = 0;
	if (f != NULL)
		for (int i = 0; i < Size; i++)
			for (int j = i + 1; j < Size; j++)
				Data[i][j] = (*f)(i, j);
}


#endif

