#ifndef _LogVector_h_
#define _LogVector_h_

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <algorithm>

#include "LogMatrix.h"
#include "Constants.h"

template<typename FunctionTypeName>
class LogVector
{
public:
	int Size;	//Attention, Size est la taille du signal + 1
	double *Data;
	bool Column;
	LogVector();
	LogVector(int N, bool Col);
	LogVector(int N, FunctionTypeName *f, bool Col);
	LogVector(LogVector &V);
	LogVector(LogMatrix<FunctionTypeName> &M, bool Col);
	void Initialize(bool Col, int N=0 ,FunctionTypeName *f=NULL);
	void operator *=(const LogMatrix<FunctionTypeName> &Original);
	LogVector operator *(const LogMatrix<FunctionTypeName> &Original);

	LogVector operator =(const LogVector<FunctionTypeName> &Original)
	{
		if (this == &Original)
			return *this;
		if (Data != NULL)
			delete[] Data;
		Initialize(Original.Column, Original.Size);
		for (int i = 0; i < Size; i++)
				Data[i] = Original.Data[i];
		return *this;
	}

	~LogVector()
	{
		if (Data != NULL)
			delete[] Data;
	}


};


template<typename FunctionTypeName>
LogVector<FunctionTypeName>::LogVector()
{
	Initialize(false);
}

template<typename FunctionTypeName>
LogVector<FunctionTypeName>::LogVector(int N, bool Col)
{
	Initialize(Col, N);
}

template<typename FunctionTypeName>
LogVector<FunctionTypeName>::LogVector(LogVector &V)
{
	Initialize(V.Column, V.Size);
	for (int i = 0; i < Size; i++)
		Data[i] = V.Data[i];
	Column = V.Column;
}

template<typename FunctionTypeName>
LogVector<FunctionTypeName>::LogVector(LogMatrix<FunctionTypeName> &M, bool Col)
{
	Initialize(Col, M.Size);
	if (Col)
		for (int i = 0; i < Size; i++)
			Data[i] = M.Data[i][Size-1];
	else
		for (int i = 0; i < Size; i++)
			Data[i] = M.Data[0][i];
}

template<typename FunctionTypeName>
void LogVector<FunctionTypeName>::Initialize(bool Col, int N, FunctionTypeName *f)
{
	Size = N;
	Column = Col;
	Data = new double [Size];
	for (int i = 0; i < Size; i++)
		Data[i] = MINUS_INFINITY;
	if (f != NULL)
	{
	  if (Col)
		for (int i = 0; i < Size; i++)
			Data[i] = (*f)(i, Size);
	  else
		for (int i = 0; i < Size; i++)
			Data[i] = (*f)(0,i);
	}
}

template<typename FunctionTypeName>
LogVector<FunctionTypeName> LogVector<FunctionTypeName>::operator *(const LogMatrix<FunctionTypeName> &Original)
{
	LogVector<FunctionTypeName> Res(Size, Column);
	if (Column)
		for (int i=0; i<Size; i++)
			for (int p=0; p<Size; p++)
				Res.Data[i]= lsum(Res.Data[i],sumoflogs(Original.Data[i][p],Data[p]));
	else
		for (int i=0; i<Size; i++)
			for (int p=0; p<Size; p++)
				Res.Data[i]= lsum(Res.Data[i],sumoflogs(Original.Data[p][i],Data[p]));
	return Res;
}

template<typename FunctionTypeName>
void LogVector<FunctionTypeName>::operator *=(const LogMatrix<FunctionTypeName> &Original)
{
	LogVector<FunctionTypeName> Aux(Original.Size, Column);
	if (Column)
		for (int i=0; i<Size; i++)
			for (int p=0; p<Size; p++)
				Aux.Data[i]= lsum(Aux.Data[i],sumoflogs(Original.Data[i][p],Data[p]));
	else
		for (int i=0; i<Size; i++)
			for (int p=0; p<Size; p++)
				Aux.Data[i]= lsum(Aux.Data[i],sumoflogs(Original.Data[p][i],Data[p]));
	Aux.Column = Column;
	(*this) = Aux;
}
#endif

