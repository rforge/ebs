/*
 *  Observations.h
 *  DistributionRuptures
 *
 *  Created by Alice Cleynen on 24/03/11.
 *  Copyright 2011 INRA, INA. All rights reserved.
 *
 */



#ifndef _Observations_h_
#define _Observations_h_

#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "MyVector.h"
#include "GeneralFunctions.h"

template<typename T>
class Observations
{
public:
  MyVector<T> y;
  Observations();
  Observations(MyVector<T> &v);
  Observations(std::string FileName, bool Binaire = false);
  Observations(std::string FileName, int FirstIndex, int LastIndex, bool Binaire = false);
  T SumInSegment(int a, int b);
  double MeanInSegment(int a, int b);
  double VarInSegment(int a, int b);
  double LogFactorialInSegment(int a, int b);
  std::ostream &operator<<(std::ostream &s);
};

template<typename T>
std::ostream &operator<<(std::ostream &s, const Observations<T> &O)
{
  for (int i=0; i<O.y.size(); i++)
    s <<  O.y[i] << std::endl;
  return s;
}



template<typename T>
Observations<T>::Observations()
{
}

template<typename T>
Observations<T>::Observations(MyVector<T> &v)
{
  y = v;
}

template<typename DataTypeName>
Observations<DataTypeName>::Observations(std::string FileName, bool Binaire)
{
  std::ifstream TheFile(FileName.c_str());
  if (!TheFile.is_open())
  {
    std::cerr << "Can't open the file " << FileName << ". Getting out with errcode 200" << std::endl;
    exit(200);
  }
	if (Binaire)
	{
		TheFile.seekg(0, std::ios::end);
		int NbElements =TheFile.tellg() / sizeof(DataTypeName);
		TheFile.seekg(0, std::ios::beg);
		for (int i = 0; i < NbElements; i++)
		{
			DataTypeName CurElement;
			TheFile.read((char *) &CurElement, 1 * sizeof(DataTypeName));
			y.push_back(CurElement);
		}
	}
	else
	{
		char *Buffer;
		TheFile.seekg(0, std::ios::end);
		int FileSize = TheFile.tellg();
		Buffer = new char[FileSize];
		TheFile.seekg(0, std::ios::beg);
		TheFile.read(Buffer, FileSize * sizeof(char));
		int BuffIndex = 0;
		DataTypeName Res;
		bool SansInteret;
		while (NextNumber<DataTypeName>(Buffer, BuffIndex, FileSize, Res, SansInteret))
			y.push_back(Res);
		delete[] Buffer;
	}
	TheFile.close();
}

template<typename DataTypeName>
Observations<DataTypeName>::Observations(std::string FileName, int FirstIndex, int LastIndex, bool Binaire)
{
  std::ifstream TheFile(FileName.c_str());
  if (!TheFile.is_open())
  {
    std::cerr << "Can't open the file " << FileName << ". Getting out with errcode 200" << std::endl;
    exit(200);
  }
	if (Binaire)
	{
		int NbElements = LastIndex - FirstIndex + 1;
		int begin = (FirstIndex - 1) * sizeof(DataTypeName);
		TheFile.seekg(begin, std::ios::beg);
		for (int i = 0; i < NbElements; i++)
		{
			DataTypeName CurElement;
			TheFile.read((char *) &CurElement, 1 * sizeof(DataTypeName));
			y.push_back(CurElement);
		}
	}
	else
	{
		char *Buffer;
		TheFile.seekg(0, std::ios::end);
		int FileSize = TheFile.tellg();
		int NbElements = LastIndex - FirstIndex + 1;
		Buffer = new char[FileSize];
		TheFile.seekg(0, std::ios::beg);
		TheFile.read(Buffer, FileSize * sizeof(char));
		int BuffIndex = 0;
		DataTypeName Res;
		bool SansInteret;
		MyVector<DataTypeName> trash;
		while (trash.size() < FirstIndex - 1)
		{
			bool NumberFound = NextNumber<DataTypeName>(Buffer, BuffIndex, FileSize, Res, SansInteret);
			if (NumberFound)
				trash.push_back(Res);
		}
		while (y.size() < NbElements)
		{
			bool NumberFound = NextNumber<DataTypeName>(Buffer, BuffIndex, FileSize, Res, SansInteret);
			if (NumberFound)
				y.push_back(Res);
		}
		delete[] Buffer;
	}
	TheFile.close();
}


template<typename T>
T Observations<T>::SumInSegment(int a, int b)
{
  T S=0;
  for (int i = a; i < b; i++)
	{
    		S += y[i];
		if (y[i] < 0)
		{
			std::cerr << "Illicite data." << std::endl;
			exit(102);
		}
	}
  return S;
}

template<typename T>
double Observations<T>::MeanInSegment(int a, int b)
{
  double M=0;
  for (int i = a; i < b; i++)
	M += y[i];
  M/=(b-a);
  return M;
}

template<typename T>
double Observations<T>::VarInSegment(int a, int b)
{
  double V=0;
  double M=MeanInSegment(a,b);
  for (int i = a; i < b; i++)
	V += (y[i]-M)*(y[i]-M);
  V/=(b-a);
  return V;
}


template<typename T>
double Observations<T>::LogFactorialInSegment(int a, int b)
{
  if (a==b)
    return 0;
  double S = 0;
  for (int i = a; i < b; i++)
	for (int j=1; j < (y[i] + 1); j++)
	    S += log(j);
  return S;
}


// ##################################################

template<typename T>
class VectObservations
{
public:
  MyVector<Observations<T> > ReplicateObs;
  MyVector<T> SumInSegment(int a, int b);
  MyVector<double> LogFactorialInSegment(int a, int b);
  VectObservations();
  VectObservations(Observations<T> &);
  VectObservations(MyVector<Observations<T> > &);
  VectObservations<T> operator=(VectObservations<T> Other);
};


template<typename T>
std::ostream &operator<<(std::ostream &s, const VectObservations<T> &VO)
{
  for (int j = 0; j < VO.ReplicateObs[0].y.size(); j++)
    {
      for (int i = 0; i < VO.ReplicateObs.size(); i++)
        s <<  VO.ReplicateObs[i].y[j] << "  ";
      s << std::endl;
    }
  return s;
}

template<typename T>
VectObservations<T> VectObservations<T>::operator=(VectObservations<T> Other)
{

	if (this != &Other)
	{
		ReplicateObs.clear();
		for (int i = 0; i < Other.ReplicateObs.size(); i++)
			ReplicateObs.push_back(Other.ReplicateObs[i]);
	}
	return *this;
}


template<typename T>
VectObservations<T>::VectObservations()
{
}

template<typename T>
VectObservations<T>::VectObservations(MyVector<Observations<T> > &V)
{
	// TODO: Check this (Need a push_back for each element?)
  ReplicateObs = V;
}

template<typename T>
VectObservations<T>::VectObservations(Observations<T> &O)
{
  ReplicateObs.push_back(O);
}

template<typename T>
MyVector<T> VectObservations<T>::SumInSegment(int a, int b)
{
  MyVector<T> S;
  for (int j = 0; j < ReplicateObs.size(); j++)
    {
      T l= ReplicateObs[j].SumInSegment(a,b);
      if (l < 0)
	{
		std::cerr << "I unhappy. Getting out." << std::endl;
		exit(103);
	}
      S.push_back(l);
    }
  return S;
}

template<typename T>
MyVector<double> VectObservations<T>::LogFactorialInSegment(int a, int b)
{
  MyVector<double> S;
  for (int j = 0; j < ReplicateObs.size(); j++)
    {
      Observations<T> V=ReplicateObs[j];
      double Aux=V.LogFactorialInSegment(a,b);
      S.push_back(Aux);
    }
  return S;
}


#endif
