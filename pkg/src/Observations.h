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
  T SumInSegment(int a, int b);
  T SumSquareInSegment(int a, int b);
  double MeanInSegment(int a, int b);
  double VarInSegment(int a, int b);
  double LogFactorialInSegment(int a, int b);
  double LogGammaPhiInSegment(int a, int b, double phi);
};



template<typename T>
Observations<T>::Observations()
{
}

template<typename T>
Observations<T>::Observations(MyVector<T> &v)
{
  y = v;
}



template<typename T>
T Observations<T>::SumInSegment(int a, int b)
{
  T S=0;
  for (int i = a; i < b; i++)
    		S += y[i];
  return S;
}

template<typename T>
T Observations<T>::SumSquareInSegment(int a, int b)
{
  T S=0;
  for (int i = a; i < b; i++)
    		S += y[i]*y[i];
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
  /*
	for (int j=1; j < (y[i] + 1); j++)
	    S += log(j);
  */
    S += lgamma(y[i]+1);
  return S;
}

template<typename T>
double Observations<T>::LogGammaPhiInSegment(int a, int b, double phi)
{
  if (a==b)
    return 0;
  double S = 0;
  for (int i = a; i < b; i++)
    S += lgamma(y[i]+phi);
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
