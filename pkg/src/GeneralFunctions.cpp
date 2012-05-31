/*
 *  GeneralFunctions.cpp
 *  Segments
 *
 *  Created by Michel Koskas on 22/08/11.
 *  Copyright 2011 INRA, INA. All rights reserved.
 *
 */

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include "GeneralFunctions.h"
#include "MyVector.h"


// This function says whether the character x is a number or not
bool IsDigit(char &x)
{
  if (x < '0')
    return false;
  if (x > '9')
    return false;
  return true;
}

// returns log(exp(la)+exp(lb))
double lsum(double la, double lb)
{
  if (la>lb)
    if (lb>MINUS_INFINITY)
	return la+log1p(exp(lb-la));
    else return la;
  else
    if (la>MINUS_INFINITY)
	return lb+log1p(exp(la-lb));
    else
	return lb;
}


double sumoflogs(double la, double lb)
{
    if (lb<=MINUS_INFINITY)
	return MINUS_INFINITY;
    if (la<=MINUS_INFINITY)
	return MINUS_INFINITY;
    return (la +lb );
}

double Norma(int a, int b, int k) //log((b-a-1)choose(k-1))
{
    if (a==b)
	return 0;
    double Res = (lgamma(b-a)-lgamma(k)-lgamma(b-a-k+1));
    return Res;
}



int GetRandomNumber(int MinValue, int MaxValue)
{
  long int x = rand();
  double y = (MaxValue - MinValue + 1) * ((double) x) / RAND_MAX;
  int Res = MinValue + y;
  return Res;
}

MyVector<int> IntersectLists(const MyVector<int> &A, const MyVector<int> &B)
{
	MyVector<int> Res;
	MyVector<int>::const_iterator IA = A.begin(), IB = B.begin();
	while ((IA != A.end()) && (IB != B.end()))
		if (*IA < *IB)
			IA++;
		else if (*IB < *IA)
			IB++;
		else
		{
			Res.push_back(*IA);
			IA++;
			IB++;
		}
	return Res;
}


