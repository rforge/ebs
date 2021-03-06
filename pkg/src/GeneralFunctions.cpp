
#include <iostream>
#include <cstdlib>
#include <string>
#include <Rmath.h>
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
    double Res = (lgammafn(b-a)-lgammafn(k)-lgammafn(b-a-k+1));
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


