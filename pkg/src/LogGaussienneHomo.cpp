
#include "LogGaussienneHomo.h"
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>

LogGaussienneHomo::LogGaussienneHomo()
{
  mu0 = 0;
  n0 = 0;
  var = 1;
}

LogGaussienneHomo::LogGaussienneHomo(Observations<double> &MyObs, double v, double mu, double nn)
{
  var = v;
  mu0 = mu;
  n0 = nn;
  LesObs = MyObs;
}


double LogGaussienneHomo::operator()(int a, int b)
{
  if (a==b)
    return 0;
  double S = LesObs.SumInSegment(a,b);
  double S2 = LesObs.SumSquareInSegment(a,b);
  int n = b-a;
  double constante = log(n0/(n0+n))/2-n/2*log(2*M_PI*var);
  double Res = constante -(S2+n0*mu0*mu0-(S+mu0*n0)*(S*mu0*n0)/(n+n0))/(2*var);
  return Res;
}


LogGaussienneHomo LogGaussienneHomo::operator=(const LogGaussienneHomo &Other)
{
  if (this != &Other)
    {
      var = Other.var;
      mu0 = Other.mu0;
      n0 = Other.n0;
      LesObs = Other.LesObs;
    }
  return *this;
}


