
#include "LogGaussienne.h"
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>

LogGaussienne::LogGaussienne()
{
  nu0 = 0;
  mu0 = 0;
  s0 = 0;
  n0 = 0;
}

LogGaussienne::LogGaussienne(Observations<double> &MyObs, double nu, double mu, double s, double nn)
{
  nu0 = nu;
  mu0 = mu;
  s0 = s;
  n0 = nn;
  LesObs = MyObs;
}


double LogGaussienne::operator()(int a, int b)
{
  if (a==b)
    return 0;
  double M = LesObs.MeanInSegment(a,b);
  double V = LesObs.VarInSegment(a,b);
  int n = b-a;
  double theta = 2/(n*V+s0+n*n0*(M-mu0)*(M-mu0)/(n+n0));
  double Res = lgamma((n+nu0)/2)+(log(n0)-log(n+n0))/2+(n+nu0)/2*log(theta)+nu0/2*log(s0/2)-lgamma(nu0/2)-n/2*log(2*M_PI);
  return Res;
}


LogGaussienne LogGaussienne::operator=(const LogGaussienne &Other)
{
  if (this != &Other)
    {
      nu0 = Other.nu0;
      mu0 = Other.mu0;
      s0 = Other.s0;
      n0 = Other.n0;
      LesObs = Other.LesObs;
    }
  return *this;
}



