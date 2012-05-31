
#include "LogPoisson.h"
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>

LogPoisson::LogPoisson()
{
  alpha = 0;
  beta = 0;
}

LogPoisson::LogPoisson(Observations<int> &MyObs, double a, double b)
{
  alpha = a;
  beta = b;
  LesObs = MyObs;
}


double LogPoisson::operator()(int a, int b)
{
  if (a==b)
    return 0;
  int S = LesObs.SumInSegment(a,b);
  double L = LesObs.LogFactorialInSegment(a,b);
  int n = b-a;
  double Res = lgamma(S+alpha)-(S+alpha)*log(n+beta)-L+alpha*log(beta)-lgamma(alpha);
  return Res;
}


LogPoisson LogPoisson::operator=(const LogPoisson &Other)
{
  if (this != &Other)
    {
      alpha = Other.alpha;
      beta = Other.beta;
      LesObs = Other.LesObs;
    }
  return *this;
}



