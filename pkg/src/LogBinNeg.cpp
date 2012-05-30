
#include "LogBinNeg.h"
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>

LogBinNeg::LogBinNeg()
{
  alpha = 0;
  beta = 0;
  phi = 1;
}

LogBinNeg::LogBinNeg(Observations<int> &MyObs, double p, double a, double b)
{
  alpha = a;
  beta = b;
  phi = p;
  LesObs = MyObs;
}


double LogBinNeg::operator()(int a, int b)
{
  if (a==b)
    return 0;
  int S = LesObs.SumInSegment(a,b);
  double L = LesObs.LogFactorialInSegment(a,b);
  int n = b-a;
  double Res = lgamma(beta+n*phi)+lgamma(S+alpha)-lgamma(beta+alpha+n*phi+S)-L;
  return Res;
}


LogBinNeg LogBinNeg::operator=(const LogBinNeg &Other)
{
  if (this != &Other)
    {
      alpha = Other.alpha;
      beta = Other.beta;
      phi = Other.phi;
      LesObs = Other.LesObs;
    }
  return *this;
}



