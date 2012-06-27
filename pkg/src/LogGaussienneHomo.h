

#ifndef _LogGaussienneHomo_h_
#define _LogGaussienneHomo_h_

#include "Function.h"
#include "Observations.h"

class LogGaussienneHomo : public Function
{
public:
  double mu0;
  double n0;
  double var;
  Observations<double> LesObs;
  LogGaussienneHomo();
  LogGaussienneHomo(Observations<double> &MyObs, double v, double mu, double n);
  double operator()(int, int); // rappel: segment [a,b[
  LogGaussienneHomo operator=(const LogGaussienneHomo &Other);
};



#endif
