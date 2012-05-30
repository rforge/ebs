

#ifndef _LogGaussienne_h_
#define _LogGaussienne_h_

#include "Function.h"
#include "Observations.h"

class LogGaussienne : public Function
{
public:
  double nu0;
  double mu0;
  double s0;
  double n0;
  Observations<double> LesObs;
  LogGaussienne();
  LogGaussienne(Observations<double> &MyObs, double nu, double mu, double s, double n);
  double operator()(int, int); // rappel: segment [a,b[
  LogGaussienne operator=(const LogGaussienne &Other);
};



#endif
