

#ifndef _LogGaussienne_h_
#define _LogGaussienne_h_

#include "Function.h"
#include "Observations.h"

class LogGaussienne : public Function
{
public:

  double mu0;
  double n0;
  double nu0;
  double s0;
  Observations<double> LesObs;
  LogGaussienne();
  LogGaussienne(Observations<double> &MyObs,  double mu, double n, double nu, double s);
  double operator()(int, int); // rappel: segment [a,b[
  LogGaussienne operator=(const LogGaussienne &Other);
};



#endif
