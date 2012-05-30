

#ifndef _LogPoisson_h_
#define _LogPoisson_h_

#include "Function.h"
#include "Observations.h"

class LogPoisson : public Function
{
public:
  double alpha;
  double beta;
  Observations<int> LesObs;
  LogPoisson();
  LogPoisson(Observations<int> &MyObs, double a, double b);
  double operator()(int, int); // rappel: segment [a,b[
  LogPoisson operator=(const LogPoisson &Other);
};



#endif
