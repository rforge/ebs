

#ifndef _LogBinNeg_h_
#define _LogBinNeg_h_

#include "Function.h"
#include "Observations.h"

class LogBinNeg : public Function
{
public:
  double alpha;
  double beta;
  double phi;
  Observations<int> LesObs;
  LogBinNeg();
  LogBinNeg(Observations<int> &MyObs, double p, double a, double b);
  double operator()(int, int); // rappel: segment [a,b[
  LogBinNeg operator=(const LogBinNeg &Other);
};



#endif
