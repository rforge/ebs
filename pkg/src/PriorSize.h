

#ifndef _PriorSize_h_
#define _PriorSize_h_

#include "Function.h"
#include "Observations.h"

class PriorSize : public Function
{
public:
	double rien;
	PriorSize(double);
  double operator()(int, int); // rappel: segment [a,b[
};



#endif
