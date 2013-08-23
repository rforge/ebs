

#ifndef _PriorUnif_h_
#define _PriorUnif_h_

#include "Function.h"
#include "Observations.h"

class PriorUnif : public Function
{
public:
	double rien;
	PriorUnif(double);
  double operator()(int, int); // rappel: segment [a,b[
};



#endif
