
#include "PriorUnif.h"
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>



double PriorUnif::operator()(int a, int b)
{
  if (a==b)
    return 0;
  double Res = 1;
  return Res;
}


PriorUnif::PriorUnif(double r)
{
	rien = r;
}

