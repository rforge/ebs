
#include "PriorSize.h"
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>



double PriorSize::operator()(int a, int b)
{
  if (a==b)
    return 0;
  int n = b-a;
  double Res = log((double) n);
  return Res;
}


PriorSize::PriorSize(double r)
{
	rien = r;
}

