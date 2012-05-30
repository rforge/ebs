/*
 *  Function.h
 *  DistributionRuptures
 *
 *  Created by Alice Cleynen on 24/03/11.
 *  Copyright 2011 INRA, INA. All rights reserved.
 *
 */

#ifndef _Function_h_
#define _Function_h_


#include <vector>


class Function
{
public:
  // false by default in ALL constructors
  Function();
  ~Function();
  double operator()();
  double operator()( int, int);
};






#endif
