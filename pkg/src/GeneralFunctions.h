/*
 *  GeneralFunctions.h
 *  Segments
 *
 *  Created by Michel Koskas on 22/08/11.
 *  Copyright 2011 INRA, INA. All rights reserved.
 *
 */
#ifndef _GeneralFunctions_h_
#define _GeneralFunctions_h_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Constants.h"
#include "MyVector.h"


bool IsDigit(char &x);
double lsum(double la, double lb);
double sumoflogs(double la, double lb);
double Norma(int a, int b, int k);


template <typename T>
void Swap(T &A, T &B)
{
	T C = A;
	A = B;
	B = C;
}


#endif
