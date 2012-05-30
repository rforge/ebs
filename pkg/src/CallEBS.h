#include "MyVector.h"
#include "Observations.h"
#include "Distributions.h"

void CallEBSPoisson(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P);
void CallEBSBinNeg(int *Size, int *KMax, double* hyper, double* theta, int* Data, double* Col, double* Li, double* P);
void CallEBSGaussienne(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P);

void BreakDistrib(int *Siz, int *k, int*Kk, double* Col, double* Li, double* Dist);
void ProbaSegment(int *Siz, int*Kk, double* Col, double* Li, double *P, double* Pseg);
void GroupSegment(int *Siz, int*Kk, double* Col, double* Li, double *P, double* Pseg);
void Moyenne(int *Size, double* Data, double* P);
void ICL(int *Siz, int *Kmax, double* Col, double* Li, double *P, int* kICL);
void ComputeBIC(int *Siz, int *Kmax, double* Col, double* BIC);
void PosteriorK(int *Siz, int *Kmax, double* Col, double* BIC);
void KBIC(int *Siz, int *Kmax, double* Col, int* kBIC);
void PostMean(int *Siz, int *Kk, double * Data, double* Col, double* Li, double *P, double *Post);
