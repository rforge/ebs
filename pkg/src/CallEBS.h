#include "MyVector.h"
#include "Observations.h"
#include "Distributions.h"

void CallEBSPoisson(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P, bool* u);
void CallEBSBinNeg(int *Size, int *KMax, double* hyper, double* theta, int* Data, double* Col, double* Li, double* P, bool* u);
void CallEBSGaussienne(int *Size, int *KMax, double* hyper, double* Data, double* Col, double* Li, double* P, bool* u);
void CallEBSGaussienneHomo(int *Size, int *KMax, double* hyper, double *Var, double* Data, double* Col, double* Li, double* P, bool* u);
void CallPriorSize(int *Size, int *KMax, double* Col, double* Li, double* P, bool* u);
void CallPriorUnif(int *Size, int *KMax, double* Col, double* Li, double* P, bool* u);

void BreakDistrib(int *Siz, int *k, int*Kk, double* Col, double* Li, double* Dist);
void ProbaSegment(int *Siz, int*Kk, double* Col, double* Li, double *P, double* Pseg);
void GroupSegment(int *Siz, int*Kk, double* Col, double* Li, double *P, double* Pseg);
void Moyenne(int *Size, double* Data, double* P);
void GetICL(int *Siz, int *Kmax, double* PriorK, double* Col, double* Li, double *P, double* ICL, int* kICL, bool* u);
void PosteriorK(int *Siz, int *Kmax, double* PriorK, double* Col, double* BIC, bool* u);
void GetBIC(int *Siz, int *Kmax, double* PriorK, double* Col, double* BIC, int* kBIC, bool* u);
void PostMean(int *Siz, int *Kk, double* Data, double* Col, double* Li, double *P, double *Post);

/*
void ProfilesPoisson(int *Siz, int *I, int *Kmax, double* hyper, int* Data, double* LiI, double *ColI, double *RI);
void ProfilesBinNeg(int *Siz, int *I, int *Kmax, double* hyper, double* theta, int* Data, double* LiI, double *ColI, double *RI);
void ProfilesGaussienne(int *Siz, int *I, int *Kmax, double* hyper, double* Data, double* LiI, double *ColI, double *RI);
void ProfilesGaussienneHomo(int *Siz, int *I, int *Kmax, double* hyper, double* Var, double* Data, double* LiI, double *ColI, double *RI);
*/
void LR(int *Siz, int *I, int *ci, int *ni, int *Kmax, double* LiI, double *ColI, double *RI, double* ratio);
