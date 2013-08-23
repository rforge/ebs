#include "CallEBS.h"
#include "R.h"
#include "Rmath.h"

extern "C"
{
  void SegmentPoisson(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P, bool* u)
  {
    CallEBSPoisson(Size, KMax, hyper, Data, Col, Li, P, u);
    return;
  }

  void SegmentBinNeg(int *Size, int *KMax, double* hyper, double *theta, int* Data, double* Col, double* Li, double* P, bool* u)
  {
    CallEBSBinNeg(Size, KMax, hyper, theta, Data, Col, Li, P, u);
    return;
  }

  void SegmentGaussienne(int *Size, int *KMax, double* hyper, double* Data, double* Col, double* Li, double* P, bool* u)
  {
    CallEBSGaussienne(Size, KMax, hyper, Data, Col, Li, P, u);
    return;
  }

  void SegmentGaussienneHomo(int *Size, int *KMax, double* hyper, double *Var, double* Data, double* Col, double* Li, double* P, bool* u)
  {
    CallEBSGaussienneHomo(Size, KMax, hyper, Var, Data, Col, Li, P, u);
    return;
  }
  
  void SetPriorUnif(int *Size, int *KMax, double* Col, double* Li, double* P, bool* u)
  {
    CallPriorUnif(Size, KMax, Col, Li, P, u);
    return;
  }
  
  void SetPriorSize(int *Size, int *KMax, double* Col, double* Li, double* P, bool* u)
  {
    CallPriorSize(Size, KMax, Col, Li, P, u);
    return;
  }

  void Distribution(int *Siz, int *kk, int*KK, double* Col, double* Li, double* Dist)
  {
    BreakDistrib(Siz, kk, KK, Col, Li, Dist);
    return;
  }

  void ChooseICL(int *Siz, int *Kmax, double *PriorK, double* Col, double* Li, double *P, double *ICL, int* kICL, bool* u)
  {
    GetICL(Siz, Kmax, PriorK, Col, Li, P, ICL, kICL, u);
    return;
  }
  
  void ChooseBIC(int *Siz, int *Kmax, double *PriorK, double* Col, double* BIC, int* kBIC, bool* u)
  {
    GetBIC(Siz, Kmax, PriorK, Col, BIC, kBIC, u);
    return;
  }
  
  void PostK(int *Siz, int *Kmax, double *PriorK, double* Col, double* Post, bool* u)
  {
    PosteriorK(Siz, Kmax, PriorK, Col, Post, u);
    return;
  }

  void GetPostMean(int *Siz, int *Kk, double* Data, double* Col, double* Li, double *P, double *Post)
  {
    PostMean(Siz, Kk, Data, Col, Li, P, Post);
    return;
  }
/*
  void GetProfilesPoisson(int *Siz, int *I, int *Kmax, double* hyper, int* Data, double* LiI, double *ColI, double *RI)
  {
    ProfilesPoisson(Siz, I, Kmax, hyper, Data, LiI, ColI, RI);
    return;
  }
  
  void GetProfilesBinNeg(int *Siz, int *I, int *Kmax, double* hyper, double *theta, int* Data, double* LiI, double *ColI, double *RI)
  {
    ProfilesBinNeg(Siz, I, Kmax, hyper, theta, Data, LiI, ColI, RI);
    return;
  }
  
  void GetProfilesGaussienne(int *Siz, int *I, int *Kmax, double* hyper, double* Data, double* LiI, double *ColI, double *RI)
  {
    ProfilesGaussienne(Siz, I, Kmax, hyper, Data, LiI, ColI, RI);
    return;
  }
  
  void GetProfilesGaussienneHomo(int *Siz,int *I, int *Kmax, double* hyper, double *Var, double* Data, double* LiI, double *ColI, double *RI)
  {
    ProfilesGaussienneHomo(Siz, I, Kmax, hyper, Var, Data, LiI, ColI, RI);
    return;
  }
*/  
  void GetLR(int *Siz, int *I, int *ci, int *ni, int *Kmax, double* LiI, double *ColI, double *RI, double* ratio)
  {
  	LR(Siz, I, ci, ni, Kmax, LiI, ColI, RI, ratio);
  	return;
  }
}
