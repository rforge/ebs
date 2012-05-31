#include "CallEBS.h"
#include "R.h"
#include "Rmath.h"

extern "C"
{
  void SegmentPoisson(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P)
  {
    CallEBSPoisson(Size, KMax, hyper, Data, Col, Li, P);
    return;
  }

  void SegmentBinNeg(int *Size, int *KMax, double* hyper, double *theta, int* Data, double* Col, double* Li, double* P)
  {
    CallEBSBinNeg(Size, KMax, hyper, theta, Data, Col, Li, P);
    return;
  }

  void SegmentGaussienne(int *Size, int *KMax, double* hyper, int* Data, double* Col, double* Li, double* P)
  {
    CallEBSGaussienne(Size, KMax, hyper, Data, Col, Li, P);
    return;
  }

  void Distribution(int *Siz, int *kk, int*KK, double* Col, double* Li, double* Dist)
  {
    BreakDistrib(Siz, kk, KK, Col, Li, Dist);
    return;
  }

  void ChooseICL(int *Siz, int *Kmax, double* Col, double* Li, double *P, double *ICL, int* kICL)
  {
    GetICL(Siz, Kmax, Col, Li, P, ICL, kICL);
    return;
  }
  
  void ChooseBIC(int *Siz, int *Kmax, double* Col, double* BIC, int* kBIC)
  {
    GetBIC(Siz, Kmax, Col, BIC, kBIC);
    return;
  }
  
  void PostK(int *Siz, int *Kmax, double* Col, double* Post)
  {
    PosteriorK(Siz, Kmax, Col, Post);
    return;
  }

  void GetPostMean(int *Siz, int *Kk, double* Data, double* Col, double* Li, double *P, double *Post)
  {
    PostMean(Siz, Kk, Data, Col, Li, P, Post);
    return;
  }

}
