#ifndef TFOAM_INTEGRAND_H
#define TFOAM_INTEGRAND_H
#include "ROOT_DEF.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      Class TFOAM_INTEGRAND                                                //
//      Abstract class representing n-dimensional real integrand function    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
#include "TROOT.h"
class TFOAM_INTEGRAND : public TObject  {
#else
class TFOAM_INTEGRAND{
#endif
 public:
  TFOAM_INTEGRAND() { };
  virtual ~TFOAM_INTEGRAND() { };
  virtual double Density(int ndim, double*) = 0;
/////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TFOAM_INTEGRAND,1) //n-dimensional real positive integrand of FOAM
#endif
};
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      End of Class TFOAM_INTEGRAND                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#endif
