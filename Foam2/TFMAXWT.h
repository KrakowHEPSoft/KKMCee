#ifndef TFMAXWT_H
#define TFMAXWT_H
#include "ROOT_DEF.h"

#ifdef ROOT_DEF
#include "TROOT.h"
#endif

#include <iostream>
using namespace std;

class TH1D;
class TFHST;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      Small auxiliary class for controling MC weight.                      //
//      It provides certain measure of the maximum weight                    //
//      depending on small parameter epsilon.                                //
//      It uses 1-dimensional histograms of the TH1D class                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
class TFMAXWT : public TObject {
#else
class TFMAXWT{
#endif
 private:
  double  m_Nent;
  int     m_nBin;
  double  m_wmax;
 public:
#ifdef ROOT_DEF
  TH1D   *m_WtHst1;      // histogram of wt
  TH1D   *m_WtHst2;      // histogram of wt with wt
#else
  TFHST  *m_WtHst1;      // histogram of wt
  TFHST  *m_WtHst2;      // histogram of wt with wt
#endif
 public:
  TFMAXWT();                          // NOT IMPLEMENTED (NEVER USED)
  TFMAXWT(const double, const int);   // Principal Constructor
  TFMAXWT(TFMAXWT &From);             // Copy constructor
  ~TFMAXWT();                         // Destructor
  void Reset();                       // Reset
  TFMAXWT& operator =(TFMAXWT &);     // operator =
  void Fill(double);
  void Make(const double, double&);
  void GetMCeff(const double, double&, double&);  // get MC efficiency= <w>/wmax
////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TFMAXWT,1) //Controling of the MC weight (maximum weight)
#endif
};
#endif
