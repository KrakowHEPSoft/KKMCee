#ifndef KKsem_H
#define KKsem_H
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKsem                                              //
//                                                                          //
//   Interface (wrapper)  to KKsem package for semi-analytical programs     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include<stdlib.h>
#include<stdio.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

//------------------------------------------------------------
extern "C" void kk2f_fort_open_( const long&, char*, long);
extern "C" void kk2f_fort_close_(const long&);
//      SUBROUTINE KKsem_Initialize(xpar_input)
extern "C" void kksem_initialize_(double xpar[]);
//      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
//      SUBROUTINE KKsem_SetKFfin(KFfin)
extern "C" void kksem_setkffin_( const long& );
//      SUBROUTINE BornV_SetKF(KFferm)
extern "C" void bornv_setkf_( const long& ); // set SINGLE Final State
//      SUBROUTINE KKsem_SetKeyFoB(KeyFoB)
extern "C" void kksem_setkeyfob_( const long& );
//      SUBROUTINE KKsem_VVplot_vec(key,chak,nbin,xmin,xmax,yy)
extern "C" void kksem_vvplot_vec_(const long&, char[5], const long&, const double&, const double&, double[]);
//      SUBROUTINE KKsem_SetCrange(Cmin,Cmax)
extern "C" void kksem_setcrange_(const double&, const double&);
//      SUBROUTINE KKsem_SetKeyZet(KeyZet)
extern "C" void kksem_setkeyzet_( const long& );
//      SUBROUTINE KKsem_MakeBorn(svar,Born)
extern "C" void kksem_makeborn_(const double&, double&);

//      SUBROUTINE BornV_MakeGami(CMSene,gamiCR,gami,alfi)
extern "C" void bornv_makegami_(const double&, double&, double&, double&);

//------------------------------------------------------------


//class KKsem : public TNamed {
class KKsem{
  //double m_ypar ... to be implemented?
 public:
//------ constructors destructors -------
  KKsem(){;};
  ~KKsem(){;};
 public:
  void Initialize(TFile&);
  void VVplot( TH1 *, long , char [], long, long );
  void Cplot(  TH1 *, long , char [], long, long, double, double);
////////////////////////////////////////////////////////////////////////////
//                      ClassDef(KKsem,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
