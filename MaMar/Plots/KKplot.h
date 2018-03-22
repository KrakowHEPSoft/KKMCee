//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKplot                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// KKplot is multipurpose toolbox for KKMC testing.
//  1. Interfaces (wrappers) to KKMC and KKsem F77 subrograms
//  2. Integrand for Foam in semianalytical xcheck
//  3. A few routines for producing latex table out of histograms
//////////////////////////////////////////////////////////////////////////////

#ifndef KKplot_H
#define KKplot_H
#include<stdlib.h>
#include<stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


#include "HisNorm.h"

//------------------------------------------------------------
//  wrappers to f77 routines in KKMC and KKsem
extern "C" void kk2f_fort_open_( const int&, const char*, int);
extern "C" void kk2f_fort_close_(const int&);
//      SUBROUTINE KKsem_Initialize(xpar_input)
extern "C" void kksem_initialize_(double xpar[]);
//      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
//-----------------------
//      SUBROUTINE KKsem_SetKFfin(KFfin)
extern "C" void kksem_setkffin_( const int& );
//      SUBROUTINE BornV_SetKF(KFferm)
extern "C" void bornv_setkf_( const int& ); // set SINGLE Final State
//      SUBROUTINE KKsem_SetKeyFoB(KeyFoB)
extern "C" void kksem_setkeyfob_( const int& );
//      SUBROUTINE KKsem_VVplot_vec(key,chak,nbin,xmin,xmax,yy)
extern "C" void kksem_vvplot_vec_(const int&, char[5], const int&, const double&, const double&, double[]);
//      SUBROUTINE KKsem_SetCrange(Cmin,Cmax)
extern "C" void kksem_setcrange_(const double&, const double&);
//      SUBROUTINE KKsem_SetKeyZet(KeyZet)
extern "C" void kksem_setkeyzet_( const int& );
//      SUBROUTINE KKsem_MakeBorn(svar,Born)
extern "C" void kksem_makeborn_(const double&, double&);
//      SUBROUTINE KKsem_Afb_Calc(KeyDist,KFi,KFf,CMSene,vv,Result)
extern "C" void kksem_afb_calc_(const int&, const int&, const int&, const double& , const double& , const double&);
//      KKsem_Born_Calc(KFi,KFf, AlfRun, CMSene, xres)
extern "C" void kksem_born_calc_(const int&, const int&, const double&, const double& , double[]);
//---------------------------------
//      DOUBLE PRECISION  FUNCTION BornV_Sig0nb(CMSene)
extern "C" double bornv_sig0nb_(const double&);
//      SUBROUTINE BornV_MakeGami(CMSene,gamiCR,gami,alfi)
extern "C" void bornv_makegami_(const double&, double&, double&, double&);
//      DOUBLE PRECISION  FUNCTION BornV_Simple(KFi,KFf,svar,costhe)
extern "C" double bornv_simple_( const int&,  const int&, const double&, const double&);
//------------------------------------------------------------
//      SUBROUTINE BornV_InterpoGSW(KFf,svar,CosThe)
extern "C" double bornv_interpogsw_( const int&,  const double&, const double&);
//      DOUBLE PRECISION FUNCTION BornV_Dizet(Mode,KFi,KFf,svar,CosThe,eps1,eps2,ta,tb)
extern "C" double bornv_dizet_(const int&, const int&, const int&,
		const double&, const double&, const double&, const double&, const double&, const double& );
//      SUBROUTINE BornV_GetGammZ(GammZ)
extern "C" void bornv_getgammz_( const double& );
//      SUBROUTINE BornV_GetMZ(MZ)
extern "C" void bornv_getmz_( const double& );
//      SUBROUTINE BornV_GetQEDcor(QEDcor)
extern "C" void bornv_getqedcor_( const double& );
//------------------------------------------------------------
//      SUBROUTINE GPS_BornF(KFi,KFf,PX,CosThe,p1,m1,p2,m2,p3,m3,p4,m4,Xborn)
extern "C" void gps_bornf_(const int&, const int&, double[], const double&,
		    double[], const double&, double[], const double&, double[], const double&, double[], const double&,
		    const double&);
//------------------------------------------------------------
//      DOUBLE PRECISION  FUNCTION GPS_MakeRhoFoam(XNorm)
extern "C" double gps_makerhofoam_(const double&);

class KKplot{
// Interface and extensions to KKsem toolbox
 public:
    int       m_jmax;
    double    m_ypar[10000];   // input parameters of KKMC
//
 public:
 	double m_CMSene;
 	double m_Mmin;
 	double m_Mmax;
 	double m_vvmax;
 	//
 	double m_gnanob;
 	double m_pi;
 	double m_ceuler;
 	double m_alfinv;
 	double m_alfpi;
 	//
 	double m_beam;
 	double m_chini;
 	//
 	double m_fin;
 	double m_chfin;
 	int   m_KFini;  // electron
 	int   m_KFf;    // muon
 	//
 	int    m_kDim;
 	int    m_nCells;
 	int    m_nSampl;
 	int    m_KeyISR;
 	int    m_KeyFSR;
//
//
//------ constructors destructors -------
 public:
  KKplot(){;}
  ~KKplot(){;}
  KKplot(const char* Name);
public:
  // Interfaces to KKsem integration routines using Gauss method
  void Initialize(TFile&);
  void Initialize(double[]);
  void VVplot( TH1D *&, int , char [], int, int );
  void VVmake( TH1D *&v, TH1D *&, int, char[], int, int, double);
  void Cplot(  TH1 *, int , char [], int, int, double, double);
  void ReaData(const char*, int, double[]);
  void Ord1fill( TH1D*&, int );

////////////////////////////////////////////////////////////////////////////
};// KKplot



#endif
