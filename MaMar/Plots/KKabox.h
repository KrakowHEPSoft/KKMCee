//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKabox                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// KKabox is multipurpose toolbox for KKMC testing.
//  1. Interfaces (wrappers) to KKMC and KKsem F77 subrograms
//  2. Integrand for Foam in semianalytical xcheck
//  3. A few routines for producing latex table out of histograms
//////////////////////////////////////////////////////////////////////////////

#ifndef KKabox_H
#define KKabox_H
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

#include "TRandom3.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"

#include "HisNorm.h"

//------------------------------------------------------------
//  wrappers to f77 routines ik KKMC
extern "C" void kk2f_fort_open_( const long&, char*, long);
extern "C" void kk2f_fort_close_(const long&);
//      SUBROUTINE KKsem_Initialize(xpar_input)
extern "C" void kksem_initialize_(double xpar[]);
//      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
//-----------------------
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
//---------------------------------
//      DOUBLE PRECISION  FUNCTION BornV_Sig0nb(CMSene)
extern "C" double bornv_sig0nb_(const double&);
//      SUBROUTINE BornV_MakeGami(CMSene,gamiCR,gami,alfi)
extern "C" void bornv_makegami_(const double&, double&, double&, double&);
//      DOUBLE PRECISION  FUNCTION BornV_Simple(KFi,KFf,svar,costhe)
extern "C" double bornv_simple_( const long&,  const long&, const double&, const double&);
//------------------------------------------------------------
//      SUBROUTINE BornV_InterpoGSW(KFf,svar,CosThe)
extern "C" double bornv_interpogsw_( const long&,  const double&, const double&);
//      DOUBLE PRECISION FUNCTION BornV_Dizet(Mode,KFi,KFf,svar,CosThe,eps1,eps2,ta,tb)
extern "C" double bornv_dizet_(const long&, const long&, const long&,
		const double&, const double&, const double&, const double&, const double&, const double& );
//------------------------------------------------------------
//      SUBROUTINE GPS_BornF(KFi,KFf,PX,CosThe,p1,m1,p2,m2,p3,m3,p4,m4,Xborn)
extern "C" void gps_bornf_(const long&, const long&, double[], const double&,
		    double[], const double&, double[], const double&, double[], const double&, double[], const double&,
		    const double&);
//------------------------------------------------------------
//      SUBROUTINE GPS_BornFoam(Mode,KFi,KFf,CMSene,CosThe,Xborn)
extern "C" void gps_bornfoam_(const long&,   const long&,   const long&,
		                      const double&, const double&, const double&);
//      DOUBLE PRECISION  FUNCTION GPS_MakeRhoFoam(XNorm)
extern "C" double gps_makerhofoam_(const double&);

class KKabox: public TFoamIntegrand{
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
 	long   m_KFini;  // electron
 	long   m_KFf;    // muon
 	//
 	int    m_kDim;
 	int    m_nCells;
 	int    m_nSampl;
 	int    m_KeyISR;
 	int    m_KeyFSR;
//
 	int    m_Mode;   // operation mode for Density
 	double m_del;
//******** MC EVENT ********
 	double m_CosTheta;
 	double m_vv;  // ISR
 	double m_uu;  // FSR
 	double m_r1;  // IFI
 	double m_r2;  // IFI
 	double m_xx;  // total
 	//
 	double m_Mka;
 	//
 	double m_p1[4];
 	double m_p2[4];
 	double m_p3[4];
 	double m_p4[4];
 	//
 	long   m_count;
//
//[[[// Latex output from ROOT histograms
//[[[	int       m_lint;         // type of input
//------ constructors destructors -------
 public:
  KKabox(){;}
  ~KKabox(){;}
  KKabox(const char* Name);
public:
  // Interfaces to KKsem integration routines using Gauss method
  void Initialize(TFile&);
  void VVplot( TH1 *, long , char [], long, long );
  void Cplot(  TH1 *, long , char [], long, long, double, double);

  // Latex output from ROOT histograms
  void PlInitialize(FILE *, int );
  void PlEnd(FILE *);
  void PlTable2(int, TH1D* iHst[], FILE*, Char_t *Capt[], Char_t[] , const char* , int,int,int);

  // Foam integrand
  double Fyfs( double );
  double gamISR( double );
  double gamFSR( double );
  double gamIFI( double );
  double Rho_isr(double, double );
  double Rho_fsr(double, double );
  double Rho_ifi(double, double , double );
  void MapIFI1( double, double, double, double &, double &);
  void MapIFI2( double, double, double, double &, double &);
  void MapIFI(  double, double, double, double &, double &);
  Double_t Density(int, Double_t*);
  Double_t Density3(int, Double_t*);
  Double_t Density5(int, Double_t*);

////////////////////////////////////////////////////////////////////////////
};// KKabox



#endif
