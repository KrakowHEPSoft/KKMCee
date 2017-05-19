//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   TMCgenFOAM                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// TMCgenFOAM is multipurpose toolbox for KKMC testing.
//  1. Interfaces (wrappers) to KKMC and KKsem F77 subrograms
//  2. Integrand for Foam in semianalytical xcheck
//  3. A few routines for producing latex table out of histograms
//////////////////////////////////////////////////////////////////////////////

#ifndef TMCgenFOAM_H
#define TMCgenFOAM_H
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

#include "TMCgen.h"

#include "TRandom3.h"
//#include "TFoam.h"
//#include "TFoamIntegrand.h"

//------------------------------------------------------------
//  wrappers to f77 routines in KKMC and KKsem
extern "C" void kk2f_fort_open_( const long&, char*, long);
extern "C" void kk2f_fort_close_(const long&);
//      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
//-----------------------
//      SUBROUTINE BornV_SetKF(KFferm)
extern "C" void bornv_setkf_( const long& ); // set SINGLE Final State
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

///////////////////////////////////////////////////////////////////
class TMCgenFOAM: public TMCgen{
// Interface and extensions to KKsem toolbox
//------ constructors destructors -------
 public:
    TMCgenFOAM();                // explicit default constructor for streamer
    TMCgenFOAM(const char*);     // user constructor
    ~TMCgenFOAM();               // explicit destructor
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
/// member functions
public:
  // Interfaces to KKsem integration routines using Gauss method
  void Initialize( double ypar[10000]);
//  void VVplot( TH1 *, long , char [], long, long );
//  void Cplot(  TH1 *, long , char [], long, long, double, double);

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
  ClassDef(TMCgenFOAM,2); // This is for ROOT persistency
////////////////////////////////////////////////////////////////////////////
};// TMCgenFOAM


#endif
