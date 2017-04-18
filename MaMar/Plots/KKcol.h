#ifndef KKcol_H
#define KKcol_H
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKcol                                              //
//                                                                          //
//   Interface (wrapper)  to KKsem package for semi-analytical programs     //
//   Some extras on toop of ROOT histogramming package                      //
//   and more...
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
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
//      DOUBLE PRECISION  FUNCTION BornV_Sig0nb(CMSene)
extern "C" double bornv_sig0nb_(const double&);
//      SUBROUTINE BornV_MakeGami(CMSene,gamiCR,gami,alfi)
extern "C" void bornv_makegami_(const double&, double&, double&, double&);
//      DOUBLE PRECISION  FUNCTION BornV_Simple(KFi,KFf,svar,costhe)
extern "C" double bornv_simple_( const long&,  const long&, const double&, const double&);
//------------------------------------------------------------


class KKcol: public TFoamIntegrand{
 public:
    int       m_jmax;
    double    m_ypar[10000];   // input parameters of KKMC
// Latex output from ROOT histograms
	int       m_lint;         // type of input
//------ constructors destructors -------
  KKcol(){;};
  ~KKcol(){;};
 public:
  // Interfaces to KKsem integration routines using Gauss method
  void Initialize(TFile&);
  void VVplot( TH1 *, long , char [], long, long );
  void Cplot(  TH1 *, long , char [], long, long, double, double);
  // Latex output from ROOT histograms
  void PlInitialize(FILE *, int );
  void PlEnd(FILE *);
  void PlTable2(int, TH1D* iHst[], FILE*, Char_t *Capt[], Char_t[] , const char* , int,int,int);
  Double_t Density(int nDim, Double_t *Xarg){;};

////////////////////////////////////////////////////////////////////////////
};// KKsem

////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
class Rho4Foam: public TFoamIntegrand{
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

	double m_beam;
	double m_chini;

	double m_fin;

	int    m_kDim;
	int    m_nCells;
	int    m_nSampl;
	int    m_KeyISR;
	int    m_KeyFSR;

///******** MC EVENT ********
	double m_x1;
	double m_x2;
	double m_vv;
	double m_uu;
	double m_Mll;
	double m_Mka;

////////////////////////////////////////////////////////////////////////////////
//
// Class defining/providing integrand for Foam
//
////////////////////////////////////////////////////////////////////////////////
public:
	///_____________________________________________________________
Rho4Foam(const char* Name)
	{
	//! all defaults defined here can be changed by the user
	//! before calling TMCgen::Initialize

	  m_gnanob  = 389.37966e3;
	  m_pi      = 3.1415926535897932;
	  m_ceuler  = 0.57721566;
	  m_alfinv  = 137.035;
	  m_alfpi   = 1/m_alfinv/m_pi;
	  m_vvmax   = 0.20;
//{{{{{{{{{{{
//	  m_beam    = 0.510999e-3;  // electron
//	  m_chini   = 1.0;          // electron

//	  m_beam    = 0.005;        // u quark
//	  m_chini   = 2./3.;        // u quark

	  m_beam    = 0.010;        // d quark
	  m_chini   = 1./3.;        // d quark

	  m_fin     = 0.105;        // final ferm. muon
//	  m_fin     = 1.777;        // final ferm. tau

	  m_kDim    =    3;         // No. of dim. for Foam, =2,3 Machine energy spread OFF/ON
	  m_nCells  = 2000;         // No. of cells, optional, default=2000
	  m_nSampl  =  200;         // No. of MC evts/cell in exploration, default=200

	  m_KeyISR  = 2;            // Type of ISR/QED switch, 0,1,2
//}}}}}}}}}}}
	cout<< "----> Rho4Foam USER Constructor "<<endl;
	}///

//
//double sqr( const double_t x ){ return x*x;};
///------------------------------------------------------------------------
double gamISR( double svar){
	  return  sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
}

///------------------------------------------------------------------------
double gamFSR( double svar){
	  return              2*m_alfpi*( log(svar/sqr(m_fin)) -1);
}

///------------------------------------------------------------------------
double Rho_isr(double svar, double vv){
/// ISR rho-function for ISR

  double alf1   = m_alfpi;
  double gami   = gamISR(svar);
  //gami = sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
///
  double gamfac = exp(-m_ceuler*gami)/TMath::Gamma(1+gami);
  double delb   = gami/4 +alf1*(-0.5  +sqr(m_pi)/3.0);
  double ffact  = gamfac*exp(delb);

  double rho,dels,delh;
  if(       m_KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = 0;
	rho  = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 1){
/// first  order
	dels = gami/2;   /// NLO part =0 as for vector boson???
    delh = vv*(-1 +vv/2);
    rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 2){
/// second order without NLO part
    dels = gami/2 +sqr(gami)/8;
    delh = vv*(-1+vv/2.0)
          +gami*0.5*(-0.25*(4.0-6.0*vv+3.0*vv*vv)*log(1-vv)-vv);
    rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenH::Rho4Foam: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
  }
///
  return rho;
}//Rho_isr



///------------------------------------------------------------------------
double Rho_fsr(double svar, double uu){
/// ISR+FSR rho-function

  double alf1   = m_alfpi;
  double gamf   = gamFSR(svar*(1-uu));
///
  /////delb   = Chf2* alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
  /////delb   = delb -betf/2 *dlog(1-uu)
  double gamfac = exp(-m_ceuler*gamf)/TMath::Gamma(1+gamf);
  double delb   = gamf/4 +alf1*(-0.5  +sqr(m_pi)/3.0)
		         -gamf/2 *log(1-uu);
  double ffact  = gamfac*exp(delb);

  double rho,dels,delh;
  if(       m_KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = 0;
	rho  = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 1){
/// first  order
	dels = gamf/2;   /// NLO part =0 as for vector boson???
    delh = uu*(-1 +uu/2);
    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 2){
/// second order without NLO part
    /////dels  = betf/2d0 +betf**2/8d0
    /////delh  = uu*(-1d0+uu/2d0)
    /////$        +betf*(-0.5d0*uu-0.25d0*uu*(-1d0+0.5d0*uu)*log(1d0-uu))
    dels = gamf/2 +sqr(gamf)/8;
    delh = uu*(-1+uu/2.0)
          +gamf*0.5*( -0.5*uu -0.25*uu*(-1.0 +0.5*uu)*log(1-uu));
    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenH::Rho4Foam: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
  }
///
  return rho;
}//Rho_fsr

///////////////////////////////////////////////////////////////
Double_t Density(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	Double_t Dist=1;
	double xmin=0.000001;
	double xmax=0.999999;

	double svar = sqr(m_CMSene);
	double svarCum = svar;
// ******** mapping for PDFs *******
	m_x1 = xmin+(xmax-xmin)*Xarg[0];
	m_x2 = xmin+(xmax-xmin)*Xarg[1];
	Dist *= sqr(xmax-xmin);

	// Valence 2*x*u(x):   XUPV = 2.18   * X**0.5D0    * (1.D0-X)**3.D0
	// Valence x*d(x):     XDNV = 1.23   * X**0.5D0    * (1.D0-X)**4.D0
	// Sea     x*s(x):     XSEA = 0.6733 * X**(-0.2D0) * (1.D0-X)**7.D0
	//                Valence UP and UP-bar of sea
    //SF12 = 2.18 *m_x1**3.D0 *z1**0.5D0   *0.6733 *m_x2**7.D0 *z2**(-0.2D0)/z1/z2
	double PDFu1   = 2.18   *exp(3.0 *log(1-m_x1))  *exp( 0.5*log(m_x1))/m_x1;
	double PDFd1   = 1.23   *exp(4.0 *log(1-m_x1))  *exp( 0.5*log(m_x1))/m_x1;
    double PDFsea1 = 0.6733 *exp(7.0 *log(1-m_x1))  *exp(-0.2*log(m_x1))/m_x1;
    double PDFsea2 = 0.6733 *exp(7.0 *log(1-m_x2))  *exp(-0.2*log(m_x2))/m_x2;
    //{{{{
    //double SF12  = 2*(PDFu1+ PDFsea1/6.) * (PDFsea2/6.);  //  u-ubar
    double SF12  = 2*(PDFd1+ PDFsea1/6.) * (PDFsea2/6.);  //  d-dbar
    //}}}}
    Dist *= SF12;  // (u-ubar)+(ubar-u)
	double svar1= svar*m_x1*m_x2;
	svarCum *=m_x1*m_x2;

// ******** mapping for ISR *******
	double gamiCR,gami,alfi;
	double CMSene1= sqrt(svar1);
	bornv_makegami_( CMSene1, gamiCR,gami,alfi);   // from KKMC
	//[[[[ debug
	//gami = gamISR(CMSene1);
	//]]]]
	// cout<<" CMSene1,gami= "<< CMSene1 <<"  "<< gami <<endl;
	double R= Xarg[2];
	m_vv = exp(1.0/gami *log(R)) *m_vvmax; // mapping
	Dist *= m_vv/R/gami ;                  // Jacobian
	//[[[[ primitive for debug
	//m_vv = R*m_vvmax;
	//Dist *= m_vvmax;
	//cout << "Density: svar1, gami = "<< svar1 <<"  "<<gami<<endl;
	//double Rho1 = gami*exp(gami*log(m_vv))/m_vv; // ISR distribution
	//]]]]
	if( gami < 0 )      return 0.0;    // temporary fix
    if( m_vv < 1e-200 ) return 0.0;    // temporary fix
	double Rho2 = Rho_isr(svar1,m_vv);               // remember take care of m_mbeam!!!
	Dist *= Rho2;
	svarCum *= (1-m_vv);
	double svar2 = svar1*(1-m_vv);

// ******** mapping for FSR *******
    if( m_KeyFSR > 0){
      double rr= Xarg[3];
      double gamf   = gamFSR(svar2);
      m_uu = exp(1.0/gamf *log(rr));     // mapping
      Dist *= m_uu/rr/gamf;              // Jacobian

  	  double Rho3 = Rho_fsr(svar2,m_uu);           // remember take care of m_mbeam!!!
  	  Dist *= Rho3;
      svarCum *= (1-m_uu);

    }// m_KeyFSR

// ******** finally Born factor *******
	long KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
	kksem_setkeyfob_( KeyFob );
	double xBorn;
	kksem_makeborn_( svar2, xBorn);
	xBorn *= 1./3.;  // colour factor corrected by hand
//[[[[[ debug
//  double Sig0nb= bornv_sig0nb_(m_CMSene);
//	xBorn = bornv_simple_( 2, 13,m_Mll*m_Mll,0e0);
//	xBorn *= 1/(1e0-m_vv)/(m_x1*m_x2)*Sig0nb;
//]]]]
	Dist *= xBorn;

// ******* inline cutoff for better efficiency *********
	m_Mll = sqrt(svarCum); // final
	m_Mka = sqrt(svar2);   // after ISR
	if( m_Mka < m_Mmin ) Dist=0;
	//if( m_Mll < m_Mmin ) Dist=0;
	//if( m_Mka < m_Mmin || m_Mka >m_Mmax ) Dist=0;
	//if( m_Mll < m_Mmin || m_Mll >m_Mmax ) Dist=0;

	return Dist;
}// Density
};// class Rho4Foam


#endif
