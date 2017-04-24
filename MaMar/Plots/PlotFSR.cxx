//////////////////////////////////////////////////////////////////////
//    make PlotFSR
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include <math.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"

#include "TFoam.h"
#include "TFile.h"
#include "TApplication.h"
#include "TRandom3.h"
#include "TFoamIntegrand.h"

#include "KKMC.h"
#include "KKabox.h"
#include "HisNorm.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
//{{{{ additional input in the code //}}}}
TFile DiskFileA("../workFSR/rmain_ddBeams_FSRon_100M.root");            // ddbar->mu-mu+
//TFile DiskFileA("../workFSR/rmain_uuBeams_FSRon_1G.root");            // uubar->mu-mu+
//TFile DiskFileA("../workFSR/rmain_uuBeams_FSRon_tau_1G.root"); // uubar->tau-tau+
//TFile DiskFileA("../workFSR/rmain.root");
TFile DiskFileB("PlotFSR.root","RECREATE","histograms");
//=============================================================================

//double sqr( const double_t x ){ return x*x;};
// Auxiliary procedures for plotting
//#include "Marker.h"

KKabox LibSem;

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

////////////////////////////////////////////////////////////////////////////////
//
//  Quick Monte Carlo exercise based on Foam
//  Two PDFs, QED and Born from KKMC
//
////////////////////////////////////////////////////////////////////////////////
//______________________________________________________________________________
void ISRgener()
{
  cout<<"--- ISRgener started ---"<<endl;

  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  double vvmax   = HST_KKMC_NORMA->GetBinContent(17);

  DiskFileB.cd();
  TH1D *hst_Mll_main  = (TH1D*)DiskFileA.Get("hst_Mll_main"); // for cloning only
  TH1D *HST_Mka = (TH1D*)hst_Mll_main->Clone("HST_Mka");
  HST_Mka->Reset();
  TH1D *HST_Mll = (TH1D*)hst_Mll_main->Clone("HST_Mll");
  HST_Mll->Reset();

  Rho4Foam *Rho1= new Rho4Foam("Rho4Foam");
  double Mmin= 60;
  double Mmax=160;
  Rho1->m_Mmin   = Mmin;
  Rho1->m_Mmax   = Mmax;
  Rho1->m_CMSene = CMSene;
  Rho1->m_vvmax  = vvmax;
  Rho1->m_KeyISR = 2;
  cout<<"%%%%%   ISRgener:  vvmax= "<< vvmax << " from KKMC" <<endl;
//------------------------------------------
  TRandom  *PseRan   = new TRandom3();  // Create random number generator
  PseRan->SetSeed(4357);
//------------------------------------------------------------------
// WARNIBG!!! MC_isr->Density() is common for both FOAM objects
// hence MC_isr->m_KeyFSR has to be adjusted before every Foam call!
  //----------------------------------------------------------------
// Setting up Foam objects !st FOAM for ISR only
  TFoam   *MC_isr    = new TFoam("MC_isr");   // Create Simulator
  MC_isr->SetkDim(3);         // No. of dimensions, obligatory!
  MC_isr->SetnCells(  5000);  // No. of cells, can be omitted, default=2000
  MC_isr->SetnSampl(100000);  // No. of MC evts/cell in exploration, default=200
  MC_isr->SetRho(Rho1);       // Set 2-dim distribution, included above
  MC_isr->SetPseRan(PseRan);  // Set random number generator, mandatory!
  MC_isr->SetOptRej(0);       // wted events (=0), default wt=1 events (=1)
  Rho1->m_KeyFSR = 0;
  MC_isr->Initialize();       // Initialize simulator, may take time...
// 2nd FOAM for ISR+FSR
  TFoam   *MC_fisr    = new TFoam("MC_fisr");   // Create Simulator
  MC_fisr->SetkDim(4);         // No. of dimensions, obligatory!
  MC_fisr->SetnCells( 10000);   // No. of cells, can be omitted, default=2000
  MC_fisr->SetnSampl(100000);   // No. of MC evts/cell in exploration, default=200
  MC_fisr->SetRho(Rho1);       // Set 2-dim distribution, included above
  MC_fisr->SetPseRan(PseRan);  // Set random number generator, mandatory!
  MC_fisr->SetOptRej(0);       // wted events (=0), default wt=1 events (=1)
  Rho1->m_KeyFSR = 2;
  MC_fisr->Initialize();       // Initialize simulator, may take time...

// loop over MC events
  double wt,Mll, wt2, Mka;
  long NevTot = 2000000;  // 2M
  NevTot =      8000000;  // 8M
  //NevTot =    100000000;  // 100M
  for(long loop=0; loop<NevTot; loop++)
  {
	//----------------------------------------------------------
	/// Generate ISR event
	Rho1->m_KeyFSR = 0;  ///!!!
    MC_isr->MakeEvent();            // generate MC event
    MC_isr->GetMCwt(wt);
    Mka=Rho1->m_Mka;
    HST_Mka->Fill(Mka,wt);

	//----------------------------------------------------------
	/// Generate ISR+FSR event
	Rho1->m_KeyFSR = 2;  ///!!!
    MC_fisr->MakeEvent();            // generate MC event
    MC_fisr->GetMCwt(wt2);
    Mll = Rho1->m_Mll;
    HST_Mll->Fill(Mll,wt2);
    if( 1000000*(loop/1000000) == loop) cout<<" Nev ="<< loop<< endl;
  }// loop
//renormalizing histo
  double Xsav, dXsav;
  MC_isr->GetIntNorm(Xsav,dXsav);
  HisNorm0( NevTot, Xsav, HST_Mka);
  //
  MC_fisr->GetIntNorm(Xsav,dXsav);
  HisNorm0( NevTot, Xsav, HST_Mll);
//
  cout<<"--- ISRgener ended ---"<<endl;
}//ISRgener


///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueMain") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vAlepCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vXGenCeex2") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_s1Ceex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_svk") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Mll_main") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Mll_eex0") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Mll_eex2") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Mka_eex0") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Mka_eex2") );
}

///////////////////////////////////////////////////////////////////////////////////
void KKsemMake(){
// Here we produce semianalytical plots using KKsem program, No plotting
//------------------------------------------------------------------------  
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"xxxxxxxxxxxxxxxx KKsemMake  BEGIN xxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  // initilalization of KKcol
  //KKcol LibSem;
  LibSem.Initialize(DiskFileA);
  //{{{{{{
  long KF=13; // muon ?????
  //}}}}}}
  long KeyDis, KeyFob;
  char chak[5];
  //KeyDis = 302;   // ISR O(alf2)
  //KeyDis = 304;   // ISR O(alf3) GribovLL
  //KeyDis = 303;   // ISR O(alf3)
  //KeyDis = 305;   // ISR O(alf3) GribovLL +NLL
  //KeyDis = 662;  // Unexp ????
  //KeyDis = 302302;   // ISR*FSR O(alf3)
  //
  KeyFob=   10; // BornV_Dizet, with EW and without integration ???
  KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
//  KeyFob=  -10; // KKsem_BornV, NO EW, NO integration OK!
//  KeyFob= -100; // KKsem_BornV, NO EW, WITH integration, OK
//  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!

  kksem_setkeyfob_( KeyFob );
  double svar= 500*500;
  double xBorn;
  kksem_makeborn_( svar, xBorn);
  cout<< "xBorn [nb]= "<<xBorn<<endl;

  cout<<"xxxxxxxxxxxxxxxx KKsemMake END xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
//------------------------------------------------------------------------  
//------------------------------------------------------------------------  
}//  KKsemMake


///////////////////////////////////////////////////////////////////////////////////
void FigCalib()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCalib =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  char captEne[100];
  //{{{{{{
  sprintf(captEne,"#sqrt{s} =%4.0fGeV, d#bar{d} -> #mu^{+}#mu^{-}", CMSene);
  //sprintf(captEne,"#sqrt{s} =%4.0fGeV, u#bar{u} -> #mu^{+}#mu^{-}", CMSene);
  //sprintf(captEne,"#sqrt{s} =%4.0fGeV, u#bar{u} -> #tau^{+}#tau^{-}", CMSene);
  //}}}}}}
  Long_t   Nevt = HST_KKMC_NORMA->GetEntries();
  Double_t Xsav = HST_KKMC_NORMA->GetBinContent(0)/Nevt; // NANOBARNS

  TH1D *hst_Mll_main   = (TH1D*)DiskFileA.Get("hst_Mll_main"); // From KKMC
  TH1D *hst_Mll_eex0   = (TH1D*)DiskFileA.Get("hst_Mll_eex0"); // From KKMC
  TH1D *hst_Mll_eex2   = (TH1D*)DiskFileA.Get("hst_Mll_eex2"); // From KKMC
  TH1D *hst_Mka_eex0   = (TH1D*)DiskFileA.Get("hst_Mka_eex0"); // From KKMC ISR only
  TH1D *hst_Mka_eex2   = (TH1D*)DiskFileA.Get("hst_Mka_eex2"); // From KKMC ISR only

  TH1D *HST_Mll        = (TH1D*)DiskFileB.Get("HST_Mll");      // From FOAM (cloned)
  TH1D *HST_Mka        = (TH1D*)DiskFileB.Get("HST_Mka");      // From FOAM (cloned)

  //***** Misccelaneous xchecks
  double Sig0nb= bornv_sig0nb_(CMSene);
  cout<< "%%%%%%%%%%%%%%%%%%%FigCalib%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
  cout<< "%%%%% Pointlike Sig0nb(CMSene)  = " << Sig0nb << endl;
  cout<< "%%%%%%%%%% Crude from KKMC "<<endl;
  cout<< "%%%%% Xsav = "<< 1000*Xsav <<" [pb],   Nevt =  "<<Nevt<<endl;
  //cout<< "%%%%% Xsav/8 [pb] = " << 1000*Xsav/8  <<endl;
  cout<< "%%%%% Xsav/Sig0nb   [R] = " << Xsav/Sig0nb << endl;
  //cout<< "%%%%% Xsav/Sig0nb/8 [R] = " << Xsav/Sig0nb/8. << endl;

  //****** integrate over bins histo from KKMC
  TH1D *Hst1, *Hst2, *Hst3, *Hst4, *HST1, *HST3;
  Hst2 =hst_Mll_main;
  //Hst1 = hst_Mka_eex0; // ISR only
  Hst1 = hst_Mka_eex2;   // KKMC ISR only
  Hst3 = hst_Mll_eex2;   // KKMC ISR+FSR
  HST1 = HST_Mka;        // FOAM ISR only
  HST3 = HST_Mll;        // FOAM ISR only

  /// Getting integrated xsection out of Hst1 (KKMC, ISR only)
  int      nbX  = Hst1->GetNbinsX();
  Double_t Xmax = Hst1->GetXaxis()->GetXmax();
  Double_t Xmin = Hst1->GetXaxis()->GetXmin();
  double dx = (Xmax-Xmin)/nbX;
  double xsum1 = 0, xssum1=0;
  for(int ix=1; ix <= nbX; ix++){
	 xsum1   += Hst1->GetBinContent(  ix ) *dx;
	 xssum1  +=  sqr(Hst1->GetBinError(  ix ) *dx );
//	 cout<< "ix="<< ix <<"  xsum1="<< xsum<<endl;
  }//ix
  char capt1[100];
  xssum1 = sqrt(xssum1);
  cout<< "%%%%%%%%%  KKMC from histo, presently EEX0, %%%%%%%%%%%%%%%%"<<endl;
  cout<< "%%%%% SigKKMC [pb]= " << 1000*xsum1 <<" +- "<< xssum1 << endl;
  cout<< "%%%%% SigKKMC [R] = " << xsum1/Sig0nb << endl;
  sprintf(capt1,"#sigma_{KKMC} =%9.3f +- %7.3f [pb]", 1000*xsum1, 1000*xssum1);

  /// Getting integrated xsection out of Hst0 (FOAM, ISR only)
  nbX  = HST1->GetNbinsX();
  Xmax = HST1->GetXaxis()->GetXmax();
  Xmin = HST1->GetXaxis()->GetXmin();
  dx = (Xmax-Xmin)/nbX;
  double xsum2 = 0, xssum2 = 0;
  for(int ix=1; ix <= nbX; ix++){
	 xsum2  += HST1->GetBinContent(  ix ) *dx;
	 xssum2 += sqr(HST1->GetBinError(  ix ) *dx );
  }//ix
  xssum2 = sqrt(xssum2);
  char capt2[100];
  sprintf(capt2,"#sigma_{FOAM} =%9.3f +- %7.3f [pb]", 1000*xsum2, 1000*xssum2);
  cout<< "%%%%%%%%%   From external Foam %%%%%%%%%%%%%%%%"<<endl;
  cout<< "%%%%% SigFoam [pb]= " << 1000*xsum2 << " +- " << 1000*xssum2 << endl;
  cout<< "%%%%% SigFoam [R] = " << xsum2/Sig0nb << endl;
  cout<< "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
//
  TH1D *HST_One  = new TH1D("HST_one" , "One",1, Xmin, Xmax );
  HST_One->SetBinContent(1,1);

//------------------------------------------------------------------------
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCalib = new TCanvas("cFigCalib","FigCalib: general info ", 50, 100,    1000,  800);
  //                                      Name    Title               xoff,yoff, WidPix,HeiPix
  cFigCalib->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigCalib->Divide( 2,  2);
  //cFigCalib->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //=====================plot1=========================
  cFigCalib->cd(1);
  gPad->SetLogy(); // !!!!!!

  //Hst1->SetStats(0);
  Hst1->SetTitle(0);

  Hst1->SetLineColor(kBlue);
  Hst1->DrawCopy("h");

  HST1->SetLineColor(kRed);
  HST1->DrawCopy("hsame");

  CaptT->DrawLatex(0.05,0.95,"d#sigma/dM [nb/GeV] KKMC and FOAM (red); ISR only");
  CaptT->DrawLatex(0.45,0.70, captEne);
  CaptT->DrawLatex(0.45,0.60, capt1);
  CaptT->DrawLatex(0.45,0.65, capt2);

  //=====================plot2=========================
  cFigCalib->cd(2);
  //  ISR only EEX
  TH1D *Hst10_Rat =(TH1D*)Hst1->Clone("Hst10_Rat");
  Hst10_Rat->Divide(HST1);
  Hst10_Rat->SetMinimum(0.90);
  Hst10_Rat->SetMaximum(1.10);
  Hst10_Rat->DrawCopy("h");
  HST_One->DrawCopy("hsame");

  CaptT->DrawLatex(0.10,0.95,"Ratio:  KKMC/FOAM; ISR only ");

  //=====================plot3=========================
  cFigCalib->cd(3);
  gPad->SetLogy(); // !!!!!!

  //Hst3->SetStats(0);
  Hst3->SetTitle(0);

  Hst3->SetLineColor(kBlue);
  Hst3->DrawCopy("h");

  HST3->SetLineColor(kRed);
  HST3->DrawCopy("hsame");

  CaptT->DrawLatex(0.05,0.95,"d#sigma/dM [nb/GeV] KKMC and FOAM(red), ISR+FSR");

  //=====================plot4=========================
  cFigCalib->cd(4);

  TH1D *Hst30_Rat =(TH1D*)Hst3->Clone("Hst30_Rat");
  Hst30_Rat->Divide(HST3);
  Hst30_Rat->SetMinimum(0.90);
  Hst30_Rat->SetMaximum(1.10);
  Hst30_Rat->DrawCopy("h");
  HST_One->DrawCopy("hsame");

  CaptT->DrawLatex(0.10,0.95,"Ratio:  KKMC/FOAM, ISR+FSR");

//----------------------------
  cFigCalib->cd();
  cout<<" ========================= FigCalib END ======================= "<<endl;
}// FigCalib


///////////////////////////////////////////////////////////////////////////////////
void FigFSR()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigFSR =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR

  TH1D *hst_Mll_eex2   = (TH1D*)DiskFileA.Get("hst_Mll_eex2"); // From KKMC ISR+FSR
  TH1D *hst_Mka_eex2   = (TH1D*)DiskFileA.Get("hst_Mka_eex2"); // From KKMC ISR only
  TH1D *Hst1= hst_Mll_eex2;
  TH1D *Hst2= hst_Mka_eex2;

  long   nbX  = Hst1->GetNbinsX();
  double Xmax = Hst1->GetXaxis()->GetXmax();
  double Xmin = Hst1->GetXaxis()->GetXmin();
  //
  TH1D *HST_One2  = new TH1D("HST_one2" , "One2",1, Xmin, Xmax );
  HST_One2->SetBinContent(1,1);

//------------------------------------------------------------------------
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigFSR = new TCanvas("cFigFSR","FigFSR: general info ", 100, 150,    1000,  800);
  //                                Name    Title                  xoff,yoff, WidPix,HeiPix
  cFigFSR->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigFSR->Divide( 2,  2);
  //cFigFSR->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);

  //==========plot1==============
  cFigFSR->cd(1);
  gPad->SetLogy(); // !!!!!!

  hst_Mll_eex2->DrawCopy("h");
  hst_Mka_eex2->SetLineColor(kRed);
  hst_Mka_eex2->DrawCopy("hsame");
  CaptT->DrawLatex(0.05,0.95,"d#sigma/dM [nb/GeV] KKMC: FSR off (red) and on (blue)");

  //==========plot2==============
  cFigFSR->cd(2);

  TH1D *Hst12_Rat =(TH1D*)Hst1->Clone("Hst12_Rat");
  Hst12_Rat->Divide(Hst2);
  Hst12_Rat->SetMinimum(0.8);
  Hst12_Rat->SetMaximum(2.6);
  Hst12_Rat->DrawCopy("h");
  HST_One2->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"KKMC 2nd Order: (ISR+FSR)/ISR");

  //==========plot3==============
  cFigFSR->cd(3);
  gPad->SetLogy(); // !!!!!!


  //==========plot4==============
  cFigFSR->cd(4);

//----------------------------
  cFigFSR->cd();
}// FigFSR



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  KKsemMake();        //
  //========== PLOTTING ==========

  ISRgener();

  FigCalib();

  FigFSR();

  //++++++++++++++++++++++++++++++++++++++++
  cout<< "###############  root file of MC       ###################"<<endl;
  DiskFileA.ls();
  cout<< "###############  root file of analysis ###################"<<endl;
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
