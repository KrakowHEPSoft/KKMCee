//////////////////////////////////////////////////////////////////////
//    make Plot1
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


#include "KKsem.h"

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
TFile DiskFileA("../workFSR/rmain.root");
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
//=============================================================================

double sqr( const double_t x ){ return x*x;};
// Auxiliary procedures for plotting
#include "HisNorm.h"
#include "Marker.h"

//_____________________________________________________________________________
class RhoISR: public TFoamIntegrand{
public:
	double m_CMSene;
	double m_Mmin;
	double m_Mmax;
	double m_Mll;
	double m_x1;
	double m_x2;
public:
Double_t Density(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	Double_t Dist=1;
	double xmin=0.000001;
	double xmax=0.999999;

	m_x1 = xmin+(xmax-xmin)*Xarg[0];
	m_x2 = xmin+(xmax-xmin)*Xarg[1];
	Dist *= sqr(xmax-xmin);

	// Valence 2*x*u(x):   XUPV = 2.18   * X**0.5D0    * (1.D0-X)**3.D0
	// Valence x*d(x):     XDNV = 1.23   * X**0.5D0    * (1.D0-X)**4.D0
	// Sea     x*s(x):     XSEA = 0.6733 * X**(-0.2D0) * (1.D0-X)**7.D0
	//                Valence UP and UP-bar of sea
    //SF12 = 2.18 *m_x1**3.D0 *z1**0.5D0   *0.6733 *m_x2**7.D0 *z2**(-0.2D0)/z1/z2
	double PDFu1   = 2.18   *exp(3.0 *log(1-m_x1))  *exp( 0.5*log(m_x1))/m_x1;
    double PDFsea1 = 0.6733 *exp(7.0 *log(1-m_x1))  *exp(-0.2*log(m_x1))/m_x1;
    double PDFsea2 = 0.6733 *exp(7.0 *log(1-m_x2))  *exp(-0.2*log(m_x2))/m_x2;
    Dist *= (PDFu1+ PDFsea1/6.) * (PDFsea2/6.);  //  u-ubar
    Dist *= 2;  // (u-ubar)+(ubar-u)

	double svar= sqr(m_CMSene);
	double shat= svar*m_x1*m_x2;

//    CALL BornV_MakeGami(m_CMSene,gamiCR,gami,alfi)           ! make gamiCR at CMSene
//    m_vv  = R**(1d0/gami)*m_vvmax
//    Rho   = Rho* m_vv/R/gami*m_vvmax

	m_Mll = sqrt(shat);

	if( m_Mll < 60.0 || m_Mll >160.0 ) Dist=0;

	long KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
	kksem_setkeyfob_( KeyFob );
	double xBorn;
	kksem_makeborn_( shat, xBorn);
	xBorn *= 1./3.;  // colour factor corrected by hand
	// xBorn *= 1000.;  // switching to picobarns
	Dist *= xBorn;

	return Dist;
}// Density
public:
};//

//______________________________________________________________________________
void ISRgener()
{
  cout<<"--- demo_small started ---"<<endl;
  DiskFileB.cd();

  double Mmin= 60;
  double Mmax=160;
  TH1D  *hst_Mll = new TH1D("hst_Mll" ,  "Mass distr.", 50,Mmin,Mmax);
  hst_Mll->Sumw2();

  //TFoamIntegrand *Rho1= new RhoISR();
  RhoISR *Rho1= new RhoISR();
  Rho1->m_Mmin = Mmin;
  Rho1->m_Mmax = Mmax;
  Rho1->m_CMSene = 8000;
// Setting up Foam object
  TRandom  *PseRan   = new TRandom3();  // Create random number generator
  PseRan->SetSeed(4357);
  TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
  FoamX->SetkDim(2);         // No. of dimensions, obligatory!
  FoamX->SetnCells( 2000);   // No. of cells, can be omitted, default=2000
  FoamX->SetnSampl( 2000);   // No. of MC evts/cell in exploration, default=200
  FoamX->SetRho(Rho1);       // Set 2-dim distribution, included above
  FoamX->SetPseRan(PseRan);  // Set random number generator, mandatory!
  FoamX->SetOptRej(0);       // wted events (=0), default wt=1 events (=1)
  FoamX->Initialize();       // Initialize simulator, may take time...
// loop over MC events
  double wt,Mll;
  long   NevTot = 1000000;
  for(long loop=0; loop<NevTot; loop++)
  {
    FoamX->MakeEvent();            // generate MC event
    FoamX->GetMCwt(wt);
    //Rho1->GetMll(Mll);
    Mll = Rho1->m_Mll;
    //cout<<"Mll =  "<< Mll <<endl;
    hst_Mll->Fill(Mll,wt);
  }// loop
//renormalizing histo
  double Xsav, dXsav;
  FoamX->GetIntNorm(Xsav,dXsav);
  HisNorm0( NevTot, Xsav, hst_Mll);
// plotting result
/*
  TCanvas *cMass = new TCanvas("cMass","Canvas for plotting",600,600);
  cMass->cd();
  gPad->SetLogy(); // !!!!!!
  hst_Mll->DrawCopy("h");  // final plot
  cMass->Update();
*/
  //
}//ISRgener

//_____________________________________________________________________________
class TFDISTR: public TFoamIntegrand{
public:
  TFDISTR(){};
  Double_t Density(int nDim, Double_t *Xarg)
  { // 2-dimensional distribution for Foam, normalized to one (within 1e-5)
    Double_t x=Xarg[0];
    Double_t y=Xarg[1];
    Double_t GamSq= sqr(0.100e0);
    Double_t Dist= 0;
    Dist +=exp(-(sqr(x-1./3) +sqr(y-1./3))/GamSq)/GamSq/TMath::Pi();
    Dist +=exp(-(sqr(x-2./3) +sqr(y-2./3))/GamSq)/GamSq/TMath::Pi();
    return 0.5*Dist;
  }//
};
//______________________________________________________________________________
void demo_small()
{
  cout<<"--- demo_small started ---"<<endl;
  TH2D     *hst_xy = new TH2D("hst_xy" ,  "x-y plot", 50,0,1.0, 50,0,1.0);
  Double_t *MCvect = new Double_t[2]; // 2-dim vector generated in the MC run
  TFoamIntegrand *Camel2= new TFDISTR();
  //
  TRandom  *PseRan   = new TRandom3();  // Create random number generator
  PseRan->SetSeed(4357);
  TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
  FoamX->SetkDim(2);         // No. of dimensions, obligatory!
  FoamX->SetnCells(500);     // No. of cells, can be omitted, default=2000
  FoamX->SetRho(Camel2);     // Set 2-dim distribution, included above
  FoamX->SetPseRan(PseRan);  // Set random number generator, mandatory!
  FoamX->Initialize();       // Initialize simulator, may take time...
  //
  // From now on FoamX is ready to generate events according to Camel2(x,y)
  // now hst_xy will be ploted, visualising generated distribution
  TCanvas *cKanwa = new TCanvas("cKanwa","Canvas for plotting",600,600);
  cKanwa->cd();
  int nshow=5000;
  for(long loop=0; loop<100000; loop++)
  {
    FoamX->MakeEvent();            // generate MC event
    FoamX->GetMCvect( MCvect);     // get generated vector (x,y)
    Double_t x=MCvect[0];
    Double_t y=MCvect[1];
    if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
    hst_xy->Fill(x,y);
  }// loop
  //
  hst_xy->DrawCopy("lego2");  // final plot
  cKanwa->Update();
  //
  Double_t MCresult, MCerror;
  FoamX->GetIntegMC( MCresult, MCerror);  // get MC integral, should be one
  cout << " MCresult= " << MCresult << " +- " << MCerror <<endl;
  cout<<"--- demo_small ended ---"<<endl;
}//demo_small



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
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_M100mu") );
}

///////////////////////////////////////////////////////////////////////////////////
void KKsemMakeHisto(){
// Here we produce semianalytical plots using KKsem program, No plotting
//------------------------------------------------------------------------  
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"xxxxxxxxxxxxxxxx KKsemMakeHisto  BEGIN xxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  // initilalization of KKsem
  KKsem LibSem;
  LibSem.Initialize(DiskFileA);
  //
  long KF=13; // muon ?????
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

//------------------------------------------------------------------------
//   MuMu  dsigma/dv
//------------------------------------------------------------------------  
//  TH1D *hstVtemplate = (TH1D*)DiskFileA.Get("hst_vTrueMain");
//  TH1D *hstCtemplate = (TH1D*)DiskFileA.Get("hst_Cost1Ceex2");
  // ISR*FSR
//  KeyDis = 302302;        // ISR*FSR O(alf2)
//  sprintf(chak,"XRHO2");  // ISR*FSR Mff
//  TH1D *vdis_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vdis_ISR2_FSR2");
//  LibSem.VVplot(vdis_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
  // ISR only
//  KeyDis = 303;           // ISR O(alf3)
//  sprintf(chak,"VRHO2");  // ISR only
//  TH1D *vdis_ISR2 =(TH1D*)hstVtemplate->Clone("vdis_ISR2");
//  LibSem.VVplot(vdis_ISR2, KF, chak, KeyDis, KeyFob);


  cout<<"xxxxxxxxxxxxxxxx KKsemMakeHisto END xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
//------------------------------------------------------------------------  
//------------------------------------------------------------------------  
}//  KKsemMakeHisto



///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  TH1D *hst_nPhAll     = (TH1D*)DiskFileA.Get("hst_nPhAll");
  TH1D *hst_nPhVis     = (TH1D*)DiskFileA.Get("hst_nPhVis");
  TH1D *hst_weight     = (TH1D*)DiskFileA.Get("hst_weight");

  TH1D *hst_vTrueCeex2 = (TH1D*)DiskFileA.Get("hst_vTrueCeex2");
  TH1D *hst_vAlepCeex2 = (TH1D*)DiskFileA.Get("hst_vAlepCeex2");
  TH1D *hst_vXGenCeex2 = (TH1D*)DiskFileA.Get("hst_vXGenCeex2");

  TH1D *hst_s1Ceex2  = (TH1D*)DiskFileA.Get("hst_s1Ceex2");
  TH1D *hst_svk      = (TH1D*)DiskFileA.Get("hst_svk");

  TH1D *hst_M100mu      = (TH1D*)DiskFileA.Get("hst_M100mu");

//------------------------------------------------------------------------  
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo: general info ", 50, 80,    1000,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  2);
  //cFigInfo->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //==========plot1==============
  cFigInfo->cd(1);
  hst_nPhVis->DrawCopy("h");
  hst_nPhAll->SetLineColor(2);
  hst_nPhAll->DrawCopy("hsame");
  //==========plot2==============
  cFigInfo->cd(2);
  hst_weight->DrawCopy("h");
  //==========plot3==============
  cFigInfo->cd(3);
  gPad->SetLogy(); // !!!!!!
  hst_vTrueCeex2->SetStats(0);
  hst_vTrueCeex2->SetTitle(0);
  hst_vTrueCeex2->DrawCopy("h");
  //
  hst_vAlepCeex2->SetLineColor(2);
  hst_vAlepCeex2->DrawCopy("hsame");
  //
  hst_vXGenCeex2->SetLineColor(4);
  hst_vXGenCeex2->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/dv (Ceex2); Black=Bare, Red=Aleph, Blue=Gener");
  //==========plot4==============
  cFigInfo->cd(4);

  //hst_svk->SetLineColor(4);
  //hst_svk->DrawCopy("h");
  hst_M100mu->SetLineColor(4);
  hst_M100mu->DrawCopy("h");

  //hst_s1Ceex2->DrawCopy("hsame");
//----------------------------
  cFigInfo->cd();
}// FigInfo

///////////////////////////////////////////////////////////////////////////////////
void FigVtest()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVtest =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  TH1D *hst_vTrueCeex2 = (TH1D*)DiskFileA.Get("hst_vTrueCeex2");
  TH1D *hst_vXGenCeex2 = (TH1D*)DiskFileA.Get("hst_vXGenCeex2");
  //
  TH1D *vdis_ISR2      = (TH1D*)DiskFileB.Get("vdis_ISR2");
  TH1D *vdis_ISR2_FSR2 = (TH1D*)DiskFileB.Get("vdis_ISR2_FSR2");
  //
  //****************************************************************************************
  //****************************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVtest = new TCanvas("cFigVtest","FigVtest: photonic2", 50, 50,    1000, 800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigVtest->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVtest->Divide( 2,  2);
  //cFigVtest->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.03);
  //==========plot1==============
  cFigVtest->cd(1);
  gPad->SetLogy(); // !!!!!!
  hst_vTrueCeex2->SetStats(0);
  hst_vTrueCeex2->SetTitle(0);
  hst_vTrueCeex2->DrawCopy("h");  // black
  //
  //vdis_ISR2->SetLineColor(kBlue); // blue
  //vdis_ISR2->DrawCopy("hsame");
  //
  vdis_ISR2_FSR2->SetLineColor(kMagenta); // magenta
  vdis_ISR2_FSR2->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv (Ceex2); Black=Bare, Red=Gener, Blue=ISR, Mag=ISR+FSR");
  //==========plot2==============
  cFigVtest->cd(2);
  hst_vTrueCeex2->Divide(vdis_ISR2_FSR2);
  hst_vTrueCeex2->SetStats(0);
  hst_vTrueCeex2->SetTitle(0);
  hst_vTrueCeex2->SetMinimum(0.85);
  hst_vTrueCeex2->SetMaximum(1.15);
  hst_vTrueCeex2->DrawCopy("h");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv (Ceex2); red: Gener/KKsemISR+FSR");
  //==========plot3==============
  cFigVtest->cd(3);
  gPad->SetLogy(); // !!!!!!
  hst_vXGenCeex2->SetStats(0);
  hst_vXGenCeex2->SetTitle(0);
  hst_vXGenCeex2->SetLineColor(kRed); // red
  hst_vXGenCeex2->DrawCopy("h");
  //
  vdis_ISR2->SetLineColor(kBlue); // blue
  vdis_ISR2->DrawCopy("hsame");
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv (Ceex2); Red=Gener, Blue=ISR");
  //==========plot4==============
  cFigVtest->cd(4);
  hst_vXGenCeex2->Divide(vdis_ISR2);
  hst_vXGenCeex2->SetStats(0);
  hst_vXGenCeex2->SetTitle(0);
  hst_vXGenCeex2->SetMinimum(0.85);
  hst_vXGenCeex2->SetMaximum(1.15);
  hst_vXGenCeex2->DrawCopy("h");  // black
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv (Ceex2); red: Gener/KKsemISR");
  //----------------------------
  cFigVtest->cd();
  //================================================
}//FigVtest



///////////////////////////////////////////////////////////////////////////////////
void FigMass()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigMass =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  char captEne[100];
  sprintf(captEne,"#sqrt{s} =%4.0fGeV, u-ubar", CMSene);

  // From FOAM
  TH1D *hst_Mll      = (TH1D*)DiskFileB.Get("hst_Mll");
  // From KKMC
  TH1D *hst_M100mu   = (TH1D*)DiskFileA.Get("hst_M100mu");

  // integrate over bins
  TH1D *Hst1 =hst_M100mu;
  TH1D *Hst2 =hst_Mll;
  int      nbX  = Hst1->GetNbinsX();
  Double_t Xmax = Hst1->GetXaxis()->GetXmax();
  Double_t Xmin = Hst1->GetXaxis()->GetXmin();
  double dx = (Xmax-Xmin)/nbX;
  double xsum1 = 0;
  for(int ix=1; ix <= nbX; ix++){
	 xsum1  += Hst1->GetBinContent(  ix ) *dx;
//	 cout<< "ix="<< ix <<"  xsum1="<< xsum<<endl;
  }//ix
  nbX  = Hst2->GetNbinsX();
  Xmax = Hst2->GetXaxis()->GetXmax();
  Xmin = Hst2->GetXaxis()->GetXmin();
  dx = (Xmax-Xmin)/nbX;
  double xsum2 = 0;
  for(int ix=1; ix <= nbX; ix++){
	 xsum2  += Hst2->GetBinContent(  ix ) *dx;
  }//ix

  char capt1[100];
  sprintf(capt1,"#sigma_{KKMC} =%9.6f [pb]", 1000*xsum1);
  char capt2[100];
   sprintf(capt2,"#sigma_{FOAM} =%9.6f [pb]", 1000*xsum2);
   //

//------------------------------------------------------------------------
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigMass = new TCanvas("cFigMass","FigMass: general info ", 50, 80,    1000,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigMass->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigMass->Divide( 2,  2);
  //cFigMass->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //==========plot1==============
  cFigMass->cd(1);

  gPad->SetLogy(); // !!!!!!

  //Hst1->SetStats(0);
  Hst1->SetTitle(0);

  Hst1->SetLineColor(kBlue);
  Hst1->DrawCopy("h");

  hst_Mll->SetLineColor(kRed);
  hst_Mll->DrawCopy("hsame");

  CaptT->DrawLatex(0.10,0.95,"d#sigma/dM [nb/GeV] KKMC and FOAM (red)");
  CaptT->DrawLatex(0.40,0.85, captEne);
  CaptT->DrawLatex(0.40,0.75, capt1);

  //==========plot2==============
  cFigMass->cd(2);
  gPad->SetLogy(); // !!!!!!

  Hst2->SetTitle(0);
  Hst2->DrawCopy("h");

  CaptT->DrawLatex(0.10,0.95,"d#sigma/dM [nb/GeV] FOAM");
  CaptT->DrawLatex(0.40,0.75, capt2);

  //==========plot3==============
  cFigMass->cd(3);
  gPad->SetLogy(); // !!!!!!
  //==========plot4==============
  cFigMass->cd(4);

//----------------------------
  cFigMass->cd();
}// FigMass




///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  KKsemMakeHisto();        //
  //========== PLOTTING ==========
  //FigInfo();
  //FigVtest(); //***

  ISRgener();
  FigMass();

  //demo_small();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
