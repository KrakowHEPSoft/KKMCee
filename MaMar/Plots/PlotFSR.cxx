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

Double_t sqr( const Double_t x ){ return x*x;};
// Auxiliary procedures for plotting
#include "HisNorm.h"
#include "Marker.h"

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
  char capt1[100];
  sprintf(capt1,"#sqrt{s} =%4.0fGeV, u-ubar", CMSene);

  //
  TH1D *hst_M100mu      = (TH1D*)DiskFileA.Get("hst_M100mu");
  // integrate over bins
  TH1D *Hst1 =hst_M100mu;
  int      nbX  = Hst1->GetNbinsX();
  Double_t Xmax = Hst1->GetXaxis()->GetXmax();
  Double_t Xmin = Hst1->GetXaxis()->GetXmin();
  double dx = (Xmax-Xmin)/nbX;
  double xsum = 0;
  for(int ix=1; ix <= nbX; ix++){
	 xsum  += Hst1->GetBinContent(  ix ) *dx;
//	 cout<< "ix="<< ix <<"  xsum="<< xsum<<endl;
  }//ix
  xsum *= 1./3.; // colour factor by hand
  char capt2[100];
  sprintf(capt2,"#sigma =%9.6f [pb]", 1000*xsum);

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

  CaptT->DrawLatex(0.10,0.95,"d#sigma/dM [nb/GeV]");
  CaptT->DrawLatex(0.40,0.85, capt1);
  CaptT->DrawLatex(0.40,0.75, capt2);

  //==========plot2==============
  cFigMass->cd(2);
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
  FigMass();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
