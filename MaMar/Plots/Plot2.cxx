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


// ROOT headers
#include "TROOT.h"
#include "TFile.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
TFile DiskFileA("../test2/rmain.root");
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
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_LnThPhAll") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_LnThPhVis") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotMain") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueMain") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vXGenCeex2") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPLCeex2") );
  //
}


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

  TH1D *hst_LnThPhAll  = (TH1D*)DiskFileA.Get("hst_LnThPhAll");
  TH1D *hst_LnThPhVis  = (TH1D*)DiskFileA.Get("hst_LnThPhVis");

  TH1D *hst_vTrueMain  = (TH1D*)DiskFileA.Get("hst_vTrueMain");
  TH1D *hst_vXGenCeex2 = (TH1D*)DiskFileA.Get("hst_vXGenCeex2");
  //
  TH1D *hst_vPhotMain  = (TH1D*)DiskFileA.Get("hst_vPhotMain");

  TH1D *hst_CosPLCeex2 = (TH1D*)DiskFileA.Get("hst_CosPLCeex2");
//------------------------------------------------------------------------  
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo: general info ", 50, 80,    1000,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  2);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //==========plot1==============
  cFigInfo->cd(1);
  gPad->SetLogy(); // !!!!!!
  hst_nPhAll->SetMinimum( 1e-3*hst_nPhAll->GetMaximum());
  hst_nPhAll->SetTitle(0);
  hst_nPhAll->SetLineColor(kBlue);
  hst_nPhAll->DrawCopy("h");
  hst_nPhVis->SetLineColor(kRed);
  hst_nPhVis->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"Photon multiplicity. All (blue) and visible (red). 161GeV");
  //==========plot2==============
  cFigInfo->cd(2);
  hst_weight->DrawCopy("h");
  //==========plot3==============
  cFigInfo->cd(3);
  gPad->SetLogy(); // !!!!!!
  hst_vTrueMain->SetStats(0);
  hst_vTrueMain->SetTitle(0);
  hst_vTrueMain->SetLineColor(kBlue);
  hst_vTrueMain->DrawCopy("h");
  ///
  hst_vPhotMain->SetLineColor(kRed);
  hst_vPhotMain->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/dv, Photon energy. All (blue) and visible (red)");
  //==========plot4==============
  cFigInfo->cd(4);
  //-----------------------------
  //hst_CosPLCeex2->SetTitle(0);
  //hst_CosPLCeex2->SetStats(0);
  //hst_CosPLCeex2->SetLineColor(4);
  //hst_CosPLCeex2->DrawCopy("h");  
  //CaptT->DrawLatex(0.10,0.95,"d#sigma/dc (Ceex2); Black=   #theta_{1}, Red=PRD, Blue=PL");
//  gPad->SetLogy(); // !!!!!!
//  hst_LnThPhAll->SetMinimum( 1e-4*hst_LnThPhAll->GetMaximum());
  hst_LnThPhAll->SetTitle(0);
  hst_LnThPhAll->SetStats(0);
  hst_LnThPhAll->SetMinimum(0);
  hst_LnThPhAll->DrawCopy("h");
  hst_LnThPhVis->SetLineColor(kRed);
  hst_LnThPhVis->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/d ln_{10}(sin(#theta)),   all and visible photons");
  //----------------------------
  cFigInfo->cd();
}


///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  //========== PLOTTING ==========
  FigInfo();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
