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
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueMain") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotMain") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueCeex2") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_mPhotCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_mPhotCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_mPhotCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotNuel") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotNumu") );
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
  TH1D *hst_vPhotCeex2 = (TH1D*)DiskFileA.Get("hst_vPhotCeex2");
  //
  TH1D *hst_vPhotMain  = (TH1D*)DiskFileA.Get("hst_vPhotMain");

  TH1D *hst_CosPLCeex2 = (TH1D*)DiskFileA.Get("hst_CosPLCeex2");
//------------------------------------------------------------------------  
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo: general info ", 50, 50,    1000,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  2);
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
  hst_LnThPhAll->SetLineColor(kBlue);
  hst_LnThPhAll->DrawCopy("h");
  hst_LnThPhVis->SetLineColor(kRed);
  hst_LnThPhVis->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/d ln_{10}(sin(#theta)),   all (blue) and visible (red) photons");
  //----------------------------
  cFigInfo->cd();
}//FigInfo

///////////////////////////////////////////////////////////////////////////////////
void FigCEEX21()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCEEX21 =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  ///
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  ///
  TH1D *hst_vPhotCeex1     = (TH1D*)DiskFileA.Get("hst_vPhotCeex1");
  TH1D *hst_vPhotCeex2     = (TH1D*)DiskFileA.Get("hst_vPhotCeex2");
  TH1D *hst_vPhotCeex12    = (TH1D*)DiskFileA.Get("hst_vPhotCeex12");
  ///
  TH1D *Hst1  = hst_vPhotCeex1;
  TH1D *Hst2  = hst_vPhotCeex2;
//!////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  TH1D *H_Vline0  = new TH1D("H_Vline0","one",  1, 0.0, 1.0);
  H_Vline0->SetBinContent(1,0);
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  int ndiv=2;
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  400*ndiv;
  TCanvas *cCeex21 = new TCanvas("cCeex21","cCeex21",100,100, WidPix,HeiPix);
  cCeex21->SetFillColor(10);
  cCeex21->Draw();
  cCeex21->SetFillColor(10);
  cCeex21->Divide(0, ndiv);
 //!-----------------------------
  cCeex21->cd(1);
  //gPad->SetLogy(); // !!!!!!
  Hst1->SetStats(0);
  Hst1->SetTitle(0);
  Hst1->SetLineColor(kBlue);
  Hst1->DrawCopy("h");
  ///
  Hst2->SetLineColor(kRed);
  Hst2->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/dv, Photon energy. CEEX1 (blue) and CEEX2 (red), neutrinos");
 //!-----------------------------
  cCeex21->cd(2);
  H_Vline0->SetStats(0);
  //H_Vline0->SetTitle(0);
  H_Vline0->GetYaxis()->SetTitleSize(0.06);
  H_Vline0->GetYaxis()->SetTitleOffset(1.2);
  H_Vline0->GetYaxis()->CenterTitle();
  H_Vline0->GetXaxis()->SetTitleOffset(0.6);
  H_Vline0->GetXaxis()->SetTitleSize(0.07);
  H_Vline0->GetXaxis()->SetTitle("v=Ephot/Ebeam");
  H_Vline0->SetMaximum( +0.1);
  H_Vline0->SetMinimum( -0.1);
  H_Vline0->GetXaxis()->SetTitleOffset(0.6);
  H_Vline0->GetXaxis()->SetTitleSize(0.07);
  H_Vline0->SetTitle("(CEEX1-CEEX2)/CEEX2");
  H_Vline0->DrawCopy("h");
  ///(TH1D*)hst_vPhotCeex1->Clone("RAT1");
  TH1D *RAT_hvCeex12 = (TH1D*)hst_vPhotCeex12->Clone("RAT_hvCeex12") ;  ///  CEEX1-CEEX2
  RAT_hvCeex12->Divide(Hst2);             /// divide over CEEX2
  RAT_hvCeex12->SetLineColor(kBlue);
  RAT_hvCeex12->DrawCopy("hsame");
  //!-----------------------------
  cCeex21->Update();
  cCeex21->cd();
  //!-----------------------------
}//FigCEEX21


//!/////////////////////////////////////////////////////////////////////////////////
void FigCEEX21mu()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCEEX21mu =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  ///
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  ///
  TH1D *hst_mPhotCeex1     = (TH1D*)DiskFileA.Get("hst_mPhotCeex1");
  TH1D *hst_mPhotCeex2     = (TH1D*)DiskFileA.Get("hst_mPhotCeex2");
  TH1D *hst_mPhotCeex12    = (TH1D*)DiskFileA.Get("hst_mPhotCeex12");
  ///
  TH1D *Hst1  = hst_mPhotCeex1;
  TH1D *Hst2  = hst_mPhotCeex2;
//!////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  TH1D *H_Vline10  = new TH1D("H_Vline10","one",  1, 0.0, 1.0);
  H_Vline10->SetBinContent(1,0);
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  int ndiv=2;
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  400*ndiv;
  TCanvas *cCeex21mu = new TCanvas("cCeex21mu","cCeex21mu",150,150, WidPix,HeiPix);
  cCeex21mu->SetFillColor(10);
  cCeex21mu->Draw();
  cCeex21mu->SetFillColor(10);
  cCeex21mu->Divide(0, ndiv);
 //!-----------------------------
  cCeex21mu->cd(1);
  //gPad->SetLogy(); // !!!!!!
  Hst1->SetStats(0);
  Hst1->SetTitle(0);
  Hst1->SetLineColor(kBlue);
  Hst1->DrawCopy("h");
  ///
  Hst2->SetLineColor(kRed);
  Hst2->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/dv, Photon energy. CEEX1 (blue) and CEEX2 (red), muon-pairs");
 //!-----------------------------
  cCeex21mu->cd(2);
  H_Vline10->SetStats(0);
  //H_Vline10->SetTitle(0);
  H_Vline10->GetYaxis()->SetTitleSize(0.06);
  H_Vline10->GetYaxis()->SetTitleOffset(1.2);
  H_Vline10->GetYaxis()->CenterTitle();
  H_Vline10->GetXaxis()->SetTitleOffset(0.6);
  H_Vline10->GetXaxis()->SetTitleSize(0.07);
  H_Vline10->GetXaxis()->SetTitle("v=Ephot/Ebeam");
  H_Vline10->SetMaximum( +0.1);
  H_Vline10->SetMinimum( -0.1);
  H_Vline10->GetXaxis()->SetTitleOffset(0.6);
  H_Vline10->GetXaxis()->SetTitleSize(0.07);
  H_Vline10->SetTitle("(CEEX1-CEEX2)/CEEX2");
  H_Vline10->DrawCopy("h");
  ///
  TH1D *RAT_hvCeex12 = (TH1D*)hst_mPhotCeex12->Clone("RAT_hvCeex12") ;  ///  CEEX1-CEEX2
  RAT_hvCeex12->Divide(Hst2);             /// divide over CEEX2
  RAT_hvCeex12->SetLineColor(kBlue);
  RAT_hvCeex12->DrawCopy("hsame");
  //!-----------------------------
  cCeex21mu->Update();
  cCeex21mu->cd();
  //!-----------------------------
}///FigCEEX21mu



///////////////////////////////////////////////////////////////////////////////////
void FigNuDiff()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNuDiff =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  ///
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  ///
  TH1D *hst_vPhotNuel     = (TH1D*)DiskFileA.Get("hst_vPhotNuel");
  TH1D *hst_vPhotNumu     = (TH1D*)DiskFileA.Get("hst_vPhotNumu");
  ///
  TH1D *Hst1  = hst_vPhotNuel;
  TH1D *Hst2  = hst_vPhotNumu;
//!////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  TH1D *H_Vline1  = new TH1D("H_Vline1","one",  1, 0.0, 1.0);
  H_Vline1->SetBinContent(1,1);
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  int ndiv=2;
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  400*ndiv;
  TCanvas *cNuDiff = new TCanvas("cNuDiff","cNuDiff",200,200, WidPix,HeiPix);
  cNuDiff->SetFillColor(10);
  cNuDiff->Draw();
  cNuDiff->SetFillColor(10);
  cNuDiff->Divide(0, ndiv);
 //!-----------------------------
  cNuDiff->cd(1);
  Hst1->SetStats(0);
  Hst1->SetTitle(0);
  Hst1->SetLineColor(kBlue);
  Hst1->DrawCopy("h");
  ///
  Hst2->SetLineColor(kRed);
  Hst2->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/dv, Photon energy. NuElectron (blue) and NuMuon (red)");

  cNuDiff->cd(2);
  H_Vline1->SetStats(0);
  //H_Vline1->SetTitle(0);
  H_Vline1->GetYaxis()->SetTitleSize(0.06);
  H_Vline1->GetYaxis()->SetTitleOffset(1.2);
  H_Vline1->GetYaxis()->CenterTitle();
  H_Vline1->GetXaxis()->SetTitleOffset(0.6);
  H_Vline1->GetXaxis()->SetTitleSize(0.07);
  H_Vline1->GetXaxis()->SetTitle("v=Ephot/Ebeam");
  H_Vline1->SetMaximum( 3);
  H_Vline1->SetMinimum( 0);
  H_Vline1->GetXaxis()->SetTitleOffset(0.6);
  H_Vline1->GetXaxis()->SetTitleSize(0.07);
  H_Vline1->SetTitle("NuElectron/NuMuon");
  H_Vline1->DrawCopy("h");
  ///
  TH1D *RAT_NuelNumu = (TH1D*)Hst1->Clone("RAT_NuelNumu");
  RAT_NuelNumu->Divide(Hst2);             /// divide over Numu
  RAT_NuelNumu->SetLineColor(kBlue);
  RAT_NuelNumu->DrawCopy("hsame");
 //!-----------------------------
  cNuDiff->Update();
  cNuDiff->cd();
}//FigNuDiff

//!/////////////////////////////////////////////////////////////////////////////////
void FigCEEX21rat()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCEEX21rat =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  ///
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  ///
  TH1D *hst_vPhotCeex1     = (TH1D*)DiskFileA.Get("hst_vPhotCeex1");  /// 3 neutrinos
  TH1D *hst_vPhotCeex2     = (TH1D*)DiskFileA.Get("hst_vPhotCeex2");  /// 3 neutrinos
  TH1D *hst_vPhotCeex12    = (TH1D*)DiskFileA.Get("hst_vPhotCeex12");
  TH1D *hst_mPhotCeex1     = (TH1D*)DiskFileA.Get("hst_mPhotCeex1");  /// muon
  TH1D *hst_mPhotCeex2     = (TH1D*)DiskFileA.Get("hst_mPhotCeex2");  /// muon
  TH1D *hst_mPhotCeex12    = (TH1D*)DiskFileA.Get("hst_mPhotCeex12");
  ///
//!////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  TH1D *H_Vline12  = new TH1D("H_Vline12","one",  1, 0.0, 1.0);
  //H_Vline12->SetBinContent(1,1);
  H_Vline12->SetBinContent(1,0);
  /// nu/mu ratios, CEEX1 and CEEX2
  TH1D *RAT1 = (TH1D*)hst_vPhotCeex1->Clone("RAT1");
  RAT1->Divide(hst_mPhotCeex1);
  TH1D *RAT2 = (TH1D*)hst_vPhotCeex2->Clone("RAT2");
  RAT2->Divide(hst_mPhotCeex2);
  /// (CEEX1-CEEX2(/CEEX2 ratios for nu and mu
  TH1D *DELnu = (TH1D*)hst_vPhotCeex12->Clone("DELnu");
  DELnu->Divide(hst_vPhotCeex2);
  TH1D *DELmu = (TH1D*)hst_mPhotCeex12->Clone("DELmu");
  DELmu->Divide(hst_mPhotCeex2);
  DELmu->Scale(-1e0);
  DELnu->Add(DELmu);
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  int ndiv=2;
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  400*ndiv;
  TCanvas *cCeex21rat = new TCanvas("cCeex21rat","cCeex21rat",250,250, WidPix,HeiPix);
  cCeex21rat->SetFillColor(10);
  cCeex21rat->Draw();
  cCeex21rat->SetFillColor(10);
  cCeex21rat->Divide(0, ndiv);
 //!-----------------------------
  cCeex21rat->cd(1);
  //gPad->SetLogy(); // !!!!!!
  
  RAT1->SetStats(0);
  RAT1->SetTitle(0);
  RAT1->SetLineColor(kBlue);
  RAT1->DrawCopy("h");
  RAT2->SetLineColor(kRed);
  RAT2->DrawCopy("hsame");
  ///
  CaptT->DrawLatex(0.10,0.95,"R_{i}=#sigma_{#nu#nu#gamma}/#sigma_{#mu#mu#gamma}, R_{2} for CEEX2 (red), R_{1} for CEEX1 (blue),");
 //!-----------------------------
  cCeex21rat->cd(2);
  H_Vline12->SetStats(0);
  //H_Vline12->SetTitle(0);
  H_Vline12->GetYaxis()->SetTitleSize(0.06);
  H_Vline12->GetYaxis()->SetTitleOffset(1.2);
  H_Vline12->GetYaxis()->CenterTitle();
  H_Vline12->GetXaxis()->SetTitleOffset(0.6);
  H_Vline12->GetXaxis()->SetTitleSize(0.07);
  H_Vline12->GetXaxis()->SetTitle("v=Ephot/Ebeam");
  //H_Vline12->SetMaximum( 1+0.02);
  //H_Vline12->SetMinimum( 1-0.02);
  H_Vline12->SetMaximum(  +0.004);
  H_Vline12->SetMinimum(  -0.004);
  H_Vline12->GetXaxis()->SetTitleOffset(0.6);
  H_Vline12->GetXaxis()->SetTitleSize(0.07);
  H_Vline12->SetTitle(" R_{1}/R_{2}-1, (CEEX1-CEEX2)/CEEX2");
  H_Vline12->DrawCopy("h");
  ///
  TH1D *RAT_21 = (TH1D*)RAT2->Clone("RAT2");
  RAT_21->Divide(RAT1);             /// divide over CEEX2
  RAT_21->SetLineColor(kBlue);
  //RAT_21->DrawCopy("hsame");
  
  DELnu->SetLineColor(kBlue);
  DELnu->DrawCopy("hsame");
  
  
  //!-----------------------------
  cCeex21rat->Update();
  cCeex21rat->cd();
  //!-----------------------------
}///FigCEEX21mu




///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  //========== PLOTTING ==========
  FigInfo();
  FigCEEX21();
  FigCEEX21mu();
  FigNuDiff();
  FigCEEX21rat();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
