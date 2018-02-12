//////////////////////////////////////////////////////////////////////
//    ( make Plot2; ./Plot2 )
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
#include "TFile.h"

#include "HisNorm.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
TFile DiskFileA("../workZinv/rmain.root");
//=============================================================================

//Double_t sqr( const Double_t x ){ return x*x;};
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot; // from KKMC run
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
//
float  gXcanv = 0, gYcanv = 0;
///////////////////////////////////////////////////////////////////////////////////



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
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueMu") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_mPhotCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_mPhotCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_mPhotCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_nPhotCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_nPhotCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_nPhotCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_lPhotCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_lPhotCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_lPhotCeex12") );
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
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  ///
  Double_t CMSene = HST_KKMC_NORMA->GetBinContent(1)/HST_KKMC_NORMA->GetBinContent(511);
  char capt1[100];
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  //
  TH1D *hst_nPhAll     = (TH1D*)DiskFileA.Get("hst_nPhAll");
  TH1D *hst_nPhVis     = (TH1D*)DiskFileA.Get("hst_nPhVis");
  TH1D *hst_weight     = (TH1D*)DiskFileA.Get("hst_weight");

  TH1D *hst_LnThPhAll  = (TH1D*)DiskFileA.Get("hst_LnThPhAll");
  TH1D *hst_LnThPhVis  = (TH1D*)DiskFileA.Get("hst_LnThPhVis");

  TH1D *hst_vPhotCeex2 = (TH1D*)DiskFileA.Get("hst_vPhotCeex2");
  TH1D *hst_vTrueMain  = (TH1D*)DiskFileA.Get("hst_vTrueMain");
  TH1D *hst_vTrueMu    = (TH1D*)DiskFileA.Get("hst_vTrueMu");
  //
  TH1D *hst_vPhotMain  = (TH1D*)DiskFileA.Get("hst_vPhotMain");

  TH1D *hst_CosPLCeex2 = (TH1D*)DiskFileA.Get("hst_CosPLCeex2");
//------------------------------------------------------------------------  
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.060);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo: general info ", 50, 50,    1000,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  2);
  //==========plot1==============
  cFigInfo->cd(1);
  gPad->SetLogy(); // !!!!!!
  hst_nPhAll->SetStats(0);
  hst_nPhAll->SetMinimum( 1e-3*hst_nPhAll->GetMaximum());
  hst_nPhAll->SetTitle(0);
  hst_nPhAll->SetLineColor(kBlue);
  hst_nPhAll->GetXaxis()->SetLabelSize(0.06);
  hst_nPhAll->DrawCopy("h");
  //
  hst_nPhVis->SetLineColor(kRed);
  hst_nPhVis->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.93,"No. of #gamma's. tagged (red) and All (blue)");
  CaptE->DrawLatex(0.70,0.85, capt1);
  ///==========plot2==============
  cFigInfo->cd(2);
  //-----------------------------
  //hst_CosPLCeex2->SetTitle(0);
  //hst_CosPLCeex2->SetStats(0);
  //hst_CosPLCeex2->SetLineColor(4);
  //hst_CosPLCeex2->DrawCopy("h");  
  //CaptT->DrawLatex(0.10,0.95,"d#sigma/dc (Ceex2); Black=   #theta_{1}, Red=PRD, Blue=PL");
///
//  gPad->SetLogy(); // !!!!!!
  hst_LnThPhAll->SetTitle(0);
  hst_LnThPhAll->SetStats(0);
  hst_LnThPhAll->SetMinimum(0);
  hst_LnThPhAll->SetLineColor(kBlue);
  hst_LnThPhAll->GetXaxis()->SetLabelSize(0.06);
  hst_LnThPhAll->DrawCopy("h");
  //
  hst_LnThPhVis->SetLineColor(kRed);
  hst_LnThPhVis->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.93,"d#sigma/d(ln_{10}sin #theta),  tagged #gamma's and All");
  CaptE->DrawLatex(0.70,0.85, capt1);
  ///==========plot3==============
  cFigInfo->cd(3);
  gPad->SetLogy(); // !!!!!!
  hst_vTrueMain->SetStats(0);
  hst_vTrueMain->SetTitle(0);
  hst_vTrueMain->GetXaxis()->SetLabelSize(0.05);
  hst_vTrueMain->SetLineColor(kBlue);
  hst_vTrueMain->DrawCopy("h");
  ///
  hst_vPhotMain->SetLineColor(kRed);
  hst_vPhotMain->DrawCopy("hsame");
  ///
  CaptT->DrawLatex(0.10,0.93,
    "d#sigma/dv,  v=1-M^{2}_{#nu#bar{#nu}}/s,  tagged #gamma's and All");
  CaptE->DrawLatex(0.70,0.85, capt1);
 ///==========plot4==============
  cFigInfo->cd(4);
  ///
  hst_weight->GetXaxis()->SetLabelSize(0.05);
  hst_weight->DrawCopy("h");
  ///
  gPad->SetLogy(); // !!!!!!
  hst_vTrueMain->DrawCopy("h");
  ///
  hst_vTrueMu->Scale(3);
  hst_vTrueMu->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.93, 
    "d#sigma/dv, #nu#bar{#nu} pairs and 3x(#mu^{-}#mu^{+}), untagged #gamma's");
  CaptE->DrawLatex(0.70,0.85, capt1);
  //
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
  CMSene  = HST_KKMC_NORMA->GetBinContent(1)/HST_KKMC_NORMA->GetBinContent(511);
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
  double vmin = hst_vPhotCeex12->GetXaxis()->GetXmin();
  double vmax = hst_vPhotCeex12->GetXaxis()->GetXmax();
  TH1D *H_Vline0  = new TH1D("H_Vline0","one",  1, vmin, vmax);
  H_Vline0->SetBinContent(1,0);
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  int ndiv=2;
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  400*ndiv;
  TCanvas *cCeex21 = new TCanvas("cCeex21","cCeex21",75,75, WidPix,HeiPix);
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
  H_Vline0->SetMaximum( +0.04);
  H_Vline0->SetMinimum( -0.04);
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
  CMSene  = HST_KKMC_NORMA->GetBinContent(1)/HST_KKMC_NORMA->GetBinContent(511);
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
  double vmin = hst_mPhotCeex12->GetXaxis()->GetXmin();
  double vmax = hst_mPhotCeex12->GetXaxis()->GetXmax();
  TH1D *H_Vline10  = new TH1D("H_Vline10","one",  1, vmin, vmax);
  H_Vline10->SetBinContent(1,0);
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  int ndiv=2;
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  400*ndiv;
  TCanvas *cCeex21mu = new TCanvas("cCeex21mu","cCeex21mu",100,100, WidPix,HeiPix);
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
  H_Vline10->SetMaximum( +0.04);
  H_Vline10->SetMinimum( -0.04);
  H_Vline10->GetXaxis()->SetTitleOffset(0.6);
  H_Vline10->GetXaxis()->SetTitleSize(0.08);
  H_Vline10->SetTitleSize(0.08);
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
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  ///
  Double_t CMSene  = HST_KKMC_NORMA->GetBinContent(1)/HST_KKMC_NORMA->GetBinContent(511);
  char capt1[100];
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  ///
  TH1D *hst_vPhotNuel     = (TH1D*)DiskFileA.Get("hst_vPhotNuel");
  TH1D *hst_vPhotNumu     = (TH1D*)DiskFileA.Get("hst_vPhotNumu");
  ///
  TH1D *Hel  = hst_vPhotNuel;
  TH1D *Hmu  = hst_vPhotNumu;
//!////////////////////////////////////////////
  double SigEl= Hel->Integral();
  double SigMu= Hmu->Integral();
  cout<< "@@@@@@@@@ SigEl="<<SigEl<<endl;
  cout<< "@@@@@@@@@ SigMu="<<SigMu<<endl;
  double DelTch=(SigEl-SigMu)/(3*SigMu); /// t_chanel contrib.
  cout<< "@@@@@@@@@ R="<< DelTch <<endl;
  char capt2[300];
  sprintf(capt2,"Integrated R_{t} = %2.4f ", DelTch);
//!////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
  CaptE->SetTextSize(0.070);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.06);
  double vmin = hst_vPhotNuel->GetXaxis()->GetXmin();
  double vmax = hst_vPhotNuel->GetXaxis()->GetXmax();
  TH1D *H_Vline1  = new TH1D("H_Vline1","one",  1, vmin, vmax);
  H_Vline1->SetBinContent(1,0);
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  int ndiv=2;
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  400*ndiv;
  TCanvas *cNuDiff = new TCanvas("cNuDiff","cNuDiff",125,125, WidPix,HeiPix);
  cNuDiff->SetFillColor(10);
  cNuDiff->Draw();
  cNuDiff->SetFillColor(10);
  cNuDiff->Divide(0, ndiv);
 //!-----------------------------
  cNuDiff->cd(1);
  Hel->SetStats(0);
  Hel->SetTitle(0);
  Hel->SetLineColor(kBlue);
  Hel->SetLabelSize(0.05);
  Hel->DrawCopy("h");
  ///
  Hmu->SetLineColor(kRed);
  Hmu->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,
   "d#sigma/dv,       v=E_{#gamma}/E_{beam}.      #nu_{el} (blue) and #nu_{#mu} (red)");
  CaptE->DrawLatex(0.80,0.85, capt1);

  cNuDiff->cd(2);
  H_Vline1->SetStats(0);
  H_Vline1->SetTitle(0);
  H_Vline1->GetYaxis()->SetTitleSize(0.06);
  H_Vline1->GetYaxis()->SetTitleOffset(1.2);
  H_Vline1->GetYaxis()->CenterTitle();
  H_Vline1->GetXaxis()->SetTitleOffset(0.6);
  H_Vline1->GetXaxis()->SetTitleSize(0.07);
  H_Vline1->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  H_Vline1->SetLabelSize(0.05);
  H_Vline1->SetMaximum( 1.15);
  H_Vline1->SetMinimum( 0.85);
  H_Vline1->SetMaximum( +0.10);
  H_Vline1->SetMinimum( -0.10);
  H_Vline1->GetXaxis()->SetTitleOffset(0.6);
  H_Vline1->GetXaxis()->SetTitleSize(0.07);
  H_Vline1->GetYaxis()->SetLabelSize(0.05);
  H_Vline1->DrawCopy("h");
  ///
  TH1D *RAT_NuelNumu = (TH1D*)Hmu->Clone("RAT_NuelNumu"); // numu
  RAT_NuelNumu->Scale(-1);    // -numu
  RAT_NuelNumu->Add(Hel);    //  nuel-numu
  RAT_NuelNumu->Divide(Hel); //  (nuel-numu)/numu
  RAT_NuelNumu->Scale(0.333333); // (nuel-numu)/(3*numu)
  RAT_NuelNumu->SetLineColor(kBlue);
  RAT_NuelNumu->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"t-channel W contrib.    R_{t}(v)=(#nu_{el}-#nu_{#mu})/(3 #nu_{#mu})");
  CaptE->DrawLatex(0.80,0.85, capt1);
  CaptE->DrawLatex(0.50,0.75, capt2);
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
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  ///
  Double_t CMSene = HST_KKMC_NORMA->GetBinContent(1)/HST_KKMC_NORMA->GetBinContent(511);
  char capt1[100];
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  ///
  TH1D *hst_vPhotCeex1     = (TH1D*)DiskFileA.Get("hst_vPhotCeex1");  /// 3 neutrinos
  TH1D *hst_vPhotCeex2     = (TH1D*)DiskFileA.Get("hst_vPhotCeex2");  /// 3 neutrinos
  TH1D *hst_vPhotCeex12    = (TH1D*)DiskFileA.Get("hst_vPhotCeex12");
  TH1D *hst_mPhotCeex1     = (TH1D*)DiskFileA.Get("hst_mPhotCeex1");  /// muon
  TH1D *hst_mPhotCeex2     = (TH1D*)DiskFileA.Get("hst_mPhotCeex2");  /// muon
  TH1D *hst_mPhotCeex12    = (TH1D*)DiskFileA.Get("hst_mPhotCeex12");
  ///
  TH1D *hst_nPhotCeex1     = (TH1D*)DiskFileA.Get("hst_nPhotCeex1");  /// 3 neutrinos
  TH1D *hst_nPhotCeex2     = (TH1D*)DiskFileA.Get("hst_nPhotCeex2");  /// 3 neutrinos
  TH1D *hst_nPhotCeex12    = (TH1D*)DiskFileA.Get("hst_nPhotCeex12");
  TH1D *hst_lPhotCeex1     = (TH1D*)DiskFileA.Get("hst_lPhotCeex1");  /// muon
  TH1D *hst_lPhotCeex2     = (TH1D*)DiskFileA.Get("hst_lPhotCeex2");  /// muon
  TH1D *hst_lPhotCeex12    = (TH1D*)DiskFileA.Get("hst_lPhotCeex12");
  ///
//!////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
  CaptE->SetTextSize(0.070);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.060);
//!////////////////////////////////////////////

  /// nu/mu ratios, CEEX1 and CEEX2
  TH1D *RAT1 = (TH1D*)hst_vPhotCeex1->Clone("RAT1");
  RAT1->Divide(hst_mPhotCeex1);
  TH1D *RAT2 = (TH1D*)hst_vPhotCeex2->Clone("RAT2");
  RAT2->Divide(hst_mPhotCeex2);

  /// (CEEX1-CEEX2)/CEEX2 ratios for nu and mu
  TH1D *DELnu = (TH1D*)hst_vPhotCeex12->Clone("DELnu");
  DELnu->Divide(hst_vPhotCeex2);
  TH1D *DELmu = (TH1D*)hst_mPhotCeex12->Clone("DELmu");
  DELmu->Divide(hst_mPhotCeex2);
  DELmu->Scale(-1e0);
  DELnu->Add(DELmu);

  /// another (CEEX1-CEEX2(/CEEX2 ratios for nu and mu
  TH1D *DELnu2 = (TH1D*)hst_nPhotCeex12->Clone("DELnu");
  DELnu2->Divide(hst_nPhotCeex2);
  TH1D *DELmu2 = (TH1D*)hst_lPhotCeex12->Clone("DELmu");
  DELmu2->Divide(hst_lPhotCeex2);
  DELmu2->Scale(-1e0);
  DELnu2->Add(DELmu2);

  TH1D *hZero = (TH1D*)hst_nPhotCeex12->Clone("hZero"); // zero line
  hZero->Reset();

//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  int ndiv=2;
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  400*ndiv;
  TCanvas *cCeex21rat = new TCanvas("cCeex21rat","cCeex21rat",150,150, WidPix,HeiPix);
  cCeex21rat->SetFillColor(10);
  cCeex21rat->Draw();
  cCeex21rat->SetFillColor(10);
  cCeex21rat->Divide(0, ndiv);
 //!-----------------------------
  cCeex21rat->cd(1);
  //gPad->SetLogy(); // !!!!!!
  
  RAT1->SetStats(0);
  RAT1->SetTitle(0);
  RAT1->GetYaxis()->SetLabelSize(0.05);
  RAT1->SetLabelSize(0.05);
  RAT1->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  RAT1->GetXaxis()->CenterTitle();
  RAT1->GetXaxis()->SetTitleSize(0.06);
  RAT1->GetXaxis()->SetTitleOffset(-0.6);
  RAT1->SetMaximum(7.0);
  RAT1->SetMinimum(5.5);
  RAT1->SetLineColor(kBlue);
  RAT1->DrawCopy("h");
  RAT2->SetLineColor(kRed);
  RAT2->DrawCopy("hsame");
  ///
  CaptT->DrawLatex(0.10,0.92,
   "     R=#sigma_{#nu#nu#gamma}/#sigma_{#mu#mu#gamma} for CEEX2 (red) and CEEX1 (blue),");
  CaptE->DrawLatex(0.80,0.85, capt1);
 //!-----------------------------
  cCeex21rat->cd(2);
  TH1D *HST2 = DELnu2;
  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetYaxis()->SetTitleSize(0.06);
  HST2->GetYaxis()->SetTitleOffset(1.2);
  HST2->GetYaxis()->CenterTitle();
  HST2->GetYaxis()->SetLabelSize(0.05);
  HST2->GetXaxis()->SetLabelSize(0.05);
  HST2->GetXaxis()->CenterTitle();
  HST2->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
  HST2->GetXaxis()->SetTitleSize(0.07);
  HST2->GetXaxis()->SetTitleOffset(-0.6); // does not work???
  //HST2->SetMaximum( 1+0.01); HST2->SetMinimum( 1-0.01);
  //HST2->SetMaximum(  +0.001);   HST2->SetMinimum(  -0.001);
  //HST2->SetMaximum(  +0.0002);  HST2->SetMinimum(  -0.0002);  // !!!
  HST2->GetXaxis()->SetTitleOffset(0.6);
  HST2->GetXaxis()->SetTitleSize(0.07);

  HST2->SetLineColor(kBlue);
  HST2->DrawCopy("h");

  hZero->DrawCopy("hsame");

  /// OBSOLETE errors are huge
  //TH1D *RAT_21 = (TH1D*)RAT2->Clone("RAT2");
  //RAT_21->Divide(RAT1);             /// divide over CEEX2
  //RAT_21->SetLineColor(kRed);
  //RAT_21->DrawCopy("hsame");  /// errors are huge
  /// errors ok
  //DELnu->SetLineColor(kBlue);
  //DELnu->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.10,0.92,"         h.o. QED effect (#Delta R)/R,     (CEEX1-CEEX2)/CEEX2");
  CaptE->DrawLatex(0.80,0.85, capt1);
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

  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  int Nodes;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);

  HistNormalize();     // Renormalization of MC histograms
  //========== PLOTTING ==========
  FigInfo();
  FigCEEX21();
  FigCEEX21mu();
  FigNuDiff();    /// t-chanel cotrib.
  FigCEEX21rat(); /// h.o. QED
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
//
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
