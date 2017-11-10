//    make PlotAFB2-run
//    Plot for FCC week 2016 in Rome and "Beyond precision..." paper

// This is only AFB study for several energies,
// renomalizing scattergrams is no necesssary!

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

#include "KKplot.h"
#include "HisNorm.h"

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
// Latest from /workKKMC
TFile DiskFileA("../workKKMC/histo.root_95GeV_26G");   // oct.2017
TFile DiskFileA2("../workKKMC/histo.root_88GeV_2.5G"); // oct.2017
TFile DiskFileA3("../workKKMC/histo.root_91GeV_3.5G"); // oct.2017
//TFile DiskFileA4("../workKKMC/histo.root_10GeV_5.8G"); // oct.2017
TFile DiskFileA4("../workKKMC/histo.root"); // oct.2017

// Archive from /workAFB, as for Rome
//TFile DiskFileA("../workAFB/rmain.root_95GeV_100M");
//TFile DiskFileA2("../workAFB/rmain.root_88GeV_100M");
//TFile DiskFileA3("../workAFB/rmain.root_91GeV_48M");
//TFile DiskFileA4("../workAFB/rmain.root_10GeV_30M");
/// Current
//TFile DiskFileA("../workAFB/rmain.root");
//TFile DiskFileA2("../test0/rmain.root");
TFile DiskFileB("RhoAFB.root","RECREATE","histograms");
//=============================================================================

//Double_t sqr( const Double_t x ){ return x*x;};
// Auxiliary procedures for plotting
//#include "HisNorm.h"

///////////////////////////////////////////////////////////////////////////////////
void PlotSame2(TH1D *HST, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
  double xcapt = 0.50;
  CaptT->DrawLatex(xcapt,ycapt, opis);
  CaptT->DrawLatex(xcapt-0.05,ycapt, label);
  //
  TLatex *CaptS = new TLatex();
  CaptS->SetTextSize(0.040);
  CaptS->SetTextAlign(21);
  CaptS->SetTextColor(kolor);
  int ib = HST->FindBin(xx);
  double yy= HST->GetBinContent(ib);
  CaptS->DrawLatex(xx,yy,label);
}// PlotSame2




///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  // renomalizing scattergrams is not necesssary!
  //  BIG scatergrams
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2") );
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2n") );
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2") );
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vKcPL_Ceex2") );
}

///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	// Here we produce semianalytical plots using KKsem program, No plotting
	// also some MC histos are preprocessed
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR

  //****************************************************************************************
  // Pure MC reprocessing part
  //
  TH2D *sct9_vAcPR_Ceex2  = (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2");
  TH2D *sct9_vAcPR_Ceex2n = (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2n");
  //
  TH2D *sct8_vAcPR_Ceex2  = (TH2D*)DiskFileA2.Get("sct_vAcPR_Ceex2");
  TH2D *sct8_vAcPR_Ceex2n = (TH2D*)DiskFileA2.Get("sct_vAcPR_Ceex2n");
  //
  TH2D *sctZ_vAcPR_Ceex2  = (TH2D*)DiskFileA3.Get("sct_vAcPR_Ceex2");
  TH2D *sctZ_vAcPR_Ceex2n = (TH2D*)DiskFileA3.Get("sct_vAcPR_Ceex2n");
  ///
  TH2D *sct1_vAcPR_Ceex2  = (TH2D*)DiskFileA4.Get("sct_vAcPR_Ceex2");
  TH2D *sct1_vAcPR_Ceex2n = (TH2D*)DiskFileA4.Get("sct_vAcPR_Ceex2n");
  //!!!!!
  TH2D *sct9_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2");
  TH2D *sct9_vTcPL_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2n");
  //
  TH2D *sct8_vTcPL_Ceex2  = (TH2D*)DiskFileA2.Get("sct_vTcPL_Ceex2");
  TH2D *sct8_vTcPL_Ceex2n = (TH2D*)DiskFileA2.Get("sct_vTcPL_Ceex2n");
  //
  TH2D *sctZ_vTcPL_Ceex2  = (TH2D*)DiskFileA3.Get("sct_vTcPL_Ceex2");
  TH2D *sctZ_vTcPL_Ceex2n = (TH2D*)DiskFileA3.Get("sct_vTcPL_Ceex2n");
  ///
  TH2D *sct1_vTcPL_Ceex2  = (TH2D*)DiskFileA4.Get("sct_vTcPL_Ceex2");
  TH2D *sct1_vTcPL_Ceex2n = (TH2D*)DiskFileA4.Get("sct_vTcPL_Ceex2n");

  // ****************************************************************************************
  /// Distributions of v with limited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  int nbMax=0;   // cosThetaMax = 1.0, no cut
  nbMax=50;      // cosThetaMax = 50/50=1.00
  nbMax=45;      // cosThetaMax = 45/50=0.90
  // ---------------------- 95GeV ----------------------------------
  TH1D  *Hsig9_vAcPR_Ceex2  = HstProjV("Hsig9_vAcPR_Ceex2", sct9_vAcPR_Ceex2, nbMax);
  TH1D  *Hafb9_vAcPR_Ceex2  = HstProjA("Hafb9_vAcPR_Ceex2", sct9_vAcPR_Ceex2, nbMax);
  //******** IFI off
  TH1D  *Hsig9_vAcPR_Ceex2n  = HstProjV("Hsig9_vAcPR_Ceex2n", sct9_vAcPR_Ceex2n, nbMax);
  TH1D  *Hafb9_vAcPR_Ceex2n  = HstProjA("Hafb9_vAcPR_Ceex2n", sct9_vAcPR_Ceex2n, nbMax);
  //
  // ---------------------- 88GeV ----------------------------------
  TH1D  *Hsig8_vAcPR_Ceex2  = HstProjV("Hsig8_vAcPR_Ceex2", sct8_vAcPR_Ceex2, nbMax);
  TH1D  *Hafb8_vAcPR_Ceex2  = HstProjA("Hafb8_vAcPR_Ceex2", sct8_vAcPR_Ceex2, nbMax);
  Hafb8_vAcPR_Ceex2->Scale(-1.0);
  //******** IFI off
  TH1D  *Hsig8_vAcPR_Ceex2n  = HstProjV("Hsig8_vAcPR_Ceex2n", sct8_vAcPR_Ceex2n, nbMax);
  TH1D  *Hafb8_vAcPR_Ceex2n  = HstProjA("Hafb8_vAcPR_Ceex2n", sct8_vAcPR_Ceex2n, nbMax);
  Hafb8_vAcPR_Ceex2n->Scale(-1.0);
  //
  // ---------------------- 91.2GeV ----------------------------------
  TH1D  *HsigZ_vAcPR_Ceex2  = HstProjV("HsigZ_vAcPR_Ceex2", sctZ_vAcPR_Ceex2, nbMax);
  TH1D  *HafbZ_vAcPR_Ceex2  = HstProjA("HafbZ_vAcPR_Ceex2", sctZ_vAcPR_Ceex2, nbMax);
  //******** IFI off
  TH1D  *HsigZ_vAcPR_Ceex2n  = HstProjV("HsigZ_vAcPR_Ceex2n", sctZ_vAcPR_Ceex2n, nbMax);
  TH1D  *HafbZ_vAcPR_Ceex2n  = HstProjA("HafbZ_vAcPR_Ceex2n", sctZ_vAcPR_Ceex2n, nbMax);
  // ---------------------- 10GeV ----------------------------------
  TH1D  *Hsig1_vAcPR_Ceex2  = HstProjV("Hsig1_vAcPR_Ceex2", sct1_vAcPR_Ceex2, nbMax);
  TH1D  *Hafb1_vAcPR_Ceex2  = HstProjA("Hafb1_vAcPR_Ceex2", sct1_vAcPR_Ceex2, nbMax);
  //******** IFI off
  TH1D  *Hsig1_vAcPR_Ceex2n  = HstProjV("Hsig1_vAcPR_Ceex2n", sct1_vAcPR_Ceex2n, nbMax);
  TH1D  *Hafb1_vAcPR_Ceex2n  = HstProjA("Hafb1_vAcPR_Ceex2n", sct1_vAcPR_Ceex2n, nbMax);

  // Warning! scattergrams for vTcPL noIFI were missing in Rome version
  // ---------------------- 95GeV ----------------------------------
  TH1D  *Hsig9_vTcPL_Ceex2  = HstProjV("Hsig9_vTcPL_Ceex2", sct9_vTcPL_Ceex2, nbMax);
  TH1D  *Hafb9_vTcPL_Ceex2  = HstProjA("Hafb9_vTcPL_Ceex2", sct9_vTcPL_Ceex2, nbMax);
  //******** IFI off
  TH1D  *Hsig9_vTcPL_Ceex2n  = HstProjV("Hsig9_vTcPL_Ceex2n", sct9_vTcPL_Ceex2n, nbMax); //???? not present
  TH1D  *Hafb9_vTcPL_Ceex2n  = HstProjA("Hafb9_vTcPL_Ceex2n", sct9_vTcPL_Ceex2n, nbMax);
  //
  // ---------------------- 88GeV ----------------------------------
  TH1D  *Hsig8_vTcPL_Ceex2  = HstProjV("Hsig8_vTcPL_Ceex2", sct8_vTcPL_Ceex2, nbMax);
  TH1D  *Hafb8_vTcPL_Ceex2  = HstProjA("Hafb8_vTcPL_Ceex2", sct8_vTcPL_Ceex2, nbMax);
  Hafb8_vTcPL_Ceex2->Scale(-1.0);
  //******** IFI off
  TH1D  *Hsig8_vTcPL_Ceex2n  = HstProjV("Hsig8_vTcPL_Ceex2n", sct8_vTcPL_Ceex2n, nbMax);
  TH1D  *Hafb8_vTcPL_Ceex2n  = HstProjA("Hafb8_vTcPL_Ceex2n", sct8_vTcPL_Ceex2n, nbMax);
  Hafb8_vTcPL_Ceex2n->Scale(-1.0);
  //
  // ---------------------- 91.2GeV ----------------------------------
  TH1D  *HsigZ_vTcPL_Ceex2  = HstProjV("HsigZ_vTcPL_Ceex2", sctZ_vTcPL_Ceex2, nbMax);
  TH1D  *HafbZ_vTcPL_Ceex2  = HstProjA("HafbZ_vTcPL_Ceex2", sctZ_vTcPL_Ceex2, nbMax);
  //******** IFI off
  TH1D  *HsigZ_vTcPL_Ceex2n  = HstProjV("HsigZ_vTcPL_Ceex2n", sctZ_vTcPL_Ceex2n, nbMax);
  TH1D  *HafbZ_vTcPL_Ceex2n  = HstProjA("HafbZ_vTcPL_Ceex2n", sctZ_vTcPL_Ceex2n, nbMax);
  // ---------------------- 10GeV ----------------------------------
  TH1D  *Hsig1_vTcPL_Ceex2  = HstProjV("Hsig1_vTcPL_Ceex2", sct1_vTcPL_Ceex2, nbMax);
  TH1D  *Hafb1_vTcPL_Ceex2  = HstProjA("Hafb1_vTcPL_Ceex2", sct1_vTcPL_Ceex2, nbMax);
  //******** IFI off
  TH1D  *Hsig1_vTcPL_Ceex2n  = HstProjV("Hsig1_vTcPL_Ceex2n", sct1_vTcPL_Ceex2n, nbMax);
  TH1D  *Hafb1_vTcPL_Ceex2n  = HstProjA("Hafb1_vTcPL_Ceex2n", sct1_vTcPL_Ceex2n, nbMax);


}//ReMakeMChisto


///////////////////////////////////////////////////////////////////////////////////
void FigTempl()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTempl =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cTempl = new TCanvas("cTempl","cTempl", 25,  25,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  cTempl->SetFillColor(10);
  cTempl->Divide( 2,  0);
  //cTempl->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cTempl->cd(1);
  CaptT->DrawLatex(0.12,0.95,"A_{FB}(v_{max}), ????");
  //-------------------------------------
  cTempl->cd(2);
  //
  cTempl->cd();
//
}// FigTempl



///////////////////////////////////////////////////////////////////////////////////
void FigIFIvAa()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigIFIvAa =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb9_vAcPR_Ceex2");
  TH1D *Hafb9_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb9_vAcPR_Ceex2n");
  //
  TH1D *Hafb8_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb8_vAcPR_Ceex2");
  TH1D *Hafb8_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb8_vAcPR_Ceex2n");
  //
  TH1D *HafbZ_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("HafbZ_vAcPR_Ceex2");
  TH1D *HafbZ_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("HafbZ_vAcPR_Ceex2n");
  //
  TH1D *Hafb1_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb1_vAcPR_Ceex2");
  TH1D *Hafb1_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb1_vAcPR_Ceex2n");
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigIFIvAa = new TCanvas("cFigIFIvAa","cFigIFIvAa", 25,  25,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  cFigIFIvAa->SetFillColor(10);
  cFigIFIvAa->Divide( 2,  1);
//*****************************************************************************
  cFigIFIvAa->cd(1);
  Hafb9_vAcPR_Ceex2->SetTitle(0);
  Hafb9_vAcPR_Ceex2->SetStats(0);
  Hafb9_vAcPR_Ceex2->GetXaxis()->SetTitle("v_{max}");
  //Hafb9_vAcPR_Ceex2->GetYaxis()->SetTitle("A_{FB}(v_{max})");
  Hafb9_vAcPR_Ceex2->SetLineColor(kBlue);
  Hafb9_vAcPR_Ceex2->SetMaximum( 0.33);
  Hafb9_vAcPR_Ceex2->SetMinimum( 0.15);
  Hafb9_vAcPR_Ceex2->DrawCopy("h");
  //
  Hafb9_vAcPR_Ceex2n->SetLineColor(kBlack);
  Hafb9_vAcPR_Ceex2n->DrawCopy("hsame");
  //
  Hafb8_vAcPR_Ceex2->SetLineColor(kBlue);
  Hafb8_vAcPR_Ceex2->DrawCopy("hsame");
  Hafb8_vAcPR_Ceex2n->SetLineColor(kBlack);
  Hafb8_vAcPR_Ceex2n->DrawCopy("hsame");

  CaptT->DrawLatex(0.02,0.95," Black=IFIoff,  Blue=IFIon, v=v_{ALEPH}");
  CaptT->DrawLatex(0.50,0.20,"  A_{FB}(v_{max}), #sqrt{s}=94.3GeV ");
  CaptT->DrawLatex(0.50,0.83," -A_{FB}(v_{max}), #sqrt{s}=87.9GeV ");
  //*****************************************************************************
  cFigIFIvAa->cd(2);
  TH1D *hZero = (TH1D*)Hafb8_vAcPR_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  TH1D *Hafb9_vAcPR_IFIdiff = HstDiff("Hafb9_vAcPR_IFIdiff",  Hafb9_vAcPR_Ceex2, Hafb9_vAcPR_Ceex2n,  kBlack);
  TH1D *Hafb8_vAcPR_IFIdiff = HstDiff("Hafb8_vAcPR_IFIdiff",  Hafb8_vAcPR_Ceex2, Hafb8_vAcPR_Ceex2n,  kBlue);
  Hafb8_vAcPR_IFIdiff->Scale(-1.0); // undoing sign change
  TH1D *HafbZ_vAcPR_IFIdiff = HstDiff("HafbZ_vAcPR_IFIdiff",  HafbZ_vAcPR_Ceex2, HafbZ_vAcPR_Ceex2n,  kMagenta);
  TH1D *Hafb1_vAcPR_IFIdiff = HstDiff("Hafb1_vAcPR_IFIdiff",  Hafb1_vAcPR_Ceex2, Hafb1_vAcPR_Ceex2n,  kGreen);
  TH1D *HDifSum             = HstDiff("HDifSum",            Hafb9_vAcPR_IFIdiff, Hafb8_vAcPR_IFIdiff,  kRed);
  //HDifSum->SetLineWidth(2);

//  TH1D *Ddiff = Hafb9_vAcPR_IFIdiff;
  TH1D *Ddiff = Hafb1_vAcPR_IFIdiff;
  Ddiff->SetTitle(0);
  Ddiff->SetStats(0);
  //Ddiff->GetYaxis()->SetTitle("#Delta A^{IFI}_{FB}(v_{max})");

  Ddiff->SetMaximum( 0.05); Ddiff->SetMinimum(-0.03);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, " A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  double ycapt =0.33;
  PlotSame2(Hafb9_vAcPR_IFIdiff, ycapt, kBlack,   0.020, "(a)", "#sqrt{s}=94.3GeV");
  PlotSame2(Hafb8_vAcPR_IFIdiff, ycapt, kBlue,    0.045, "(b)", "#sqrt{s}=87.9GeV");
  PlotSame2(HafbZ_vAcPR_IFIdiff, ycapt, kMagenta, 0.060, "(c)", "#sqrt{s}=M_{Z}");
  PlotSame2(Hafb1_vAcPR_IFIdiff, ycapt, kGreen,   0.090, "(d)", "#sqrt{s}=10GeV");
  PlotSame2(HDifSum,             ycapt, kRed,     0.030, "(e)", "= (a) - (b) ");
  //
  hZero->DrawCopy("hsame");
  cFigIFIvAa->cd();
//
}// FigIFIvAa


///////////////////////////////////////////////////////////////////////////////////
void FigIFIvAb()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigIFIvAb =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb9_vAcPR_Ceex2");
  TH1D *Hafb9_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb9_vAcPR_Ceex2n");
  //
  TH1D *Hafb8_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb8_vAcPR_Ceex2");
  TH1D *Hafb8_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb8_vAcPR_Ceex2n");
  //
  TH1D *HafbZ_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("HafbZ_vAcPR_Ceex2");
  TH1D *HafbZ_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("HafbZ_vAcPR_Ceex2n");
  //
  TH1D *Hafb1_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb1_vAcPR_Ceex2");
  TH1D *Hafb1_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb1_vAcPR_Ceex2n");

  TH1D *Hafb9_vAcPR_IFIdiff = HstDiff("Hafb9_vAcPR_IFIdiff",  Hafb9_vAcPR_Ceex2, Hafb9_vAcPR_Ceex2n,  kBlack);
  TH1D *Hafb8_vAcPR_IFIdiff = HstDiff("Hafb8_vAcPR_IFIdiff",  Hafb8_vAcPR_Ceex2, Hafb8_vAcPR_Ceex2n,  kBlue);
  Hafb8_vAcPR_IFIdiff->Scale(-1.0); // undoing sign change
  TH1D *HafbZ_vAcPR_IFIdiff = HstDiff("HafbZ_vAcPR_IFIdiff",  HafbZ_vAcPR_Ceex2, HafbZ_vAcPR_Ceex2n,  kMagenta);
  TH1D *Hafb1_vAcPR_IFIdiff = HstDiff("Hafb1_vAcPR_IFIdiff",  Hafb1_vAcPR_Ceex2, Hafb1_vAcPR_Ceex2n,  kGreen);

  TH1D *hZero = (TH1D*)Hafb8_vAcPR_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigIFIvAb = new TCanvas("cFigIFIvAb","cFigIFIvAb", 50,  50,   600,  600);
  //                                   Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  cFigIFIvAb->SetFillColor(10);

  TH1D *Ddiff = Hafb9_vAcPR_IFIdiff;
  Ddiff->SetTitle(0);
  Ddiff->SetStats(0);

  Ddiff->SetMaximum( 0.05); Ddiff->SetMinimum(-0.03);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, " A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  double ycapt =0.33;
  PlotSame2(Hafb9_vAcPR_IFIdiff, ycapt, kBlack,   0.020, "(a)", "#sqrt{s}=94.3GeV");
  PlotSame2(Hafb8_vAcPR_IFIdiff, ycapt, kBlue,    0.045, "(b)", "#sqrt{s}=87.9GeV");
  PlotSame2(HafbZ_vAcPR_IFIdiff, ycapt, kMagenta, 0.060, "(c)", "#sqrt{s}=M_{Z}");
  PlotSame2(Hafb1_vAcPR_IFIdiff, ycapt, kGreen,   0.090, "(d)", "#sqrt{s}=10GeV");
  hZero->DrawCopy("hsame");

  cFigIFIvAb->SaveAs("cFigIFIvAb.pdf");
  //
  }// FigIFIvAb


///////////////////////////////////////////////////////////////////////////////////
void FigIFIvTa()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigIFIvTa =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex2");
  TH1D *Hafb9_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex2n");
  //
  TH1D *Hafb8_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex2");
  TH1D *Hafb8_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex2n");
  //
  TH1D *HafbZ_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("HafbZ_vTcPL_Ceex2");
  TH1D *HafbZ_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("HafbZ_vTcPL_Ceex2n");
  //
  TH1D *Hafb1_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb1_vTcPL_Ceex2");
  TH1D *Hafb1_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb1_vTcPL_Ceex2n");
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigIFIvTa = new TCanvas("cFigIFIvTa","cFigIFIvTa", 25,  25,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  cFigIFIvTa->SetFillColor(10);
  cFigIFIvTa->Divide( 2,  1);
//*****************************************************************************
  cFigIFIvTa->cd(1);
  Hafb9_vTcPL_Ceex2->SetTitle(0);
  Hafb9_vTcPL_Ceex2->SetStats(0);
  Hafb9_vTcPL_Ceex2->GetXaxis()->SetTitle("v_{max}");
  //Hafb9_vTcPL_Ceex2->GetYaxis()->SetTitle("A_{FB}(v_{max})");
  Hafb9_vTcPL_Ceex2->SetLineColor(kBlue);
  Hafb9_vTcPL_Ceex2->SetMaximum( 0.33);
  Hafb9_vTcPL_Ceex2->SetMinimum( 0.15);
  Hafb9_vTcPL_Ceex2->DrawCopy("h");
  //
  Hafb9_vTcPL_Ceex2n->SetLineColor(kBlack);
  Hafb9_vTcPL_Ceex2n->DrawCopy("hsame");
  //
  Hafb8_vTcPL_Ceex2->SetLineColor(kBlue);
  Hafb8_vTcPL_Ceex2->DrawCopy("hsame");
  Hafb8_vTcPL_Ceex2n->SetLineColor(kBlack);
  Hafb8_vTcPL_Ceex2n->DrawCopy("hsame");

  CaptT->DrawLatex(0.02,0.95," Black=IFIoff,  Blue=IFIon, v=v_{True}");
  CaptT->DrawLatex(0.50,0.20,"  A_{FB}(v_{max}), #sqrt{s}=94.3GeV ");
  CaptT->DrawLatex(0.50,0.83," -A_{FB}(v_{max}), #sqrt{s}=87.9GeV ");
  //*****************************************************************************
  cFigIFIvTa->cd(2);
  TH1D *hZero = (TH1D*)Hafb8_vTcPL_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  TH1D *Hafb9_vTcPL_IFIdiff = HstDiff("Hafb9_vTcPL_IFIdiff",  Hafb9_vTcPL_Ceex2, Hafb9_vTcPL_Ceex2n,  kBlack);
  TH1D *Hafb8_vTcPL_IFIdiff = HstDiff("Hafb8_vTcPL_IFIdiff",  Hafb8_vTcPL_Ceex2, Hafb8_vTcPL_Ceex2n,  kBlue);
  Hafb8_vTcPL_IFIdiff->Scale(-1.0); // undoing sign change
  TH1D *HafbZ_vTcPL_IFIdiff = HstDiff("HafbZ_vTcPL_IFIdiff",  HafbZ_vTcPL_Ceex2, HafbZ_vTcPL_Ceex2n,  kMagenta);
  TH1D *Hafb1_vTcPL_IFIdiff = HstDiff("Hafb1_vTcPL_IFIdiff",  Hafb1_vTcPL_Ceex2, Hafb1_vTcPL_Ceex2n,  kGreen);
  TH1D *HDifSum             = HstDiff("HDifSum",            Hafb9_vTcPL_IFIdiff, Hafb8_vTcPL_IFIdiff,  kRed);
  //HDifSum->SetLineWidth(2);

  TH1D *Ddiff = Hafb9_vTcPL_IFIdiff;
  Ddiff->SetTitle(0);
  Ddiff->SetStats(0);
  //Ddiff->GetYaxis()->SetTitle("#Delta A^{IFI}_{FB}(v_{max})");

  Ddiff->SetMaximum( 0.05); Ddiff->SetMinimum(-0.03);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, " A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  double ycapt =0.33;
  PlotSame2(Hafb9_vTcPL_IFIdiff, ycapt, kBlack,   0.020, "(a)", "#sqrt{s}=94.3GeV");
  PlotSame2(Hafb8_vTcPL_IFIdiff, ycapt, kBlue,    0.045, "(b)", "#sqrt{s}=87.9GeV");
  PlotSame2(HafbZ_vTcPL_IFIdiff, ycapt, kMagenta, 0.060, "(c)", "#sqrt{s}=M_{Z}");
  PlotSame2(Hafb1_vTcPL_IFIdiff, ycapt, kGreen,   0.090, "(d)", "#sqrt{s}=10GeV");
  PlotSame2(HDifSum,             ycapt, kRed,     0.030, "(e)", "= (a) - (b) ");
  //
  hZero->DrawCopy("hsame");
  cFigIFIvTa->cd();
//
}// FigIFIvTa




///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  //KKsemMakeHisto();    // prepare histos for plotting
  ReMakeMChisto();     // reprocessing MC histos
  //========== PLOTTING ==========
  // Template empty canvas  with 2 figures
  //FigTempl();
  FigIFIvAa();
  FigIFIvAb();
  //
  FigIFIvTa();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
