//    make PlotAFB2
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
// Archive
// Archive
TFile DiskFileA("../workAFB/rmain.root_95GeV_100M");
TFile DiskFileA2("../workAFB/rmain.root_88GeV_100M");
TFile DiskFileA3("../workAFB/rmain.root_91GeV_48M");
TFile DiskFileA4("../workAFB/rmain.root_10GeV_30M");
/// Current
//TFile DiskFileA("../workAFB/rmain.root");
//TFile DiskFileA2("../test0/rmain.root");
TFile DiskFileB("RhoAFB.root","RECREATE","histograms");
//=============================================================================

//Double_t sqr( const Double_t x ){ return x*x;};
// Auxiliary procedures for plotting
//#include "HisNorm.h"

///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  // renomalizing scattergrams is no necesssary!
  //  BIG scatergrams
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vKcPL_Ceex2") );
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
  TH2D *sct9_vKcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vKcPL_Ceex2");

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
void FigCdifIFI()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCdifIFI =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
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
  TCanvas *cCdifIFI = new TCanvas("cCdifIFI","cCdifIFI", 75,  75,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  cCdifIFI->SetFillColor(10);
  cCdifIFI->Divide( 2,  1);
  //cCdifIFI->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //*****************************************************************************
  cCdifIFI->cd(1);
  Hafb9_vAcPR_Ceex2->SetTitle(0);
  Hafb9_vAcPR_Ceex2->SetStats(0);
  Hafb9_vAcPR_Ceex2->GetXaxis()->SetTitle("v_{max,ALEPH}");
  Hafb9_vAcPR_Ceex2->GetYaxis()->SetTitle("A_{FB}(v_{max})");
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

  CaptT->DrawLatex(0.12,0.95,"A_{FB}(v_{max}), Black=IFIoff,  Blue=IFIon,    v=v_{ALEPH}");
  CaptT->DrawLatex(0.65,0.20," #sqrt{s}=94.3GeV ");
  CaptT->DrawLatex(0.65,0.83," #sqrt{s}=87.9GeV ");
  //*****************************************************************************
  cCdifIFI->cd(2);
  TH1D *hZero = (TH1D*)Hafb8_vAcPR_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  TH1D *Hafb9_vAcPR_IFIdiff= (TH1D*)Hafb9_vAcPR_Ceex2->Clone("Hafb9_vAcPR_IFIdiff");
  Hafb9_vAcPR_IFIdiff->Add(Hafb9_vAcPR_IFIdiff,Hafb9_vAcPR_Ceex2n,1.0,-1.0);

  Hafb9_vAcPR_IFIdiff->SetTitle(0);
  Hafb9_vAcPR_IFIdiff->SetStats(0);
  Hafb9_vAcPR_IFIdiff->GetYaxis()->SetTitle("#Delta A^{IFI}_{FB}(v_{max})");
  Hafb9_vAcPR_IFIdiff->SetMaximum( 0.05);
  Hafb9_vAcPR_IFIdiff->SetMinimum(-0.05);
  Hafb9_vAcPR_IFIdiff->SetLineColor(kBlack);
  Hafb9_vAcPR_IFIdiff->DrawCopy("h");
  //
  TH1D *Hafb8_vAcPR_IFIdiff= (TH1D*)Hafb8_vAcPR_Ceex2->Clone("Hafb8_vAcPR_IFIdiff");
  Hafb8_vAcPR_IFIdiff->Add(Hafb8_vAcPR_IFIdiff,Hafb8_vAcPR_Ceex2n,1.0,-1.0);
  Hafb8_vAcPR_IFIdiff->SetLineColor(kBlue);
  Hafb8_vAcPR_IFIdiff->DrawCopy("hsame");
  //
  TH1D *HafbZ_vAcPR_IFIdiff= (TH1D*)HafbZ_vAcPR_Ceex2->Clone("HafbZ_vAcPR_IFIdiff");
  HafbZ_vAcPR_IFIdiff->Add(HafbZ_vAcPR_IFIdiff,HafbZ_vAcPR_Ceex2n,1.0,-1.0);
  HafbZ_vAcPR_IFIdiff->SetLineColor(kMagenta);
  HafbZ_vAcPR_IFIdiff->DrawCopy("hsame");
  //
  TH1D *Hafb1_vAcPR_IFIdiff= (TH1D*)Hafb1_vAcPR_Ceex2->Clone("Hafb1_vAcPR_IFIdiff");
  Hafb1_vAcPR_IFIdiff->Add(Hafb1_vAcPR_IFIdiff,Hafb1_vAcPR_Ceex2n,1.0,-1.0);
  Hafb1_vAcPR_IFIdiff->SetLineColor(kMagenta);
  Hafb1_vAcPR_IFIdiff->DrawCopy("hsame");
  //
  TH1D *HDifSum = (TH1D*)Hafb8_vAcPR_IFIdiff->Clone("HDifSum");
  HDifSum->Add(  Hafb8_vAcPR_IFIdiff, Hafb9_vAcPR_IFIdiff, 1.0, 1.0);
  HDifSum->SetLineColor(kRed);
  HDifSum->SetLineWidth(2);
  HDifSum->DrawCopy("hsame");
  //
  hZero->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.18,0.95,"#Delta A^{IFI}_{FB} =A^{IFIon}_{FB}(v_{max})-A^{IFIoff}_{FB}(v_{max}),     v=v_{ALEPH}");
  CaptT->DrawLatex(0.20,0.85,"Black: #sqrt{s_{+}}=94.3GeV,    Blue: #sqrt{s_{-}}=87.9GeV");
  CaptT->DrawLatex(0.25,0.80,"Magenta: #sqrt{s}=91.2GeV, 10GeV");
  CaptT->DrawLatex(0.35,0.75,
     "Red: #Delta A^{IFI}_{FB}(s_{+}) - #Delta A^{IFI}_{FB}(s_{-})");
  cCdifIFI->cd();
  cCdifIFI->SaveAs("cCdifIFI.jpg");
//
}// FigCdifIFI


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
  FigCdifIFI();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
