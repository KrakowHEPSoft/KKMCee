//    make PlotAFB-run

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
// New MainKKMC
TFile DiskFileA("../workKKMC/histo.root");  // current
//
//TFile DiskFileA("../workKKMC/histo.root_88GeV.new");  // current
//TFile DiskFileA("../workKKMC/histo.root_95GeV.new");  // current
// Archive
//TFile DiskFileA("../workKKMC/histo.root_95GeV_26G");  // last
//TFile DiskFileA("../workKKMC/histo.root_95GeV.4G"); // newest!
//TFile DiskFileA("../workKKMC/histo.root_95GeV_1200M");
///////////////////////////
// Old Prod0
//TFile DiskFileA("../test0/rmain.root");
//TFile DiskFileA("../workAFB/rmain.root");
// Archive
//TFile DiskFileA("../workAFB/rmain.root_95GeV_100M"); //
//TFile DiskFileA("../workAFB/rmain.root_189GeV_100M"); //
//TFile DiskFileA("../test0/rmain.root_88GeV_100M"); // lacks sct_vAcPL_Ceex2
//TFile DiskFileA("../workAFB/rmain.root_91GeV_48M"); //
//TFile DiskFileA("../workAFB/rmain.root_10GeV_30M"); //
// special ISR only
//TFile DiskFileA("../workAFB/rmain.root_95GeV_ISRonly_1M"); //
//TFile DiskFileA("../workAFB/rmain.root_95GeV_ISR-EEX_1M"); //
//TFile DiskFileA("../workAFB/rmain.root_95GeV_IFIoff_1M"); //
//TFile DiskFileA("../test0/rmain.root_88GeV_ISRonly_1M"); //
//TFile DiskFileA("../test0/rmain.root_88GeV_ISR-EEX_1M"); //
//TFile DiskFileA("../test0/rmain.root_88GeV_IFIoff_1M"); //
TFile DiskFileB("RhoAFB.root","RECREATE","histograms");
//=============================================================================
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
//
float  gXcanv = 50, gYcanv = 50;

//Double_t sqr( const Double_t x ){ return x*x;};

///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  // 1-dim histos
  /*
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueMain") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vACeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vACeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vACeex21F") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vACeex21B") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vACeex21") );
  */
  //
  //  BIG scatergrams
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vKcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vAcPL_Ceex2") );
}

///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
// Some MC histos from KKMC are preprocessed
//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  // CMSene=xpar(1) stored in normalization histo
  CMSene  = HST_KKMC_NORMA->GetBinContent(1);
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  //****************************************************************************************
  // Pure MC reprocessing part
  //
  //
  TH2D *sct_vAcPR_Ceex2  = (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2");
  TH2D *sct_vAcPR_Ceex2n = (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2n");
  TH2D *sct_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2");
  TH2D *sct_vKcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vKcPL_Ceex2");
  TH2D *sct_vAcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vAcPL_Ceex2");

  // ****************************************************************************************
  /// Distributions of v with limited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  int nbMax=0;        // cosThetaMax = 1.0, nno cut
  nbMax=50;           // cosThetaMax = 50/50=1.00
  nbMax=45;           // cosThetaMax = 45/50=0.90
  TH1D                    *Hsig_vAcPR_Ceex2, *Hafb_vAcPR_Ceex2;
  ProjV( sct_vAcPR_Ceex2,  Hsig_vAcPR_Ceex2,  Hafb_vAcPR_Ceex2, nbMax);
  Hsig_vAcPR_Ceex2->SetName("Hsig_vAcPR_Ceex2");
  Hafb_vAcPR_Ceex2->SetName("Hafb_vAcPR_Ceex2");
  //if(CMSene < 91.0 ) Hafb_vAcPR_Ceex2->Scale(-1.0);
  // IFI off
  TH1D                     *Hsig_vAcPR_Ceex2n, *Hafb_vAcPR_Ceex2n;
  ProjV( sct_vAcPR_Ceex2n,  Hsig_vAcPR_Ceex2n,  Hafb_vAcPR_Ceex2n, nbMax);
  Hsig_vAcPR_Ceex2n->SetName("Hsig_vAcPR_Ceex2n");
  Hafb_vAcPR_Ceex2n->SetName("Hafb_vAcPR_Ceex2n");
  //if(CMSene < 91.0 ) Hafb_vAcPR_Ceex2n->Scale(-1.0);
  // Different variable: vTrue, bare muons
  TH1D                    *Hsig_vTcPL_Ceex2,  *Hafb_vTcPL_Ceex2;
  ProjV( sct_vTcPL_Ceex2,  Hsig_vTcPL_Ceex2,   Hafb_vTcPL_Ceex2, nbMax);
  Hsig_vTcPL_Ceex2->SetName("Hsig_vTcPL_Ceex2");
  Hafb_vTcPL_Ceex2->SetName("Hafb_vTcPL_Ceex2");
  //if(CMSene < 91.0 ) Hafb_vTcPL_Ceex2->Scale(-1.0);
  // Different variable: vKarlud, pure ISR (unphysical)
  TH1D                    *Hsig_vKcPL_Ceex2,  *Hafb_vKcPL_Ceex2;
  ProjV( sct_vKcPL_Ceex2,  Hsig_vKcPL_Ceex2,   Hafb_vKcPL_Ceex2, nbMax);
  Hsig_vKcPL_Ceex2->SetName("Hsig_vKcPL_Ceex2");
  Hafb_vKcPL_Ceex2->SetName("Hafb_vKcPL_Ceex2");
  //if(CMSene < 91.0 ) Hafb_vKcPL_Ceex2->Scale(-1.0);
  // Different cost(theta) of PL type
  TH1D                    *Hsig_vAcPL_Ceex2,  *Hafb_vAcPL_Ceex2;
  ProjV( sct_vAcPL_Ceex2,  Hsig_vAcPL_Ceex2,   Hafb_vAcPL_Ceex2, nbMax);
  Hsig_vAcPL_Ceex2->SetName("Hsig_vAcPL_Ceex2");
  Hafb_vAcPL_Ceex2->SetName("Hafb_vAcPL_Ceex2");
  //if(CMSene < 91.0 ) Hafb_vAcPL_Ceex2->Scale(-1.0);

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
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cTempl = new TCanvas("cTempl","cTempl", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
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
void FigScatA()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigScat =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  //
  TH2D *sct_vAcPR_Ceex2  = (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2");
  TH2D *sct_vAcPR_Ceex2n = (TH2D*)DiskFileA.Get("sct_vAcPR_Ceex2n");
  TH2D *sct_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2");
  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cScatA = new TCanvas("cScatA","cScatA", gXcanv,  gXcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cScatA->SetFillColor(10);
  cScatA->Divide( 2,  0);
  //cScatA->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  TString OptSurf;
  //OptSurf="      "; // 2D scatergram, points
  //OptSurf="col   "; // 2D histogram, color
  //OptSurf="colz  "; // 2D kolorowe paski, ze skala
  //OptSurf="surf1 "; // 3D surface color
  OptSurf="lego2 "; // 3D histogram color
  //OptSurf="surf3 "; // 3D histogram, z plotem "na dachu"
  //OptSurf="surf2z"; // 3D kolorowe paski, ze skala
  //OptSurf="surf2 "; // 3D kolorowe paski bez skali
  //OptSurf="surf4 "; // 3D gladka powierchnia
  //-------------------------------------
  TH2D *SCT1, *SCT2, *SCT3, *SCT4;
  cScatA->cd(1);
  SCT1=sct_vAcPR_Ceex2;
  //SCT1=sct_vAcPR_Ceex2n;
  gPad->SetTheta(25);
  gPad->SetPhi( -38);
  //
  gPad->SetLogz(); // !!!!!!
  double zmax = SCT1->GetMaximum();
  SCT1->SetMaximum(zmax*1.2);
  SCT1->SetMinimum(zmax*1e-3);
  SCT1->Draw(OptSurf);
  //-------------------------------------
  cScatA->cd(2);
  SCT2=sct_vTcPL_Ceex2;
  gPad->SetTheta(25);
  gPad->SetPhi( -38);
  gPad->SetLogz(); // !!!!!!
  SCT2->SetMaximum(zmax*1.2);
  SCT2->SetMinimum(zmax*1e-3);
  SCT2->Draw(OptSurf);
  //
  cScatA->cd();
} //FigScatA




///////////////////////////////////////////////////////////////////////////////////
void FigVsig()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVsig =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hsig_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hsig_vAcPR_Ceex2");
  TH1D *Hsig_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hsig_vAcPR_Ceex2n");
  TH1D *Hsig_vKcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hsig_vKcPL_Ceex2");
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cVsig = new TCanvas("cVsig","cVsig", gXcanv,  gYcanv,   1000,  400);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cVsig->SetFillColor(10);
  cVsig->Divide( 2,  0);
  //cVsig->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cVsig->cd(1);
  //gPad->SetLogy(); // !!!!!!
  Hsig_vAcPR_Ceex2->SetTitle(0);
  Hsig_vAcPR_Ceex2->SetStats(0);
  Hsig_vAcPR_Ceex2->SetMinimum(0);
  Hsig_vAcPR_Ceex2->GetXaxis()->SetTitle("v_{max,ALEPH}");
  Hsig_vAcPR_Ceex2->SetLineColor(kBlack);
  Hsig_vAcPR_Ceex2->DrawCopy("h");
  //
  Hsig_vAcPR_Ceex2n->SetLineColor(kRed);
  Hsig_vAcPR_Ceex2n->DrawCopy("hsame");
  CaptT->DrawLatex(0.12,0.95,"#sigma(v_{ALEPH,max}),  Red=IFIon,  Black=IFIoff,  ");
  CaptT->DrawLatex(0.60,0.75, "KKMC ISR+FSR");
  //-------------------------------------
  cVsig->cd(2);
  TH1D *hRat1 = (TH1D*)Hsig_vAcPR_Ceex2n->Clone("hRat1"); // IFI off/on
  hRat1->Divide(Hsig_vAcPR_Ceex2);
  hRat1->SetMaximum(1.20);
  hRat1->SetMinimum(0.90);
  hRat1->SetTitle(0);
  hRat1->SetStats(0);
  hRat1->GetXaxis()->SetTitle("v_{max}");
  hRat1->DrawCopy("h");
  TH1D *hRat2 = (TH1D*)Hsig_vKcPL_Ceex2->Clone("hRat2"); //
  hRat2->Divide(Hsig_vAcPR_Ceex2);
  hRat2->SetLineColor(kBlue);
  hRat2->DrawCopy("hsame");
  CaptT->DrawLatex(0.12,0.95,
    "#sigma_{IFIon}/#sigma_{IFIoff}(v_{max}),  Red: v=v_{ALEPH},  Blue: v=v_{Bare}");
  CaptT->DrawLatex(0.60,0.75,TextEne);
  //
  cVsig->cd();
//
}// FigVsig


///////////////////////////////////////////////////////////////////////////////////
void FigAfbIFI()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfbIFI =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);

  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb_vAcPR_Ceex2");
  TH1D *Hafb_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb_vAcPR_Ceex2n");
  //TH1D *Hafb_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb_vTcPL_Ceex2");
  //TH1D *Hafb_vKcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb_vKcPL_Ceex2");
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbIFI = new TCanvas("cAfbIFI","cAfbIFI", gXcanv,  gYcanv,   1200,  500);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbIFI->SetFillColor(10);
  cAfbIFI->Divide( 2,  1);
  //cAfbIFI->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cAfbIFI->cd(1);
  Hafb_vAcPR_Ceex2->SetTitle(0);
  Hafb_vAcPR_Ceex2->SetStats(0);
  Hafb_vAcPR_Ceex2->GetXaxis()->SetTitle("v_{max,ALEPH}");
  Hafb_vAcPR_Ceex2->SetLineColor(kBlue);
  Hafb_vAcPR_Ceex2->DrawCopy("h");
  //
  Hafb_vAcPR_Ceex2n->SetLineColor(kBlack);
  Hafb_vAcPR_Ceex2n->DrawCopy("hsame");
  CaptT->DrawLatex(0.12,0.95,"A_{FB}(v_{max}), Black=IFIoff,  Blue=IFIon,    v=v_{ALEPH}");
  CaptT->DrawLatex(0.60,0.70,"KKMC ISR+FSR");
  //-------------------------------------
  cAfbIFI->cd(2);
  TH1D *Hafb_vAcPR_IFIdiff= (TH1D*)Hafb_vAcPR_Ceex2->Clone("Hafb_vAcPR_IFIdiff");
  Hafb_vAcPR_IFIdiff->Add(Hafb_vAcPR_IFIdiff,Hafb_vAcPR_Ceex2n,1.0,-1.0);
  Hafb_vAcPR_IFIdiff->SetTitle(0);
  Hafb_vAcPR_IFIdiff->SetStats(0);
  Hafb_vAcPR_IFIdiff->SetLineColor(kBlue);
  Hafb_vAcPR_IFIdiff->DrawCopy("h");
  //
  CaptT->DrawLatex(0.12,0.95,"#Delta A^{IFI}_{FB}(v_{ALEPH,max}) = A^{IFIon}_{FB} - A^{IFIoff}_{FB}");
  CaptT->DrawLatex(0.60,0.70,TextEne);
  //-------------------------------------
  cAfbIFI->cd();
//
}// FigAfbIFI



///////////////////////////////////////////////////////////////////////////////////
void FigAfbKin()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfbKin =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);

  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb_vAcPR_Ceex2");
  TH1D *Hafb_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb_vAcPR_Ceex2n");
  TH1D *Hafb_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb_vTcPL_Ceex2");
  TH1D *Hafb_vKcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb_vKcPL_Ceex2");
  TH1D *Hafb_vAcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb_vAcPL_Ceex2");
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbKin = new TCanvas("cAfbKin","cAfbKin", gXcanv,  gYcanv,   1200,  500);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbKin->SetFillColor(10);
  cAfbKin->Divide( 2,  1);
  //cAfbKin->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //-------------------------------------
  cAfbKin->cd(1);
  Hafb_vAcPR_Ceex2->SetLineColor(kBlue);
  Hafb_vAcPR_Ceex2->DrawCopy("h");
  //
  Hafb_vTcPL_Ceex2->SetLineColor(kRed);
  Hafb_vTcPL_Ceex2->DrawCopy("hsame");
  //
  Hafb_vKcPL_Ceex2->SetLineColor(kBlack);
  Hafb_vKcPL_Ceex2->DrawCopy("hsame");
  // cannot be distinguished
  //Hafb_vAcPL_Ceex2->SetLineColor(kMagenta);
  //Hafb_vAcPL_Ceex2->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.12,0.95,"A_{FB}(v_{max}), Blue=v_{ALEPH}, Red=v_{Bare}, Black=v_{ISR}");
  CaptT->DrawLatex(0.50,0.70,"KKMC ISR+FSR, IFI on");
  //-------------------------------------
  cAfbKin->cd(2);
  TH1D *Hafb_vAcPR_Vdiff= (TH1D*)Hafb_vKcPL_Ceex2->Clone("Hafb_vAcPR_Vdiff");
  Hafb_vAcPR_Vdiff->Add(Hafb_vAcPR_Ceex2, Hafb_vKcPL_Ceex2,1.0,-1.0);
  Hafb_vAcPR_Vdiff->SetLineColor(kBlue);
  //
  TH1D *Hafb_vTcPR_Vdiff= (TH1D*)Hafb_vKcPL_Ceex2->Clone("Hafb_vTcPR_Vdiff");
  Hafb_vTcPR_Vdiff->Add(Hafb_vTcPL_Ceex2, Hafb_vKcPL_Ceex2,1.0,-1.0);
  Hafb_vTcPR_Vdiff->SetLineColor(kRed);
  //
  Hafb_vAcPR_Vdiff->SetTitle(0);
  Hafb_vAcPR_Vdiff->SetStats(0);
  Hafb_vAcPR_Vdiff->SetMaximum( 0.020);
  Hafb_vAcPR_Vdiff->SetMinimum(-0.020);
  Hafb_vAcPR_Vdiff->DrawCopy("h");
  Hafb_vTcPR_Vdiff->DrawCopy("hsame");
  //
  TH1D *Hafb_vAcPL_Ceex2diff= (TH1D*)Hafb_vAcPL_Ceex2->Clone("Hafb_vAcPL_Ceex2diff");
  Hafb_vAcPL_Ceex2diff->Add( Hafb_vAcPL_Ceex2diff , Hafb_vAcPR_Ceex2 , 1.0,-1.0);
  Hafb_vAcPL_Ceex2diff->SetLineColor(kMagenta);
  Hafb_vAcPL_Ceex2diff->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"#DeltaA_{FB}(v_{max}),  Blue=ALEPH-ISR,  Red=Bare-ISR");
  CaptT->DrawLatex(0.60,0.70,TextEne);
  //
  cAfbKin->cd();
//
}// FigAfbKin



///////////////////////////////////////////////////////////////////////////////////
void Fig_vA_Ceex21()
{
//------------------------------------------------------------------------
  cout<<" ========================= Fig_vA_Ceex21 =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);

// AFB with v_Aleph and ThetaPRD
  TH1D *hst_vACeex1    = (TH1D*)DiskFileA.Get("hst_vACeex1");   // total CEEX1
  TH1D *hst_vACeex2    = (TH1D*)DiskFileA.Get("hst_vACeex2");   // total CEEX2
  TH1D *hst_vACeex2F   = (TH1D*)DiskFileA.Get("hst_vACeex2F");  // total CEEX2 Forward
  TH1D *hst_vACeex21   = (TH1D*)DiskFileA.Get("hst_vACeex21");  // total CEEX2
  TH1D *hst_vACeex21F  = (TH1D*)DiskFileA.Get("hst_vACeex21F"); // CEEX2-CEEX1 Forward
// crudely approximate
  //TH1D *HAfb_vACeex21 = HstAFB3( "HAfb_vACeex21", hst_vACeex21F, hst_vACeex21, hst_vACeex2 );
// well approximate
  TH1D *HAfb_vACeex21 = HstAFB4( "HAfb_vACeex21", hst_vACeex21F, hst_vACeex21, hst_vACeex2F, hst_vACeex2 );

  TH1D *hZero      = (TH1D*)hst_vACeex2->Clone("hZero");  // zero line
  TH1D *hZeroPlus  = (TH1D*)hst_vACeex2->Clone("hZeroPlus");  //
  TH1D *hZeroMinus = (TH1D*)hst_vACeex2->Clone("hZeroMinus");  //
  for(int i=1; i <= hZero->GetNbinsX() ; i++) {
    hZero->SetBinContent(i, 0);          hZero->SetBinError(i, 0);
    hZeroPlus->SetBinContent(i,  3e-5);  hZeroPlus->SetBinError(i, 0);
    hZeroMinus->SetBinContent(i,-3e-5);  hZeroMinus->SetBinError(i, 0);
    }

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cDifCeex = new TCanvas("cDifCeex","cDifCeex",gXcanv, gXcanv,   600,  500);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cDifCeex->SetFillColor(10);
  cDifCeex->Divide( 1,  0);
  //cDifCeex->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cDifCeex->cd(1);
  //gPad->SetLogy(); // !!!!!!

  TH1D *Hst = HAfb_vACeex21;
  Hst->SetStats(0);
  Hst->SetTitle(0);
//  Hst->SetMaximum( 0.3e-3); Hst->SetMinimum(-0.3e-3);
//  Hst->SetMaximum( 0.0015); Hst->SetMinimum(-0.0015);
  Hst->GetXaxis()->SetTitleSize(0.05);
  //Hst->GetXaxis()->SetTitle("v_{max}");
  Hst->GetXaxis()->SetTitle("v_{max,ALEPH}");
  Hst->DrawCopy("h");

  hZero->SetLineColor(kBlack);
  hZero->DrawCopy("hsame");
  hZeroPlus->SetLineColor(kRed);
  hZeroPlus->DrawCopy("hsame");
  hZeroMinus->SetLineColor(kRed);
  hZeroMinus->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"#Delta A_{FB}(v_{max}),    CEEX2-CEEX1,    KKMC ISR+FSR+IFI");
//  CaptT->DrawLatex(0.12,0.95,"#Delta A_{FB}(v_{max}),    CEEX2-CEEX1,    KKMC ISR only");
//  CaptT->DrawLatex(0.12,0.95,"#Delta A_{FB}(v_{max}),    CEEX2-CEEX1,    KKMC ISR+FSR, IFI off");
//  CaptT->DrawLatex(0.12,0.95,"#Delta A_{FB}(v_{max}),    EEX3-EEX2,    KKMC ISR only");
  CaptT->DrawLatex(0.25,0.85,TextEne);
  CaptT->DrawLatex(0.12,0.55," #delta#alpha/#alpha = 10^{-4}");
  //-------------------------------------
  //cDifCeex->cd(2);
  //gPad->SetLogy(); // !!!!!!
  //
  cDifCeex->cd();
//
}// Fig_vA_Ceex21



///////////////////////////////////////////////////////////////////////////////////
void Fig_vT_Ceex21()
{
//------------------------------------------------------------------------
  cout<<" ========================= Fig_vT_Ceex21 =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);

// AFB with v_true and costhetaPL
  TH1D *hst_vT_Ceex1    = (TH1D*)DiskFileA.Get("hst_vT_Ceex1");    // total CEEX1
  TH1D *hst_vT_Ceex2    = (TH1D*)DiskFileA.Get("hst_vT_Ceex2");    // total CEEX2
  TH1D *hst_vT_Ceex2_F  = (TH1D*)DiskFileA.Get("hst_vT_Ceex2_F");  // total CEEX2 Forward
  TH1D *hst_vT_Ceex21   = (TH1D*)DiskFileA.Get("hst_vT_Ceex21");   // total CEEX2
  TH1D *hst_vT_Ceex21_F = (TH1D*)DiskFileA.Get("hst_vT_Ceex21_F"); // CEEX2-CEEX1 Forward
// Exact formula for AFB from weight differences
  TH1D *HAfb_vT_Ceex21 = HstAFB4( "HAfb_vT_Ceex21", hst_vT_Ceex21_F, hst_vT_Ceex21, hst_vT_Ceex2_F, hst_vT_Ceex2 );
//
  TH1D *hZero      = (TH1D*)hst_vT_Ceex2->Clone("hZero");  // zero line
  TH1D *hZeroPlus  = (TH1D*)hst_vT_Ceex2->Clone("hZeroPlus");  //
  TH1D *hZeroMinus = (TH1D*)hst_vT_Ceex2->Clone("hZeroMinus");  //
  for(int i=1; i <= hZero->GetNbinsX() ; i++) {
    hZero->SetBinContent(i, 0);          hZero->SetBinError(i, 0);
    hZeroPlus->SetBinContent(i,  3e-5);  hZeroPlus->SetBinError(i, 0);
    hZeroMinus->SetBinContent(i,-3e-5);  hZeroMinus->SetBinError(i, 0);
    }

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cDifCeex = new TCanvas("cDifCeex","cDifCeex",gXcanv, gXcanv,   600,  500);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cDifCeex->SetFillColor(10);
  cDifCeex->Divide( 1,  0);
  //cDifCeex->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cDifCeex->cd(1);
  //gPad->SetLogy(); // !!!!!!

  TH1D *Hst = HAfb_vT_Ceex21;
//  Hst->SetStats(0);
  Hst->SetTitle(0);
//  Hst->SetMaximum( 0.3e-3); Hst->SetMinimum(-0.3e-3);
  Hst->SetMaximum( 1.0e-4); Hst->SetMinimum(-1.0e-4);
  Hst->GetXaxis()->SetTitleSize(0.05);
  //Hst->GetXaxis()->SetTitle("v_{max}");
  Hst->GetXaxis()->SetTitle("v_{max}");
  Hst->DrawCopy("h");

  hZero->SetLineColor(kBlack);
  hZero->DrawCopy("hsame");
  hZeroPlus->SetLineColor(kRed);
  hZeroPlus->DrawCopy("hsame");
  hZeroMinus->SetLineColor(kRed);
  hZeroMinus->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"#Delta A_{FB}(v_{max}),    CEEX2-CEEX1,    KKMC ISR+FSR+IFI");
  CaptT->DrawLatex(0.25,0.85,TextEne);
  CaptT->DrawLatex(0.12,0.55," #delta#alpha/#alpha = 10^{-4}");
  //-------------------------------------
  //cDifCeex->cd(2);
  //gPad->SetLogy(); // !!!!!!
  //
  cDifCeex->cd();
//
}// Fig_vT_Ceex21


///////////////////////////////////////////////////////////////////////////////////
void FigPhi()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigPhi =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted

  TH1D *hst_phi1_Ceex2_pol    = (TH1D*)DiskFileA.Get("hst_phi1_Ceex2_pol");

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPhi = new TCanvas("cPhi","cPhi", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cPhi->SetFillColor(10);
//  cPhi->Divide( 2,  0);
  //cPhi->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cPhi->cd(1);

  TH1D *Hst = hst_phi1_Ceex2_pol;
//  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitleSize(0.05);
  Hst->GetXaxis()->SetTitle(" #phi ");
  Hst->DrawCopy("h");

  CaptT->DrawLatex(0.12,0.95,"d#sigma/d#phi:         P_{y}^{1}=+0.10,    P_{y}^{2}=-0.10 ");
  //-------------------------------------
  cPhi->cd(2);

  cPhi->SaveAs("cPhi.pdf");

  //
  cPhi->cd();
//
}// FigPhi



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  ReMakeMChisto();     // reprocessing MC histos
  //========== PLOTTING ==========
  // Template empty canvas  with 2 figures
  //FigTempl();
  FigPhi();
  // Raw MC double distr. of v and cost(theta)
  //FigScatA();
  //FigVsig();
  //FigAfbIFI();
  //FigAfbKin();
  //
  //Fig_vT_Ceex21();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
