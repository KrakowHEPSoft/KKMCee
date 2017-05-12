//////////////////////////////////////////////////////////////////////
//    make PlotFoam1-run
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
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
#include "TObjString.h"
#include "TFile.h"

#include "HisNorm.h"
//#include "KKabox.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
//
TFile DiskFileA("../workAFB/rmain.root_95GeV_100M");
// current
TFile DiskFileF("../workFOAM/rmain.root");
//
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");

// Interface to KKabox and some extra plotting facilities
//KKabox LibSem;

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
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Cost1Ceex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPLCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPRCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPREex2") );
  //
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2") );
  //
}


///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA5");

  TH1D *HST_xx_Ceex2n = (TH1D*)DiskFileF.Get("HST_xx_Ceex2n");  // FOAM
  TH2D *SCA_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCA_xc_Ceex2n");  // FOAM
  TH2D *SCA_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCA_xc_Ceex2");   // FOAM


  HisNorm1(HST_FOAM_NORMA3, HST_xx_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA3, SCA_xc_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA5, SCA_xc_Ceex2 );   // normalizing

  // sigma(vmax) direct histogramming
  TH1D *HST_xmax_Ceex2n;
  MakeCumul(HST_xx_Ceex2n,HST_xmax_Ceex2n);
  HST_xmax_Ceex2n->SetName("HST_xmax_Ceex2n");

  // sigma(vmax) and AFB(vmax) from sacttergram
  int nbMax=0;   // <--- CosThetaMax = 1.0
  TH1D                 *Htot_xmax_Ceex2n, *Hafb_xmax_Ceex2n;
  ProjV( SCA_xc_Ceex2n, Htot_xmax_Ceex2n,  Hafb_xmax_Ceex2n, nbMax);  //!!!!
  Htot_xmax_Ceex2n->SetName("Htot_xmax_Ceex2n");
  Hafb_xmax_Ceex2n->SetName("Hafb_xmax_Ceex2n");
  //
  TH1D                *Htot_xmax_Ceex2, *Hafb_xmax_Ceex2;
  ProjV( SCA_xc_Ceex2, Htot_xmax_Ceex2,  Hafb_xmax_Ceex2, nbMax);  //!!!!
  Htot_xmax_Ceex2->SetName("Htot_xmax_Ceex2");
  Hafb_xmax_Ceex2->SetName("Hafb_xmax_Ceex2");

////////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;

  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR

  //****************************************************************************************
  // Pure MC reprocessing part
  //
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n");

  ///****************************************************************************************
  /// Distributions of v=vTrue with unlimited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  nbMax=0;   // cosThetaMax = 1.0
  TH1D                    *HTot_vTcPR_Eex2, *HAfb_vTcPR_Eex2;
  ProjV( sca_vTcPR_Eex2,  HTot_vTcPR_Eex2,  HAfb_vTcPR_Eex2, nbMax);  //!!!!
  HTot_vTcPR_Eex2->SetName("HTot_vTcPR_Eex2");
  HAfb_vTcPR_Eex2->SetName("HAfb_vTcPR_Eex2");
  nbMax=0;   // cosThetaMax = 1.0
  TH1D                    *HTot_vTcPR_Ceex2, *HAfb_vTcPR_Ceex2;
  ProjV( sca_vTcPR_Ceex2,  HTot_vTcPR_Ceex2,  HAfb_vTcPR_Ceex2, nbMax);  //!!!!
  HTot_vTcPR_Ceex2->SetName("HTot_vTcPR_Ceex2");
  HAfb_vTcPR_Ceex2->SetName("HAfb_vTcPR_Ceex2");
  //
  // IFI off
  nbMax=0;   // cosThetaMax = 1.0
  TH1D                    *HTot_vTcPR_Ceex2n, *HAfb_vTcPR_Ceex2n;
  ProjV( sca_vTcPR_Ceex2n,  HTot_vTcPR_Ceex2n,  HAfb_vTcPR_Ceex2n, nbMax);  //!!!!
  HTot_vTcPR_Ceex2n->SetName("HTot_vTcPR_Ceex2n");
  HAfb_vTcPR_Ceex2n->SetName("HAfb_vTcPR_Ceex2n");
  //
  ///****************************************************************************************
  //  dsigma/dv unlimited cos(theta)
  TH1D *Hpro_vT_Ceex2;
  ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2);
  Hpro_vT_Ceex2->SetName("Hpro_vT_Ceex2");

  //  dsigma/dv unlimited cos(theta)
  TH1D *Hpro_vT_Ceex2n;
  ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2n);
  Hpro_vT_Ceex2->SetName("Hpro_vT_Ceex2n");

  cout<<"================ ReMakeMChisto ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto


///////////////////////////////////////////////////////////////////////////////////
void FigVdist()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVdist =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  TH1D *HTot_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2n");  // KKMC sigma(vmax) from scat.
  //
  TH1D *Hpro_vT_Ceex2n    = (TH1D*)DiskFileB.Get("Hpro_vT_Ceex2n");     // KKMC dsigma/dv IFI off, from scat.
  //
  TH1D *HST_xx_Ceex2n     = (TH1D*)DiskFileF.Get("HST_xx_Ceex2n");      // Foam dsigma/d(v) direct
  TH1D *HST_xmax_Ceex2n   = (TH1D*)DiskFileB.Get("HST_xmax_Ceex2n");    // Foam sigma(vmax) direct
  //
  TH1D *Htot_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Htot_xmax_Ceex2n");    //Foam sigma(vmax) scatt.ISR+FSR
  TH1D *Htot_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Htot_xmax_Ceex2");     //Foam sigma(vmax) scatt.ISR+FSR+IFI
  //
//  TH1D *vcum_ISR2_FSR2    = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");     // KKsem  sigma(vmax)
//  TH1D *vdis_ISR2_FSR2    = (TH1D*)DiskFileB.Get("vdis_ISR2_FSR2");     // KKsem  dsigma/d(v)
  //
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVdist = new TCanvas("cFigVdist","FigVdist: NEW", 50, 50,    1000, 800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigVdist->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVdist->Divide( 2,  2);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //====================plot1========================
  //                sigma(vmax)
  cFigVdist->cd(1);
  //gPad->SetLogy(); // !!!!!!
  // MC v-true direct
  TH1D *Hst1 = HTot_vTcPR_Ceex2n;  //  KKMC sigma(vmax) from scat.
  //
  Hst1->SetStats(0);
  Hst1->SetTitle(0);
  //Hst1->SetMinimum(1e-3*Hst1->GetMaximum());
  Hst1->DrawCopy("h");
  //
  HST_xmax_Ceex2n->SetLineColor(kRed);     // red
  HST_xmax_Ceex2n->DrawCopy("hsame");      // Foam sigma(vmax) direct.
  //
  Htot_xmax_Ceex2n->SetLineColor(kBlue);   // blue ISR+FSR
  Htot_xmax_Ceex2n->DrawCopy("hsame");     // Foam sigma(vmax) scatt.
  //
  Htot_xmax_Ceex2->SetLineColor(kGreen);   // green ISR+FSR+IFI
  Htot_xmax_Ceex2->DrawCopy("hsame");      // Foam sigma(vmax) scatt.
  //
  //vcum_ISR2_FSR2->SetLineColor(kBlue);   // blue
  //vcum_ISR2_FSR2->DrawCopy("hsame");     // KKsem sigma(vmax)
  //
  CaptT->DrawLatex(0.02,0.95, "d#sigma/dv(ISR+FSR) Black KKMC_CEEX2, Blue FOAM");
  //====================plot2========================
  cFigVdist->cd(2);
  TH1D *Hst1_ratio =(TH1D*)Hst1->Clone("Hst1_ratio");
  //Hst1_ratio->Divide(vcum_ISR2_FSR2);   // divide by KKsem
  //Hst1_ratio->Divide(HST_xmax_Ceex2n);  // divide by Foam direct
  Hst1_ratio->Divide(Htot_xmax_Ceex2n);   // divide by Foam scatt.
  Hst1_ratio->SetMinimum(0.95);
  Hst1_ratio->SetMaximum(1.10);
  Hst1_ratio->SetLineColor(kBlue);
  Hst1_ratio->DrawCopy("h");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR+FSR); Ratio KKMC/FOAM");
  //====================plot3========================
  //                 dsigma/d(v)
  cFigVdist->cd(3);
  gPad->SetLogy(); // !!!!!!
  TH1D *Hst3 = Hpro_vT_Ceex2n;           // KKMC dsigma/dv from scat.
  Hst3->SetStats(0);
  Hst3->SetTitle(0);
  Hst3->SetLineColor(kRed); // red
  Hst3->DrawCopy("h");
  //
  // KKsem ISR+FSR
  //vdis_ISR2_FSR2->SetLineColor(kMagenta); // magenta
  //vdis_ISR2_FSR2->DrawCopy("hsame");      // KKsem dsigma/d(v) ISR+FSR

  HST_xx_Ceex2n->SetLineColor(kBlue);     // blue
  HST_xx_Ceex2n->DrawCopy("hsame");       // Foam dsigma/d(v)

  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR+FSR),  Red KKMC_CEEX2, Blue FOAM");
  //====================plot4========================
  cFigVdist->cd(4);

  TH1D *Hst3_ratio =(TH1D*)Hst3->Clone("Hst3_ratio");
  //Hst3_ratio->Divide(vdis_ISR2_FSR2);
  Hst3_ratio->Divide(HST_xx_Ceex2n);  // direct

  Hst3_ratio->SetStats(0);
  Hst3_ratio->SetTitle(0);
//  Hst3_ratio->SetMinimum(0.85);
//  Hst3_ratio->SetMaximum(1.15);
//  Hst3_ratio->SetLineColor(kRed);
  Hst3_ratio->DrawCopy("h");  // black

  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR+FSR ); Ratio KKMC/FOAM");
  //----------------------------
  cFigVdist->cd();
  //================================================
}//FigVdist


///////////////////////////////////////////////////////////////////////////////////
void FigAfb()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR

  TH1D *HAfb_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2");  // KKMC
  TH1D *HAfb_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2n"); // KKMC

//  TH1D *afbv_ISR2_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");    // KKsem

  TH1D *Hafb_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb_xmax_Ceex2n");  // FOAM scatt.
  TH1D *Hafb_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb_xmax_Ceex2");   // FOAM scatt.
 //
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigAfb = new TCanvas("cFigAfb","FigAfb: NEW", 20, 300,   1000, 550);
  //                            Name    Title                   xoff,yoff, WidPix,HeiPix
  cFigAfb->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigAfb->Divide( 2,  0);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //====================plot1========================
  //                AFB(vmax)
  cFigAfb->cd(1);
  //gPad->SetLogy(); // !!!!!!
  // MC v-true direct
  TH1D *Hst1 = HAfb_vTcPR_Ceex2n;  //  KKMC AFB(vmax) from scat. IFI off
  TH1D *Hst2 = HAfb_vTcPR_Ceex2;   //  KKMC AFB(vmax) from scat. IFI on
  //
  Hst2->SetStats(0);
  Hst2->SetTitle(0);
//  Hst2->SetMinimum(-0.02);  // 10GeV
//  Hst2->SetMaximum( 0.08);  // 10GeV
  Hst2->SetLineColor(kMagenta);            // magenta
  Hst2->DrawCopy("h");                     // KKMC AFB(vmax) from scat. IFI on
  //
  Hst1->SetLineColor(kBlack);              // black
  Hst1->DrawCopy("hsame");                 // KKMC AFB(vmax) from scat. IFI off
  //
//  afbv_ISR2_FSR2->SetLineColor(kRed);      // red
//  afbv_ISR2_FSR2->DrawCopy("hsame");       // KKsem AFB(vmax) direct. IFI off
  //
  Hafb_xmax_Ceex2n->SetLineColor(kBlue);   // blue FOAM ISR+FSR IFI off
  Hafb_xmax_Ceex2n->DrawCopy("hsame");     // Foam AFB(vmax) scatt.
  //
  Hafb_xmax_Ceex2->SetLineColor(kGreen);   // green FOAM ISR+FSR+IFI, IFI on
  Hafb_xmax_Ceex2->DrawCopy("hsame");      // FOAM ISR+FSR+IFI
  //
  CaptT->DrawLatex(0.02,0.95, "A_{FB}(v_{max}) (ISR+FSR) Black KKMC, Blue Foam, Red KKsem");
  //====================plot2========================
  cFigAfb->cd(2);
//  TH1D *Hst1_diff1 =(TH1D*)Hst1->Clone("Hst1_diff1");
  TH1D *Hst1_diff2 =(TH1D*)Hst1->Clone("Hst1_diff2");
  TH1D *Hst2_diff1 =(TH1D*)Hst2->Clone("Hst2_diff1");
  TH1D *Hst1_diff4 =(TH1D*)Hafb_xmax_Ceex2->Clone("Hst1_diff4");
  TH1D *Hst2_diff2 =(TH1D*)Hst2->Clone("Hst2_diff2");

  Hst2_diff1->Add(Hst2_diff1, Hst1,             1.0, -1.0); // KKMC_IFI   minus KKMC  noIFI black
//  Hst1_diff1->Add(Hst1_diff1, afbv_ISR2_FSR2,   1.0, -1.0); // KKMC_noIFI minus KKsem_noIFI red
  Hst1_diff2->Add(Hst1_diff2, Hafb_xmax_Ceex2n, 1.0, -1.0); // KKMC_noIFI minus FOAM  noIFI blue
  Hst1_diff4->Add(Hst1_diff4, Hafb_xmax_Ceex2n, 1.0, -1.0); // Foam_IFI   minus FOAM  noIFI magenta
  Hst2_diff2->Add(Hst2_diff2, Hafb_xmax_Ceex2,  1.0, -1.0); // KKMC_IFI   minus FOAM  IFI   green

//  Hst2_diff1->SetMinimum(-0.02);  // 189GeV, 10GeV
//  Hst2_diff1->SetMaximum( 0.06);  // 189GeV, 10GeV
//  Hst2_diff1->SetMinimum(-0.002);  // 91GeV
//  Hst2_diff1->SetMaximum( 0.015);  // 91GeV
  Hst2_diff1->SetMinimum(-0.010);  // 95GeV, 88FeV
  Hst2_diff1->SetMaximum( 0.025);  // 95GeV, 88FeV

  Hst2_diff1->SetLineColor(kBlack);
  Hst2_diff1->DrawCopy("h");

//  Hst1_diff1->SetLineColor(kRed);
//  Hst1_diff1->DrawCopy("hsame");

  Hst1_diff2->SetLineColor(kBlue);
  Hst1_diff2->DrawCopy("hsame");

  Hst1_diff4->SetLineColor(kMagenta);
  Hst1_diff4->DrawCopy("hsame");

  Hst2_diff2->SetLineColor(kGreen);
  Hst2_diff2->DrawCopy("hsame");
  //
  //CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR+FSR); Ratio KKMC/FOAM");
  CaptT->DrawLatex(0.12,0.85,"A^{KKMC}_{FB}-A^{Foam}_{FB}, (ISR+FSR) ");

  cFigAfb->cd();
  //================================================
}//FigAfb




///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  cout<< " FigInfo:  CMSene="<<CMSene<<endl;

  TH1D *hst_weight3  = (TH1D*)DiskFileF.Get("hst_weight3"); // Foam3
  TH1D *hst_weight5  = (TH1D*)DiskFileF.Get("hst_weight5"); // Foam5
 //
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo ", 140, 500,   1000, 550);
  //                            Name    Title                     xoff,yoff, WidPix,HeiPix
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  0);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //====================plot1========================
  //                sigma(vmax)
  cFigInfo->cd(1);
  //gPad->SetLogy(); // !!!!!!
  // MC v-true direct
  TH1D *Hst1 = hst_weight3;  //  weight of Foam
  //
  //Hst1->SetStats(0);
  Hst1->SetTitle(0);
  Hst1->DrawCopy("h");

  CaptT->DrawLatex(0.02,0.95, " Weight distribution Foam");
  //====================plot2========================
  cFigInfo->cd(2);

  Hst1 = hst_weight5;
  Hst1->SetTitle(0);
  Hst1->DrawCopy("h");

  cFigInfo->cd();
  //================================================
}//FigInfo


///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  //LibSem.Initialize(DiskFileA);

  DiskFileB.cd();
  HistNormalize();     // Renormalization of MC histograms
//  KKsemMakeHisto();    // prepare histos from KKsem
  ReMakeMChisto();     // reprocessing MC histos from KKC and Foam
//========== PLOTTING ==========
  FigVdist();  // sigma(v) and sigma(vmax) KKMC/Foam
  FigAfb();    // AFB(vmax) KKMC/Foam
  FigInfo();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();

  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


