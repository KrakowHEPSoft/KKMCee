//////////////////////////////////////////////////////////////////////
//    make PlotFoam0-run
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
#include "KKplot.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
////  *** KKMC
//TFile DiskFileA("../workKKMC/histo.root");
//
TFile DiskFileA("../workKKMC/histo.root_95GeV_2.5G"); // obsolete
//TFile DiskFileA("../workKKMC/histo.root_95GeV_1200M"); // obsolete
//TString XparFile="../workKKMC/workKKMC_95GeV.input";

//TFile DiskFileA("../workAFB/rmain.root");
//TFile DiskFileA("../workAFB/rmain_95GeV.root"); // 100M new
//TFile DiskFileA("../workAFB/rmain.root_189GeV_100M"); // obsolete

////  *** FOAM
//TFile DiskFileF("../workFOAM/histo.root"); // current
TFile DiskFileF("../workFOAM/histo.root_95GeV_36G");
//TFile DiskFileF("../workFOAM/histo_95GeV_241M.root");

//  *** older FOAM
//TFile DiskFileF("../workFoam0/rmain.root");
//TFile DiskFileF("../workFoam0rmain_95GeV_64M.root");
//
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");

// Interface to KKabox and some extra plotting facilities
//KKabox LibSem;

///////////////////////////////////////////////////////////////////////////////////
KKplot LibSem("KKplot");

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

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n") );
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
  TH2D *SCA_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCA_xc_Ceex2n");  // FOAM big range
  TH2D *SCA_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCA_xc_Ceex2");   // FOAM big range
  TH2D *SCT_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2n");  // FOAM small range
  TH2D *SCT_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2");   // FOAM small range

  HisNorm1(HST_FOAM_NORMA3, HST_xx_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA3, SCA_xc_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA5, SCA_xc_Ceex2 );   // normalizing
  HisNorm2(HST_FOAM_NORMA3, SCT_xc_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA5, SCT_xc_Ceex2 );   // normalizing

  // sigma(vmax) direct histogramming
  TH1D *HST_xmax_Ceex2n;
  MakeCumul(HST_xx_Ceex2n,HST_xmax_Ceex2n);
  HST_xmax_Ceex2n->SetName("HST_xmax_Ceex2n");

  // sigma(vmax) and AFB(vmax) from scattergram vmax<1
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

  // sigma(vmax) and AFB(vmax) from scattergram vmax<0.2
  nbMax=0;   // <--- CosThetaMax = 1.0
  TH1D                 *Htot2_xmax_Ceex2n, *Hafb2_xmax_Ceex2n;
  ProjV( SCT_xc_Ceex2n, Htot2_xmax_Ceex2n,  Hafb2_xmax_Ceex2n, nbMax);  //!!!!
  Htot2_xmax_Ceex2n->SetName("Htot2_xmax_Ceex2n");
  Hafb2_xmax_Ceex2n->SetName("Hafb2_xmax_Ceex2n");
  //
  TH1D                *Htot2_xmax_Ceex2, *Hafb2_xmax_Ceex2;
  ProjV( SCT_xc_Ceex2, Htot2_xmax_Ceex2,  Hafb2_xmax_Ceex2, nbMax);  //!!!!
  Htot2_xmax_Ceex2->SetName("Htot2_xmax_Ceex2");
  Hafb2_xmax_Ceex2->SetName("Hafb2_xmax_Ceex2");

  ////////////////////////////////////////////////////////////////////
//********************************************************************
// Pure KKMC reprocessing part
  cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;

  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene        /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  // Wide range, vmax<1.
  TH2D *sca_vTcPR_Eex2   = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  TH2D *sca_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n");
  ///****************************************************************************************
  /// Distributions of v=vTrue with unlimited c=cos(theta)
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
  // IFI off
  nbMax=0;   // cosThetaMax = 1.0
  TH1D                    *HTot_vTcPR_Ceex2n, *HAfb_vTcPR_Ceex2n;
  ProjV( sca_vTcPR_Ceex2n, HTot_vTcPR_Ceex2n,  HAfb_vTcPR_Ceex2n, nbMax);  //!!!!
  HTot_vTcPR_Ceex2n->SetName("HTot_vTcPR_Ceex2n");
  HAfb_vTcPR_Ceex2n->SetName("HAfb_vTcPR_Ceex2n");
  ///****************************************************************************************
  //  dsigma/dv unlimited cos(theta)
  TH1D *Hpro_vT_Ceex2;
  ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2);
  Hpro_vT_Ceex2->SetName("Hpro_vT_Ceex2");
  //  dsigma/dv unlimited cos(theta)
  TH1D *Hpro_vT_Ceex2n;
  ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2n);
  Hpro_vT_Ceex2n->SetName("Hpro_vT_Ceex2n");

  cout<<"================ ReMakeMChisto ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto


///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto2(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto2  BEGIN ============================"<<endl;
////////////////////////////////////////////////////////////////////
// Pure KKMC reprocessing part
// from bigger scattergram and restricted vmax<0.2
  //////////////////////////////////////////////////////////////////
    cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;

    TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
    cout<<"ReMakeMChisto2 [1]"<<endl;
   // Wide range, vmax<1.
    TH2D *sct_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2");
    TH2D *sct_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n");
    cout<<"ReMakeMChisto2 [2]"<<endl;
    ///****************************************************************************************
    ///****************************************************************************************
    /// Distributions of v=vTrue with unlimited c=cos(theta)
    int nbMax=0;   // cosThetaMax = 1.0
    TH1D                    *HTot2_vTcPR_Ceex2, *HAfb2_vTcPR_Ceex2;
    ProjV( sct_vTcPR_Ceex2,  HTot2_vTcPR_Ceex2,  HAfb2_vTcPR_Ceex2, nbMax);  //!!!!
    HTot2_vTcPR_Ceex2->SetName("HTot2_vTcPR_Ceex2");
    HAfb2_vTcPR_Ceex2->SetName("HAfb2_vTcPR_Ceex2");
    // IFI off
    nbMax=0;   // cosThetaMax = 1.0
    TH1D                    *HTot2_vTcPR_Ceex2n, *HAfb2_vTcPR_Ceex2n;
    ProjV( sct_vTcPR_Ceex2n, HTot2_vTcPR_Ceex2n,  HAfb2_vTcPR_Ceex2n, nbMax);  //!!!!
    HTot2_vTcPR_Ceex2n->SetName("HTot2_vTcPR_Ceex2n");
    HAfb2_vTcPR_Ceex2n->SetName("HAfb2_vTcPR_Ceex2n");
    ///****************************************************************************************

  cout<<"================ ReMakeMChisto2 ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto2



///////////////////////////////////////////////////////////////////////////////////
void KKsemMakeHisto(){
  // Here we produce semianalytical plots using KKsem program, No plotting
  //------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem MakeHisto  BEGIN ============================"<<endl;
  //
  long KF=13; // muon
  long KeyDis, KeyFob;
  char chak[5];
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
//------------------------------------------------------------------------
  TH1D *hstVtemplate = (TH1D*)DiskFileB.Get("HTot2_vTcPR_Ceex2");
//------------------------------------------------------------------------
//   MuMu  Sigma(vmax) with ulimited c=cos(theta)
//------------------------------------------------------------------------
// ISR*FSR
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  LibSem.VVplot(vcum_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
//-------------------------------------------------
//    AFB(vmax) for unlimited c=cos(theta)
//-------------------------------------------------
  kksem_setcrange_(0, 25.0/25); // forward cos(theta)
  TH1D *afbv_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("afbv_ISR2_FSR2");
  LibSem.VVplot(afbv_ISR2_FSR2, KF, chak, KeyDis, KeyFob);// Forward
  afbv_ISR2_FSR2->Add(afbv_ISR2_FSR2, vcum_ISR2_FSR2, 2.0, -1.0) ; // numerator F-B = 2F-(F+B)
  afbv_ISR2_FSR2->Divide(vcum_ISR2_FSR2);                          // finally (F-B)(F+B)
  kksem_setcrange_(-1.0, 1.0); // undoing forward,
  //

  cout<<"================ KKsem MakeHisto ENDs ============================="<<endl;
  cout<<"==================================================================="<<endl;
//------------------------------------------------------------------------
//------------------------------------------------------------------------
}//  KKsemMakeHisto


///////////////////////////////////////////////////////////////////////////////////
void FigVdist()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVdist =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
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
  TCanvas *cFigVdist = new TCanvas("cFigVdist","FigVdist", 50, 50,    1000, 800);
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
  //Hst1_ratio->SetMinimum(0.95);
  //Hst1_ratio->SetMaximum(1.10);
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
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);


  TH1D *HAfb_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2");  // KKMC
  TH1D *HAfb_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2n"); // KKMC

//  TH1D *afbv_ISR2_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");    // KKsem

  TH1D *Hafb_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb_xmax_Ceex2n");  // FOAM scatt.
  TH1D *Hafb_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb_xmax_Ceex2");   // FOAM scatt.
 //
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigAfb = new TCanvas("cFigAfb","FigAfb", 20, 300,   1000, 550);
  //                            Name    Title                   xoff,yoff, WidPix,HeiPix
  cFigAfb->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigAfb->Divide( 2,  0);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
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
  Hafb_xmax_Ceex2->SetLineColor(kGreen);   // green FOAM IFI on
  Hafb_xmax_Ceex2->DrawCopy("hsame");      // FOAM ISR+FSR+IFI
  //
  Hafb_xmax_Ceex2n->SetLineColor(kBlue);   // blue FOAM IFI off
  Hafb_xmax_Ceex2n->DrawCopy("hsame");     // Foam AFB(vmax) scatt.
  //
  CaptT->DrawLatex(0.12,0.95, "A_{FB}^{IF on}(v_{max}):  KKMC=magenta, Foam=green");
  CaptT->DrawLatex(0.12,0.85, "A_{FB}^{IFI off}(v_{max}): KKMC=black,  Foam=blue");
  CaptT->DrawLatex(0.60,0.75,TextEne);
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
//  Hst2_diff1->SetMinimum(-0.002);  // zoom
//  Hst2_diff1->SetMaximum( 0.002);  // zoom

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
  CaptT->DrawLatex(0.12,0.95,"A_{FB}^{IFI}(v_{max}): Black KKMC, Magenta=FOAM");
  CaptT->DrawLatex(0.12,0.85,"A^{KKMC}_{FB}-A^{FOAM}: Green=IFI, Blue=NOIFI");

  cFigAfb->cd();
  //================================================
}//FigAfb



///////////////////////////////////////////////////////////////////////////////////
void FigAfb2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb2 =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);

  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n");  //
  TH1D *HAfb2_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2");  //

  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // FOAM scatt.
  TH1D *Hafb2_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2");   // FOAM scatt.

  TH1D *afbv_ISR2_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");    // KKsem
//
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigAfb2 = new TCanvas("cFigAfb2","FigAfb2", 70, 350,   1000, 550);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  cFigAfb2->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigAfb2->Divide( 2,  0);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  //====================plot1========================
  //                AFB(vmax)
  cFigAfb2->cd(1);
  //gPad->SetLogy(); // !!!!!!
  TH1D *Hst1 = HAfb2_vTcPR_Ceex2n;         // KKMC AFB(vmax) from scat. IFI off
  TH1D *Hst2 = HAfb2_vTcPR_Ceex2;          // KKMC AFB(vmax) from scat. IFI on
  //
  Hst2->SetStats(0);
  Hst2->SetTitle(0);
  Hst2->SetLineColor(kMagenta);            // magenta
  Hst2->DrawCopy("h");                     // KKMC AFB(vmax) from scat. IFI on
  //
  Hst1->SetLineColor(kBlack);              // black
  Hst1->DrawCopy("hsame");                 // KKMC AFB(vmax) from scat. IFI off
  //
  Hafb2_xmax_Ceex2n->SetLineColor(kBlue);
  Hafb2_xmax_Ceex2n->DrawCopy("hsame");   // Foam IFI OFF
  Hafb2_xmax_Ceex2->SetLineColor(kGreen);
  Hafb2_xmax_Ceex2->DrawCopy("hsame");    // Foam IFI ON

  afbv_ISR2_FSR2->SetLineColor(kRed);     // KKsem  RED!!!
  afbv_ISR2_FSR2->DrawCopy("hsame");

  //
  CaptT->DrawLatex(0.12,0.95, "A_{FB}^{IFI on}(v_{max})  KKMC=magenta, Foam=green");
  CaptT->DrawLatex(0.12,0.85, "A_{FB}^{IFI off}(v_{max}) KKMC=black,  Foam=blue");

  CaptT->DrawLatex(0.60,0.75,TextEne);
  //====================plot2========================
  cFigAfb2->cd(2);

  TH1D *Hst21_diff =(TH1D*)Hst2->Clone("Hst21_diff");
  Hst21_diff->Add(Hst21_diff, Hst1,  1.0, -1.0); // KKMC_IFI

  TH1D *HST21_diff =(TH1D*)Hafb2_xmax_Ceex2->Clone("HST21_diff");
  HST21_diff->Add(HST21_diff, Hafb2_xmax_Ceex2n,  1.0, -1.0); // FOAMC_IFI

  Hst21_diff->SetMinimum(-0.004);  // zoom
  Hst21_diff->SetMaximum( 0.004);  // zoom

  Hst21_diff->SetLineColor(kBlack);              // blue, KKMC
  Hst21_diff->DrawCopy("h");
  HST21_diff->SetLineColor(kMagenta);               // red, FOAM
  HST21_diff->DrawCopy("hsame");

  TH1D *HstKF_diff =(TH1D*)Hst2->Clone("HstKF_diff");
  HstKF_diff->Add(HstKF_diff, Hafb2_xmax_Ceex2,  1.0, -1.0); // KKMC-Foam IFIon

  TH1D *HstKFn_diff =(TH1D*)Hst1->Clone("HstKFn_diff");
  HstKFn_diff->Add(HstKFn_diff, Hafb2_xmax_Ceex2n,  1.0, -1.0); // KKMC-Foam IFIoff

  HstKF_diff->SetLineColor(kGreen);               // KKMC-Foam IFIon
  HstKF_diff->DrawCopy("hsame");

  HstKFn_diff->SetLineColor(kBlue);              // KKMC-Foam IFIon
  HstKFn_diff->DrawCopy("hsame");

// zero line
  TH1D *hZero = (TH1D*)Hst1->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  hZero->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"A_{FB}^{IFI}(v_{max}): Black KKMC, Magenta=FOAM");
  CaptT->DrawLatex(0.12,0.85,"A^{KKMC}_{FB}-A^{FOAM}: Green=IFI, Blue=NOIFI");

  cFigAfb2->cd();
  //================================================
}//FigAfb2




///////////////////////////////////////////////////////////////////////////////////
void FigTech()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTech =========================== "<<endl;

// AFB(vmax)
  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // AFB FOAM scatt.
  TH1D *afbv_ISR2_FSR2     = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");     // AFB KKsem
  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n"); // KKMC AFB(vmax) from scat
// sigma(vmax)
  TH1D *Htot2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Htot2_xmax_Ceex2n");  // FOAM scatt.
  TH1D *vcum_ISR2_FSR2     = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");     // KKsem
  TH1D *HTot2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HTot2_vTcPR_Ceex2n"); // KKMC sigma(vmax) from scat.

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigTech = new TCanvas("cFigTech","FigTech", 100, 400,   1000, 550);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  cFigTech->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigTech->Divide( 2,  0);
  //====================plot1========================
  //                AFB(vmax)
  cFigTech->cd(1);

  TH1D *HstTech_ratio =(TH1D*)Htot2_xmax_Ceex2n->Clone("HstTech_ratio");
  HstTech_ratio->Divide(HstTech_ratio, vcum_ISR2_FSR2); // sigma(xmax), KKsem/Foam IFIoff

  TH1D *HstTech_ratio2 =(TH1D*)HTot2_vTcPR_Ceex2n->Clone("HstTech_ratio2");
  HstTech_ratio2->Divide(HstTech_ratio2, vcum_ISR2_FSR2); // sigma(xmax), KKMC/Foam IFIoff

  HstTech_ratio->SetStats(0);
  HstTech_ratio->SetTitle(0);


  HstTech_ratio->SetMinimum(1 -0.0005);  // zoom
  HstTech_ratio->SetMaximum(1 +0.0005);  // zoom
  HstTech_ratio->SetLineColor(kGreen);
  HstTech_ratio->DrawCopy("h");
  HstTech_ratio2->SetLineColor(kBlue);
  HstTech_ratio2->DrawCopy("hsame");


  CaptT->DrawLatex(0.12,0.95,"#sigma^{IFIoff}(v_{max}): FOAM/KKsem(green), KKMC/KKsem");

  //====================plot2========================
  cFigTech->cd(2);

  TH1D *HstTech_diff =(TH1D*)Hafb2_xmax_Ceex2n->Clone("HstTech_diff");
  HstTech_diff->Add(HstTech_diff, afbv_ISR2_FSR2,  1.0, -1.0); // AFB, KKMC-Foam IFIoff

  TH1D *HstTech_diff2 =(TH1D*)HAfb2_vTcPR_Ceex2n->Clone("HstTech_diff2");
  HstTech_diff2->Add(HstTech_diff2, afbv_ISR2_FSR2,  1.0, -1.0); // AFB, KKMC-Foam IFIoff

  HstTech_diff->SetStats(0);
  HstTech_diff->SetTitle(0);
  HstTech_diff->SetMinimum(-0.0005);  // zoom
  HstTech_diff->SetMaximum( 0.0005);  // zoom
  HstTech_diff->SetLineColor(kGreen);
  HstTech_diff->DrawCopy("h");
  HstTech_diff2->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"A_{FB}^{IFIoff}(v_{max}): FOAM-KKsem(green), KKMC-KKsem");

  cFigTech->cd();
  //================================================

}//FigTech


///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
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
  /////////////////////////////////////////////////////////
  LibSem.Initialize(DiskFileA);  // for non-farm case
  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  /*
  double xpar[10001];
  int jmax = 10000;
  LibSem.ReaData("../../.KK2f_defaults",     jmax, xpar);  // numbering as in input!!!
  char dname[100];  sprintf(dname,XparFile);
  LibSem.ReaData(dname, -jmax, xpar);  // jmax<0 for append mode
  double CMSene  = xpar[ 1];
  cout<< "////// Main: CMSene = "<<  CMSene  <<endl;
  LibSem.Initialize(xpar);  // for non-farm case
  */
  //////////////////////////////////////////////////////////////////////////
  DiskFileB.cd();
  HistNormalize();     // Renormalization of MC histograms
  ReMakeMChisto();     // reprocessing MC histos from KKC and Foam
  ReMakeMChisto2();    // reprocessing MC histos from KKC and Foam
  KKsemMakeHisto();    // prepare histos from KKsem
//========== PLOTTING ==========
// vmax=1
  FigVdist();  // sigma(v) and sigma(vmax) KKMC/Foam
  FigAfb();    // AFB(vmax) KKMC/Foam
// vmax =0.2
  FigAfb2();
  FigTech();
// weight distribution
  FigInfo();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //cout<<"------------------------------A.ls----------------------------------"<<endl;
  //DiskFileA.ls();
  //cout<<"------------------------------F.ls----------------------------------"<<endl;
  //DiskFileF.ls();
  //cout<<"------------------------A.GetListOfKeys-----------------------------"<<endl;
  //DiskFileA.GetListOfKeys()->Print();
  cout<<"------------------------F.GetListOfKeys-----------------------------"<<endl;
  DiskFileF.GetListOfKeys()->Print();
  //cout<<"------------------------------end---------------------------------"<<endl;
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


