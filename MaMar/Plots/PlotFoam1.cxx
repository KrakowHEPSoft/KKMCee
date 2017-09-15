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
// August2017 runs
//TFile DiskFileA("../workKKMC/histo.root_10GeV_1G"); //
//TFile DiskFileA("../workKKMC/histo.root_88GeV_2.1G"); //
TFile DiskFileA("../workKKMC/histo.root_95GeV_16G");
//TFile DiskFileA("../workKKMC/histo.root_91GeV_9G"); ///????

////  *** FOAM
TFile DiskFileF("../workFOAM1/histo.root"); // current

TFile DiskFileB("RhoSemi.root","RECREATE","histograms");

//////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot, gNevTot2; // from KKMC and KKfoam MC runs (histograms)
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int  kGold=92, kPine=46;
//
//int    gNbMax=50;         // gCosTheta = 45/50=0.90
//double gCosTheta=1.00;    // to be synchronized with gNbMax
//
KKplot LibSem("KKplot");
///////////////////////////////////////////////////////////////////////////////////



void HistNormKKMC(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  //DiskFileA.ls("");
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
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPL_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPL_Eex2") );
  //
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_EEX2") );
  //
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2n") );

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2") );
  //
}

///////////////////////////////////////////////////////////////////////////////////
void ReMakeKKMC(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeKKMC  BEGIN ============================"<<endl;
////////////////////////////////////////////////////////////////////
// Pure KKMC reprocessing part
// from bigger scattergram and restricted vmax<0.2
//////////////////////////////////////////////////////////////////
    cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;

    TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

    ///****************************************************************************************
    // Pure KKMC reprocessing part
    int nbMax;
    cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;
    // Wide range, vmax<1.
    TH2D *sca_vTcPR_Eex2   = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
    TH2D *sca_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
    TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n");

    //****************************************************************************************
    //  dsigma/dv unlimited cos(theta)
    //****************************************************************************************
    TH1D *Hpro_vT_Ceex2;
    ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2);
    Hpro_vT_Ceex2->SetName("Hpro_vT_Ceex2");
    //  dsigma/dv unlimited cos(theta)
    TH1D *Hpro_vT_Ceex2n;
    ProjX1(sca_vTcPR_Ceex2n, Hpro_vT_Ceex2n);
    Hpro_vT_Ceex2n->SetName("Hpro_vT_Ceex2n");

    ///****************************************************************************************
    // More Wide range, vmax<1.
    TH2D *sca_vTcPL_Eex2   = (TH2D*)DiskFileA.Get("sca_vTcPL_Eex2");
    TH2D *sca_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sca_vTcPL_Ceex2");
    TH2D *sca_vTcPL_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPL_Ceex2n");
    //
    TH1D                    *HTot_vTcPL_Ceex2, *HAfb_vTcPL_Ceex2;
    ProjV( sca_vTcPL_Ceex2,  HTot_vTcPL_Ceex2,  HAfb_vTcPL_Ceex2, nbMax);  //!!!!
    HTot_vTcPL_Ceex2->SetName("HTot_vTcPL_Ceex2");
    HAfb_vTcPL_Ceex2->SetName("HAfb_vTcPL_Ceex2");
    //
    TH1D                    *HTot_vTcPL_Ceex2n, *HAfb_vTcPL_Ceex2n;
    ProjV( sca_vTcPL_Ceex2n,  HTot_vTcPL_Ceex2n,  HAfb_vTcPL_Ceex2n, nbMax);  //!!!!
    HTot_vTcPL_Ceex2n->SetName("HTot_vTcPL_Ceex2n");
    HAfb_vTcPL_Ceex2n->SetName("HAfb_vTcPL_Ceex2n");


  cout<<"================ ReMakeKKMC ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//ReMakeKKMC



///////////////////////////////////////////////////////////////////////////////////
void ReMakeFoam1(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeFoam1  BEGIN    ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA1 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA1");

  TH1D *HST_xx_Ord1  = (TH1D*)DiskFileF.Get("HST_xx_Ord1");  // FOAM
  TH1D *HST_xx_Crd1  = (TH1D*)DiskFileF.Get("HST_xx_Crd1");  // FOAM
  TH1D *HST_xx_Ord1n = (TH1D*)DiskFileF.Get("HST_xx_Ord1n");  // FOAM
  TH1D *HST_xx_Crd1n = (TH1D*)DiskFileF.Get("HST_xx_Crd1n");  // FOAM
  TH1D *HST_xx_Hrd1  = (TH1D*)DiskFileF.Get("HST_xx_Hrd1");  // FOAM
  TH1D *HST_xx_Srd1  = (TH1D*)DiskFileF.Get("HST_xx_Srd1");  // FOAM
  TH1D *HST_xx_Ird1  = (TH1D*)DiskFileF.Get("HST_xx_Ird1");  // FOAM

  HisNorm1(HST_FOAM_NORMA1, HST_xx_Ord1 );   // normalizing
  HisNorm1(HST_FOAM_NORMA1, HST_xx_Crd1 );   // normalizing
  HisNorm1(HST_FOAM_NORMA1, HST_xx_Ord1n );  // normalizing
  HisNorm1(HST_FOAM_NORMA1, HST_xx_Crd1n );  // normalizing
  HisNorm1(HST_FOAM_NORMA1, HST_xx_Hrd1 );   // normalizing
  HisNorm1(HST_FOAM_NORMA1, HST_xx_Srd1 );   // normalizing
  HisNorm1(HST_FOAM_NORMA1, HST_xx_Ird1 );   // normalizing

  //---------------------------------------------------
  // sigma(vmax) direct histogramming, (0,1) range
  TH1D *HST_xxcum_Ord1;
  MakeCumul(HST_xx_Ord1,HST_xxcum_Ord1);
  HST_xxcum_Ord1->SetName("HST_xxcum_Ord1");

  TH1D *HST_xxcum_Ord1n;
  MakeCumul(HST_xx_Ord1n,HST_xxcum_Ord1n);
  HST_xxcum_Ord1n->SetName("HST_xxcum_Ord1n");

  //---------------------------------------------------
  // sigma^*(vmax) direct histogramming, (0,1) range
  TH1D *HST_xxcum_Crd1;
  MakeCumul(HST_xx_Crd1,HST_xxcum_Crd1);
  HST_xxcum_Crd1->SetName("HST_xxcum_Crd1");

  TH1D *HST_xxcum_Crd1n;
  MakeCumul(HST_xx_Crd1n,HST_xxcum_Crd1n);
  HST_xxcum_Crd1n->SetName("HST_xxcum_Crd1n");

  TH1D *HST_xxcum_Hrd1;
  MakeCumul(HST_xx_Hrd1,HST_xxcum_Hrd1);
  HST_xxcum_Hrd1->SetName("HST_xxcum_Hrd1");

  TH1D *HST_xxcum_Srd1;
  MakeCumul(HST_xx_Srd1,HST_xxcum_Srd1);
  HST_xxcum_Srd1->SetName("HST_xxcum_Srd1");

  TH1D *HST_xxcum_Ird1;
  MakeCumul(HST_xx_Ird1,HST_xxcum_Ird1);
  HST_xxcum_Ird1->SetName("HST_xxcum_Ird1");

  //---------------------------------------------------
  //  and finally AFB from ratio sigma^*(vmax)/sigma(vmax)
  TH1D *HST_xxAfb_Ord1 = (TH1D*)HST_xxcum_Crd1->Clone("HST_xxAfb_Ord1");
  HST_xxAfb_Ord1->Divide(HST_xxcum_Ord1n);   // <-- divided by NO IFI !!!!
  HST_xxAfb_Ord1->Scale(3.0/4.0);

  TH1D *HST_xxAfb_Ord1n = (TH1D*)HST_xxcum_Crd1n->Clone("HST_xxAfb_Ord1n");
  HST_xxAfb_Ord1n->Divide(HST_xxcum_Ord1n);
  HST_xxAfb_Ord1n->Scale(3.0/4.0);

  TH1D *HST_xxAfb_Hrd1 = (TH1D*)HST_xxcum_Hrd1->Clone("HST_xxAfb_Hrd1");
  HST_xxAfb_Hrd1->Divide(HST_xxcum_Ord1n);
  HST_xxAfb_Hrd1->Scale(3.0/4.0);

  TH1D *HST_xxAfb_Srd1 = (TH1D*)HST_xxcum_Srd1->Clone("HST_xxAfb_Srd1");
  HST_xxAfb_Srd1->Divide(HST_xxcum_Ord1n);
  HST_xxAfb_Srd1->Scale(3.0/4.0);

  TH1D *HST_xxAfb_Ird1 = (TH1D*)HST_xxcum_Ird1->Clone("HST_xxAfb_Ird1");
  HST_xxAfb_Ird1->Divide(HST_xxcum_Ord1n);
  HST_xxAfb_Ird1->Scale(3.0/4.0);

  cout<<"================ ReMakeFoam1 ENDs     ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//ReMakeFoam1


///////////////////////////////////////////////////////////////////////////////////
void KKsemMakeHisto(){
  // Here we produce semianalytical plots using KKsem program, No plotting
  //------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem MakeHisto  BEGIN ============================"<<endl;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  //
  int KF=13; // muon
  int KeyDis, KeyFob;
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
  KeyFob=  -10; // KKsem_BornV, NO EW, NO integration OK!
  KeyFob= -100; // KKsem_BornV, NO EW, WITH integration, OK
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!

//  KeyFob=  -10; // KKsem_BornV, NO EW, NO integration OK!
//------------------------------------------------------------------------
  TH1D *hstVtemplate = (TH1D*)DiskFileA.Get("hst_vTrueCeex2");
  TH1D *hstCtemplate = (TH1D*)DiskFileA.Get("hst_Cost1Ceex2");
//------------------------------------------------------------------------
  cout<<"  MuMu  Sigma(vmax) with ulimited c=cos(theta) "<<endl;
//------------------------------------------------------------------------
// ISR*FSR
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  LibSem.VVplot(vcum_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
//-------------------------------------------------
  cout<<"  AFB(vmax) for unlimited c=cos(theta) "<<endl;
//-------------------------------------------------
  kksem_setcrange_(0, 25.0/25); // forward cos(theta)
  TH1D *afbv_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("afbv_ISR2_FSR2");
  LibSem.VVplot(afbv_ISR2_FSR2, KF, chak, KeyDis, KeyFob);// Forward
  afbv_ISR2_FSR2->Add(afbv_ISR2_FSR2, vcum_ISR2_FSR2, 2.0, -1.0) ; // numerator F-B = 2F-(F+B)
  afbv_ISR2_FSR2->Divide(vcum_ISR2_FSR2);                          // finally (F-B)(F+B)
  kksem_setcrange_(-1.0, 1.0); // undoing forward,
  //
  //------------------------------------------------------------------------
  cout<<"  MuMu  dsigma/dv, unlimited cos(theta)"<<endl;
  //------------------------------------------------------------------------
  // ISR*FSR
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XRHO2");  // ISR*FSR Mff
  TH1D *vdis_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vdis_ISR2_FSR2");
  LibSem.VVplot(vdis_ISR2_FSR2, KF, chak, KeyDis, KeyFob);


  cout<<"================ KKsem MakeHisto ENDs ============================="<<endl;
  cout<<"==================================================================="<<endl;
//------------------------------------------------------------------------
//------------------------------------------------------------------------
}//  KKsemMakeHisto



///////////////////////////////////////////////////////////////////////////////////
void FigVV()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVV =========================== "<<endl;
 //
  TH1D *Hpro_vT_Ceex2n    = (TH1D*)DiskFileB.Get("Hpro_vT_Ceex2n");     // KKMC dsigma/dv IFI off, from scat.
  //HST_xx_Ord1
  TH1D *HST_xx_Ord1       = (TH1D*)DiskFileF.Get("HST_xx_Ord1");        // Foam1
  TH1D *HST_xx_Ord1n      = (TH1D*)DiskFileF.Get("HST_xx_Ord1n");        // Foam1
  //
  TH1D *vdis_ISR2_FSR2    = (TH1D*)DiskFileB.Get("vdis_ISR2_FSR2");     // KKsem  dsigma/d(v)

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVV = new TCanvas("cFigVV","FigVV", 50, 50,    1000, 550);
  //                            Name    Title   xoff,yoff, WidPix,HeiPix
  cFigVV->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVV->Divide( 2,  1);
  //====================plot1========================
  //                 dsigma/d(v)
  cFigVV->cd(1);
  gPad->SetLogy(); // !!!!!!
  TH1D *Hst1 = HST_xx_Ord1n;
  Hst1->SetStats(0);
  Hst1->SetTitle(0);
  Hst1->DrawCopy("h");

  Hpro_vT_Ceex2n->SetLineColor(kRed);   // red
  Hpro_vT_Ceex2n->DrawCopy("hsame");    // KKMC dsigma/dv from scat.
  //
  vdis_ISR2_FSR2->SetLineColor(kGreen);  //
  vdis_ISR2_FSR2->DrawCopy("hsame");     // KKsem
  //
  HST_xx_Ord1n->SetLineColor(kBlue);
  HST_xx_Ord1n->DrawCopy("hsame");        // TMCgenFoam1, IFI off
  //
  HST_xx_Ord1->SetLineColor(kYellow);
  HST_xx_Ord1->DrawCopy("hsame");         // TMCgenFoam1, IFI on

  //====================plot2========================
  cFigVV->cd(2);

  TH1D *Hst1_ratio =(TH1D*)HST_xx_Ord1n->Clone("Hst1_ratio");
  Hst1_ratio->Divide(Hpro_vT_Ceex2n);

  Hst1_ratio->SetStats(0);
  Hst1_ratio->SetTitle(0);
  Hst1_ratio->SetMinimum(0);
  Hst1_ratio->SetMaximum(2);
  Hst1_ratio->DrawCopy("h");
  cout<<" ========================= FigVV end======================== "<<endl;
  //================================================
}//FigVV



///////////////////////////////////////////////////////////////////////////////////
void FigVcum()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVcum =========================== "<<endl;
 //
  TH1D *HTot_vTcPL_Ceex2n = (TH1D*)DiskFileB.Get("HTot_vTcPL_Ceex2n");  // KKMC sigma(vmax) from scat.
  TH1D *HTot_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HTot_vTcPL_Ceex2");  // KKMC sigma(vmax) from scat.

  TH1D *vcum_ISR2_FSR2    = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");     // KKsem  sigma(vmax)

  TH1D *HST_xxcum_Ord1    = (TH1D*)DiskFileB.Get("HST_xxcum_Ord1");      // Foam1
  TH1D *HST_xxcum_Ord1n   = (TH1D*)DiskFileB.Get("HST_xxcum_Ord1n");      // Foam1

  cout<<"FigVcum calculating PRD xsection "<<endl;

  TH1D *HST_SigPRD =(TH1D*)vcum_ISR2_FSR2->Clone("HST_SigPRD");
  LibSem.Ord1fill(HST_SigPRD,301);
  HST_SigPRD->SetLineColor(kCyan);

  cout<<"FigVcum calculating PRD xsection END "<<endl;


  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVcum = new TCanvas("cFigVcum","FigVcum", 100, 100,    1000, 550);
  //                            Name    Title   xoff,yoff, WidPix,HeiPix
  cFigVcum->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVcum->Divide( 2,  1);

  //====================plot1========================
  //                 sigma(vmax)
    cFigVcum->cd(1);

    TH1D *Hst3 = HTot_vTcPL_Ceex2n;
    Hst3->SetStats(0);
    Hst3->SetTitle(0);

    Hst3->SetMinimum(0.32);     // 95GeV
    Hst3->SetMaximum(0.41);    // 95GeV

    Hst3->DrawCopy("h");                     // TO BE REPEATED BELOW

    HTot_vTcPL_Ceex2n->SetLineColor(kRed);   // red
    HTot_vTcPL_Ceex2n->DrawCopy("hsame");    // KKMC  IFI OFF

    HTot_vTcPL_Ceex2->SetLineColor(kBlack);  // black
    HTot_vTcPL_Ceex2->DrawCopy("hsame");     // KKMC  IFI on

    vcum_ISR2_FSR2->SetLineColor(kGreen);    // Green
    vcum_ISR2_FSR2->DrawCopy("hsame");       // KKsem IFI OFF

    HST_xxcum_Ord1n->SetLineColor(kBlue);     // blue
    HST_xxcum_Ord1n->DrawCopy("hsame");       // Foam1 IFI off

    HST_xxcum_Ord1->SetLineColor(kGold);     // Gold
    HST_xxcum_Ord1->DrawCopy("hsame");       // Foam1 IFI on

    HST_SigPRD->SetLineColor(kMagenta);     // Magenta
    HST_SigPRD->DrawCopy("hsame");          // PRD41 IFI off

    float ycapt = 0.50;
    CaptT->SetTextColor(kBlack); ycapt += -0.04;
    CaptT->DrawLatex(0.50,ycapt,gTextEne);

    CaptT->SetTextColor(kRed); ycapt += -0.04;
    CaptT->DrawLatex(0.50,ycapt, "KKMC  IFI off ");

    CaptT->SetTextColor(kBlack); ycapt += -0.04;
    CaptT->DrawLatex(0.50,ycapt, "KKMC  IFI on ");

    CaptT->SetTextColor(kGreen); ycapt += -0.04;
    CaptT->DrawLatex(0.50,ycapt, "KKsem  IFI off ");

    CaptT->SetTextColor(kBlue); ycapt += -0.04;
    CaptT->DrawLatex(0.50,ycapt, "Foam1  IFI off ");

    CaptT->SetTextColor(kGold); ycapt += -0.04;
    CaptT->DrawLatex(0.50,ycapt, "Foam1  IFI on ");

    CaptT->SetTextColor(kMagenta); ycapt += -0.04;
    CaptT->DrawLatex(0.50,ycapt, "PRD41 IFI off ");

    //====================plot1========================
    cFigVcum->cd(2);

    TH1D *Hst3_ratio =(TH1D*)HST_xxcum_Ord1n->Clone("Hst3_ratio");
    Hst3_ratio->Divide(HST_SigPRD);  // direct

    Hst3_ratio->SetStats(0);
    Hst3_ratio->SetTitle(0);

    Hst3_ratio->SetMinimum(1-0.02);  //  95GeV
    Hst3_ratio->SetMaximum(1+0.02);  //  95GeV

    Hst3_ratio->SetLineColor(kBlue);
    Hst3_ratio->DrawCopy("h");

    ycapt = 0.90;
    CaptT->SetTextColor(kBlack); ycapt += -0.04;
    CaptT->DrawLatex(0.40,ycapt,gTextEne);

    CaptT->SetTextColor(kBlue); ycapt += -0.04;
    CaptT->DrawLatex(0.40,ycapt, "Foam1/PRD41,  IFI off ");

    //----------------------------
    cFigVcum->cd();

  cout<<" ========================= FigVcum end======================== "<<endl;
  //================================================
}//FigVcum


///////////////////////////////////////////////////////////////////////////////////
void FigAfb()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb =========================== "<<endl;
  //
  TH1D *HAfb_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPL_Ceex2");  // KKMC
  TH1D *HAfb_vTcPL_Ceex2n = (TH1D*)DiskFileB.Get("HAfb_vTcPL_Ceex2n"); // KKMC

  TH1D *HST_xxAfb_Ord1     = (TH1D*)DiskFileB.Get("HST_xxAfb_Ord1");   // FOAM1 IFI on
  TH1D *HST_xxAfb_Ord1n    = (TH1D*)DiskFileB.Get("HST_xxAfb_Ord1n");  // FOAM1 IFI off

  TH1D *HST_xxAfb_Hrd1     = (TH1D*)DiskFileB.Get("HST_xxAfb_Hrd1");   // FOAM1 IFI hard non-IR
  TH1D *HST_xxAfb_Srd1     = (TH1D*)DiskFileB.Get("HST_xxAfb_Srd1");   // FOAM1 IFI soft non-IR
  TH1D *HST_xxAfb_Ird1     = (TH1D*)DiskFileB.Get("HST_xxAfb_Ird1");   // FOAM1 IFI alone

  TH1D *afbv_ISR2_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");    // KKsem

  TH1D *HST_AfbPRD =(TH1D*)HAfb_vTcPL_Ceex2->Clone("HST_AfbPRD");
  LibSem.Ord1fill(HST_AfbPRD,101);
  HST_AfbPRD->SetLineColor(kRed);

  TH1D *HST_AfbPRD_PL =(TH1D*)HAfb_vTcPL_Ceex2->Clone("HST_AfbPRD_PL");
  LibSem.Ord1fill(HST_AfbPRD_PL,100);
  HST_AfbPRD_PL->SetLineColor(kGold);

  TH1D *HST_AfbPL =(TH1D*)HAfb_vTcPL_Ceex2->Clone("HST_AfbPL");
  LibSem.Ord1fill(HST_AfbPL,102);
  HST_AfbPL->SetLineColor(kMagenta);

  TH1D *HST_AfbPLhard =(TH1D*)HAfb_vTcPL_Ceex2->Clone("HST_AfbPLhard");
  LibSem.Ord1fill(HST_AfbPLhard,105);
  HST_AfbPLhard->SetLineColor(kBlue);

  TH1D *HST_AfbPLsoft =(TH1D*)HAfb_vTcPL_Ceex2->Clone("HST_AfbPLsoft");
  LibSem.Ord1fill(HST_AfbPLsoft,106);
  HST_AfbPLsoft->SetLineColor(kGold);


  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  float ycapt = 0.90;

  // zero line
  TH1D *hZero0 = (TH1D*)HST_xxAfb_Ord1->Clone("hZero0");  // zero line
  hZero0->SetLineColor(kBlack);
  for(int i=1; i <= hZero0->GetNbinsX() ; i++) { hZero0->SetBinContent(i, 0); hZero0->SetBinError(i, 0);}

  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigAfb = new TCanvas("cFigAfb","FigAfb", 150, 150,   1000, 550);
  //                            Name    Title                   xoff,yoff, WidPix,HeiPix
  cFigAfb->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigAfb->Divide( 2,  0);
  //
  //====================plot1========================
  cFigAfb->cd(1); // AFB(vmax)
  //
  TH1D *Hst1 = HAfb_vTcPL_Ceex2n;  //  KKMC AFB(vmax) from scat. IFI off
  TH1D *Hst2 = HAfb_vTcPL_Ceex2;   //  KKMC AFB(vmax) from scat. IFI on
  //

  Hst2->SetStats(0);
  Hst2->SetTitle(0);

  if( fabs(gCMSene -10.0) < 1.0) {
    Hst2->SetMinimum(-0.02); // 10GeV
    Hst2->SetMaximum(+0.06); // 10GeV
  }else if( fabs(gCMSene -95.0) < 1.0) {
    Hst2->SetMinimum( 0.15); // 95GeV
    Hst2->SetMaximum(+0.20); // 95GeV
  } else if( fabs(gCMSene -88.0) < 1.0) {
    Hst2->SetMinimum(-0.32); // 88GeV
    Hst2->SetMaximum(-0.20); // 88GeV
  } else if( fabs(gCMSene -91.0) < 1.0) {
    Hst2->SetMinimum(-0.000); // 91GeV
    Hst2->SetMaximum(+0.025); // 91GeV
  }
  Hst2->SetLineColor(kMagenta);            // magenta
  Hst2->DrawCopy("h");                     // KKMC  IFI on AFB(vmax) from scat.
  //
  //hZero0->DrawCopy("hsame");               // zero line
  //
  Hst1->SetLineColor(kBlack);              // black
//  Hst1->DrawCopy("hsame");                 // KKMC IFI off AFB(vmax) from scat.

  HST_xxAfb_Ord1n->SetLineColor(kCyan);
  HST_xxAfb_Ord1n->DrawCopy("hsame");        // cyan, Foam1 MC IFI off

  HST_xxAfb_Ord1->SetLineColor(kPine);
//  HST_xxAfb_Ord1->DrawCopy("hsame");       // Pine, Foam1 MC IFI on

  HST_AfbPRD->SetLineColor(kRed);
  HST_AfbPRD->DrawCopy("hsame");            // red, PRD41 formula, IFI off

  HST_AfbPRD_PL->SetLineColor(kGold);
//  HST_AfbPRD_PL->DrawCopy("hsame");        // red, PRD41+PL219 IFI on

  afbv_ISR2_FSR2->SetLineColor(kBlue);
//  afbv_ISR2_FSR2->DrawCopy("hsame");        // KKsem

  ycapt = 0.90;
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.50,ycapt,gTextEne);

  CaptT->SetTextColor(kMagenta); ycapt += -0.04;
  CaptT->DrawLatex(0.50,ycapt, "KKMC  IFI on ");

  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.50,ycapt, "KKMC  IFI off ");

  CaptT->SetTextColor(kBlue); ycapt += -0.04;
  CaptT->DrawLatex(0.50,ycapt, "KKsem  IFI off ");

  CaptT->SetTextColor(kCyan); ycapt += -0.04;
  CaptT->DrawLatex(0.50,ycapt, "Foam1  IFI off ");

  CaptT->SetTextColor(46); ycapt += -0.04;
  CaptT->DrawLatex(0.50,ycapt, "Foam1  IFI on ");

  CaptT->SetTextColor(kRed); ycapt += -0.04;
  CaptT->DrawLatex(0.50,ycapt, "PRD41  IFI off");

  CaptT->SetTextColor(kGold); ycapt += -0.04;
  CaptT->DrawLatex(0.50,ycapt, "PRD41+PLB219  IFIon");


  //====================plot2========================
  cFigAfb->cd(2);

  TH1D *Hst2_diff1 =(TH1D*)Hst2->Clone("Hst2_diff1");
  Hst2_diff1->Add(Hst2_diff1, Hst1,    1.0, -1.0); // KKMC_IFI   minus KKMC  noIFI black
  Hst2_diff1->SetLineColor(kBlack);

  TH1D *HstPRD_diff =(TH1D*)HST_xxAfb_Ord1n->Clone("HstPRD_diff");
  HstPRD_diff->Add(HstPRD_diff,  HST_AfbPRD,    1.0, -1.0); //
  HstPRD_diff->SetLineColor(kCyan);

  TH1D *Hst_Foam1_diff =(TH1D*)HST_xxAfb_Ord1->Clone("Hst_Foam1_diff");
  Hst_Foam1_diff->Add(Hst_Foam1_diff, HST_xxAfb_Ord1n,    1.0, -1.0); //
  Hst_Foam1_diff->SetLineColor(kRed);

  if( fabs(gCMSene -10.0) < 1.0) {
    Hst2_diff1->SetMinimum(-0.02);   // 189GeV, 10GeV
    Hst2_diff1->SetMaximum( 0.06);   // 189GeV, 10GeV
  }else if( fabs(gCMSene -95.0) < 1.0) {
    Hst2_diff1->SetMinimum(-0.004);  // 95GeV
    Hst2_diff1->SetMaximum( 0.004);  // 95GeV
  } else if( fabs(gCMSene -88.0) < 1.0) {
    Hst2_diff1->SetMinimum(-0.010);  // 88GeV
    Hst2_diff1->SetMaximum( 0.025);  // 88GeV
  }else if( fabs(gCMSene -91.0) < 1.0) {
    Hst2_diff1->SetMinimum(-0.002);  // 91GeV
    Hst2_diff1->SetMaximum( 0.006);  // 91GeV
  }
  Hst2_diff1->SetLineColor(kBlack);
  Hst2_diff1->DrawCopy("h");

  HstPRD_diff->DrawCopy("hsame");

  HST_xxAfb_Hrd1->SetLineColor(46);  // Brown
  HST_xxAfb_Hrd1->DrawCopy("hsame");

  HST_xxAfb_Srd1->SetLineColor(71);  // PineGreen
  HST_xxAfb_Srd1->DrawCopy("hsame");

  HST_xxAfb_Ird1->SetLineColor(kRed);  // Red
  HST_xxAfb_Ird1->DrawCopy("hsame");

  HST_AfbPLhard->DrawCopy("hsame");  // Blue

  HST_AfbPLsoft->DrawCopy("hsame");  // Gold

  HST_AfbPL->DrawCopy("hsame");      // Magenta

//  Hst_Foam1_diff->DrawCopy("hsame");  // Red

  hZero0->DrawCopy("hsame");

  ycapt = 0.90;
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);

  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, "KKMC_IFIon - KKMC_IFIoff ");

  CaptT->SetTextColor(kRed); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, "Foam1: IFI only");

  CaptT->SetTextColor(kCyan); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, "IFI off: Foam1-PRD41 ");

  CaptT->SetTextColor(kPine); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, "Foam1: non-IR IFI hard ");

  CaptT->SetTextColor(kMagenta); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, "PLB219: IFI only ");

  CaptT->SetTextColor(kBlue); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, "PLB219: non-IR IFI hard ");

  CaptT->SetTextColor(kGold); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, "PLB219:   IFI soft, IR subt. ");

  CaptT->SetTextColor(71); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, "Foam1: IFI soft, IR subt. ");

  cFigAfb->cd();
  //================================================

  cout<<" ========================= FigAfb end======================== "<<endl;

}//FigAfb


///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;

  TH1D *hst_weight1  = (TH1D*)DiskFileF.Get("hst_weight1"); // Foam3
//
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo ", 200, 200,   1000, 550);
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
  TH1D *Hst1 = hst_weight1;  //  weight of Foam
  //
  Hst1->SetStats(0);
  Hst1->SetTitle(0);
  Hst1->DrawCopy("h");

  CaptT->DrawLatex(0.02,0.95, " Weight distribution Foam");
  //====================plot2========================
  cFigInfo->cd(2);
  gPad->SetLogy(); // !!!!!!

  //Hst2->SetStats(0);
  //Hst2->SetTitle(0);
  //Hst2->DrawCopy("h");

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

  int Nodes, Nodes2;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
  //
  TH1D *HST_FOAM_NORMA1 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA1");
  Nodes2   =  HST_FOAM_NORMA1->GetBinContent(511);    // No of farm nodes (trick)
  double  CMSeneF  = HST_FOAM_NORMA1->GetBinContent(1)/Nodes2; // CMSene=xpar(1)
  if( fabs(gCMSene/CMSeneF-1) >1e-4 ){
	  cout<<" +++++ Wrong input files !!!! KKMC "<< gCMSene <<"GeV and  FOAM "<< CMSeneF<<"GeV"<<endl;
	  exit(19);
  }
  gNevTot2  = HST_FOAM_NORMA1->GetEntries();       // MC statistics from FOAM
  sprintf(gTextNev2,"FOAM:%10.2e events", gNevTot2);


//////////////////////////////////////////////////////////////////////////
// ========= Preparing plots ==========
  DiskFileB.cd();
  HistNormKKMC();     // Renormalization of MC histograms
  ReMakeKKMC();       // reprocessing MC histos from KKC
  ReMakeFoam1();      // reprocessing MC histos from Foam1
  KKsemMakeHisto();   // prepare histos from KKsem
//========== PLOTTING ==========
// vmax=1
  FigVV();     // dsigma/dv KKMC/Foam
  FigVcum();   // sigma(vmax) KKMC/Foam
  FigAfb();    // AFB(vmax) KKMC/Foam
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
  cout<<"------------------------A.GetListOfKeys-----------------------------"<<endl;
//  DiskFileA.GetListOfKeys()->Print();
  cout<<"------------------------F.GetListOfKeys-----------------------------"<<endl;
//  DiskFileF.GetListOfKeys()->Print();
  //
  cout<< "KKMC: CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "FOAM: CMSene[GeV] = "<< CMSeneF<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  cout<< "FOAM: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<<endl;
  //cout<<"------------------------------end---------------------------------"<<endl;
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


