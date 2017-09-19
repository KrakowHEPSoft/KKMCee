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
//TFile DiskFileA("../workKKMC/histo.root_91GeV_9G"); ????
// July2017 runs
//TFile DiskFileA("../workKKMC/histo.root_91GeV_6G"); //
//TFile DiskFileA("../workKKMC/histo.root_88GeV_4G"); //
//TFile DiskFileA("../workKKMC/histo.root_10GeV_5.7G"); //
//TFile DiskFileA("../workKKMC/histo.root_95GeV.4G");   //

////  *** FOAM
//TFile DiskFileF("../workFOAM/histo.root"); // current
//TFile DiskFileF("../workFOAM/histo.root_10GeV_37G_vmax0.2");
//TFile DiskFileF("../workFOAM/histo.root_88GeV_16G");
//TFile DiskFileF("../workFOAM/histo.root_91GeV_45G");
TFile DiskFileF("../workFOAM/histo.root_95GeV_10G");
//TFile DiskFileF("../workFOAM/histo.root_10GeV_32G");
//
//TFile DiskFileF("../workFOAM/histo.root_91GeV_35G");
//TFile DiskFileF("../workFOAM/histo.root_88GeV_32G");
//TFile DiskFileF("../workFOAM/histo.root_10GeV_15G");
//TFile DiskFileF("../workFOAM/histo.root_95GeV_15G");

//  ****** older FOAM and KKMC progs ******
//TFile DiskFileF("../workFoam0/rmain.root");
//TFile DiskFileF("../workFoam0rmain_95GeV_64M.root");
//
//TString XparFile="../workKKMC/workKKMC_95GeV.input"; // obsolete?
//
//TFile DiskFileA("../workAFB/rmain.root");
//TFile DiskFileA("../workAFB/rmain_95GeV.root"); // 100M new
//TFile DiskFileA("../workAFB/rmain.root_189GeV_100M"); // obsolete

TFile DiskFileB("RhoSemi.root","RECREATE","histograms");

// Interface to KKabox and some extra plotting facilities


///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot, gNevTot2; // from KKMC and KKfoam MC runs (histograms)
char   gTextEne[100], gTextNev[100], gTextNev2[100];
//
//int    gNbMax=45;         // gCosTheta = 45/50=0.90
//double gCosTheta=0.90;    // to be synchronized with gNbMax
//
int    gNbMax=50;         // gCosTheta = 45/50=0.90
double gCosTheta=1.00;    // to be synchronized with gNbMax
//
KKplot LibSem("KKplot");
///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////
void PlotSame(TH1D *HST, double &ycapt, Int_t kolor, TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
  //double yy=ycapt;
  CaptT->DrawLatex(0.40,ycapt, opis);
}// PlotSame



void TestNorm(){
// testing/debugging normalization
cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
cout<<"<<<<<<<<<<<<<<<<<<<<<<< TestNorm >>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;

TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
TH1D *HST_tx_Ceex2n   = (TH1D*)DiskFileF.Get("HST_tx_Ceex2n");  // FOAM testing norm.

//int  NevTot = 4.5e6;
int NevTot  = HST_FOAM_NORMA3->GetEntries();
//NevTot = 4.5e6;
double Xsav = HST_FOAM_NORMA3->GetBinContent(0)/NevTot; // NANOBARNS
cout<<"TestNorm: Xsav = "<<Xsav<<"  NevTot =  "<< NevTot<<endl;
// first 3 bins only for vmax=0.006
double SWT =0;
double SSWT =0;
//for( int i=0; i <= 3; i++){ // 3 bins v<0.006
//for( int i=0; i <= 3; i++){ // 100 bins v<0.200
for( int i=0; i <= 3; i++){ // 100 bins v<0.020
  SWT  += HST_tx_Ceex2n->GetBinContent(i);
  SSWT += sqr(HST_tx_Ceex2n->GetBinError(i));
}//for

double AWT = SWT/NevTot;
double Sigma = sqrt(SSWT/NevTot-sqr(AWT));
double Error = Sigma/sqrt(NevTot);
double Errel = Error/AWT;

double Sigm0 = sqrt(SSWT/NevTot);       // approximation used in ROOT
double Erre0 = Sigm0/sqrt(NevTot)/AWT;  // approximation used in ROOT

cout<<" TestNorm:  AWT="<< AWT << endl;
cout<<" TestNorm:  Errel="<< Errel << endl;
cout<<" TestNorm:  Erre0="<< Erre0 << " ROOT approx. " <<endl;
cout<<" TestNorm:  ratio="<< Erre0/Errel<<endl;
cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;

}//TestNorm


void HistNormKKMC(){
  //
  cout<<"----------------------------- HistNormKKMC ------------------------------------"<<endl;
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

   // Wide range, vmax<1.
    TH2D *sct_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2");
    TH2D *sct_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n");
    TH2D *sct_vTcPR_EEX2   = (TH2D*)DiskFileA.Get("sct_vTcPR_EEX2");
    TH2D *sct_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2");
    TH2D *sct_vTcPL_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2n");
    cout<<"ReMakeKKMC [2]"<<endl;
    //****************************************************************************************
    // Distributions of v=vTrue<vmax<0.20, c=cos(theta) with 100 bins
    //****************************************************************************************
    //gNbMax=45;         // cosThetaMax = 45/50=0.90 Now global variable
    // IFI on
    TH1D                    *HTot2_vTcPR_Ceex2, *HAfb2_vTcPR_Ceex2;
    ProjV( sct_vTcPR_Ceex2,  HTot2_vTcPR_Ceex2,  HAfb2_vTcPR_Ceex2, gNbMax);  //!!!!
    HTot2_vTcPR_Ceex2->SetName("HTot2_vTcPR_Ceex2");
    HAfb2_vTcPR_Ceex2->SetName("HAfb2_vTcPR_Ceex2");
    // IFI off
    TH1D                    *HTot2_vTcPR_Ceex2n, *HAfb2_vTcPR_Ceex2n;
    ProjV( sct_vTcPR_Ceex2n, HTot2_vTcPR_Ceex2n,  HAfb2_vTcPR_Ceex2n, gNbMax);  //!!!!
    HTot2_vTcPR_Ceex2n->SetName("HTot2_vTcPR_Ceex2n");
    HAfb2_vTcPR_Ceex2n->SetName("HAfb2_vTcPR_Ceex2n");
    // IFI off
    TH1D                    *HTot2_vTcPR_EEX2, *HAfb2_vTcPR_EEX2;
    ProjV( sct_vTcPR_EEX2, HTot2_vTcPR_EEX2,  HAfb2_vTcPR_EEX2, gNbMax);  //!!!!
    HTot2_vTcPR_EEX2->SetName("HTot2_vTcPR_EEX2");
    HAfb2_vTcPR_EEX2->SetName("HAfb2_vTcPR_EEX2");
    // IFI on
    TH1D                    *HTot2_vTcPL_Ceex2, *HAfb2_vTcPL_Ceex2;
    ProjV( sct_vTcPL_Ceex2,  HTot2_vTcPL_Ceex2,  HAfb2_vTcPL_Ceex2, gNbMax);  //!!!!
    HTot2_vTcPL_Ceex2->SetName("HTot2_vTcPL_Ceex2");
    HAfb2_vTcPL_Ceex2->SetName("HAfb2_vTcPL_Ceex2");
    // IFI off
    TH1D                    *HTot2_vTcPL_Ceex2n, *HAfb2_vTcPL_Ceex2n;
    ProjV( sct_vTcPL_Ceex2n, HTot2_vTcPL_Ceex2n,  HAfb2_vTcPL_Ceex2n, gNbMax);  //!!!!
    HTot2_vTcPL_Ceex2n->SetName("HTot2_vTcPL_Ceex2n");
    HAfb2_vTcPL_Ceex2n->SetName("HAfb2_vTcPL_Ceex2n");

  ///****************************************************************************************
  // Pure KKMC reprocessing part
    int nbMax;
    cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;
    // Wide range, vmax<1.
    TH2D *sca_vTcPR_Eex2   = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
    TH2D *sca_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
    TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n");
    //****************************************************************************************
    // Distributions of v=vTrue<1.0 unlimited c=cos(theta), 50 bins
    //****************************************************************************************
    nbMax=0;            // cosThetaMax = 1.0
    TH1D                    *HTot_vTcPR_Eex2, *HAfb_vTcPR_Eex2;
    ProjV( sca_vTcPR_Eex2,  HTot_vTcPR_Eex2,  HAfb_vTcPR_Eex2, nbMax);  //!!!!
    HTot_vTcPR_Eex2->SetName("HTot_vTcPR_Eex2");
    HAfb_vTcPR_Eex2->SetName("HAfb_vTcPR_Eex2");
    //
    TH1D                    *HTot_vTcPR_Ceex2, *HAfb_vTcPR_Ceex2;
    ProjV( sca_vTcPR_Ceex2,  HTot_vTcPR_Ceex2,  HAfb_vTcPR_Ceex2, nbMax);  //!!!!
    HTot_vTcPR_Ceex2->SetName("HTot_vTcPR_Ceex2");
    HAfb_vTcPR_Ceex2->SetName("HAfb_vTcPR_Ceex2");
    // IFI off
    TH1D                    *HTot_vTcPR_Ceex2n, *HAfb_vTcPR_Ceex2n;
    ProjV( sca_vTcPR_Ceex2n, HTot_vTcPR_Ceex2n,  HAfb_vTcPR_Ceex2n, nbMax);  //!!!!
    HTot_vTcPR_Ceex2n->SetName("HTot_vTcPR_Ceex2n");
    HAfb_vTcPR_Ceex2n->SetName("HAfb_vTcPR_Ceex2n");
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

  cout<<"================ ReMakeKKMC ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//ReMakeKKMC



///////////////////////////////////////////////////////////////////////////////////
void ReMakeFoam35(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeFoam35  BEGIN  ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA5");

  TH1D *HST_xx_Ceex2n = (TH1D*)DiskFileF.Get("HST_xx_Ceex2n");  // FOAM
  TH1D *HST_tx_Ceex2n = (TH1D*)DiskFileF.Get("HST_tx_Ceex2n");  // FOAM testing norm.
  TH2D *SCA_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCA_xc_Ceex2n");  // FOAM big   range x<0.99
  TH2D *SCT_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2n");  // FOAM small range x<0.20
  TH2D *SCT_xc_EEX2   = (TH2D*)DiskFileF.Get("SCT_xc_EEX2");    // FOAM small range x<0.20

  TH2D *SCA_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCA_xc_Ceex2");   // FOAM big   range x<0.99
  TH2D *SCT_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2");   // FOAM small range x<0.20

  HisNorm1(HST_FOAM_NORMA3, HST_xx_Ceex2n );  // normalizing
  HisNorm1(HST_FOAM_NORMA3, HST_tx_Ceex2n );  // normalizing testing norm.
  HisNorm2(HST_FOAM_NORMA3, SCA_xc_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA3, SCT_xc_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA3, SCT_xc_EEX2 );    // normalizing

  HisNorm2(HST_FOAM_NORMA5, SCA_xc_Ceex2 );   // normalizing
  HisNorm2(HST_FOAM_NORMA5, SCT_xc_Ceex2 );   // normalizing

  // sigma(vmax) direct histogramming (0,1) range
  TH1D *HST_xmax_Ceex2n;
  MakeCumul(HST_xx_Ceex2n,HST_xmax_Ceex2n);
  HST_xmax_Ceex2n->SetName("HST_xmax_Ceex2n");
  // sigma(vmax) direct histogramming, (0,0.2) range
  TH1D *HST_txmax_Ceex2n;
  MakeCumul(HST_tx_Ceex2n,HST_txmax_Ceex2n);
  HST_txmax_Ceex2n->SetName("HST_txmax_Ceex2n");

  // sigma(vmax) and AFB(vmax) from KKfoam, vmax<1, 50 bins in costheta
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

  // sigma(vmax) and AFB(vmax) from KKfoam scat. vmax<0.2, 100 bins in ctheta
  //gNbMax=45;           // cosThetaMax = 45/50=0.90 Now global variable
  TH1D                *Htot2_xmax_Ceex2n, *Hafb2_xmax_Ceex2n;
  ProjV( SCT_xc_Ceex2n, Htot2_xmax_Ceex2n,  Hafb2_xmax_Ceex2n, gNbMax);  //!!!!
  Htot2_xmax_Ceex2n->SetName("Htot2_xmax_Ceex2n");
  Hafb2_xmax_Ceex2n->SetName("Hafb2_xmax_Ceex2n");
  //
  TH1D                *Htot2_xmax_EEX2, *Hafb2_xmax_EEX2;
  ProjV( SCT_xc_EEX2, Htot2_xmax_EEX2,  Hafb2_xmax_EEX2, gNbMax);  //!!!!
  Htot2_xmax_EEX2->SetName("Htot2_xmax_EEX2");
  Hafb2_xmax_EEX2->SetName("Hafb2_xmax_EEX2");
  //
  TH1D                *Htot2_xmax_Ceex2, *Hafb2_xmax_Ceex2;
  ProjV( SCT_xc_Ceex2, Htot2_xmax_Ceex2,  Hafb2_xmax_Ceex2, gNbMax);  //!!!!
  Htot2_xmax_Ceex2->SetName("Htot2_xmax_Ceex2");
  Hafb2_xmax_Ceex2->SetName("Hafb2_xmax_Ceex2");


  cout<<"================ ReMakeFoam35 ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//ReMakeFoam35


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
  //double gCosTheta = 0.9; // now glogal variable
  kksem_setcrange_(-gCosTheta, gCosTheta); // F+B cos(theta) range
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  LibSem.VVplot(vcum_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
//-------------------------------------------------
//    AFB(vmax) for limited/unlimited c=cos(theta)
//-------------------------------------------------
  kksem_setcrange_(0, gCosTheta); // forward cos(theta)
  TH1D *afbv_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("afbv_ISR2_FSR2");
  LibSem.VVplot(afbv_ISR2_FSR2, KF, chak, KeyDis, KeyFob);// Forward
  afbv_ISR2_FSR2->Add(afbv_ISR2_FSR2, vcum_ISR2_FSR2, 2.0, -1.0) ; // numerator F-B = 2F-(F+B)
  afbv_ISR2_FSR2->Divide(vcum_ISR2_FSR2);                          // finally (F-B)(F+B)
  kksem_setcrange_(-1.0, 1.0); // back to default full range
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
 //
  TH1D *Hpro_vT_Ceex2n    = (TH1D*)DiskFileB.Get("Hpro_vT_Ceex2n");     // KKMC dsigma/dv IFI off, from scat.
  TH1D *HTot_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2n");  // KKMC sigma(vmax) from scat.
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
  CaptT->DrawLatex(0.02,0.95, "#sigma(v_{max}) ISR+FSR, Black KKMC_CEEX2, Blue FOAM");
  CaptT->DrawLatex(0.60,0.75,gTextEne);

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
  CaptT->DrawLatex(0.02,0.95,"#sigma(v_{max}) ISR+FSR, Ratio KKMC/FOAM");
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
  CaptT->DrawLatex(0.60,0.75,gTextEne);

  //----------------------------
  cFigVdist->cd();
  //================================================
}//FigVdist


///////////////////////////////////////////////////////////////////////////////////
void FigAfb()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb =========================== "<<endl;

  TH1D *HAfb_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2");  // KKMC
  TH1D *HAfb_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2n"); // KKMC
  //
  TH1D *HAfb_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPL_Ceex2");  // KKMC
  TH1D *HAfb_vTcPL_Ceex2n = (TH1D*)DiskFileB.Get("HAfb_vTcPL_Ceex2n"); // KKMC

//  TH1D *afbv_ISR2_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");    // KKsem

  TH1D *Hafb_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb_xmax_Ceex2n");  // FOAM scatt.
  TH1D *Hafb_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb_xmax_Ceex2");   // FOAM scatt.
 //
  TH1D *HST_PLBZ =(TH1D*)HAfb_vTcPR_Ceex2->Clone("HST_PLBZ");
  LibSem.Ord1fill(HST_PLBZ,100);
  HST_PLBZ->SetLineColor(kCyan);
  //
  TH1D *HST_IFI5 =(TH1D*)HAfb_vTcPL_Ceex2->Clone("HST_IFI5");
  LibSem.Ord1fill(HST_IFI5,105);
  HST_IFI5->SetLineColor(kCyan);
 //

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
  //TH1D *Hst1 = HAfb_vTcPR_Ceex2n;  //  KKMC AFB(vmax) from scat. IFI off
  //TH1D *Hst2 = HAfb_vTcPR_Ceex2;   //  KKMC AFB(vmax) from scat. IFI on
  //
  TH1D *Hst1 = HAfb_vTcPL_Ceex2n;  //  KKMC AFB(vmax) from scat. IFI off
  TH1D *Hst2 = HAfb_vTcPL_Ceex2;   //  KKMC AFB(vmax) from scat. IFI on
  // zero line
  TH1D *hZero0 = (TH1D*)Hst2->Clone("hZero0");  // zero line
  hZero0->SetLineColor(kBlack);
  for(int i=1; i <= hZero0->GetNbinsX() ; i++) { hZero0->SetBinContent(i, 0); hZero0->SetBinError(i, 0);}
  //
  Hst2->SetStats(0);
  Hst2->SetTitle(0);
  if( fabs(gCMSene-10e0) <0.01 ) {
    Hst2->SetMinimum(-0.02);  // 10GeV
    Hst2->SetMaximum( 0.08);  // 10GeV
  }
  Hst2->SetLineColor(kMagenta);            // magenta
  Hst2->DrawCopy("h");                     // KKMC AFB(vmax) from scat. IFI on
  //
  hZero0->DrawCopy("hsame");               // zero line
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
  HST_PLBZ->DrawCopy("hsame");             // analytical formula
  //
  CaptT->DrawLatex(0.12,0.95, "A_{FB}^{IF on}(v_{max}):  KKMC=magenta, Foam=green");
  CaptT->DrawLatex(0.12,0.85, "A_{FB}^{IFI off}(v_{max}): KKMC=black,  Foam=blue");
  CaptT->DrawLatex(0.60,0.75,gTextEne);
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

  HST_IFI5->DrawCopy("hsame");    // analytical
  //
  CaptT->DrawLatex(0.12,0.95,"A_{FB}^{IFI}(v_{max}): Black KKMC, Magenta=FOAM");
  CaptT->DrawLatex(0.12,0.85,"A^{KKMC}_{FB}-A^{FOAM}: Green=IFI, Blue=NOIFI");
  CaptT->DrawLatex(0.60,0.75,gTextEne);

  cFigAfb->cd();
  //================================================
}//FigAfb



///////////////////////////////////////////////////////////////////////////////////
void FigAfb2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb2 =========================== "<<endl;

  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n");  //
  TH1D *HAfb2_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2");  //
  //
  TH1D *HAfb2_vTcPL_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2n");  //
  TH1D *HAfb2_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2");  //
  //
  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // FOAM scatt.
  TH1D *Hafb2_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2");   // FOAM scatt.

  TH1D *afbv_ISR2_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");    // KKsem


  TH1D *HST_PLBZ2 =(TH1D*)HAfb2_vTcPR_Ceex2->Clone("HST_PLBZ2");
  LibSem.Ord1fill(HST_PLBZ2,101);
  HST_PLBZ2->SetLineColor(kCyan);
/*
  // A_FB from PLB219,p103 ]]]
  double alfinv  = 137.035989;
  double alfpi   = 1/alfinv/3.1415926535;
  TH1D *HST_PL =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_PL");
  HST_PL->SetLineColor(kMagenta);
  int Nbin    = HST_PL->GetNbinsX();
  double vmax = HST_PL->GetXaxis()->GetXmax();
  for(int i=1; i <= Nbin ; i++) {
	  double vv = (i*vmax)/Nbin;
	  double afb = 3.0/2.0 *alfpi*( 3*vv+log(1-vv/2) ); // only gamma
	  cout<< "FigAfb2: ib, vv, afb ="<<i<< "  "<< vv << "   "<<afb<<endl;
	  HST_PL->SetBinContent(i, afb);
	  HST_PL->SetBinError(i, 0);
  }// i
*/
  //
   TH1D *HST_IFI4 =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_IFI4");
   LibSem.Ord1fill(HST_IFI4,105);
   HST_IFI4->SetLineColor(kCyan);
 //


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
  //TH1D *Hst1 = HAfb2_vTcPR_Ceex2n;         // KKMC AFB(vmax) from scat. IFI off
  //TH1D *Hst2 = HAfb2_vTcPR_Ceex2;          // KKMC AFB(vmax) from scat. IFI on
  TH1D *Hst1 = HAfb2_vTcPL_Ceex2n;         // KKMC AFB(vmax) from scat. IFI off
  TH1D *Hst2 = HAfb2_vTcPL_Ceex2;          // KKMC AFB(vmax) from scat. IFI on
  //
  Hst2->SetStats(0);
  Hst2->SetTitle(0);
  Hst2->SetLineColor(kMagenta);            // magenta

  Hst2->SetMinimum(-1);  //
  Hst2->SetMaximum( 1);  //

  if( fabs(gCMSene-10e0) <0.01 ) {
    Hst2->SetMinimum(-0.02);  // 10GeV
    Hst2->SetMaximum( 0.08);  // 10GeV
  }
  Hst2->DrawCopy("h");                     // KKMC AFB(vmax) from scat. IFI on
  //
  Hst1->SetLineColor(kBlack);              // black
  Hst1->DrawCopy("hsame");                 // KKMC AFB(vmax) from scat. IFI off
  //
  Hafb2_xmax_Ceex2n->SetLineColor(kBlue);
  Hafb2_xmax_Ceex2n->DrawCopy("hsame");   // Foam IFI OFF
  //
  Hafb2_xmax_Ceex2->SetLineColor(kGreen);
  Hafb2_xmax_Ceex2->DrawCopy("hsame");    // Foam IFI ON

  afbv_ISR2_FSR2->SetLineColor(kRed);     // KKsem  RED!!!
  //afbv_ISR2_FSR2->DrawCopy("hsame");

  HST_PLBZ2->DrawCopy("hsame");     // analytical

  //
  CaptT->DrawLatex(0.12,0.95, "A_{FB}^{IFI on}(v_{max})  KKMC=magenta, Foam=green");
  CaptT->DrawLatex(0.12,0.85, "A_{FB}^{IFI off}(v_{max}) KKMC=black,  Foam=blue");
  CaptT->DrawLatex(0.60,0.75,gTextEne);
  //====================plot2========================
  cFigAfb2->cd(2);

  TH1D *Hst21_diff =(TH1D*)Hst2->Clone("Hst21_diff");
  Hst21_diff->Add(Hst21_diff, Hst1,  1.0, -1.0); // KKMC_IFI
  Hst21_diff->SetLineColor(kBlack);              // blue, KKMC

  TH1D *HST21_diff =(TH1D*)Hafb2_xmax_Ceex2->Clone("HST21_diff");
  HST21_diff->Add(HST21_diff, Hafb2_xmax_Ceex2n,  1.0, -1.0); // FOAMC_IFI
  HST21_diff->SetLineColor(kMagenta);

  TH1D *HstKF_diff =(TH1D*)Hst2->Clone("HstKF_diff");
  HstKF_diff->Add(HstKF_diff, Hafb2_xmax_Ceex2,  1.0, -1.0); // KKMC-Foam IFIon
  HstKF_diff->SetLineColor(kGreen);                          // KKMC-Foam IFIon

  TH1D *HstKFn_diff =(TH1D*)Hst1->Clone("HstKFn_diff");
  HstKFn_diff->Add(HstKFn_diff, Hafb2_xmax_Ceex2n,  1.0, -1.0); // KKMC-Foam IFIoff
  HstKFn_diff->SetLineColor(kBlue);                             // KKMC-Foam IFIoff

  TH1D *HstPL_diff =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HstPL_diff");
  HstPL_diff->Add(HstPL_diff, Hafb2_xmax_Ceex2,  1.0, -1.0);    // KKMC-Foam IFIon
  HstPL_diff->SetLineColor(kRed);

  Hst21_diff->SetMinimum(-0.004);  // zoom
  Hst21_diff->SetMaximum( 0.004);  // zoom

  Hst21_diff->DrawCopy("h");
  HST21_diff->DrawCopy("hsame");
  HstPL_diff->DrawCopy("hsame"); //!!! cosThetaPL !!!
  HstKF_diff->DrawCopy("hsame");
  HstKFn_diff->DrawCopy("hsame");
  //
  ////HST_PL->DrawCopy("hsame"); // !!!???
  HST_IFI4->DrawCopy("hsame"); // !!!???

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
  cout<<" ========================= FigTech IFI off =========================== "<<endl;

// sigma(vmax)
  TH1D *HTot2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HTot2_vTcPR_Ceex2n"); // KKMC sigma(vmax) from scat.
  TH1D *HTot2_vTcPR_EEX2   = (TH1D*)DiskFileB.Get("HTot2_vTcPR_EEX2");   // KKMC sigma(vmax) from scat.
  TH1D *Htot2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Htot2_xmax_Ceex2n");  // FOAM3 scatt. GPS Born
  TH1D *Htot2_xmax_EEX2    = (TH1D*)DiskFileB.Get("Htot2_xmax_EEX2");    // FOAM3 scatt. EEX Born
  TH1D *vcum_ISR2_FSR2     = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");     // KKsem
// AFB(vmax)
  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n"); // KKMC AFB(vmax) from scat
  TH1D *HAfb2_vTcPR_EEX2   = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_EEX2");   // KKMC AFB(vmax) from scat
  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // AFB FOAM3 scatt. GPS Born
  TH1D *Hafb2_xmax_EEX2    = (TH1D*)DiskFileB.Get("Hafb2_xmax_EEX2");    // AFB FOAM3 scatt. EEX Born
  TH1D *afbv_ISR2_FSR2     = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");     // AFB KKsem

  TH1D *HST_txmax_Ceex2n   = (TH1D*)DiskFileB.Get("HST_txmax_Ceex2n");   // FOAM direct

  /*
  // correcting missing phase space beta factor in GPS Born
  double Mmu = 0.105, bin,vv, beta;
  TH1D *HST_bad = Htot2_xmax_EEX2;
  int Nbin    = HST_bad->GetNbinsX();
  double vmax = HST_bad->GetXaxis()->GetXmax();
  for(int i=1; i <= Nbin ; i++) {
	  bin= HST_bad->GetBinContent(i);
	  vv = (i*vmax)/Nbin;
	  beta = sqrt(1-4*sqr(Mmu/gCMSene)/(1-vv) );
	  cout<< " vv, beta ="<< vv << "   "<<beta<<endl;
	  HST_bad->SetBinContent(i, bin*beta);
  }*/

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
  TH1D *HstTech_ratio0  = HstRatio("HstTech_ratio0",  Htot2_xmax_EEX2,  vcum_ISR2_FSR2, kGreen);
  TH1D *HstTech_ratio1  = HstRatio("HstTech_ratio1",  Htot2_xmax_Ceex2n,  vcum_ISR2_FSR2, kGreen);
  TH1D *HstTech_ratio3  = HstRatio("HstTech_ratio3",  HST_txmax_Ceex2n,  vcum_ISR2_FSR2, kRed);
  TH1D *HstTech_ratio2  = HstRatio("HstTech_ratio2",  HTot2_vTcPR_Ceex2n,  vcum_ISR2_FSR2, kBlack);
  TH1D *HstTech_ratio4  = HstRatio("HstTech_ratio4",   HTot2_vTcPR_EEX2,    vcum_ISR2_FSR2, kMagenta);

  TH1D *HstTech_ratio = HstTech_ratio0;   // FoamEEX2/KKsem IFIoff magenta

  HstTech_ratio->SetStats(0);
  HstTech_ratio->SetTitle(0);
  HstTech_ratio->SetMinimum(1 -0.0007); HstTech_ratio->SetMaximum(1 +0.0007);  // zoom
  HstTech_ratio->GetXaxis()->SetTitle("v_{max}");
  HstTech_ratio->DrawCopy("h");

  HstTech_ratio1->DrawCopy("hsame");      // FoamGPS/KKsem IFIoff   green
  //HstTech_ratio3->DrawCopy("hsame");    // testing norm. Foam/KKsem
  HstTech_ratio4->DrawCopy("hsame");      // KKMCeex/KKsem IFIoff   magenta
  HstTech_ratio2->DrawCopy("hsame");      // KKMCceexn/KKsem IFIoff black

  double ycapt = 0.40; // starting value, to be decremented below
  PlotSame(HstTech_ratio0,    ycapt, kBlue,       "FoamEEX2/KKsem IFIoff ");
  PlotSame(HstTech_ratio1,    ycapt, kGreen,       "FoamGPS/KKsem  IFIoff ");
//  PlotSame(HstTech_ratio3,    ycapt, kRed,         "testing norm. Foam/KKsem ");
  PlotSame(HstTech_ratio4,    ycapt, kMagenta,     "KKMCeex/KKsem   IFIoff ");
  PlotSame(HstTech_ratio2,    ycapt, kBlack,       "KKMCceexn/KKsem IFIoff");

  TH1D *hOne = (TH1D*)HstTech_ratio->Clone("hOne");  // unity line
  for(int i=1; i <= hOne->GetNbinsX() ; i++) { hOne->SetBinContent(i, 1); hOne->SetBinError(i, 0);}
  hOne->SetLineColor(kBlack);
  hOne->DrawCopy("hsame");


  CaptT->DrawLatex(0.12,0.95,"#sigma^{IFIoff}(v_{max}) ");

  //====================plot2========================
  cFigTech->cd(2);
  TH1D *HstTech_diff0  = HstDiff("HstTech_diff0",   Hafb2_xmax_EEX2,    afbv_ISR2_FSR2, kMagenta);
  TH1D *HstTech_diff1  = HstDiff("HstTech_diff1",   Hafb2_xmax_Ceex2n,  afbv_ISR2_FSR2, kGreen);
  TH1D *HstTech_diff2  = HstDiff("HstTech_diff2",   HAfb2_vTcPR_Ceex2n, afbv_ISR2_FSR2, kBlack);
  TH1D *HstTech_diff4  = HstDiff("HstTech_diff4",   HAfb2_vTcPR_EEX2,   afbv_ISR2_FSR2, kBlue);

  TH1D *HstTech_diff= HstTech_diff4;
  HstTech_diff->SetStats(0); HstTech_diff->SetTitle(0);
  HstTech_diff->SetMinimum(-0.0005); HstTech_diff->SetMaximum( 0.0005);  // zoom

  HstTech_diff->GetXaxis()->SetTitle("v_{max}");
  HstTech_diff->DrawCopy("h");

  ycapt = 0.90; // starting value, to be decremented below
  PlotSame(HstTech_diff4,    ycapt, kBlue,       "KKMCeex  -KKsem IFIoff ");
  PlotSame(HstTech_diff2,    ycapt, kBlack,      "KKMCceexn-KKsem IFIoff ");
  PlotSame(HstTech_diff0,    ycapt, kMagenta,    "Foam3EEX -KKsem IFIoff ");
  PlotSame(HstTech_diff1,    ycapt, kGreen,      "Foam3GPS -KKsem IFIoff ");

  TH1D *hZero = (TH1D*)HstTech_diff->Clone("hZero");  // unity line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}
  hZero->SetLineColor(kRed);
  hZero->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"A_{FB}^{IFIoff}(v_{max})");
  ycapt =0.30;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev2); ycapt += -0.04;

  cFigTech->cd();
  //================================================

}//FigTech


///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;

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
  Hst1->SetStats(0);
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
  int Nodes, Nodes2;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
  //
  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  Nodes2   =  HST_FOAM_NORMA3->GetBinContent(511);    // No of farm nodes (trick)
  double  CMSeneF  = HST_FOAM_NORMA3->GetBinContent(1)/Nodes2; // CMSene=xpar(1)
  if( fabs(gCMSene/CMSeneF-1) >1e-4 ){
	  cout<<" +++++ Wrong input files !!!! KKMC "<< gCMSene <<"GeV and  FOAM "<< CMSeneF<<"GeV"<<endl;
	  exit(19);
  }
  gNevTot2  = HST_FOAM_NORMA3->GetEntries();       // MC statistics from KKMC
  sprintf(gTextNev2,"FOAM:%10.2e events", gNevTot2);

//////////////////////////////////////////////////////////////////////////
// ========= Preparing plots ==========
  DiskFileB.cd();
  TestNorm();          // special test of normalizaion
  HistNormKKMC();      // Renormalization of MC histograms
  ReMakeFoam35();      // reprocessing MC histos from KKC and Foam
  ReMakeKKMC();        // reprocessing MC histos from KKC and Foam
  KKsemMakeHisto();    // prepare histos from KKsem
//========== PLOTTING ==========
// vmax=1
  FigVdist();  // sigma(v) and sigma(vmax) KKMC/Foam
  FigAfb();    // AFB(vmax) KKMC/Foam
// vmax =0.2
  FigAfb2();
  FigTech();
// weight distribution
  //igInfo();
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
  DiskFileA.GetListOfKeys()->Print();
  cout<<"------------------------F.GetListOfKeys-----------------------------"<<endl;
  DiskFileF.GetListOfKeys()->Print();
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  cout<< "FOAM: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<<endl;
  //cout<<"------------------------------end---------------------------------"<<endl;
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


