//////////////////////////////////////////////////////////////////////
//    make PlotFoam2-run
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
//TFile DiskFileA("../workKKMC/histo.root"); // current
//TFile DiskFileA("../workKKMC/histo.root_91GeV_6G"); //
//TFile DiskFileA("../workKKMC/histo.root_88GeV_4G"); //
TFile DiskFileA("../workKKMC/histo.root_10GeV_5.7G"); //
//TFile DiskFileA("../workKKMC/histo.root_95GeV.4G");   //

////  *** FOAM
//TFile DiskFileF("../workFOAM/histo.root"); // current
TFile DiskFileF("../workFOAM/histo.root_10GeV_37G_vmax0.2");
//TFile DiskFileF("../workFOAM/histo.root_88GeV_16G");
//TFile DiskFileF("../workFOAM/histo.root_91GeV_45G");
//TFile DiskFileF("../workFOAM/histo.root_95GeV_10G");
//TFile DiskFileF("../workFOAM/histo.root_10GeV_32G");


TFile DiskFileB("RhoSemi.root","RECREATE","histograms");

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot, gNevTot2; // from KKMC and KKfoam MC runs (histograms)
char   gTextEne[100], gTextNev[100], gTextNev2[100];

//
int    gNbMax=50;         // gCosTheta = 45/50=0.90
double gCosTheta=1.00;    // to be synchronized with gNbMax
//
KKplot LibSem("KKplot");
///////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA5");

  TH2D *SCT_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2n");  // FOAM small range x<0.20
  TH2D *SCT_xc_EEX2   = (TH2D*)DiskFileF.Get("SCT_xc_EEX2");    // FOAM small range x<0.20
  TH2D *SCT_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2");   // FOAM small range x<0.20

  // sigma(vmax) and AFB(vmax) from KKfoam scat. vmax<0.2, 100 bins in ctheta
  //gNbMax=45;           // cosThetaMax = 45/50=0.90 Now global variable
  TH1D                 *Htot2_xmax_Ceex2n, *Hafb2_xmax_Ceex2n;
  ProjV( SCT_xc_Ceex2n, Htot2_xmax_Ceex2n,  Hafb2_xmax_Ceex2n, gNbMax);  //!!!!
  Htot2_xmax_Ceex2n->SetName("Htot2_xmax_Ceex2n");
  Hafb2_xmax_Ceex2n->SetName("Hafb2_xmax_Ceex2n");
  //
  TH1D               *Htot2_xmax_EEX2, *Hafb2_xmax_EEX2;
  ProjV( SCT_xc_EEX2, Htot2_xmax_EEX2,  Hafb2_xmax_EEX2, gNbMax);  //!!!!
  Htot2_xmax_EEX2->SetName("Htot2_xmax_EEX2");
  Hafb2_xmax_EEX2->SetName("Hafb2_xmax_EEX2");
  //
  TH1D                *Htot2_xmax_Ceex2, *Hafb2_xmax_Ceex2;
  ProjV( SCT_xc_Ceex2, Htot2_xmax_Ceex2,  Hafb2_xmax_Ceex2, gNbMax);  //!!!!
  Htot2_xmax_Ceex2->SetName("Htot2_xmax_Ceex2");
  Hafb2_xmax_Ceex2->SetName("Hafb2_xmax_Ceex2");

  cout<<"================ ReMakeMChisto ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto





///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto2(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto2  BEGIN ============================"<<endl;

  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_EEX2") );

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2") );

////////////////////////////////////////////////////////////////////
// Pure KKMC reprocessing part
// from bigger scattergram and restricted vmax<0.2
  //////////////////////////////////////////////////////////////////
   // Wide range, vmax<1.
    TH2D *sct_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2");
    TH2D *sct_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n");
    TH2D *sct_vTcPR_EEX2   = (TH2D*)DiskFileA.Get("sct_vTcPR_EEX2");
    TH2D *sct_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2");
    cout<<"ReMakeMChisto2 [2]"<<endl;
    ///****************************************************************************************
    ///****************************************************************************************
    /// Distributions of v=vTrue<vmax<0.20, c=cos(theta) with 100 bins
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


    ///****************************************************************************************

  cout<<"================ ReMakeMChisto2 ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto2





///////////////////////////////////////////////////////////////////////////////////
void FigAfb2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb2 =========================== "<<endl;

  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n");  //
  TH1D *HAfb2_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2");  //
  //
  TH1D *HAfb2_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2");  //

  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // FOAM scatt.
  TH1D *Hafb2_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2");   // FOAM scatt.


  // A_FB from PLB219,p103 ]]]
  double alfinv  = 137.035989;
  double alfpi   = 1/alfinv/3.1415926535;
  TH1D *HST_PL =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_PL");
  HST_PL->SetLineColor(kMagenta);
  int Nbin    = HST_PL->GetNbinsX();
  double vmax = HST_PL->GetXaxis()->GetXmax();
  for(int i=1; i <= Nbin ; i++) {
	  double vv = (i*vmax)/Nbin;
	  double afb = alfpi*( 3*vv+log(1-vv/2) ); // only gamma
	  cout<< " vv, afb ="<< vv << "   "<<afb<<endl;
	  HST_PL->SetBinContent(i, afb);
	  HST_PL->SetBinError(i, 0);
  }// i

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
  //
  Hafb2_xmax_Ceex2->SetLineColor(kGreen);
  Hafb2_xmax_Ceex2->DrawCopy("hsame");    // Foam IFI ON

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
  HST_PL->DrawCopy("hsame"); // !!!???

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
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  /////////////////////////////////////////////////////////
  LibSem.Initialize(DiskFileA);  // for non-farm case
  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
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
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  cout<< "FOAM: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<<endl;
//////////////////////////////////////////////////////////////////////////
// ========= Preparing plots ==========
  DiskFileB.cd();
  ReMakeMChisto();     // reprocessing MC histos from KKC and Foam
  ReMakeMChisto2();    // reprocessing MC histos from KKC and Foam
//========== PLOTTING ==========
  FigAfb2();
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


