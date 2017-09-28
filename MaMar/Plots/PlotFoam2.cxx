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

// Sept. 2017 runs
TFile DiskFileA("../workKKMC/histo.root_95GeV_55M"); //

// August2017 runs
//TFile DiskFileA("../workKKMC/histo.root_10GeV_5.7G"); //
//TFile DiskFileA("../workKKMC/histo.root_88GeV_2.1G"); //
//TFile DiskFileA("../workKKMC/histo.root_95GeV_16G");
//TFile DiskFileA("../workKKMC/histo.root_91GeV_9G"); ///????

////  *** FOAM 5dim
TFile DiskFileF("../workFOAM/histo.root"); // current
//TFile DiskFileF("../workFOAM/histo.root_10GeV_37G_vmax0.2");
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
int  kGold=92, kBrune=46, kPine=71;

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
  double xcapt = 0.40;
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
}// PlotSame




///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA5");

  HisNorm2(HST_FOAM_NORMA3, (TH2D*)DiskFileF.Get("SCT_xc_Ceex2n") );
  HisNorm2(HST_FOAM_NORMA3, (TH2D*)DiskFileF.Get("SCT_xc_EEX2") );
  HisNorm2(HST_FOAM_NORMA5, (TH2D*)DiskFileF.Get("SCT_xc_Ceex2") );

  TH2D *SCT_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2n");  // FOAM small range x<0.20
  TH2D *SCT_xc_EEX2   = (TH2D*)DiskFileF.Get("SCT_xc_EEX2");    // FOAM small range x<0.20
  TH2D *SCT_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2");   // FOAM small range x<0.20

  // KKFoam1: sigma(vmax) and AFB(vmax) from scat. vmax<0.2, 100 bins in ctheta
  // gNbMax=45;           // cosThetaMax = 45/50=0.90 Now global variable

  TH1D  *Htot2_xmax_Ceex2n   = HstProjV("Htot2_xmax_Ceex2n",SCT_xc_Ceex2n,gNbMax);
  TH1D  *Hafb2_xmax_Ceex2n   = HstProjA("Hafb2_xmax_Ceex2n",SCT_xc_Ceex2n,gNbMax);

  TH1D  *Htot2_xmax_Ceex2   = HstProjV("Htot2_xmax_Ceex2",SCT_xc_Ceex2,gNbMax);
  TH1D  *Hafb2_xmax_Ceex2   = HstProjA("Hafb2_xmax_Ceex2",SCT_xc_Ceex2,gNbMax);

  TH1D  *Htot2_xmax_EEX2    = HstProjV("Htot2_xmax_EEX2", SCT_xc_EEX2,gNbMax);
  TH1D  *Hafb2_xmax_EEX2    = HstProjA("Hafb2_xmax_EEX2", SCT_xc_EEX2,gNbMax);

  ///****************************************************************************************
  //     Calculating AFB from <cos(theta)>
  TH1D *Afb5T_Ceex2  = HstTildeAFB("Afb5T_Ceex2", (TH1D*)DiskFileF.Get("HST5_xc_Ceex2"),
		                                          (TH1D*)DiskFileF.Get("HST5_xx_Ceex2") );

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
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2n") );

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex0") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex0n") );

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("hst_vT_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("hst_vTcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("hst_vT_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("hst_vTcPL_Ceex2n") );

////////////////////////////////////////////////////////////////////
// Pure KKMC reprocessing from scattergram with restricted vmax<0.2
////////////////////////////////////////////////////////////////////

// Wide range, vmax<1.
  TH2D *sct_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2");
  TH2D *sct_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n");
  TH2D *sct_vTcPR_EEX2   = (TH2D*)DiskFileA.Get("sct_vTcPR_EEX2");
  TH2D *sct_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2");
  TH2D *sct_vTcPL_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2n");

  TH2D *sct_vTcPL_Ceex0  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex0");
  TH2D *sct_vTcPL_Ceex0n = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex0n");

  cout<<"ReMakeMChisto2 [2]"<<endl;
///****************************************************************************************
/// Distributions of v=vTrue<vmax<0.20, c=cos(theta) with 100 bins
//gNbMax=45;         // cosThetaMax = 45/50=0.90 Now global variable
// KKMC IFI on
  TH1D  *HTot2_vTcPR_EEX2   = HstProjV("HTot2_vTcPR_EEX2",sct_vTcPR_EEX2,gNbMax);
  TH1D  *HAfb2_vTcPR_EEX2   = HstProjA("HAfb2_vTcPR_EEX2",sct_vTcPR_EEX2,gNbMax);
// KKMC IFI off
  TH1D  *HTot2_vTcPR_Ceex2  = HstProjV("HTot2_vTcPR_Ceex2",sct_vTcPR_Ceex2,gNbMax);
  TH1D  *HAfb2_vTcPR_Ceex2  = HstProjA("HAfb2_vTcPR_Ceex2",sct_vTcPR_Ceex2,gNbMax);
// KKMC IFI off
  TH1D  *HTot2_vTcPR_Ceex2n = HstProjV("HTot2_vTcPR_Ceex2n",sct_vTcPR_Ceex2n,gNbMax);
  TH1D  *HAfb2_vTcPR_Ceex2n = HstProjA("HAfb2_vTcPR_Ceex2n",sct_vTcPR_Ceex2n,gNbMax);
// KKMC IFI off
  TH1D  *HTot2_vTcPL_Ceex2  = HstProjV("HTot2_vTcPL_Ceex2",sct_vTcPL_Ceex2,gNbMax);
  TH1D  *HAfb2_vTcPL_Ceex2  = HstProjA("HAfb2_vTcPL_Ceex2",sct_vTcPL_Ceex2,gNbMax);
// KKMC IFI off
  TH1D  *HTot2_vTcPL_Ceex2n = HstProjV("HTot2_vTcPL_Ceex2n",sct_vTcPL_Ceex2n,gNbMax);
  TH1D  *HAfb2_vTcPL_Ceex2n = HstProjA("HAfb2_vTcPL_Ceex2n",sct_vTcPL_Ceex2n,gNbMax);

// KKMC IFI on
  TH1D  *HTot2_vTcPL_Ceex0  = HstProjV("HTot2_vTcPL_Ceex0",sct_vTcPL_Ceex0,gNbMax);
  TH1D  *HAfb2_vTcPL_Ceex0  = HstProjA("HAfb2_vTcPL_Ceex0",sct_vTcPL_Ceex0,gNbMax);
// KKMC IFI on
  TH1D  *HTot2_vTcPL_Ceex0n = HstProjV("HTot2_vTcPL_Ceex0n",sct_vTcPL_Ceex0n,gNbMax);
  TH1D  *HAfb2_vTcPL_Ceex0n = HstProjA("HAfb2_vTcPL_Ceex0n",sct_vTcPL_Ceex0n,gNbMax);


///****************************************************************************************
//      AFB fron <cos(theta)>
  TH1D *AfbT_Ceex2  = HstTildeAFB("AfbT_Ceex2", (TH1D*)DiskFileA.Get("hst_vTcPL_Ceex2"),
		                                        (TH1D*)DiskFileA.Get("hst_vT_Ceex2")   );
  TH1D *AfbT_Ceex2n = HstTildeAFB("AfbT_Ceex2n",(TH1D*)DiskFileA.Get("hst_vTcPL_Ceex2n"),
		                                        (TH1D*)DiskFileA.Get("hst_vT_Ceex2n")  );
  cout<<"================ ReMakeMChisto2 ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto2


///////////////////////////////////////////////////////////////////////////////////
void FigAfb3a()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb3a =========================== "<<endl;

  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n"); // KKMC[PR]
  TH1D *HAfb2_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2");  // KKMC[PR]
  //
  TH1D *HAfb2_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2");  // KKMC[PL]

  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // KKFoam scat.
  TH1D *Hafb2_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2");   // KKFoam scat.

  TH1D *HST_PLBhard =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_PLBhard");
  LibSem.Ord1fill(HST_PLBhard,105);                           // IFI only, hard part, PLB219
  HST_PLBhard->SetLineColor(kCyan);

  // This HST_PL is temporrary insert
  double alfinv  = 137.035989;
  double alfpi   = 1/alfinv/3.1415926535897932;
  TH1D *HST_PL =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_PL");
  HST_PL->SetLineColor(kYellow);
  int Nbin    = HST_PL->GetNbinsX();
  double vmax = HST_PL->GetXaxis()->GetXmax();
  for(int i=1; i <= Nbin ; i++) {
	  double vv = (i*vmax)/Nbin;
	  // A_FB from PLB219,p103, pure gamma exch.
	  double afb = 3e0/2e0* alfpi*( 3*vv+log(1-vv/2) ); // only gamma
      double KFi=11; int KFf=13;
	  HST_PL->SetBinContent(i, afb);
	  HST_PL->SetBinError(i, 0);
  }// i

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
//
  //*****************************************************************************
  TCanvas *cFigAfb3a = new TCanvas("cFigAfb3a","FigAfb3a", 50, 50,   600, 600);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  cFigAfb3a->SetFillColor(10);
  cFigAfb3a->cd();

  TH1D *Hst21_diff   = HstDiff("Hst21_diff",    HAfb2_vTcPR_Ceex2,  HAfb2_vTcPR_Ceex2n, kBlack);
  TH1D *HST21_diff   = HstDiff("HST21_diff",    Hafb2_xmax_Ceex2,   Hafb2_xmax_Ceex2n,  kMagenta);
  TH1D *HstKF_diff   = HstDiff("HstKF_diff",    HAfb2_vTcPR_Ceex2,  Hafb2_xmax_Ceex2,   kGreen);
  TH1D *HstKFn_diff  = HstDiff("HstKFn_diff",   HAfb2_vTcPR_Ceex2n, Hafb2_xmax_Ceex2n,  kBlue);
  TH1D *HstPL_diff   = HstDiff("HstPL_diff",    HAfb2_vTcPL_Ceex2,  Hafb2_xmax_Ceex2,   kRed);

  Hst21_diff->SetStats(0);
  Hst21_diff->SetTitle(0);
  if( fabs(gCMSene -95.0) < 1.0) { Hst21_diff->SetMinimum(-0.004);  Hst21_diff->SetMaximum( 0.006);}  // 95GeV
  if( fabs(gCMSene -88.0) < 1.0) { Hst21_diff->SetMinimum(-0.002);  Hst21_diff->SetMaximum( 0.008);}  // 88GeV
  if( fabs(gCMSene -91.0) < 1.0) { Hst21_diff->SetMinimum(-0.0005); Hst21_diff->SetMaximum( 0.0035);} // 91GeV
  if( fabs(gCMSene -10.0) < 1.0) { Hst21_diff->SetMinimum(-0.002);  Hst21_diff->SetMaximum( 0.006);}  // 10GeV
  Hst21_diff->GetXaxis()->SetTitle("v_{max}");
  Hst21_diff->DrawCopy("h");

  double ycapt = 0.90; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);
  CaptT->DrawLatex(0.06,0.95, "A_{FB}(v_{max}) ");

  PlotSame2(Hst21_diff,   ycapt, kBlack,     0.120, "(a)", "KKMC[vPR]: IFIon-IFIoff ");
  PlotSame2(HST21_diff,   ycapt, kMagenta,   0.140, "(b)", "KKFoam:    IFIon-IFIoff ");
  PlotSame2(HstKFn_diff,  ycapt, kBlue,      0.160, "(c)", "KKMC[vPR]-KKFoam:  IFIoff ");
  PlotSame2(HstKF_diff,   ycapt, kGreen,     0.120, "(d)", "KKMC[vPR]-KKFoam:  IFIon");
  PlotSame2(HstPL_diff,   ycapt, kRed,       0.140, "(e)", "KKMC[PL] -KKFoam:  IFIon");
  PlotSame2(HST_PLBhard,  ycapt, kCyan,      0.030, "(f)", "PLB219, IFI only,  hard");
  PlotSame2(HST_PL,       ycapt, kGold,      0.030, "(f)", "PLB219, IFI, gamma-gamma");

// zero line
  TH1D *hZero = (TH1D*)HAfb2_vTcPR_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  hZero->DrawCopy("hsame");

}//FigAfb3a



///////////////////////////////////////////////////////////////////////////////////
void FigAfb3b()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb3b =========================== "<<endl;
  //
  TH1D *HAfb2_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2");  // KKMC[PL] CCEX2
  TH1D *HAfb2_vTcPL_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2n"); // KKMC[PL] CCEX2

  TH1D *HAfb2_vTcPL_Ceex0  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex0");  // KKMC[PL] CCEX0
  TH1D *HAfb2_vTcPL_Ceex0n = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex0n"); // KKMC[PL] CCEX0

  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // KKFoam scat.
  TH1D *Hafb2_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2");   // KKFoam scat.

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
//
  //*****************************************************************************
  TCanvas *cFigAfb3b = new TCanvas("cFigAfb3b","FigAfb3b", 150, 100,   600, 600);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  cFigAfb3b->SetFillColor(10);
  cFigAfb3b->cd();

  TH1D *HstC2_diff = HstDiff("HstC2_diff",    HAfb2_vTcPL_Ceex2,  HAfb2_vTcPL_Ceex2n, kBlack);  // KK_Ceex2:IFIon-IFIoff
  TH1D *HstC0_diff = HstDiff("HstC0_diff",    HAfb2_vTcPL_Ceex0,  HAfb2_vTcPL_Ceex0n, kBlack);  // KK_Ceex0:IFIon-IFIoff
  TH1D *HST21_diff = HstDiff("HST21_diff",    Hafb2_xmax_Ceex2,   Hafb2_xmax_Ceex2n,  kMagenta);// KKFoam5: IFIon-IFIoff
  TH1D *HstPL_diff = HstDiff("HstPL_diff",    HAfb2_vTcPL_Ceex2,  Hafb2_xmax_Ceex2,   kRed);    // KK_Ceex2-KKFoam:IFIon
  TH1D *HstXX_diff = HstDiff("HstXX_diff",    HstC0_diff,         HST21_diff,         kBlue);   //

  HST21_diff->SetStats(0);
  HST21_diff->SetTitle(0);
  if( fabs(gCMSene -95.0) < 1.0) { HST21_diff->SetMinimum(-0.004);  HST21_diff->SetMaximum( 0.006);}  // 95GeV
  if( fabs(gCMSene -88.0) < 1.0) { HST21_diff->SetMinimum(-0.002);  HST21_diff->SetMaximum( 0.008);}  // 88GeV
  if( fabs(gCMSene -91.0) < 1.0) { HST21_diff->SetMinimum(-0.0005); HST21_diff->SetMaximum( 0.0035);} // 91GeV
  if( fabs(gCMSene -10.0) < 1.0) { HST21_diff->SetMinimum(-0.002);  HST21_diff->SetMaximum( 0.006);}  // 10GeV
  HST21_diff->GetXaxis()->SetTitle("v_{max}");
  HST21_diff->DrawCopy("h");

  double ycapt = 0.90; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);
  CaptT->DrawLatex(0.06,0.95, "A_{FB}(v_{max}) ");

  PlotSame2(HstC2_diff,   ycapt, kBlack,     0.100, "(a)", "KK_Ceex2: IFIon-IFIoff ");
  PlotSame2(HstC0_diff,   ycapt, kGold,      0.100, "(b)", "KK_Ceex0: IFIon-IFIoff ");
  PlotSame2(HST21_diff,   ycapt, kMagenta,   0.120, "(c)", "KKFoam5:  IFIon-IFIoff ");
  PlotSame2(HstPL_diff,   ycapt, kRed,       0.140, "(d)", "KK_Ceex2 - KKFoam:  IFIon");
  PlotSame2(HstXX_diff,   ycapt, kBlue,      0.040, "(e)", "Special: (b)-(c)");

  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev2); ycapt += -0.04;

// zero line
  TH1D *hZero = (TH1D*)HAfb2_vTcPL_Ceex2->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  hZero->DrawCopy("hsame");

}//FigAfb3b




///////////////////////////////////////////////////////////////////////////////////
void FigAfb4()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb4 =========================== "<<endl;

  TH1D *HAfb2_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2");   // KKMC[PL]
  TH1D *HAfb2_vTcPL_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2n");  // KKMC[PL]

  TH1D *Hafb2_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2"); // KKFoam5 scat.
  TH1D *Afb5T_Ceex2        = (TH1D*)DiskFileB.Get("Afb5T_Ceex2");   // KKFoam5, AFB form <cos(theta)>


  TH1D *AfbT_Ceex2          =(TH1D*)DiskFileB.Get("AfbT_Ceex2");
  TH1D *AfbT_Ceex2n         =(TH1D*)DiskFileB.Get("AfbT_Ceex2n");

// zero line
  TH1D *hZero2 = (TH1D*)HAfb2_vTcPL_Ceex2->Clone("hZero2");  // zero line
  for(int i=1; i <= hZero2->GetNbinsX() ; i++) { hZero2->SetBinContent(i, 0); hZero2->SetBinError(i, 0);}

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  //*****************************************************************************
  TCanvas *cFigAfb4 = new TCanvas("cFigAfb4","FigAfb4", 250, 150,   600, 600);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  cFigAfb4->SetFillColor(10);
  cFigAfb4->cd();
  //=======================================================

  TH1D *HstPL0_diff = HstDiff("HstPL0_diff",  AfbT_Ceex2,   HAfb2_vTcPL_Ceex2,  kBlack);
  TH1D *HstPL1_diff = HstDiff("HstPL1_diff",  AfbT_Ceex2n,  HAfb2_vTcPL_Ceex2n, kBlue);
  TH1D *HstKF1_diff = HstDiff("HstKF1_diff",  Afb5T_Ceex2,  Hafb2_xmax_Ceex2,   kMagenta);
  TH1D *HstKFX_diff = HstDiff("HstKFX_diff",  AfbT_Ceex2,   Afb5T_Ceex2,        kRed);

//  HstPL1_diff->SetMinimum(-0.0004);  HstPL1_diff->SetMaximum( 0.0004);  // zoom
  HstPL1_diff->SetMinimum(-0.004);  HstPL1_diff->SetMaximum( 0.006);  // zoom
  HstPL1_diff->SetStats(0);
  HstPL1_diff->SetTitle(0);
  HstPL1_diff->GetXaxis()->SetTitle("v_{max}");
  HstPL1_diff->DrawCopy("h");

  double ycapt = 0.90; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);
  CaptT->DrawLatex(0.06,0.95, "A_{FB}(v_{max}) ");

  PlotSame2(HstPL1_diff,   ycapt, kBlue,      0.020, "(a)", "KKMC:    #tilde{A}_{FB} - A_{FB}, IFIoff ");
  PlotSame2(HstKF1_diff,   ycapt, kMagenta,   0.030, "(b)", "KKFoam:  #tilde{A}_{FB} - A_{FB}, IFIon ");
  PlotSame2(HstPL0_diff,   ycapt, kBlack,     0.040, "(c)", "KKMC:    #tilde{A}_{FB} - A_{FB}, IFIon ");
  PlotSame2(HstKFX_diff,   ycapt, kRed,       0.100, "(d)", "KKMC-KKFoam:    #tilde{A}_{FB},   IFIon ");

  hZero2->DrawCopy("hsame");

  cFigAfb4->cd();
  //================================================
}//FigAfb4


///////////////////////////////////////////////////////////////////////////////////
void FigAfb5()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb5 =========================== "<<endl;

  TH1D *HAfb2_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2");  // KKMC[PR]
  TH1D *HAfb2_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2");  // KKMC[PL]

// zero line
  TH1D *hZero2 = (TH1D*)HAfb2_vTcPR_Ceex2->Clone("hZero2");  // zero line
  for(int i=1; i <= hZero2->GetNbinsX() ; i++) { hZero2->SetBinContent(i, 0); hZero2->SetBinError(i, 0);}

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  //*****************************************************************************
  TCanvas *cFigAfb5 = new TCanvas("cFigAfb5","FigAfb5", 350, 200,   600, 600);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  cFigAfb5->SetFillColor(10);
  cFigAfb5->cd();
  //=======================================================
  TH1D *HstPL2_diff = HstDiff("HstPL2_diff",  HAfb2_vTcPR_Ceex2,  HAfb2_vTcPL_Ceex2, kBlack);

  HstPL2_diff->SetMinimum(-0.0004);  HstPL2_diff->SetMaximum( 0.0004);  // zoom
  HstPL2_diff->SetStats(0);
  HstPL2_diff->SetTitle(0);
  HstPL2_diff->GetXaxis()->SetTitle("v_{max}");
  HstPL2_diff->DrawCopy("h");

  hZero2->DrawCopy("hsame");

  CaptT->DrawLatex(0.22,0.95,"KKMC: A_{FB}[#theta_{PRD}] - A_{FB}[#theta_{PL}] ");

  CaptT->DrawLatex(0.50,0.75,gTextEne);

  cFigAfb5->cd();
  //================================================
}//FigAfb5




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
  FigAfb3a();
  FigAfb3b();
  FigAfb4();
  FigAfb5();
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


