//////////////////////////////////////////////////////////////////////
//    make PlotKKsem-run
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

#include "HisReMake.h"
#include "HisNorm.h"
#include "KKplot.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
////  *** input from KKMC run
//TFile DiskFileA("../workKKMC/histo.root_88GeV_11G"); //
TFile DiskFileA("../workKKMC/histo.root_95GeV_26G");  // last

TFile DiskFileB("RhoSemi.root","RECREATE","histograms");

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot, gNevTot2; // from KKMC and KKfoam MC runs (histograms)
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3, kCyan2=kCyan+2;
//
//int    gNbMax   =45;          // for 100bins, gCosTheta = 45/50=0.90
//double gCosTheta=0.90;        // to be synchronized with gNbMax
//
int    gNbMax   =0;          // for 100bins, default=0 for gCosTheta = 1.00
double gCosTheta=1.00;       // to be synchronized with gNbMax
//
int    gNbMax2=0;            // for 50 bins, default=0 for gCosTheta = 1.00
//
float  gXcanv = 50, gYcanv = 50;
KKplot LibSem("KKplot");
///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void PlotSame2(TH1D *HST, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
//  ycapt += -0.050;
  ycapt += -0.047;
  double xcapt = 0.50;
  CaptT->DrawLatex(xcapt,ycapt, opis);
  CaptT->DrawLatex(xcapt-0.07,ycapt, label);
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
// Template is for vmax<0.20 and 100 bins in cos(theta)
  TH1D *hstVtemplate = (TH1D*)DiskFileA.Get("hst_vT_Ceex2");
//------------------------------------------------------------------------
//   MuMu  Sigma(vmax) and AFB(vmax) with ulimited c=cos(theta)
//
  KeyDis = 303302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR3_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR3_FSR2");
  TH1D *afbv_ISR3_FSR2 =(TH1D*)hstVtemplate->Clone("afbv_ISR3_FSR2");
  LibSem.VVmake( vcum_ISR3_FSR2, afbv_ISR3_FSR2, KF, chak, KeyDis, KeyFob, gCosTheta);
  //
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  TH1D *afbv_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("afbv_ISR2_FSR2");
  LibSem.VVmake( vcum_ISR2_FSR2, afbv_ISR2_FSR2, KF, chak, KeyDis, KeyFob, gCosTheta);
  //
  KeyDis = 301301;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR1_FSR1 =(TH1D*)hstVtemplate->Clone("vcum_ISR1_FSR1");
  TH1D *afbv_ISR1_FSR1 =(TH1D*)hstVtemplate->Clone("afbv_ISR1_FSR1");
  LibSem.VVmake( vcum_ISR1_FSR1, afbv_ISR1_FSR1, KF, chak, KeyDis, KeyFob, gCosTheta);
  //
  KeyDis = 300300;        // ISR*FSR O(alf0)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR0_FSR0 =(TH1D*)hstVtemplate->Clone("vcum_ISR0_FSR0");
  TH1D *afbv_ISR0_FSR0 =(TH1D*)hstVtemplate->Clone("afbv_ISR0_FSR0");
  LibSem.VVmake( vcum_ISR0_FSR0, afbv_ISR0_FSR0, KF, chak, KeyDis, KeyFob, gCosTheta);
  //
  cout<<"================ KKsem MakeHisto ENDs ============================="<<endl;
  cout<<"==================================================================="<<endl;
//------------------------------------------------------------------------
//------------------------------------------------------------------------
}//  KKsemMakeHisto


///////////////////////////////////////////////////////////////////////////////////
void Afb_eex32()
{
//------------------------------------------------------------------------
  cout<<" ========================= Afb_eex32 =========================== "<<endl;
  //--------------------------
  // AFB with v_true and costhetaPL
  TH1D *afbv_ISR3_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR3_FSR2");    // total EEX3
  TH1D *afbv_ISR2_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");    // total EEX2
  TH1D *afbv_ISR1_FSR1    = (TH1D*)DiskFileB.Get("afbv_ISR1_FSR1");    // total EEX1
  TH1D *afbv_ISR0_FSR0    = (TH1D*)DiskFileB.Get("afbv_ISR0_FSR0");    // total EEX0

  //--------------------
  TH1D *HAfb_Diff_vT_EEX32   = HstDiff("HAfb_Diff_vT_EEX32", afbv_ISR3_FSR2, afbv_ISR2_FSR2,  kBlack);
  TH1D *HAfb_Diff_vT_EEX21   = HstDiff("HAfb_Diff_vT_EEX21", afbv_ISR2_FSR2, afbv_ISR1_FSR1,  kBlack);

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);

  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfb_eex32 = new TCanvas("cAfb_eex32","cAfb_eex32", gXcanv,  gYcanv,   1200,  600);
  //                                                 Name    Title            xoff,    yoff,  WidPix, HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfb_eex32->SetFillColor(10);
  cAfb_eex32->Divide( 2,  0);
  //********************************************************************************
  cAfb_eex32->cd(1);

  TH1D *HST = afbv_ISR2_FSR2;
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->GetXaxis()->SetTitle("v_{max}");
  HST->SetMinimum(0.15);HST->SetMaximum(0.35);
  HST->DrawCopy("h");

  double ycapt =0.63;
  PlotSame2(afbv_ISR3_FSR2,       ycapt,  kBlue,    0.10, "(a)", "KKsem EEX3");
  PlotSame2(afbv_ISR2_FSR2,       ycapt,  kBlack,   0.12, "(b)", "KKsem EEX2");
  PlotSame2(afbv_ISR1_FSR1,       ycapt,  kBlack,   0.14, "(c)", "KKsem EEX1");
  PlotSame2(afbv_ISR0_FSR0,       ycapt,  kPine,    0.16, "(d)", "KKsem EEX0");

  CaptT->DrawLatex(0.20, 0.95, " KKsem EEX without IFI,  |cos(#theta)| < 1, ");
  //********************************************************************************
  cAfb_eex32->cd(2);
  HST = HAfb_Diff_vT_EEX32;

  HST->SetTitle(0); HST->SetStats(0);
  HST->GetXaxis()->SetTitle("v_{max}");
  HST->SetMaximum( 1.0e-4); HST->SetMinimum(-1.0e-4);
  HST->DrawCopy("h");

  ycapt =0.33;
  PlotSame2(HAfb_Diff_vT_EEX32,       ycapt,  kBlack,   0.18, "(a)", "EEX3-EEX2");
  PlotSame2(HAfb_Diff_vT_EEX21,       ycapt,  kBlack,   0.18, "(b)", "EEX2-EEX1");

  TH1D *hZero      = (TH1D*)HST->Clone("hZero");      // zero line
  TH1D *hZeroPlus  = (TH1D*)HST->Clone("hZeroPlus");  //
  TH1D *hZeroMinus = (TH1D*)HST->Clone("hZeroMinus"); //
  for(int i=1; i <= hZero->GetNbinsX() ; i++) {
    hZero->SetBinContent(i, 0);          hZero->SetBinError(i, 0);
    hZeroPlus->SetBinContent(i,  3e-5);  hZeroPlus->SetBinError(i, 0);
    hZeroMinus->SetBinContent(i,-3e-5);  hZeroMinus->SetBinError(i, 0);
    }// for i
  hZeroPlus->SetLineStyle(9); hZeroMinus->SetLineStyle(9);
  hZero->DrawCopy("hsame"); hZeroPlus->DrawCopy("hsame"); hZeroMinus->DrawCopy("hsame");
  CaptT->DrawLatex(0.12,0.63," #delta#alpha/#alpha = 10^{-4}");

  CaptT->DrawLatex(0.20, 0.95, "#delta A_{FB}(v_{max})");

  cAfb_eex32->SaveAs("cAfb_eex32.pdf");

  cAfb_eex32->cd();
//
}// Afb_eex32

///////////////////////////////////////////////////////////////////////////////////
void FigTempl()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTempl =========================== "<<endl;
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

//////////////////////////////////////////////////////////////////////////
// ========= Preparing plots ==========
  DiskFileB.cd();

  KKsemMakeHisto();    // prepare histos from KKsem
//========== PLOTTING ==========
  Afb_eex32();
//  FigTempl();
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
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  cout<< "FOAM: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<<endl;
  //cout<<"------------------------------end---------------------------------"<<endl;
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


