#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include <math.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TSystem.h>

#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
//#include "TH2.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"
#include "TFile.h"

#include "HisNorm.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
TFile *DiskFileA;
TFile *DiskFileF;
TFile *DiskFileB;

TString FileA= "../ProdRun/work1/histo.root";

TString FileF= "../ProdRun/workFoam/histo.root";

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot; // from KKMC run
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
//
float  gXcanv = 0, gYcanv = 0, gDcanv = 30;
///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA->ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_weight") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vvTrue") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_nPhot") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_CosTheta") );

  //[[[[[[[[[[[[[[[[[[[[[
  //cout<<"----------------DiskFileA->GetListOfKeys--------------------"<<endl;
  //DiskFileA->GetListOfKeys()->Print();
  //cout<<"----------------DiskFileA->ShowStreamerInfo-------------------------"<<endl;
  //DiskFileA->ShowStreamerInfo();
  //cout<<"----------------DiskFileF->ls--------------------"<<endl;
  //DiskFileF->ls("");
  //cout<<"----------------DiskFileF->ShowStreamerInfo-------------------------"<<endl;
  //DiskFileF->ShowStreamerInfo();
  //]]]]]]]]]]]]]]]]]]]]]]
//----------------------
  TH1D *HST_FOAM_NORMA4 = (TH1D*)DiskFileF->Get("HST_FOAM_NORMA4");
  HisNorm1(HST_FOAM_NORMA4, (TH1D*)DiskFileF->Get("HST_weight4") );
  HisNorm1(HST_FOAM_NORMA4, (TH1D*)DiskFileF->Get("HST_vv_eex3") );
//------------------------------
  TH1D *HST_FOAM_NORMA6 = (TH1D*)DiskFileF->Get("HST_FOAM_NORMA6");
  HisNorm1(HST_FOAM_NORMA6, (TH1D*)DiskFileF->Get("HST_weight6") );


}// HistNormalize

///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;
  char capt1[100];
  double CMSene=100.0;
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  //
  TH1D *hst_weight    = (TH1D*)DiskFileA->Get("hst_weight");
  TH1D *hst_nPhot     = (TH1D*)DiskFileA->Get("hst_nPhot");
  TH1D *hst_CosTheta  = (TH1D*)DiskFileA->Get("hst_CosTheta");
  //
  TH1D *HST_weight4    = (TH1D*)DiskFileF->Get("HST_weight4");
  TH1D *HST_weight6    = (TH1D*)DiskFileF->Get("HST_weight6");
  //------------------------------------------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
//  CaptT->SetTextSize(0.060);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo: general info ",  gXcanv,  gYcanv,    1000,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  2);
  ///==========plot1==============
  cFigInfo->cd(1);
  gPad->SetLogy(); // !!!!!!
  ///
  CaptT->DrawLatex(0.10,0.95,"weight distribution");
  hst_weight->GetXaxis()->SetLabelSize(0.05);
  hst_weight->DrawCopy("h");

  //==========plot2==============
  cFigInfo->cd(2);
  gPad->SetLogy(); // !!!!!!
  TH1D *HST = hst_nPhot;
  HST->SetStats(0);
  HST->SetMinimum( 1e-3*HST->GetMaximum());
  HST->SetTitle(0);
  HST->SetLineColor(kBlue);
  HST->GetXaxis()->SetLabelSize(0.06);
  HST->DrawCopy("h");
  CaptT->DrawLatex(0.10,0.93,"No. of #gamma's. tagged (red) and All (blue)");
  CaptE->DrawLatex(0.70,0.85, capt1);

  ///==========plot3==============
  cFigInfo->cd(3);
  //-----------------------------
  gPad->SetLogy(); // !!!!!!
  HST = HST_weight4;

  //HST->SetTitle(0);
  //HST->SetStats(0);
  HST->SetLineColor(4);
  //HST->SetMinimum(0);

  HST->DrawCopy("h");
  HST_weight6->DrawCopy("hsame");

  //=================
  cFigInfo->cd();

  cFigInfo->SaveAs("cFigInfo.pdf");
}//FigInfo

///////////////////////////////////////////////////////////////////////////////////
void FigVplot()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVplot =========================== "<<endl;
  char capt1[100];
  double CMSene=100.0;
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  //
  TH1D *hst_vvTrue    = (TH1D*)DiskFileA->Get("hst_vvTrue");
  TH1D *HST_vv_eex3  =  (TH1D*)DiskFileF->Get("HST_vv_eex3");
  //------------------------------------------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
//  CaptT->SetTextSize(0.060);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVplot = new TCanvas("cFigVplot","FigVplot: general info ",  gXcanv,  gYcanv,    1000,  500);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigVplot->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVplot->Divide( 2,  1);

  //==========plot1==============
  cFigVplot->cd(1);
  gPad->SetLogy(); // !!!!!!
  TH1D *HST; //
  HST = hst_vvTrue; //
  //HST = HST_vv_eex3;

  //HST->SetTitle(0);
  //HST->SetStats(0);
  HST->SetLineColor(4);
  //HST->SetMinimum(0);
  HST->DrawCopy("h");

  HST_vv_eex3->DrawCopy("hsame");

//==========plot1==============
  cFigVplot->cd(2);
//  TH1D *RAT_Hh  = HstRatio("RAT_Hh",  hst_Mll_ceex2_B,  hst_Mll_ceex2_A,  kBlack ); //
  TH1D *RAT_Hh  = HstRatio("RAT_Hh",  hst_vvTrue,  HST_vv_eex3,  kBlack ); //

  HST =RAT_Hh;
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->SetLineColor(4);
  HST->SetMinimum(1-0.05);
  HST->SetMaximum(1+0.05);
  HST->DrawCopy("h");

  TH1D *hOne = (TH1D*)RAT_Hh->Clone("hOne");  // zero line
  for(int i=1; i <= hOne->GetNbinsX() ; i++) { hOne->SetBinContent(i, 1); hOne->SetBinError(i, 0);}
  hOne->DrawCopy("hsame");

  //=================
  cFigVplot->cd();

  cFigVplot->SaveAs("cFigVplot.pdf");
}//FigVplot


///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

  gSystem->Load("../SRCee/.libs/libKKee.so");     // needed ???
  gSystem->Load("../SRCee/.libs/libKKfm.so");     // NEEDED! why?

  DiskFileA = new TFile(FileA);
  DiskFileF = new TFile(FileF);
//
  DiskFileB = new TFile("Plot1.root","RECREATE","histograms");

  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  /*
  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  int Nodes;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
 */
  HistNormalize();     // Renormalization of MC histograms
  //========== PLOTTING ==========
  FigInfo();
  FigVplot();
 //++++++++++++++++++++++++++++++++++++++++
  DiskFileA->ls();
//
//  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
//  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
//
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
