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
//#include "TH2.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"
#include "TFile.h"

#include "HisNorm.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
TFile DiskFileA("../ProdRun/work1/histo.root");
//
TFile DiskFileB("Plot1.root","RECREATE","histograms");
//=============================================================================

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
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_weight") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvTrue") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_nPhot") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvBES") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosTheta") );
  //
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_r1r2") );
  //
  cout<<"----------------DiskFileA.GetListOfKeys--------------------"<<endl;
  DiskFileA.GetListOfKeys()->Print();
  cout<<"----------------DiskFileA.ShowStreamerInfo-------------------------"<<endl;
  DiskFileA.ShowStreamerInfo();

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
  TH1D *hst_weight    = (TH1D*)DiskFileA.Get("hst_weight");
  TH1D *hst_nPhot     = (TH1D*)DiskFileA.Get("hst_nPhot");
  TH1D *hst_vvBES     = (TH1D*)DiskFileA.Get("hst_vvBES");
  TH1D *hst_vvTrue    = (TH1D*)DiskFileA.Get("hst_vvTrue");
  TH1D *hst_CosTheta  = (TH1D*)DiskFileA.Get("hst_CosTheta");
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
  //==========plot1==============
  cFigInfo->cd(1);
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

  ///==========plot2==============
  cFigInfo->cd(2);
  //-----------------------------
  gPad->SetLogy(); // !!!!!!
  //HST = hst_CosTheta;
  HST = hst_vvTrue;

  HST->SetTitle(0);
  HST->SetStats(0);
  HST->SetLineColor(4);
  //HST->SetMinimum(0);

  HST->DrawCopy("h");

  CaptT->DrawLatex(0.10,0.95,"d#sigma/dv (Ceex2) ");
  ///==========plot3==============
  cFigInfo->cd(3);
  gPad->SetLogy(); // !!!!!!
  ///
  CaptT->DrawLatex(0.10,0.95,"weight distribution");
  hst_weight->GetXaxis()->SetLabelSize(0.05);
  hst_weight->DrawCopy("h");
 ///==========plot4==============
  cFigInfo->cd(4);
  ///
  //gPad->SetLogy(); // !!!!!!
  hst_vvBES->SetTitle(0);
  hst_vvBES->DrawCopy("h");
  CaptT->DrawLatex(0.10,0.95,"Beam Energy Spread");
  CaptE->DrawLatex(0.30,0.85, capt1);
  //
  cFigInfo->cd();

  cFigInfo->SaveAs("cFigInfo.pdf");
}//FigInfo

///////////////////////////////////////////////////////////////////////////////////
void FigBES2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigBES2 =========================== "<<endl;
  TH2D *sca_r1r2    = (TH2D*)DiskFileA.Get("sca_r1r2");

  //////////////////////////////////////////////
  TString OptSurf;
  //OptSurf="      "; // 2D scatergram, points
  //OptSurf="col"; // 2D histogram, color
  //OptSurf="colz"; // 2D kolorowe paski, ze skala
  //OptSurf="surf1 "; // 3D surface color
  //OptSurf="lego2 "; // 3D histogram color
  OptSurf="surf3 "; // 3D histogram, z plotem "na dachu"
  //OptSurf="surf2z"; // 3D kolorowe paski, ze skala
  //OptSurf="surf2 "; // 3D kolorowe paski bez skali
  //OptSurf="surf4 "; // 3D gladka powierchnia
  //-------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
//  CaptT->SetTextSize(0.060);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigBES2 = new TCanvas("cFigBES2","FigBES2: general info ",  gXcanv,  gYcanv,    800,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigBES2->SetFillColor(10);

  cFigBES2->cd();
//  gPad->SetLogz(); // !!!!!!

  gPad->SetTheta(25);
  gPad->SetPhi( -70);
  TH2D *Scat1= sca_r1r2;

  Scat1->SetTitle(0);
  Scat1->SetStats(0);

  Scat1->GetYaxis()->CenterTitle();
  Scat1->GetYaxis()->SetTitleOffset(1.4);
  Scat1->GetYaxis()->SetTitleSize(0.035);
  Scat1->GetYaxis()->SetNdivisions(5);
  Scat1->GetYaxis()->SetTitle("r_{1}");
  Scat1->GetXaxis()->CenterTitle();
  Scat1->GetXaxis()->SetTitleOffset(1.4);
  Scat1->GetXaxis()->SetNdivisions(5);
  Scat1->GetXaxis()->SetTitle("r_{2}");

  sca_r1r2->DrawCopy(OptSurf);

//  CaptT->DrawLatex(0.10,0.95,"Gaussian Beam Enegy Spread");
  CaptT->DrawLatex(0.10,0.95,"Beamstrahlung Energy Spread");

  cFigBES2->SaveAs("cFigBES2.pdf");
}//FigBES2

///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  /*
  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  int Nodes;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
 */
  HistNormalize();     // Renormalization of MC histograms
  //========== PLOTTING ==========
  FigBES2();
  FigInfo();
 //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
//
//  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
//  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
//
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
