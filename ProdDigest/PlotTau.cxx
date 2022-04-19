// make PlotTau-run
//

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
//TFile DiskFileA("../ProdRun/workTau/histo.root");
//
TFile DiskFileA(   "../ProdRun/workTau/histo.root_1pi_6G");
TFile DiskFileA_P0("../ProdRun/workTau/histo.root_1pi_Phts0_1G");
TFile DiskFileB(   "../ProdRun/workTau/histo.root_1pi_Born_400M");

TFile DiskFileS0(   "../ProdRun/workTau/histo.root_1pi_Spin0_400M");
TFile DiskFileS0_P0("../ProdRun/workTau/histo.root_1pi_Spin0_Phts0_10G");
TFile DiskFileS0_B( "../ProdRun/workTau/histo.root_1pi_Spin0_Born_400M");
//
TFile DiskFileU("Plot1.root","RECREATE","histograms");
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
void PlotSame2(TH1D *HST, double &xcapt, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
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
}// PlotSame2


///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  // Normalization to [nb]
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_WtMain") );
  //HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_x1pi") );
  ////////////////////////////////////
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_x1c1") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_x1x2") );
  ////////////////////////////////////
  // Normalization to ONE
  // Pi energy, spin ON
  HisNorm0(1.0, (TH1D*)DiskFileA.Get("hst_x1pi") );     // ISR+FSR+PHTS
  HisNorm0(1.0, (TH1D*)DiskFileA_P0.Get("hst_x1pi") );  // ISR+FSR
  HisNorm0(1.0, (TH1D*)DiskFileB.Get("hst_x1pi") );     // Born PHTS=0?
  // Pi energy, spin OFF
  HisNorm0(1.0, (TH1D*)DiskFileS0.Get("hst_x1pi") );    // ISR+FSR+PHTS
  HisNorm0(1.0, (TH1D*)DiskFileS0_P0.Get("hst_x1pi") ); // ISR+FSR
  HisNorm0(1.0, (TH1D*)DiskFileS0_B.Get("hst_x1pi") );  // Born
  //
  cout<<"----------------DiskFileA.GetListOfKeys--------------------"<<endl;
  DiskFileA.GetListOfKeys()->Print();
  cout<<"----------------DiskFileA.ShowStreamerInfo-------------------------"<<endl;
  DiskFileA.ShowStreamerInfo();

}// HistNormalize


///////////////////////////////////////////////////////////////////////////////////
void FigTau1()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTau1 =========================== "<<endl;
  TH2D *sca_x1c1    = (TH2D*)DiskFileA.Get("sca_x1c1");

  //////////////////////////////////////////////
  TString OptSurf;
  //OptSurf="      "; // 2D scatergram, points
  //OptSurf="col"; // 2D histogram, color
  //OptSurf="colz"; // 2D kolorowe paski, ze skala
  //OptSurf="surf1 "; // 3D surface color
  OptSurf="lego2 "; // 3D histogram color
  //OptSurf="surf3 "; // 3D histogram, z plotem "na dachu"
  //OptSurf="surf2z"; // 3D kolorowe paski, ze skala
  //OptSurf="surf2 "; // 3D kolorowe paski bez skali
  //OptSurf="surf4 "; // 3D gladka powierchnia
  //-------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.040);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigTau1 = new TCanvas("cFigTau1","FigTau1: general info ",  gXcanv,  gYcanv,    800,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigTau1->SetFillColor(10);

  cFigTau1->cd();
//  gPad->SetLogz(); // !!!!!!

  gPad->SetTheta(12);
  gPad->SetPhi( -62);
  TH2D *Scat1= sca_x1c1;

  Scat1->SetTitle(0);
  Scat1->SetStats(0);

  Scat1->GetYaxis()->CenterTitle();
  Scat1->GetYaxis()->SetTitleOffset(1.8);
  Scat1->GetYaxis()->SetTitleSize(0.035);
  Scat1->GetYaxis()->SetNdivisions(5);
  Scat1->GetYaxis()->SetTitle("x_{1}");
  Scat1->GetXaxis()->CenterTitle();
  Scat1->GetXaxis()->SetTitleOffset(1.8);
  Scat1->GetXaxis()->SetNdivisions(5);
  Scat1->GetXaxis()->SetTitle("cos(#theta)");
  Scat1->SetMinimum(0.0);

  Scat1->DrawCopy(OptSurf);

  CaptT->DrawLatex(0.10,0.95,"#frac{d#sigma }{ dx_{1} dcos(#theta)}");

  cFigTau1->SaveAs("cFigTau1.pdf");
//  cFigTau1->SaveAs("cFigTau1.jpg");
  cFigTau1->SaveAs("cFigTau1.png");
}//FigTau1

///////////////////////////////////////////////////////////////////////////////////
void FigTau2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTau2 =========================== "<<endl;
  TH2D *sca_x1x2    = (TH2D*)DiskFileA.Get("sca_x1x2");

  //////////////////////////////////////////////
  TString OptSurf;
  //OptSurf="      "; // 2D scatergram, points
  //OptSurf="col"; // 2D histogram, color
  //OptSurf="colz"; // 2D kolorowe paski, ze skala
  //OptSurf="surf1 "; // 3D surface color
  OptSurf="lego2 "; // 3D histogram color
  //OptSurf="surf3 "; // 3D histogram, z plotem "na dachu"
  //OptSurf="surf2z"; // 3D kolorowe paski, ze skala
  //OptSurf="surf2 "; // 3D kolorowe paski bez skali
  //OptSurf="surf4 "; // 3D gladka powierchnia
  //-------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
//  CaptT->SetTextSize(0.060);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigTau2 = new TCanvas("cFigTau2","FigTau2: general info ",  gXcanv,  gYcanv,    800,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigTau2->SetFillColor(10);

  cFigTau2->cd();
//  gPad->SetLogz(); // !!!!!!

  gPad->SetTheta(35);
  gPad->SetPhi( -70);
  TH2D *Scat1= sca_x1x2;

  Scat1->SetTitle(0);
  Scat1->SetStats(0);

  Scat1->GetYaxis()->CenterTitle();
  Scat1->GetYaxis()->SetTitleOffset(1.8);
  Scat1->GetYaxis()->SetTitleSize(0.035);
  Scat1->GetYaxis()->SetNdivisions(5);
  Scat1->GetYaxis()->SetTitle("x_{1}");
  Scat1->GetXaxis()->CenterTitle();
  Scat1->GetXaxis()->SetTitleOffset(1.8);
  Scat1->GetXaxis()->SetNdivisions(5);
  Scat1->GetXaxis()->SetTitle("x_{2}");
  Scat1->SetMinimum(0.0);

  Scat1->DrawCopy(OptSurf);

  CaptT->DrawLatex(0.10,0.95,"Spin correlations in #pi energies");

  cFigTau2->SaveAs("cFigTau2.pdf");
//  cFigTau2->SaveAs("cFigTau2.jpg");
  cFigTau2->SaveAs("cFigTau2.png");
}//FigTau2


///////////////////////////////////////////////////////////////////////////////////
void FigX1pi()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigX1pi =========================== "<<endl;
  //
  TH1D *hst_x1pi       = (TH1D*)DiskFileA.Get("hst_x1pi");     // ISR+FSR+PHTS
  TH1D *hst_x1pi_phts0 = (TH1D*)DiskFileA_P0.Get("hst_x1pi");  // ISR+FSR
  TH1D *hst_x1pi_Born  = (TH1D*)DiskFileB.Get("hst_x1pi");     // Born
  //////
  TH1D *hst_x1pi_S0       = (TH1D*)DiskFileS0.Get("hst_x1pi");    // ISR+FSR+PHTS
  TH1D *hst_x1pi_S0_phts0 = (TH1D*)DiskFileS0_P0.Get("hst_x1pi"); // ISR+FSR
  TH1D *hst_x1pi_S0_Born  = (TH1D*)DiskFileS0_B.Get("hst_x1pi");  // Born
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigX1pi = new TCanvas("cFigX1pi","cFigX1pi ",  gXcanv,  gYcanv,    800,  800);
  //                                  Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigX1pi->SetFillColor(10);
//  cFigX1pi->Divide( 2,  1);
  ////////////////////////////////////////////////////////////////////////////////
  //==========plot1==============
  cFigX1pi->cd(1);
  TH1D *HST; //
  HST = hst_x1pi; //
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->GetXaxis()->SetTitle("x_{1}");
  HST->GetXaxis()->SetLabelSize(0.04);
  HST->SetMinimum(0.8);
  HST->SetMaximum(1.3);
  HST->DrawCopy("h");
  CaptT->DrawLatex(0.06,0.95, "dP/dx_{1}");
  double ycapt = 0.90; double xcapt=0.40;
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(xcapt,ycapt, "e^{+}e^{-} -> #tau^{+} #tau^{-},   #tau^{#pm} -> #pi^{#pm}#nu#bar{#nu}");
  CaptT->DrawLatex(xcapt+0.40,ycapt,gTextEne);
  PlotSame2(hst_x1pi,        xcapt, ycapt, kBlack,   0.10, "(a)", " ISR{#tau}+FSR{#tau}+FSR{#pi} ");
  hst_x1pi_phts0->SetLineStyle(2);
  PlotSame2(hst_x1pi_phts0,  xcapt, ycapt, kBlue,    0.20, "(b)", " ISR{#tau}+FSR{#tau} (dashed)");
  PlotSame2(hst_x1pi_Born,   xcapt, ycapt, kMagenta, 0.30, "(c)", " QED off ");
///////////
  PlotSame2(hst_x1pi_S0,     xcapt, ycapt, kBlack,    0.10, "(f)", " ISR{#tau}+FSR{#tau}+FSR{#pi} Spin off");
  hst_x1pi_S0_phts0->SetLineStyle(2);
  PlotSame2(hst_x1pi_S0_phts0,  xcapt, ycapt, kBlue,  0.20, "(g)", " ISR{#tau}+FSR{#tau} Spin off ");
  PlotSame2(hst_x1pi_S0_Born,   xcapt, ycapt, kRed,   0.30, "(h)", " QED off, Spin off ");
////////////////
  cFigX1pi->SaveAs("cFigX1pi.pdf");

}//FigX1pi



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.cd();
  DiskFileA.ls();
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
  FigTau1();
  FigTau2();
  FigX1pi();
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
