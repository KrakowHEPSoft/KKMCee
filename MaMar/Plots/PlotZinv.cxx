//////////////////////////////////////////////////////////////////////
//    make PlotZinv-run
//////////////////////////////////////////////////////////////////////

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
#include "TH2.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"
#include "TFile.h"

#include "HisNorm.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
TFile DiskFileA("../workZinv/rmain.root");
//
TFile DiskFileB("Plot2.root","RECREATE","histograms");
//=============================================================================

//Double_t sqr( const Double_t x ){ return x*x;};
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot; // from KKMC run
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
//
float  gXcanv = 20, gYcanv = 20, gDcanv = 30;
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
}// PlotSame2


///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_nPhAll") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_nPhVis") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_LnThPhAll") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_LnThPhVis") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueMain") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotMain") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueMu") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvNuCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvNuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvNuCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evNuCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evNuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evNuCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex1n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex2n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex12n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex2ifi") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotNuel") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotNumu") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPLCeex2") );
  //
}


///////////////////////////////////////////////////////////////////////////////////
void FigNPhont()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNPhont =========================== "<<endl;
 //

  TH1D *hst_nPhAll     = (TH1D*)DiskFileA.Get("hst_nPhAll");
  TH1D *hst_nPhVis     = (TH1D*)DiskFileA.Get("hst_nPhVis");

  TH1D *hst_LnThPhAll  = (TH1D*)DiskFileA.Get("hst_LnThPhAll");
  TH1D *hst_LnThPhVis  = (TH1D*)DiskFileA.Get("hst_LnThPhVis");

////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
 ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigNPhont = new TCanvas("cFigNPhont","cFigNPhont", gXcanv, gYcanv,    1200, 600);
  //                                      Name    Title        xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  cFigNPhont->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigNPhont->Divide( 2,  0);
  //====================plot1========================
  cFigNPhont->cd(1);
  gPad->SetLogy(); // !!!!!!
  TH1D *Hst=hst_nPhAll;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->CenterTitle();
  //Hst->GetYaxis()->SetTitleSize(0.04);
  //Hst->GetYaxis()->SetTitle("d#sigma/dcos(#theta) [nb]");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("N");

  Hst->SetLineColor(kBlue);
  Hst->SetMinimum( 1e-3*Hst->GetMaximum());
  Hst->DrawCopy("h");

  CaptT->DrawLatex(0.10,0.94,"d#sigma/N;        KKMC  e^{+}e^{-} -> #nu+#bar{#nu}+N#gamma");
  double ycapt = 0.30;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  hst_nPhAll->SetLineWidth(2);
  hst_nPhVis->SetLineWidth(2);
  PlotSame2(hst_nPhAll,  ycapt, kBlue,  +2.0, "(a)", "All photons");
  PlotSame2(hst_nPhVis,  ycapt, kRed,   +1.0, "(b)", "Tagged photons");

  //====================plot2========================
  cFigNPhont->cd(2);

  Hst=hst_LnThPhAll;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("ln_{10} sin #theta_{#gamma}");

  //Hst->SetMinimum(1.0 -0.01); Hst->SetMaximum(1.0 +0.01);
  Hst->DrawCopy("h");
  CaptT->DrawLatex(0.10,0.94,"d#sigma/d(ln_{10} sin #theta_{#gamma})        e^{+}e^{-} -> #nu+#bar{#nu}+N#gamma");
  ycapt = 0.80; // starting value, to be decremented below
  CaptT->DrawLatex(0.40,ycapt,gTextEne);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev);  ycapt += -0.01;
  //
  PlotSame2(hst_LnThPhAll,  ycapt, kBlue,    -3.0, "(a)", "All #gamma's");
  PlotSame2(hst_LnThPhVis,  ycapt, kRed,     -0.5, "(b)", "Tagged #gamma's");
  //

  //================================================
  cFigNPhont->SaveAs("cFigNPhont.pdf");

}//FigNPhont


///////////////////////////////////////////////////////////////////////////////////
void FigVPhont()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVPhont =========================== "<<endl;
 //
  TH1D *hst_vTrueMain  = (TH1D*)DiskFileA.Get("hst_vTrueMain");
  TH1D *hst_vPhotMain  = (TH1D*)DiskFileA.Get("hst_vPhotMain");
  TH1D *hst_vTrueMu    = (TH1D*)DiskFileA.Get("hst_vTrueMu");
 //
  TH1D *hst_vvMuCeex2  = (TH1D*)DiskFileA.Get("hst_vvMuCeex2");

 ////////////////////////////////////////
   TLatex *CaptT = new TLatex();
   CaptT->SetNDC(); // !!!
   CaptT->SetTextSize(0.04);
 ///////////////////////////////////////////////////////////////////////////////
   TCanvas *cFigVPhont = new TCanvas("cFigVPhont","cFigVPhont", gXcanv, gYcanv,    1200, 600);
 //                                      Name    Title        xoff,yoff, WidPix,HeiPix
   gXcanv += 25, gYcanv += 25;
   cFigVPhont->SetFillColor(10);
////////////////////////////////////////////////////////////////////////////////
   cFigVPhont->Divide( 2,  0);
   //====================plot1========================
   cFigVPhont->cd(1);
   gPad->SetLogy(); // !!!!!!
   TH1D *Hst=hst_vTrueMain;
   Hst->SetStats(0);
   Hst->SetTitle(0);
   Hst->GetXaxis()->CenterTitle();
   Hst->GetXaxis()->SetTitleSize(0.04);
   Hst->GetXaxis()->SetTitle("N");

   Hst->SetLineColor(kBlue);
//   Hst->SetMinimum( 1e-3*Hst->GetMaximum());
   Hst->DrawCopy("h");

   CaptT->DrawLatex(0.10,0.94,"d#sigma/N;        KKMC  e^{+}e^{-} -> #nu+#bar{#nu}+N#gamma");
   double ycapt = 0.80;
   CaptT->DrawLatex(0.40, ycapt,gTextEne);
//
   PlotSame2(hst_vTrueMain,  ycapt, kBlue,  +2.0, "(a)", "All photons");
   PlotSame2(hst_vPhotMain,  ycapt, kRed,   +1.0, "(b)", "Tagged photons");

   //====================plot2========================
   cFigVPhont->cd(2);

   gPad->SetLogy(); // !!!!!!
   Hst=hst_vTrueMain;
   Hst->SetStats(0);
   Hst->SetTitle(0);
   Hst->GetXaxis()->CenterTitle();
   Hst->GetXaxis()->SetTitleSize(0.04);
   Hst->GetXaxis()->SetTitle("N");

   Hst->SetLineColor(kBlue);
//   Hst->SetMinimum( 1e-3*Hst->GetMaximum());
   Hst->DrawCopy("h");

   CaptT->DrawLatex(0.10,0.94,"d#sigma/N;        KKMC  e^{+}e^{-} -> #nu+#bar{#nu}+N#gamma");
   ycapt = 0.80;
   CaptT->DrawLatex(0.40, ycapt,gTextEne);
//
   PlotSame2(hst_vTrueMain,  ycapt, kBlue,  +2.0, "(a)", "All photons");
   PlotSame2(hst_vvMuCeex2,    ycapt, kRed, +1.0, "(b)", "muons tagged");
//   PlotSame2(hst_vTrueMu,    ycapt, kRed,   +1.0, "(b)", "muons untagged");

  //================================================
  cFigVPhont->SaveAs("cFigVPhont.pdf");

}//FigVPhont



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++

  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  int Nodes;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);

  HistNormalize();     // Renormalization of MC histograms
  //========== PLOTTING ==========
  FigNPhont();
  FigVPhont();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
//
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}

