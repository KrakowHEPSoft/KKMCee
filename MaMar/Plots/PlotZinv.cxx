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
//
//TFile DiskFileA("../workZinv/rmain.root");
//  March 2019
TFile DiskFileA("../workZinv/rmain.root_E=161GeV_6G");
//TFile DiskFileA("../workZinv/rmain.root_E=105GeV_4G");
//  Febr. 2018
//TFile DiskFileA("../workZinv/rmain.root_E105GeV_cmax1"); // cost(heta)_max =1.0
//TFile DiskFileA("../workZinv/rmain.root_E161GeV_3G");
//TFile DiskFileA("../workZinv/rmain.root_E105GeV_4G");
//
// ************* FSR off, pure ISR **************
//TFile DiskFileB("../workZinv/rmain.root");
// New March 2019
TFile DiskFileB("../workZinv/rmain.root_E=161GeV_ISR_9G");
//TFile DiskFileB("../workZinv/main.root_E=105GeV_ISR_2G");
// Old Febr. 2018
//TFile DiskFileB("../workZinv/rmain.root_E161GeV_ISR_5G");
//TFile DiskFileB("../workZinv/rmain.root_E105GeV_ISR_1.5G");
//
//+++++++++++++++++++++++++++
TFile DiskFileX("Plot2.root","RECREATE","histograms");
//=============================================================================

//Double_t sqr( const Double_t x ){ return x*x;};
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot; // from KKMC run
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
int    g161GeVyes =0, g105GeVyes=0;
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
  //*************************************************************
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_nPhAll") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_nPhVis") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_LnThPhAll") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_LnThPhVis") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vtNuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vaNuCeex2") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vtMuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vaMuCeex2") );
  //============================================================
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvNuCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvNuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvNuCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex12") );
  //--------------------------------------
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex1n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex2n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex12n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuCeex2ifi") );
  //============================================================
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evNuCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evNuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evNuCeex12") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex12") );
  //--------------------------------------
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex1n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex2n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex12n") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_evMuCeex2ifi") );
  //===========================================================
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotNuel") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vPhotNumu") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPLCeex2") );
  //************************************************************
  TH1D *HST_KKMC_NORMB = (TH1D*)DiskFileB.Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvNuCeex1") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvNuCeex2") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvNuCeex12") );
//
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvMuCeex1") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvMuCeex2") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvMuCeex12") );
//
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvMuCeex1n") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvMuCeex2n") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvMuCeex12n") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_vvMuCeex2ifi") );
//===========================================================
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_evNuCeex1") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_evNuCeex2") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_evNuCeex12") );
//
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_evMuCeex1") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_evMuCeex2") );
  HisNorm1(HST_KKMC_NORMB, (TH1D*)DiskFileB.Get("hst_evMuCeex12") );
  //==========================================================
  //
}


///////////////////////////////////////////////////////////////////////////////////
void FigNPhot()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNPhot =========================== "<<endl;
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
  TCanvas *cFigNPhot = new TCanvas("cFigNPhot","cFigNPhot", gXcanv, gYcanv,    1200, 600);
  //                                      Name    Title        xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  cFigNPhot->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigNPhot->Divide( 2,  0);
  //====================plot1========================
  cFigNPhot->cd(1);
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

  CaptT->DrawLatex(0.10,0.94,"d#sigma/N;        KKMC  e^{+}e^{-} -> #nu#bar{#nu}+N#gamma");
  double ycapt = 0.40;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  hst_nPhAll->SetLineWidth(2);
  hst_nPhVis->SetLineWidth(2);
  PlotSame2(hst_nPhAll,  ycapt, kBlue,  +2.0, "(a)", "All photons");
  PlotSame2(hst_nPhVis,  ycapt, kRed,   +1.0, "(b)", "Tagged photons");

  //====================plot2========================
  cFigNPhot->cd(2);

  Hst=hst_LnThPhAll;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("ln_{10} sin #theta_{#gamma}");

  //Hst->SetMinimum(1.0 -0.01); Hst->SetMaximum(1.0 +0.01);
  Hst->DrawCopy("h");
  CaptT->DrawLatex(0.10,0.94,"d#sigma/d(ln_{10} sin #theta_{#gamma})        e^{+}e^{-} -> #nu#bar{#nu}+N#gamma");
  ycapt = 0.80; // starting value, to be decremented below
  CaptT->DrawLatex(0.40,ycapt,gTextEne);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev);  ycapt += -0.01;
  //
  PlotSame2(hst_LnThPhAll,  ycapt, kBlue,    -3.0, "(a)", "All #gamma's");
  PlotSame2(hst_LnThPhVis,  ycapt, kRed,     -0.5, "(b)", "Tagged #gamma's");
  //

  //================================================
  if( g161GeVyes) cFigNPhot->SaveAs("cFigNPhot_161GeV.pdf");
  if( g105GeVyes) cFigNPhot->SaveAs("cFigNPhot_105GeV.pdf");

}//FigNPhot


///////////////////////////////////////////////////////////////////////////////////
void FigVPhot()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVPhot =========================== "<<endl;
 //
  TH1D *hst_vtNuCeex2  = (TH1D*)DiskFileA.Get("hst_vtNuCeex2");  // Phot. untaged
  TH1D *hst_vaNuCeex2  = (TH1D*)DiskFileA.Get("hst_vaNuCeex2");  // Phot. tagged
  TH1D *hst_vtMuCeex2  = (TH1D*)DiskFileA.Get("hst_vtMuCeex2");  // Phot. untaged
  TH1D *hst_vaMuCeex2  = (TH1D*)DiskFileA.Get("hst_vaMuCeex2");  // Phot. tagged

 ////////////////////////////////////////
   TLatex *CaptT = new TLatex();
   CaptT->SetNDC(); // !!!
   CaptT->SetTextSize(0.04);
///////////////////////////////////////////////////////////////////////////////
   TCanvas *cFigVPhot = new TCanvas("cFigVPhot","cFigVPhot", gXcanv, gYcanv,    1200, 600);
//                                      Name    Title        xoff,yoff, WidPix,HeiPix
   gXcanv += 25, gYcanv += 25;
   cFigVPhot->SetFillColor(10);
////////////////////////////////////////////////////////////////////////////////
   cFigVPhot->Divide( 2,  0);
   //====================plot1========================
   cFigVPhot->cd(1);
   gPad->SetLogy(); // !!!!!!
   TH1D *Hst=hst_vtNuCeex2;
   Hst->SetStats(0);
   Hst->SetTitle(0);
   Hst->GetXaxis()->CenterTitle();
   Hst->GetXaxis()->SetTitleSize(0.04);
   Hst->GetXaxis()->SetTitle("v=1-M^{2}_{#nu#bar{#nu}}/s");

   Hst->SetLineColor(kBlue);
//   Hst->SetMinimum( 1e-3*Hst->GetMaximum());
   Hst->DrawCopy("h");

   CaptT->DrawLatex(0.05,0.94,"d#sigma/dv;     KKMC  e^{+}e^{-} -> #nu#bar{#nu}+n#gamma, #nu=#nu_{e}+#nu_{#mu}+#nu_{#tau}");
   double ycapt = 0.85;
   CaptT->DrawLatex(0.40, ycapt,gTextEne);
//
   PlotSame2(hst_vtNuCeex2,  ycapt, kBlue, 0.4, "(a)", "#gamma's untaged");
   PlotSame2(hst_vaNuCeex2,  ycapt, kRed,  0.5, "(b)", "#gamma's tagged");

   //====================plot2========================
   cFigVPhot->cd(2);

   gPad->SetLogy(); // !!!!!!
   Hst=hst_vtNuCeex2;
   Hst->SetStats(0);
   Hst->SetTitle(0);
   Hst->GetXaxis()->CenterTitle();
   Hst->GetXaxis()->SetTitleSize(0.04);
   Hst->GetXaxis()->SetTitle("v=1-M^{2}_{f #bar{f}}/s");

   Hst->SetLineColor(kBlue);
//   Hst->SetMinimum( 1e-3*Hst->GetMaximum());
   Hst->DrawCopy("h");

   CaptT->DrawLatex(0.10,0.94,"d#sigma/dv ");
   ycapt = 0.85;
   CaptT->DrawLatex(0.40, ycapt,gTextEne);
//
   PlotSame2(hst_vtNuCeex2,  ycapt, kBlue, 0.4, "(a)","e^{+}e^{-} -> #nu#bar{#nu}+N#gamma,    #gamma's untaged");
   hst_vaMuCeex2->Scale(3.0);
   PlotSame2(hst_vaMuCeex2,  ycapt, kRed,  0.5, "(b)","e^{+}e^{-} -> #mu^{+}#mu^{-}+N#gamma,  #gamma's tagged (x3)");
   hst_vaMuCeex2->Scale(0.33333333);

  //================================================
   if( g161GeVyes) cFigVPhot->SaveAs("cFigVPhot_161GeV.pdf");
   if( g105GeVyes) cFigVPhot->SaveAs("cFigVPhot_105GeV.pdf");

 ///////////////////////////////////////////////////////////////////////////////
    TCanvas *cFigVPhot1 = new TCanvas("cFigVPhot1","cFigVPhot1", gXcanv, gYcanv,    600, 600);
 //                                      Name    Title        xoff,yoff, WidPix,HeiPix
    gXcanv += 25, gYcanv += 25;
    cFigVPhot1->SetFillColor(10);
 ////////////////////////////////////////////////////////////////////////////////
    //====================plot1========================
    cFigVPhot1->cd();
    gPad->SetLogy(); // !!!!!!
    Hst=hst_vtNuCeex2;
    Hst->SetStats(0);
    Hst->SetTitle(0);
    Hst->GetXaxis()->CenterTitle();
    Hst->GetXaxis()->SetTitleSize(0.04);
    Hst->GetXaxis()->SetTitle("v=1-M^{2}_{#nu#bar{#nu}}/s");

    Hst->SetLineColor(kBlue);
    Hst->DrawCopy("h");

    CaptT->DrawLatex(0.05,0.94,"d#sigma/dv;     KKMC  e^{+}e^{-} -> #nu#bar{#nu}+n#gamma, #nu=#nu_{e}+#nu_{#mu}+#nu_{#tau}");
    ycapt = 0.85;
    CaptT->DrawLatex(0.40, ycapt,gTextEne);
 //
    PlotSame2(hst_vtNuCeex2,  ycapt, kBlue, 0.4, "(a)", "#gamma's untaged");
    PlotSame2(hst_vaNuCeex2,  ycapt, kRed,  0.5, "(b)", "#gamma's tagged");
    //================================================
    if( g161GeVyes) cFigVPhot1->SaveAs("cFigVPhot1_161GeV.pdf");
    if( g105GeVyes) cFigVPhot1->SaveAs("cFigVPhot1_105GeV.pdf");

}//FigVPhot


///////////////////////////////////////////////////////////////////////////////////
void FigNuDif1()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNuDif1 =========================== "<<endl;
  /// This is for paper with Alain
  TH1D *hst_vPhotNuel     = (TH1D*)DiskFileA.Get("hst_vPhotNuel");
  TH1D *hst_vPhotNumu     = (TH1D*)DiskFileA.Get("hst_vPhotNumu");
  //!////////////////////////////////////////////////////////////////////////
  double IntLumi = 10; // integrated luminosity [10 atobarns] at 161GeV
  if( gCMSene <160.0) IntLumi = 13;
  double vZ = 1.0 - sqr(91.1876/gCMSene);
  char CaptLum[300],CaptVZ[300];
  sprintf(CaptLum,"Integr. Lumi. = %3.1f [ab^{-1}]",IntLumi);
  sprintf(CaptVZ, "v_{Z}= 1-M_{Z}^{2}/s = %3.5f",vZ);

  IntLumi *=  1e9;     // the same in [nanobarns]
  double SigEl= hst_vPhotNuel->Integral("width");
  double SigMu= hst_vPhotNumu->Integral("width");
  cout<< "@@@@@@@@@ SigEl [pb] ="<<SigEl*1e3<<"  +- "<<SigEl*1e3/sqrt(SigMu*IntLumi) <<endl;
  cout<< "@@@@@@@@@ SigMu [pb] ="<<SigMu*1e3<<"  +- "<<SigMu*1e3/sqrt(SigMu*IntLumi) << endl;
  char CaptSigEl[300], CaptSigMu[300];
  sprintf(CaptSigEl,"#sigma_{e} = %2.5f +- %2.5f [pb]",   SigEl*1e3,SigEl/sqrt(SigMu*IntLumi)*1e3);
  sprintf(CaptSigMu,"#sigma_{#mu} = %2.5f +- %2.5f [pb]", SigMu*1e3,SigMu/sqrt(SigMu*IntLumi)*1e3);
  ///
  TH1D *Hel  = hst_vPhotNuel;
  TH1D *Hmu  = hst_vPhotNumu;
//!////////////////////////////////////////////
  cout<< "@@@@@@@@@ SigEl [pb]="<<SigEl<<" SigMu [pb] ="<<SigMu<<endl;
  double DelTch=(SigEl-SigMu)/(3*SigMu); /// t_chanel contrib.
  cout<< "@@@@@@@@@ R="<< DelTch <<endl;
  char CaptRat[300];
  sprintf(CaptRat,"Integrated R_{t} = %2.4f ", DelTch);
///////////////////////////////////////////////
  TH1D *RAT_NuelNumu = HstDiff("RAT_NuelNumu", Hel, Hmu, kBlack );
  RAT_NuelNumu->Divide(Hel);
  RAT_NuelNumu->Scale(0.333333); // (nuel-numu)/(3*numu)
//!////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

  TH1D *H_Vline0  = (TH1D*)hst_vPhotNuel->Clone("H_Vline0");  // zero line
  H_Vline0->Reset();
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cNuDif1 = new TCanvas("cNuDif1","cNuDif1", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cNuDif1->SetFillColor(10);
  cNuDif1->Draw();
  cNuDif1->Divide(2, 0);
 /////////////////////////////////////////////
  cNuDif1->cd(1);
  TH1D *Hst=Hel;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vcapt = 1.0 - sqr(91.2/gCMSene);
  CaptT->DrawLatex(0.35, ycapt,gTextEne);
  PlotSame2(Hel,  ycapt, kBlue, vcapt+0.010, "(a)"," #nu = #nu_{el}");
  PlotSame2(Hmu,  ycapt, kRed,  vcapt-0.010, "(b)"," #nu = #nu_{#mu}");
  //CaptT->DrawLatex(0.25, 0.25, CaptSigEl);
  //CaptT->DrawLatex(0.25, 0.20, CaptSigMu);
  //CaptT->DrawLatex(0.25, 0.15, CaptLum);

//
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dv [nb],   e^{+}e^{-} -> #nu#bar{#nu}+N#gamma,    #gamma's taged");
  ////////////////////////////////////////////
  cNuDif1->cd(2);
  H_Vline0->SetStats(0);
  H_Vline0->SetTitle(0);
  H_Vline0->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  H_Vline0->SetMaximum( +0.10); H_Vline0->SetMinimum( -0.10);
  H_Vline0->DrawCopy("h");
//
  RAT_NuelNumu->SetLineColor(kBlue);
  RAT_NuelNumu->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.50,0.75, CaptRat);
  CaptT->DrawLatex(0.10,0.95,"t-channel W contrib.    R_{t}(v)=(#nu_{el}-#nu_{#mu})/(3 #nu_{#mu})");
 //!-----------------------------
  cNuDif1->Update();
  cNuDif1->cd();
//
  if( g161GeVyes) cNuDif1->SaveAs("cNuDif1_161GeV.pdf");
  if( g105GeVyes) cNuDif1->SaveAs("cNuDif1_105GeV.pdf");
}//FigNuDif1


///////////////////////////////////////////////////////////////////////////////////
void FigNuDif2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNuDif2 =========================== "<<endl;
  ///
  TH1D *hst_vPhotNuel     = (TH1D*)DiskFileA.Get("hst_vPhotNuel");
  TH1D *hst_vPhotNumu     = (TH1D*)DiskFileA.Get("hst_vPhotNumu");
  //!////////////////////////////////////////////////////////////////////////
  double IntLumi = 10; // integrated luminosity [10 atobarns] at 161GeV
  if( gCMSene <160.0) IntLumi = 13;
  double vZ = 1.0 - sqr(91.1876/gCMSene);
  char CaptLum[300],CaptVZ[300];
  sprintf(CaptLum,"Integr. Lumi. = %3.1f [ab^{-1}]",IntLumi);
  sprintf(CaptVZ, "v_{Z}= 1-M_{Z}^{2}/s = %3.5f",vZ);

  IntLumi *=  1e9;     // the same in [nanobarns]
  double SigEl= hst_vPhotNuel->Integral("width");
  double SigMu= hst_vPhotNumu->Integral("width");
  cout<< "@@@@@@@@@ SigEl [pb] ="<<SigEl*1e3<<"  +- "<<SigEl*1e3/sqrt(SigMu*IntLumi) <<endl;
  cout<< "@@@@@@@@@ SigMu [pb] ="<<SigMu*1e3<<"  +- "<<SigMu*1e3/sqrt(SigMu*IntLumi) << endl;
  char CaptSigEl[300], CaptSigMu[300];
  sprintf(CaptSigEl,"#sigma_{e} = %2.5f +- %2.5f [pb]",   SigEl*1e3,SigEl/sqrt(SigMu*IntLumi)*1e3);
  sprintf(CaptSigMu,"#sigma_{#mu} = %2.5f +- %2.5f [pb]", SigMu*1e3,SigMu/sqrt(SigMu*IntLumi)*1e3);
  ///
  TH1D *Hel  = hst_vPhotNuel;
  TH1D *Hmu  = hst_vPhotNumu;
//!////////////////////////////////////////////
  cout<< "@@@@@@@@@ SigEl [pb]="<<SigEl<<" SigMu [pb] ="<<SigMu<<endl;
  double DelTch=(SigEl-SigMu)/(3*SigMu); /// t_chanel contrib.
  cout<< "@@@@@@@@@ R="<< DelTch <<endl;
  char CaptRat[300];
  sprintf(CaptRat,"Integrated R_{t} = %2.4f ", DelTch);
///////////////////////////////////////////////
  TH1D *RAT_NuelNumu = HstDiff("RAT_NuelNumu", Hel, Hmu, kBlack );
  RAT_NuelNumu->Divide(Hel);
  RAT_NuelNumu->Scale(0.333333); // (nuel-numu)/(3*numu)
//!////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

  TH1D *H_Vline0  = (TH1D*)hst_vPhotNuel->Clone("H_Vline0");  // zero line
  H_Vline0->Reset();
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cNuDif2 = new TCanvas("cNuDif2","cNuDif2", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cNuDif2->SetFillColor(10);
  cNuDif2->Draw();
  cNuDif2->Divide(2, 0);
 /////////////////////////////////////////////
  cNuDif2->cd(1);
  TH1D *Hst=Hel;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vcapt = 1.0 - sqr(91.2/gCMSene);
  CaptT->DrawLatex(0.35, ycapt,gTextEne);
  PlotSame2(Hel,  ycapt, kBlue, vcapt+0.010, "(a)"," #nu = #nu_{el}");
  PlotSame2(Hmu,  ycapt, kRed,  vcapt-0.010, "(b)"," #nu = #nu_{#mu}");
  CaptT->DrawLatex(0.25, 0.25, CaptSigEl);
  CaptT->DrawLatex(0.25, 0.20, CaptSigMu);
  CaptT->DrawLatex(0.25, 0.15, CaptLum);

//
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dv [nb],   e^{+}e^{-} -> #nu#bar{#nu}+N#gamma,    #gamma's taged");
  ////////////////////////////////////////////
  cNuDif2->cd(2);
  H_Vline0->SetStats(0);
  H_Vline0->SetTitle(0);
  H_Vline0->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  H_Vline0->SetMaximum( +0.10); H_Vline0->SetMinimum( -0.10);
  H_Vline0->DrawCopy("h");
//
  RAT_NuelNumu->SetLineColor(kBlue);
  RAT_NuelNumu->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.50,0.75, CaptRat);
  CaptT->DrawLatex(0.10,0.95,"t-channel W contrib.    R_{t}(v)=(#nu_{el}-#nu_{#mu})/(3 #nu_{#mu})");
 //!-----------------------------
  cNuDif2->Update();
  cNuDif2->cd();
//
  if( g161GeVyes) cNuDif2->SaveAs("cNuDif2_161GeV.pdf");
  if( g105GeVyes) cNuDif2->SaveAs("cNuDif2_105GeV.pdf");
}//FigNuDif2


///////////////////////////////////////////////////////////////////////////////////
// !!!!!!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!
void FigNuEle()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNuEle =========================== "<<endl;
  ///
  TH1D *hst_vPhotNuel     = (TH1D*)DiskFileA.Get("hst_vPhotNuel");
  TH1D *hst_vPhotNumu     = (TH1D*)DiskFileA.Get("hst_vPhotNumu");
  ///
  double IntLumi = 10; // integrated luminosity [10 atobarns] at 161GeV
  if( gCMSene <160.0) IntLumi = 13;
  double vZ = 1.0 - sqr(91.1876/gCMSene);
  char CaptLum[300],CaptVZ[300];
  sprintf(CaptLum,"Integr. Lumi. = %3.1f [ab^{-1}]",IntLumi);
  sprintf(CaptVZ, "v_{Z}= 1-M_{Z}^{2}/s = %3.5f",vZ);

  IntLumi *=  1e9;     // the same in [nanobarns]
  //
  TH1D *hst_vPhotNuInt    = HstDiff( "hst_vPhotNuInt", hst_vPhotNuel, hst_vPhotNumu, kBlack );
  //
  TH1D *hst_vNuTot1    = (TH1D*)hst_vPhotNuel->Clone("hst_vNuTot1");
  hst_vNuTot1->Add(hst_vPhotNumu, hst_vPhotNuInt,  3.0,  1.0);
  //
  TH1D *hst_vNuTot0    = (TH1D*)hst_vPhotNuel->Clone("hst_vNuTot0");
  hst_vNuTot0->Add(hst_vPhotNumu, hst_vPhotNuInt,  3.0,  0.0);
  //
  // translate into No of events for given integr. lumi
  //TH1D *hst_vNuel     = HstEvent("hst_vNuel", hst_vPhotNuel, IntLumi);
  //TH1D *hst_vNumu     = HstEvent("hst_vNumu", hst_vPhotNumu, IntLumi);
  //
  double Xmax = hst_vPhotNumu->GetXaxis()->GetXmax();
  double Xmin = hst_vPhotNumu->GetXaxis()->GetXmin();
  cout<< "@@@@@@@@@ Xmin, Xmax = "<<Xmin <<" "<< Xmax<<" Xmax-Xmin="<< Xmax-Xmin<<endl;

  double SigTot1= hst_vNuTot1->Integral("width");
  double SigTot0= hst_vNuTot0->Integral("width");
  //
  double SigTot1Mi = hst_vNuTot1->Integral( 1,20,"width");
  double SigTot1Pl = hst_vNuTot1->Integral(21,40,"width");
  double SigTot0Mi = hst_vNuTot0->Integral( 1,20,"width");
  double SigTot0Pl = hst_vNuTot0->Integral(21,40,"width");

  //!////////////////////////////////////////////////////////////////////////
  cout<< "@@@@@@@@@ SigTot1 [pb] ="<<SigTot1*1e3<<"  +- "<<SigTot1*1e3/sqrt(SigTot0*IntLumi) <<endl;
  cout<< "@@@@@@@@@ SigTot0 [pb] ="<<SigTot0*1e3<<"  +- "<<SigTot0*1e3/sqrt(SigTot0*IntLumi) <<endl;
  double DelTch=(SigTot1-SigTot0)/(SigTot1); /// t_chanel contrib.
  cout<< "@@@@@@@@@ A(0) ="<< DelTch <<" +- " << 1/sqrt(SigTot1*IntLumi) << endl;
  //
  cout<< "@@@@@@@@@ SigTot1Pl="<<SigTot1Pl<<"   SigTot1Mi="<<SigTot1Mi<<" sum/SigTot1="<<(SigTot1Pl+SigTot1Mi)/SigTot1<<endl;
  double ASYnuEl= (SigTot1Pl-SigTot1Mi)/(SigTot1Pl+SigTot1Mi);
  double ASYnuMu= (SigTot0Pl-SigTot0Mi)/(SigTot0Pl+SigTot0Mi);
  cout<< "@@@@@@@@@ ASYnuEl=(Pl-Mi)/Pl+Mi)="<<ASYnuEl <<" +- " << 1/sqrt(SigTot1*IntLumi) <<endl;
  cout<< "@@@@@@@@@ ASYnuMu=(Pl-Mi)/Pl+Mi)="<<ASYnuMu <<" +- " << 1/sqrt(SigTot1*IntLumi)  <<endl;
  cout<< "@@@@@@@@@   ASYnuEl - ASYnuMu =  "<< ASYnuEl - ASYnuMu <<" +- " << 1/sqrt(SigTot1*IntLumi)  <<endl;
  char CaptRat[300], CaptS0[300], CaptSigTot1[300], CaptSigTot0[300];
  sprintf(CaptRat,"(#sigma_{on}-#sigma_{off})/#sigma_{on} = %2.5f +- %2.5f ", DelTch, 1/sqrt(SigTot1*IntLumi));
  sprintf(CaptS0, "S(0)= %2.5f +- %2.5f ", ASYnuMu, 1/sqrt(SigTot1*IntLumi));
  sprintf(CaptSigTot1,"#sigma_{on } = %2.5f +- %2.5f [pb]",  SigTot1*1e3,SigTot1*1e3/sqrt(SigTot0*IntLumi));
  sprintf(CaptSigTot0,"#sigma_{off} = %2.5f +- %2.5f [pb]", SigTot0*1e3,SigTot0*1e3/sqrt(SigTot0*IntLumi));
  cout<< "@@@@@@@@@ FCCee statistics for 10atob. ="<< SigTot1*IntLumi <<endl;
  double EquEve = MCequiv(hst_vNuTot1);
  cout<< "@@@@@@@@@ No. of equiv. WT=1 MC events ="<< EquEve <<endl;
//////////////////////////////////////////////
  double etamax = 0.100;
  int nbx =200;
  TH1D *hst_DelA  = new TH1D("hst_nPhVis" , "Delta Asym",  nbx, -etamax , etamax);
  hst_DelA->Sumw2();
  double staterr= 1/sqrt(SigTot0*IntLumi);
  double dASY,eta, SigTot2Pl, SigTot2Mi, ASYnu2;
  TH1D *hst_vNuTot2    = (TH1D*)hst_vPhotNuel->Clone("hst_vNuTot2");
  for(int ib=1; ib <= nbx; ib++){
	eta = -etamax + 2*etamax*(ib-0.5)/nbx;
	dASY = eta;
	hst_vNuTot2->Add(hst_vPhotNumu, hst_vPhotNuInt,  3.0, sqrt(1+eta) );
	SigTot2Mi = hst_vNuTot2->Integral( 1,20,"width");
	SigTot2Pl = hst_vNuTot2->Integral(21,40,"width");
	ASYnu2= (SigTot2Pl-SigTot2Mi)/(SigTot2Pl+SigTot2Mi);
	dASY =ASYnu2-ASYnuEl;
	hst_DelA->SetBinContent(ib, dASY);
	hst_DelA->SetBinError(  ib, staterr);
	cout<< "%%%%%%%%% "<<ib<<"  "<< eta <<"  "<<dASY<< "  +-"<<  staterr <<endl;
  }
//////////////////////////////////////////////
  TH1D *Hel  = hst_vNuTot1;
  TH1D *Hmu  = hst_vNuTot0;
///////////////////////////////////////////////////////////////////////////
  TH1D *RAT_NuelNumu = HstDiff("RAT_NuelNumu", Hel, Hmu, kBlack );
  RAT_NuelNumu->Divide(Hel);
  RAT_NuelNumu->Scale(0.333333); // (nuel-numu)/(3*numu)
  //
  //TH1D *RAT_NuelNumu = hst_vNuInt;
  //RAT_NuelNumu->Divide(hst_vNuTot1);
//!////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

  TH1D *H_Vline0  = (TH1D*)hst_DelA->Clone("H_Vline0");  // zero line
  H_Vline0->Reset();
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cNuEle = new TCanvas("cNuEle","cNuEle", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cNuEle->SetFillColor(10);
  cNuEle->Draw();
  cNuEle->Divide(2, 0);
 /////////////////////////////////////////////
  cNuEle->cd(1);
  TH1D *Hst=Hel;
  //Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.45;
  double vcapt = vZ;
  CaptT->DrawLatex(0.35, ycapt,gTextEne);
  PlotSame2(Hel,  ycapt, kBlue, vcapt+0.010, "(a)"," W_{t} is ON");
  PlotSame2(Hmu,  ycapt, kRed,  vcapt-0.010, "(b)"," W_{t} is OFF");
  CaptT->DrawLatex(0.25, 0.30, CaptSigTot1);
  CaptT->DrawLatex(0.25, 0.25, CaptSigTot0);
  CaptT->DrawLatex(0.25, 0.20, CaptRat);
  CaptT->DrawLatex(0.25, 0.15, CaptLum);
//
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dv [nb],   e^{+}e^{-} -> 3#nu#bar{#nu}+N#gamma,    #gamma's taged");
  ////////////////////////////////////////////
  cNuEle->cd(2);
  H_Vline0->SetStats(0);
  H_Vline0->SetTitle(0);
  H_Vline0->GetXaxis()->SetTitle("#eta= 3#Gamma_{#nu_{e}}/#Gamma_{invis.}^{SM}-1");
  H_Vline0->SetMaximum( +0.0015); H_Vline0->SetMinimum( -0.0015);
  if(gCMSene < 160.0){
	  H_Vline0->SetMaximum( +0.0010); H_Vline0->SetMinimum( -0.0010);}
  H_Vline0->SetTitleSize(0.04);
  H_Vline0->SetLabelSize(0.030);
  H_Vline0->SetNdivisions(10);
  H_Vline0->DrawCopy("h");
//
  hst_DelA->SetLineColor(kBlue);
  hst_DelA->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.10,0.95, "#Delta S = S(#eta)-S(0)");
  CaptT->DrawLatex(0.40,0.85, "S = (#sigma_{+}-#sigma_{-})/(#sigma_{+}-#sigma_{-})" );
  CaptT->DrawLatex(0.40,0.80, "#sigma_{+}=#sigma(v>v_{Z}),  #sigma_{-}=#sigma(v<v_{Z})" );
  CaptT->DrawLatex(0.40,0.75, CaptVZ);
  CaptT->DrawLatex(0.40,0.70, CaptS0);
 //!-----------------------------
  cNuEle->Update();
  cNuEle->cd();
//
  if( g161GeVyes) cNuEle->SaveAs("cNuEle_161GeV.pdf");
  if( g105GeVyes) cNuEle->SaveAs("cNuEle_105GeV.pdf");
}//FigNuEle




///////////////////////////////////////////////////////////////////////////////////
void FigCeex12nu()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCeex12nu =========================== "<<endl;
  //
    TH1D *hst_evNuCeex1     = (TH1D*)DiskFileA.Get("hst_evNuCeex1");
    TH1D *hst_evNuCeex2     = (TH1D*)DiskFileA.Get("hst_evNuCeex2");
    TH1D *hst_evNuCeex12    = (TH1D*)DiskFileA.Get("hst_evNuCeex12");
  //
    TH1D *hst_vvNuCeex1     = (TH1D*)DiskFileA.Get("hst_vvNuCeex1");
    TH1D *hst_vvNuCeex2     = (TH1D*)DiskFileA.Get("hst_vvNuCeex2");
    TH1D *hst_vvNuCeex12    = (TH1D*)DiskFileA.Get("hst_vvNuCeex12");
///////////////////////////////////////////////
  TH1D *RAT_ceex12 = HstRatio("RAT_ceex12", hst_vvNuCeex12, hst_vvNuCeex2, kBlack );
//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

  TH1D *H_Eline0  = (TH1D*)hst_vvNuCeex1->Clone("H_Eline0");  // zero line
  H_Eline0->Reset();
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cCeex12nu = new TCanvas("cCeex12nu","cCeex12nu", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cCeex12nu->SetFillColor(10);
  cCeex12nu->Draw();
  cCeex12nu->Divide(2, 0);
/////////////////////////////////////////////
  cCeex12nu->cd(1);
  TH1D *Hst=hst_vvNuCeex1;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vvZ    = 1-sqr(91.2/gCMSene);
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(hst_vvNuCeex1,  ycapt, kBlue, vvZ+0.003,  "(a)"," CEEX1");
  PlotSame2(hst_vvNuCeex2,  ycapt, kRed,  vvZ-0.003,  "(b)"," CEEX2");
//
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dv;  "
		  "e^{+}e^{-} -> #nu#bar{#nu}+N#gamma,  #nu=#nu_{e}+#nu_{#mu}+#nu_{#tau}     #gamma's taged");
////////////////////////////////////////////
  cCeex12nu->cd(2);
  H_Eline0->SetStats(0);
  H_Eline0->SetTitle(0);
  H_Eline0->SetMaximum( +0.05); H_Eline0->SetMinimum( -0.05);
  H_Eline0->SetLineWidth(2);
  H_Eline0->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  H_Eline0->DrawCopy("h");
  RAT_ceex12->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dv;     (CEEX1-CEEX2)/CEEX2");
  cCeex12nu->Update();
  cCeex12nu->cd();
//
  if( g161GeVyes) cCeex12nu->SaveAs("cCeex12nu_161GeV.pdf");
  if( g105GeVyes) cCeex12nu->SaveAs("cCeex12nu_105GeV.pdf");
//
  /*
  int Nb = hst_vvNuCeex2->GetNbinsX();
  cout<<"**************************************************************************"<<endl;
  cout<<" From FigCeex12nu: gCMSene="<<gCMSene<<endl;
  cout<<" Nb=  "<<Nb<<endl;
  for(int ix=1; ix <= Nb; ix++){
	  cout<<ix<<"    "<< hst_vvNuCeex2->GetBinContent(ix)<<"    +-    "<< hst_vvNuCeex2->GetBinError(ix)<<endl;
  }
  cout<<"**************************************************************************"<<endl;
  */
//
}//FigCeex12nu


///////////////////////////////////////////////////////////////////////////////////
void FigCeex12mu()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCeex12mu =========================== "<<endl;
//
  TH1D *hst_evMuCeex1     = (TH1D*)DiskFileA.Get("hst_evMuCeex1");
  TH1D *hst_evMuCeex2     = (TH1D*)DiskFileA.Get("hst_evMuCeex2");
  TH1D *hst_evMuCeex12    = (TH1D*)DiskFileA.Get("hst_evMuCeex12");
  //
  TH1D *hst_vvMuCeex1     = (TH1D*)DiskFileA.Get("hst_vvMuCeex1");
  TH1D *hst_vvMuCeex2     = (TH1D*)DiskFileA.Get("hst_vvMuCeex2");
  TH1D *hst_vvMuCeex12    = (TH1D*)DiskFileA.Get("hst_vvMuCeex12");
///////////////////////////////////////////////
  TH1D *RAT_ceex12  = HstRatio("RAT_ceex12", hst_evMuCeex12, hst_evMuCeex2, kBlack );
  TH1D *RAT_ceex12v = HstRatio("RAT_ceex12", hst_vvMuCeex12, hst_vvMuCeex2, kBlack );
//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

////  TH1D *H_Eline0  = (TH1D*)hst_evMuCeex1->Clone("H_Eline0");  // zero line
  TH1D *H_Eline0  = (TH1D*)hst_vvMuCeex1->Clone("H_Eline0");  // zero line
  H_Eline0->Reset();
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cCeex12mu = new TCanvas("cCeex12mu","cCeex12mu", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cCeex12mu->SetFillColor(10);
  cCeex12mu->Draw();
  cCeex12mu->Divide(2, 0);
/////////////////////////////////////////////
  cCeex12mu->cd(1);
  TH1D *Hst=hst_vvMuCeex1;
  Hst->SetStats(0);
  Hst->SetTitle(0);
//  Hst->GetXaxis()->SetTitle("E_{#gamma}");
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vvZ    = 1-sqr(91.2/gCMSene);
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(hst_vvMuCeex1,  ycapt, kBlue, vvZ+0.005,  "(a)"," CEEX1");
  PlotSame2(hst_vvMuCeex2,  ycapt, kRed,  vvZ-0.005,  "(b)"," CEEX2");
//
//  CaptT->DrawLatex(0.10,0.95, "d#sigma/dE_{#gamma};  "
//		  "e^{+}e^{-} -> #mu^{+}#mu^{-}+N#gamma,       #gamma's taged");
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dv;  "
		  "e^{+}e^{-} -> #mu^{+}#mu^{-}+N#gamma,       #gamma's taged");
///////////////////////////////////////////
  cCeex12mu->cd(2);
  H_Eline0->SetStats(0);
  H_Eline0->SetTitle(0);
  H_Eline0->SetMaximum( +0.09); H_Eline0->SetMinimum( -0.09);
  H_Eline0->SetLineWidth(2);
////  H_Eline0->GetXaxis()->SetTitle("E_{#gamma}");
  H_Eline0->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  H_Eline0->DrawCopy("h");
  RAT_ceex12v->DrawCopy("hsame");
////  CaptT->DrawLatex(0.10,0.95, "d#sigma/dE_{#gamma};     (CEEX1-CEEX2)/CEEX2");
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dv;     (CEEX1-CEEX2)/CEEX2");
  cCeex12mu->Update();
  cCeex12mu->cd();
//
  if( g161GeVyes) cCeex12mu->SaveAs("cCeex12mu_161GeV.pdf");
  if( g105GeVyes) cCeex12mu->SaveAs("cCeex12mu_105GeV.pdf");
}//FigCeex12mu

///////////////////////////////////////////////////////////////////////////////////
void FigCeex12rat()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCeex12rat =========================== "<<endl;
//
  TH1D *hst_vvNuCeex1     = (TH1D*)DiskFileA.Get("hst_vvNuCeex1");
  TH1D *hst_vvNuCeex2     = (TH1D*)DiskFileA.Get("hst_vvNuCeex2");
  TH1D *hst_vvNuCeex12    = (TH1D*)DiskFileA.Get("hst_vvNuCeex12");
//
  TH1D *hst_vvMuCeex1     = (TH1D*)DiskFileA.Get("hst_vvMuCeex1");
  TH1D *hst_vvMuCeex2     = (TH1D*)DiskFileA.Get("hst_vvMuCeex2");
  TH1D *hst_vvMuCeex12    = (TH1D*)DiskFileA.Get("hst_vvMuCeex12");
//
  TH1D *hst_vvMuCeex1n     = (TH1D*)DiskFileA.Get("hst_vvMuCeex1n");
  TH1D *hst_vvMuCeex2n     = (TH1D*)DiskFileA.Get("hst_vvMuCeex2n");
  TH1D *hst_vvMuCeex12n    = (TH1D*)DiskFileA.Get("hst_vvMuCeex12n");
//
  TH1D *hst_vvMuCeex2ifi   = (TH1D*)DiskFileA.Get("hst_vvMuCeex2ifi");
///////////////////////////////////////////////
  TH1D *RAT_ceex1 = HstRatio("RAT_ceex1", hst_vvNuCeex1, hst_vvMuCeex1, kBlack );  // R=nu/mu
  TH1D *RAT_ceex2 = HstRatio("RAT_ceex2", hst_vvNuCeex2, hst_vvMuCeex2, kBlack );  // R=nu/mu
  RAT_ceex1->Scale(1.0/3.0); // R=nu/(3*mu)
  RAT_ceex2->Scale(1.0/3.0); // R=nu/(3*mu)
  // (CEEX1-CEEX2)/CEEX2
  TH1D *RAT_ceex12 = HstDiff("RAT_ceex12", RAT_ceex1, RAT_ceex2, kBlack );
  RAT_ceex12->Divide(RAT_ceex2);
  // (CEEX1-CEEX2)/CEEX2
  TH1D *DELnu12 = HstRatio("DELnu12", hst_vvNuCeex12, hst_vvNuCeex2, kBlack );  // nu
  TH1D *DELmu12 = HstRatio("DELmu12", hst_vvMuCeex12, hst_vvMuCeex2, kBlack );  // mu
  TH1D *DEL12   = HstDiff("DEL12",    DELnu12, DELmu12, kBlack );  // nu-mu
  //
  // ****** IFI off for muons, for neutrinos the same pure ISR
  TH1D *RAT_ceex1n = HstRatio("RAT_ceex1", hst_vvNuCeex1, hst_vvMuCeex1n, kBlack );  // R=nu/mu
  TH1D *RAT_ceex2n = HstRatio("RAT_ceex2", hst_vvNuCeex2, hst_vvMuCeex2n, kBlack );  // R=nu/mu
  RAT_ceex1n->Scale(1.0/3.0); // R=nu/(3*mu)
  RAT_ceex2n->Scale(1.0/3.0); // R=nu/(3*mu)
  // direct
  TH1D *RAT_ceex12n = HstDiff("RAT_ceex12n", RAT_ceex1n, RAT_ceex2n, kBlack );
  RAT_ceex12n->Divide(RAT_ceex2n);
  // shortcut
  TH1D *DELmu12n = HstRatio("DELmu12n", hst_vvMuCeex12n, hst_vvMuCeex2n, kBlack );  // mu
  TH1D *DEL12n   = HstDiff("DEL12n",    DELnu12, DELmu12n, kBlack );  // nu-mu
  // IFI
  TH1D *DEL_IFI  = HstRatio("DEL_IFI",  hst_vvMuCeex2ifi, hst_vvMuCeex2n, kBlack );  // IFI

//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

  TH1D *H_Eline0  = (TH1D*)DiskFileX.Get("H_Eline0");

//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cCeex12rat = new TCanvas("cCeex12rat","cCeex12rat", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cCeex12rat->SetFillColor(10);
  cCeex12rat->Draw();
  cCeex12rat->Divide(2, 0);
/////////////////////////////////////////////
  cCeex12rat->cd(1);
  TH1D *Hst=RAT_ceex2;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
//  Hst->SetMinimum(1.9); Hst->SetMaximum(2.5);
  if( g161GeVyes) {Hst->SetMinimum(1.6); Hst->SetMaximum(2.2);}
  if( g105GeVyes) {Hst->SetMinimum(1.6); Hst->SetMaximum(2.2);}
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.30;
  double vvZ    = 1-sqr(91.2/gCMSene);
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  RAT_ceex1->SetLineWidth(2); RAT_ceex2->SetLineWidth(2);
  PlotSame2(RAT_ceex1,  ycapt, kBlue, vvZ+0.010,"(a)"," CEEX1");
  PlotSame2(RAT_ceex2,  ycapt, kRed,  vvZ-0.015,"(b)"," CEEX2");
  //
  PlotSame2(RAT_ceex1n, ycapt, kPine, vvZ-0.015,"(c)"," CEEX1, IFIoff");
  PlotSame2(RAT_ceex2n, ycapt, kBrune,vvZ+0.010,"(d)"," CEEX2, IFIoff");
//
  CaptT->DrawLatex(0.10,0.92,
   "     R=d#sigma_{#nu#nu#gamma}/3d#sigma_{#mu#mu#gamma},     #nu=#nu_{e}+#nu_{#mu}+#nu_{#tau}");
////////////////////////////////////////////
  cCeex12rat->cd(2);
  Hst = RAT_ceex12; // direct
  Hst = DEL12;      // Small staticical errors:
  //
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->SetMaximum( +0.09); Hst->SetMinimum( -0.09);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->SetLineColor(kBlack);
  Hst->DrawCopy("h");

  // **** Direct
  //RAT_ceex12->SetLineColor(kRed);
  //RAT_ceex12->DrawCopy("hsame");
  // ***** Small staticical errors:
  //DEL12->SetLineColor(kBlack);
  //DEL12->DrawCopy("hsame");

  ycapt = 0.85;
  DEL12->SetLineWidth(2);
  PlotSame2(DEL12,   ycapt, kBlack, vvZ-0.010, "(a)"," (CEEX1-CEEX2)/CEEX2");
  PlotSame2(DEL12n,  ycapt, kRed,   vvZ+0.015, "(b)"," (CEEX1-CEEX2)/CEEX2, IFI off");
//  PlotSame2(RAT_ceex12n, ycapt,  kRed,   vvZ+0.005, "(b)"," (CEEX1-CEEX2)/CEEX2, IFI off");
  PlotSame2(DEL_IFI, ycapt, kBlue,  vvZ-0.010, "(c)"," IFI/CEEX2");

  H_Eline0->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.92,"    QED effects:  #Delta R/R ");

  cCeex12rat->cd();
//
  if( g161GeVyes) cCeex12rat->SaveAs("cCeex12rat_161GeV.pdf");
  if( g105GeVyes) cCeex12rat->SaveAs("cCeex12rat_105GeV.pdf");
}//FigCeex12rat


///////////////////////////////////////////////////////////////////////////////////
void FigCeex12isr()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCeex12isr =========================== "<<endl;
// FSR off
  TH1D *hst_vvNuCeex1b     = (TH1D*)DiskFileB.Get("hst_vvNuCeex1");
  TH1D *hst_vvNuCeex2b     = (TH1D*)DiskFileB.Get("hst_vvNuCeex2");
  TH1D *hst_vvNuCeex12b    = (TH1D*)DiskFileB.Get("hst_vvNuCeex12");
//
  TH1D *hst_vvMuCeex1b     = (TH1D*)DiskFileB.Get("hst_vvMuCeex1");
  TH1D *hst_vvMuCeex2b     = (TH1D*)DiskFileB.Get("hst_vvMuCeex2");
  TH1D *hst_vvMuCeex12b    = (TH1D*)DiskFileB.Get("hst_vvMuCeex12");

///////////////////////////////////////////////
  TH1D *RAT_ceex1b = HstRatio("RAT_ceex1b", hst_vvNuCeex1b, hst_vvMuCeex1b, kBlack );  // R=nu/mu
  TH1D *RAT_ceex2b = HstRatio("RAT_ceex2b", hst_vvNuCeex2b, hst_vvMuCeex2b, kBlack );  // R=nu/mu
  RAT_ceex1b->Scale(1.0/3.0); // R=nu/(3*mu)
  RAT_ceex2b->Scale(1.0/3.0); // R=nu/(3*mu)

// (CEEX1-CEEX2)/CEEX2
  TH1D *DELnu12b = HstRatio("DELnu12b", hst_vvNuCeex12b, hst_vvNuCeex2b, kBlack );  // nu
  TH1D *DELmu12b = HstRatio("DELmu12b", hst_vvMuCeex12b, hst_vvMuCeex2b, kBlack );  // mu
  TH1D *DEL12b   = HstDiff( "DEL12b",    DELnu12b, DELmu12b, kBlack );  // nu-mu

//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  TH1D *H_Eline0  = (TH1D*)DiskFileX.Get("H_Eline0");
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cCeex12isr = new TCanvas("cCeex12isr","cCeex12isr", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cCeex12isr->SetFillColor(10);
  cCeex12isr->Draw();
  cCeex12isr->Divide(2, 0);
/////////////////////////////////////////////
  cCeex12isr->cd(1);
  TH1D *Hst=RAT_ceex2b;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  if( g161GeVyes) {Hst->SetMinimum(1.7); Hst->SetMaximum(2.2);}
  if( g105GeVyes) {Hst->SetMinimum(1.7); Hst->SetMaximum(2.2);}
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vvZ    = 1-sqr(91.2/gCMSene);
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  RAT_ceex1b->SetLineWidth(2); RAT_ceex2b->SetLineWidth(2);
  PlotSame2(RAT_ceex1b,  ycapt, kBlue, vvZ+0.000, "(a)"," CEEX1 FSR off");
  PlotSame2(RAT_ceex2b,  ycapt, kRed,  vvZ+0.005, "(b)"," CEEX2 FSR off");
//
  CaptT->DrawLatex(0.10,0.92,
     "     R=d#sigma_{#nu#nu#gamma}/3d#sigma_{#mu#mu#gamma},     #nu=#nu_{e}+#nu_{#mu}+#nu_{#tau}");
/////////////////////////////////////////////
  cCeex12isr->cd(2);
  Hst = DEL12b;      // Small staticical errors:
  //
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->SetMaximum( +0.001); Hst->SetMinimum( -0.001);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->SetLineColor(kBlack);
  Hst->DrawCopy("h");

  H_Eline0->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.92,"    FSR off:  #Delta R/R,  (CEEX1-CEEX2)/CEEX2 ");

  cCeex12isr->cd();
//
  if( g161GeVyes) cCeex12isr->SaveAs("cCeex12isr_161GeV.pdf");
  if( g105GeVyes) cCeex12isr->SaveAs("cCeex12isr_105GeV.pdf");
}//FigCeex12rat


///////////////////////////////////////////////////////////////////////////////////
void FigCeex2fsr()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCeex2fsr =========================== "<<endl;
// FSR off
  TH1D *hst_vvNuCeex2      = (TH1D*)DiskFileA.Get("hst_vvNuCeex2"); // ISR+FSR+IFI
//
//  TH1D *hst_vvMuCeex2      = (TH1D*)DiskFileA.Get("hst_vvMuCeex2"); // ISR+FSR+IFI
  TH1D *hst_vvMuCeex2      = (TH1D*)DiskFileA.Get("hst_vvMuCeex2n"); // ISR+FSR, IFIoff!
  TH1D *hst_vvMuCeex2b     = (TH1D*)DiskFileB.Get("hst_vvMuCeex2"); // ISR only
//
///////////////////////////////////////////////
  TH1D *RAT_ceex2a = HstRatio("RAT_ceex2",  hst_vvNuCeex2, hst_vvMuCeex2,  kBlack );  // R=nu/mu
  TH1D *RAT_ceex2b = HstRatio("RAT_ceex2b", hst_vvNuCeex2, hst_vvMuCeex2b, kBlack );  // R=nu/mu
  RAT_ceex2a->Scale(1.0/3.0); // R=nu/(3*mu)
  RAT_ceex2b->Scale(1.0/3.0); // R=nu/(3*mu)

  TH1D *DELfsr   = HstDiff( "DELfsr",    hst_vvMuCeex2, hst_vvMuCeex2b, kBlack );  // nu-mu
  DELfsr->Divide(hst_vvMuCeex2b);

//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  TH1D *H_Eline0  = (TH1D*)DiskFileX.Get("H_Eline0");
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cCeex2fsr = new TCanvas("cCeex2fsr","cCeex2fsr", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cCeex2fsr->SetFillColor(10);
  cCeex2fsr->Draw();
  cCeex2fsr->Divide(2, 0);
/////////////////////////////////////////////
  cCeex2fsr->cd(1);
  TH1D *Hst= RAT_ceex2a;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  if( g161GeVyes) {Hst->SetMinimum(1.7); Hst->SetMaximum(2.2);}
  if( g105GeVyes) {Hst->SetMinimum(1.7); Hst->SetMaximum(2.2);}
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.35;
  double vvZ    = 1-sqr(92.2/gCMSene);
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  //RAT_ceex1b->SetLineWidth(2); RAT_ceex2b->SetLineWidth(2);
  PlotSame2(RAT_ceex2a,  ycapt, kBlue, vvZ+0.000,  "(a)"," CEEX2 FSR on");
  PlotSame2(RAT_ceex2b,  ycapt, kRed,  vvZ+0.005,  "(b)"," CEEX2 FSR off");
  CaptT->DrawLatex(0.10,0.92,
     "     R=d#sigma_{#nu#nu#gamma}/3d#sigma_{#mu#mu#gamma},     #nu=#nu_{e}+#nu_{#mu}+#nu_{#tau}");
///////////////////////////////////////////
  cCeex2fsr->cd(2);
  Hst = DELfsr;      // Small staticical errors:
  //
  Hst->SetStats(0);
  Hst->SetTitle(0);
  if(g105GeVyes) {Hst->SetMinimum( -0.05); Hst->SetMaximum( +0.25);}
  if(g161GeVyes) {Hst->SetMinimum( -0.10); Hst->SetMaximum( +0.10);}
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->SetLineColor(kBlack);
  Hst->DrawCopy("h");

  H_Eline0->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.92,"    (FSRon-FSRoff)/FSRoff ");
//
  cCeex2fsr->cd();
//
  if(g161GeVyes) cCeex2fsr->SaveAs("cCeex2fsr_161GeV.pdf");
  if(g105GeVyes) cCeex2fsr->SaveAs("cCeex2fsr_105GeV.pdf");
}//FigCeex12rat





///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++

  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  int Nodes= HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
  if( fabs(gCMSene-161)<0.01) g161GeVyes = 1;
  if( fabs(gCMSene-105)<0.01) g105GeVyes = 1;
  //
  TH1D *HST_KKMC_NORMB = (TH1D*)DiskFileB.Get("HST_KKMC_NORMA");
  int Nodes2       = HST_KKMC_NORMB->GetBinContent(511);       // No of farm nodes (trick)
  double gCMSene2  = HST_KKMC_NORMB->GetBinContent(1)/Nodes2;   // CMSene=xpar(1), farn adjusted
  double gNevTot2  = HST_KKMC_NORMB->GetEntries();             // MC statistics from KKMC
  //
  if( fabs(gCMSene2-gCMSene)>0.01 ){
	  cout<< "STOP: different CMSene =" << gCMSene << "GeV  "<< gCMSene2 <<"GeV  " <<endl;
	  exit(13);
  }

  HistNormalize();     // Renormalization of MC histograms
  //========== PLOTTING ==========
  //*
  FigNPhot();
  FigVPhot();
  FigCeex12nu();
  FigCeex12mu();
  FigCeex12rat();
  FigCeex12isr();
  FigCeex2fsr();
  //*/
  FigNuDif1();
  FigNuDif2();
  FigNuEle();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< "GeV"<<endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<< endl;
//
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}

