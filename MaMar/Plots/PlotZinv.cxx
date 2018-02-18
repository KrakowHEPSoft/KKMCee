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
//
//  Febr. 2018
//TFile DiskFileA("../workZinv/rmain.root_E105GeV_3G");
TFile DiskFileA("../workZinv/rmain.root_E161GeV_3G");
//
// FSR off, pure ISR
TFile DiskFileB("../workZinv/rmain.root_E161GeV_ISR_New");
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
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vtNuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vaNuCeex2") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vtMuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vaMuCeex2") );
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

  CaptT->DrawLatex(0.10,0.94,"d#sigma/N;        KKMC  e^{+}e^{-} -> #nu#bar{#nu}+N#gamma");
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
  CaptT->DrawLatex(0.10,0.94,"d#sigma/d(ln_{10} sin #theta_{#gamma})        e^{+}e^{-} -> #nu#bar{#nu}+N#gamma");
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
  TH1D *hst_vtNuCeex2  = (TH1D*)DiskFileA.Get("hst_vtNuCeex2");  // Phot. untaged
  TH1D *hst_vaNuCeex2  = (TH1D*)DiskFileA.Get("hst_vaNuCeex2");  // Phot. tagged
  TH1D *hst_vtMuCeex2  = (TH1D*)DiskFileA.Get("hst_vtMuCeex2");  // Phot. untaged
  TH1D *hst_vaMuCeex2  = (TH1D*)DiskFileA.Get("hst_vaMuCeex2");  // Phot. tagged

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
   TH1D *Hst=hst_vtNuCeex2;
   Hst->SetStats(0);
   Hst->SetTitle(0);
   Hst->GetXaxis()->CenterTitle();
   Hst->GetXaxis()->SetTitleSize(0.04);
   Hst->GetXaxis()->SetTitle("v=1-M^{2}_{#nu#bar{#nu}}/s");

   Hst->SetLineColor(kBlue);
//   Hst->SetMinimum( 1e-3*Hst->GetMaximum());
   Hst->DrawCopy("h");

   CaptT->DrawLatex(0.10,0.94,"d#sigma/dv;        KKMC  e^{+}e^{-} -> #nu#bar{#nu}+N#gamma");
   double ycapt = 0.85;
   CaptT->DrawLatex(0.40, ycapt,gTextEne);
//
   PlotSame2(hst_vtNuCeex2,  ycapt, kBlue, 0.5, "(a)", "#gamma's untaged");
   PlotSame2(hst_vaNuCeex2,  ycapt, kRed,  0.6, "(b)", "#gamma's tagged");

   //====================plot2========================
   cFigVPhont->cd(2);

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
   PlotSame2(hst_vtNuCeex2,  ycapt, kBlue, 0.5, "(a)","e^{+}e^{-} -> #nu#bar{#nu}+N#gamma,    #gamma's untaged");
   hst_vaMuCeex2->Scale(3.0);
   PlotSame2(hst_vaMuCeex2,  ycapt, kRed,  0.7, "(b)","e^{+}e^{-} -> #mu^{+}#mu^{-}+N#gamma,  #gamma's tagged (x3)");
   hst_vaMuCeex2->Scale(0.33333333);

  //================================================
  cFigVPhont->SaveAs("cFigVPhont.pdf");

}//FigVPhont


///////////////////////////////////////////////////////////////////////////////////
void FigNuDiff()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNuDiff =========================== "<<endl;
  ///
  TH1D *hst_vPhotNuel     = (TH1D*)DiskFileA.Get("hst_vPhotNuel");
  TH1D *hst_vPhotNumu     = (TH1D*)DiskFileA.Get("hst_vPhotNumu");
  ///
  TH1D *Hel  = hst_vPhotNuel;
  TH1D *Hmu  = hst_vPhotNumu;
//!////////////////////////////////////////////
  double SigEl= Hel->Integral();
  double SigMu= Hmu->Integral();
  cout<< "@@@@@@@@@ SigEl="<<SigEl<<" SigMu="<<SigMu<<endl;
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
  TCanvas *cNuDiff = new TCanvas("cNuDiff","cNuDiff", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cNuDiff->SetFillColor(10);
  cNuDiff->Draw();
  cNuDiff->Divide(2, 0);
 /////////////////////////////////////////////
  cNuDiff->cd(1);
  TH1D *Hst=Hel;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
//  double vcapt = 1.0 - sqr(89.3/gCMSene);
  double vcapt = 1.0 - sqr(90.5/gCMSene);
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(Hel,  ycapt, kBlue, vcapt,       "(a)"," #nu = #nu_{el}");
  PlotSame2(Hmu,  ycapt, kRed,  vcapt+0.003, "(b)"," #nu = #nu_{#mu}");
//
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dv,   e^{+}e^{-} -> #nu#bar{#nu}+N#gamma,    #gamma's taged");
  ////////////////////////////////////////////
  cNuDiff->cd(2);
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
  cNuDiff->Update();
  cNuDiff->cd();
//
  cNuDiff->SaveAs("cNuDiff.pdf");
}//FigNuDiff




///////////////////////////////////////////////////////////////////////////////////
void FigCeex12nu()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCeex12nu =========================== "<<endl;
//
  TH1D *hst_evNuCeex1     = (TH1D*)DiskFileA.Get("hst_evNuCeex1");
  TH1D *hst_evNuCeex2     = (TH1D*)DiskFileA.Get("hst_evNuCeex2");
  TH1D *hst_evNuCeex12    = (TH1D*)DiskFileA.Get("hst_evNuCeex12");
///////////////////////////////////////////////
  TH1D *RAT_ceex12 = HstRatio("RAT_ceex12", hst_evNuCeex12, hst_evNuCeex2, kBlack );
//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

  TH1D *H_Eline0  = (TH1D*)hst_evNuCeex1->Clone("H_Eline0");  // zero line
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
  TH1D *Hst=hst_evNuCeex1;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("E_{#gamma}");
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vvZ    = 1-sqr(91.2/gCMSene) +0.01;
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(hst_evNuCeex1,  ycapt, kBlue, vcapt,      "(a)"," CEEX1");
  PlotSame2(hst_evNuCeex2,  ycapt, kRed,  vcapt+0.25, "(b)"," CEEX2");
//
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dE_{#gamma};  "
		  "e^{+}e^{-} -> #nu#bar{#nu}+N#gamma,  #nu=#nu_{e}+#nu_{#mu}+#nu_{#tau}     #gamma's taged");
////////////////////////////////////////////
  cCeex12nu->cd(2);
  H_Eline0->SetStats(0);
  H_Eline0->SetTitle(0);
  H_Eline0->SetMaximum( +0.05); H_Eline0->SetMinimum( -0.05);
  H_Eline0->SetLineWidth(2);
  H_Eline0->GetXaxis()->SetTitle("E_{#gamma}");
  H_Eline0->DrawCopy("h");
  RAT_ceex12->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dE_{#gamma};     (CEEX1-CEEX2)/CEEX2");
  cCeex12nu->Update();
  cCeex12nu->cd();
//
  cCeex12nu->SaveAs("cCeex12nu.pdf");
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
///////////////////////////////////////////////
  TH1D *RAT_ceex12 = HstRatio("RAT_ceex12", hst_evMuCeex12, hst_evMuCeex2, kBlack );
//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);

  TH1D *H_Eline0  = (TH1D*)hst_evMuCeex1->Clone("H_Eline0");  // zero line
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
  TH1D *Hst=hst_evMuCeex1;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("E_{#gamma}");
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vvZ    = 1-sqr(91.2/gCMSene) +0.01;
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(hst_evMuCeex1,  ycapt, kBlue, vcapt,      "(a)"," CEEX1");
  PlotSame2(hst_evMuCeex2,  ycapt, kRed,  vcapt+0.25, "(b)"," CEEX2");
//
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dE_{#gamma};  "
		  "e^{+}e^{-} -> #mu^{+}#mu^{-}+N#gamma,       #gamma's taged");
////////////////////////////////////////////
  cCeex12mu->cd(2);
  H_Eline0->SetStats(0);
  H_Eline0->SetTitle(0);
  H_Eline0->SetMaximum( +0.09); H_Eline0->SetMinimum( -0.09);
  H_Eline0->SetLineWidth(2);
  H_Eline0->GetXaxis()->SetTitle("E_{#gamma}");
  H_Eline0->DrawCopy("h");
  RAT_ceex12->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95, "d#sigma/dE_{#gamma};     (CEEX1-CEEX2)/CEEX2");
  cCeex12mu->Update();
  cCeex12mu->cd();
//
  cCeex12mu->SaveAs("cCeex12mu.pdf");
}//FigCeex12mu

///////////////////////////////////////////////////////////////////////////////////
void FigCeex12rat()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCeex12rat =========================== "<<endl;
//
  TH1D *hst_evNuCeex1     = (TH1D*)DiskFileA.Get("hst_evNuCeex1");
  TH1D *hst_evNuCeex2     = (TH1D*)DiskFileA.Get("hst_evNuCeex2");
  TH1D *hst_evNuCeex12    = (TH1D*)DiskFileA.Get("hst_evNuCeex12");
//
  TH1D *hst_evMuCeex1     = (TH1D*)DiskFileA.Get("hst_evMuCeex1");
  TH1D *hst_evMuCeex2     = (TH1D*)DiskFileA.Get("hst_evMuCeex2");
  TH1D *hst_evMuCeex12    = (TH1D*)DiskFileA.Get("hst_evMuCeex12");
//
  TH1D *hst_evMuCeex1n     = (TH1D*)DiskFileA.Get("hst_evMuCeex1n");
  TH1D *hst_evMuCeex2n     = (TH1D*)DiskFileA.Get("hst_evMuCeex2n");
  TH1D *hst_evMuCeex12n    = (TH1D*)DiskFileA.Get("hst_evMuCeex12n");
//
  TH1D *hst_evMuCeex2ifi   = (TH1D*)DiskFileA.Get("hst_evMuCeex2ifi");
///////////////////////////////////////////////
  TH1D *RAT_ceex1 = HstRatio("RAT_ceex1", hst_evNuCeex1, hst_evMuCeex1, kBlack );  // R=nu/mu
  TH1D *RAT_ceex2 = HstRatio("RAT_ceex2", hst_evNuCeex2, hst_evMuCeex2, kBlack );  // R=nu/mu
  RAT_ceex1->Scale(1.0/3.0); // R=nu/(3*mu)
  RAT_ceex2->Scale(1.0/3.0); // R=nu/(3*mu)
  // (CEEX1-CEEX2)/CEEX2
  TH1D *RAT_ceex12 = HstDiff("RAT_ceex12", RAT_ceex1, RAT_ceex2, kBlack );
  RAT_ceex12->Divide(RAT_ceex2);
  // (CEEX1-CEEX2)/CEEX2
  TH1D *DELnu12 = HstRatio("DELnu12", hst_evNuCeex12, hst_evNuCeex2, kBlack );  // nu
  TH1D *DELmu12 = HstRatio("DELmu12", hst_evMuCeex12, hst_evMuCeex2, kBlack );  // mu
  TH1D *DEL12   = HstDiff("DEL12",    DELnu12, DELmu12, kBlack );  // nu-mu
  //
  // ****** IFI off for muons, for neutrinos the same pure ISR
  TH1D *RAT_ceex1n = HstRatio("RAT_ceex1", hst_evNuCeex1, hst_evMuCeex1n, kBlack );  // R=nu/mu
  TH1D *RAT_ceex2n = HstRatio("RAT_ceex2", hst_evNuCeex2, hst_evMuCeex2n, kBlack );  // R=nu/mu
  RAT_ceex1n->Scale(1.0/3.0); // R=nu/(3*mu)
  RAT_ceex2n->Scale(1.0/3.0); // R=nu/(3*mu)
  // direct
  TH1D *RAT_ceex12n = HstDiff("RAT_ceex12n", RAT_ceex1n, RAT_ceex2n, kBlack );
  RAT_ceex12n->Divide(RAT_ceex2n);
  // shortcut
  TH1D *DELmu12n = HstRatio("DELmu12n", hst_evMuCeex12n, hst_evMuCeex2n, kBlack );  // mu
  TH1D *DEL12n   = HstDiff("DEL12n",    DELnu12, DELmu12n, kBlack );  // nu-mu
  // IFI
  TH1D *DEL_IFI  = HstRatio("DEL_IFI",  hst_evMuCeex2ifi, hst_evMuCeex2n, kBlack );  // IFI

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
  Hst->GetXaxis()->SetTitle("E_{#gamma}");
  Hst->SetMinimum(1.9); Hst->SetMaximum(2.5);
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vvZ    = 1-sqr(91.2/gCMSene) +0.01;
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  RAT_ceex1->SetLineWidth(2); RAT_ceex2->SetLineWidth(2);
  PlotSame2(RAT_ceex1,  ycapt, kBlue, vcapt,      "(a)"," CEEX1");
  PlotSame2(RAT_ceex2,  ycapt, kRed,  vcapt+0.25, "(b)"," CEEX2");
  //
  PlotSame2(RAT_ceex1n, ycapt, kPine, vcapt,      "(c)"," CEEX1, IFIoff");
  PlotSame2(RAT_ceex2n, ycapt, kBrune,vcapt+0.25, "(d)"," CEEX2, IFIoff");
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
  Hst->GetXaxis()->SetTitle("E_{#gamma}");
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
  PlotSame2(DEL12,   ycapt, kBlack, vcapt,      "(a)"," (CEEX1-CEEX2)/CEEX2");
  PlotSame2(DEL12n,  ycapt, kRed,   vcapt+0.25, "(b)"," (CEEX1-CEEX2)/CEEX2, IFI off");
//  PlotSame2(RAT_ceex12n, ycapt,  kRed,   vcapt+0.25, "(b)"," (CEEX1-CEEX2)/CEEX2, IFI off");
  PlotSame2(DEL_IFI, ycapt, kBlue,  vcapt-0.25, "(c)"," IFI/CEEX2");

  H_Eline0->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.92,"    QED effects:  #Delta R/R ");

  cCeex12rat->cd();
//
  cCeex12rat->SaveAs("cCeex12rat.pdf");
}//FigCeex12rat


///////////////////////////////////////////////////////////////////////////////////
void FigCeex12isr()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCeex12isr =========================== "<<endl;
// FSR off
  TH1D *hst_evNuCeex1b     = (TH1D*)DiskFileB.Get("hst_evNuCeex1");
  TH1D *hst_evNuCeex2b     = (TH1D*)DiskFileB.Get("hst_evNuCeex2");
  TH1D *hst_evNuCeex12b    = (TH1D*)DiskFileB.Get("hst_evNuCeex12");
//
  TH1D *hst_evMuCeex1b     = (TH1D*)DiskFileB.Get("hst_evMuCeex1");
  TH1D *hst_evMuCeex2b     = (TH1D*)DiskFileB.Get("hst_evMuCeex2");
  TH1D *hst_evMuCeex12b    = (TH1D*)DiskFileB.Get("hst_evMuCeex12");

///////////////////////////////////////////////
  TH1D *RAT_ceex1b = HstRatio("RAT_ceex1b", hst_evNuCeex1b, hst_evMuCeex1b, kBlack );  // R=nu/mu
  TH1D *RAT_ceex2b = HstRatio("RAT_ceex2b", hst_evNuCeex2b, hst_evMuCeex2b, kBlack );  // R=nu/mu
  RAT_ceex1b->Scale(1.0/3.0); // R=nu/(3*mu)
  RAT_ceex2b->Scale(1.0/3.0); // R=nu/(3*mu)

// (CEEX1-CEEX2)/CEEX2
  TH1D *DELnu12b = HstRatio("DELnu12b", hst_evNuCeex12b, hst_evNuCeex2b, kBlack );  // nu
  TH1D *DELmu12b = HstRatio("DELmu12b", hst_evMuCeex12b, hst_evMuCeex2b, kBlack );  // mu
  TH1D *DEL12b   = HstDiff("DEL12b",    DELnu12b, DELmu12b, kBlack );  // nu-mu

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
  Hst->GetXaxis()->SetTitle("E_{#gamma}");
  Hst->SetMinimum(1.9); Hst->SetMaximum(2.5);
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vvZ    = 1-sqr(91.2/gCMSene) +0.01;
  double vcapt  = vvZ *gCMSene/2;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  RAT_ceex1b->SetLineWidth(2); RAT_ceex2b->SetLineWidth(2);
  PlotSame2(RAT_ceex1b,  ycapt, kBlue, vcapt,      "(a)"," CEEX1 FSR off");
  PlotSame2(RAT_ceex2b,  ycapt, kRed,  vcapt+0.25, "(b)"," CEEX2 FSR off");
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
  Hst->GetXaxis()->SetTitle("E_{#gamma}");
  Hst->SetLineColor(kBlack);
  Hst->DrawCopy("h");

  H_Eline0->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.92,"    RSR off:  #Delta R/R,  (CEEX1-CEEX2)/CEEX2 ");

  cCeex12isr->cd();
//
  cCeex12isr->SaveAs("cCeex12isr.pdf");
}//FigCeex12rat


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
  FigNuDiff();
  FigCeex12nu();
  FigCeex12mu();
  FigCeex12rat();
  FigCeex12isr();
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

