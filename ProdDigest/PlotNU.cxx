// make PlotNU-run
// make Tab2-pdf

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
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"
#include "TFile.h"

#include "HisNorm.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
TFile *DiskFileA;

//TString FileA= "../ProdRun/workNU/histo.root";

TString FileA= "../ProdRun/workNU/histo.root_161_GeV_WTed_1G";
//TString FileA= "../ProdRun/workNU/histo.root_189GeV_WT=1_1G";

//TString FileA= "../ProdRun/workNU/histo.root_105GeV_4G";
//TString FileA= "../ProdRun/workNU/histo.root_161GeV_7G";

FILE *DiskFileTeX;

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot; // from KKMC run
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
int    g161GeVyes =0, g125GeVyes=0, g105GeVyes=0;
//
float  gXcanv = 20, gYcanv = 20, gDcanv = 30;
///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////
void PlotSame2(TH1D *HST, const double &xcapt, double &ycapt, Int_t kolor, const double xx,  TString label,  TString opis)
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
  DiskFileA->ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  //
  //HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_WtMain") );
  //HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_WtFoam") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_nPhot") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_CosTheta") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_CosThOve") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vvTrue") );
//
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Eex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Eex0") );
//===========================================================
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vPhotNuel") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vPhotNumu") );

  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vtNuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vaNuCeex2") );

  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vaNuMuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vaNuElCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vaNuTaCeex2") );

  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vxNuMuCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vxNuElCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vxNuTaCeex2") );

}//HistNormalize


///////////////////////////////////////////////////////////////////////////////////
void FigWtMain()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigWtMain =========================== "<<endl;
  //
  TH1D *hst_WtMain    = (TH1D*)DiskFileA->Get("hst_WtMain");
//  TH1D *hst_WtMain    = (TH1D*)DiskFileA->Get("hst_WtMain200");
  TH1D *hst_WtFoam    = (TH1D*)DiskFileA->Get("hst_WtFoam");
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigWtMain = new TCanvas("cFigWtMain","cFigWtMain: general info ",  gXcanv,  gYcanv,    500,  500);
  //                                  Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigWtMain->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  //==========plot1==============
  gPad->SetLogy(); // !!!!!!
  TH1D *HST; //
  HST = hst_WtMain; //
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->GetXaxis()->SetTitle("WT");
  HST->GetXaxis()->SetLabelSize(0.04);
  HST->DrawCopy("h");

  CaptT->DrawLatex(0.06,0.95, "Events");
  double ycapt = 0.80; double xcapt=0.20;
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(xcapt,ycapt, "e^{+}e^{-} -> #nu#bar{#nu}(+n#gamma)");
  CaptT->DrawLatex(xcapt+0.40,ycapt,gTextEne);

  PlotSame2(hst_WtMain,    xcapt,  ycapt, kBlue,   0.30, "(A)", "  KKMCee CEEX2 ");
  PlotSame2(hst_WtFoam,    xcapt,  ycapt, kRed,    1.00, "(B)", "  Foam weight ");

  cFigWtMain->SaveAs("cFigWtMain.pdf");
}//FigWtMain()



///////////////////////////////////////////////////////////////////////////////////
void FigCosThe()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCosThe =========================== "<<endl;
  //
  TH1D *hst_CosTheta  = (TH1D*)DiskFileA->Get("hst_CosTheta");
  TH1D *hst_CosThOve  = (TH1D*)DiskFileA->Get("hst_CosThOve");
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCosThe = new TCanvas("cFigCosThe","cFigCosThe ",  gXcanv,  gYcanv,    500,  500);
  //                                  Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigCosThe->SetFillColor(10);
//  cFigCosThe->Divide( 2,  1);
  ////////////////////////////////////////////////////////////////////////////////
  //==========plot1==============
  cFigCosThe->cd(1);
  TH1D *HST; //
  HST = hst_CosTheta; //
//  HST = hst_CosThOve;
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->GetXaxis()->SetTitle("cos(#theta)");
  HST->GetXaxis()->SetLabelSize(0.04);
  HST->SetMinimum(0.0);
  HST->DrawCopy("h");
  CaptT->DrawLatex(0.06,0.95, "Events");
  double ycapt = 0.80; double xcapt=0.20;
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(xcapt,ycapt, "e^{+}e^{-} -> #nu#bar{#nu}(+n#gamma)");
  CaptT->DrawLatex(xcapt+0.40,ycapt,gTextEne);
  PlotSame2(hst_CosTheta,    xcapt, ycapt, kBlue,   0.30, "(A)", "  KKMCee CEEX2 ");
//  hst_CosThOve->Scale(100);
  PlotSame2(hst_CosThOve,    xcapt, ycapt, kRed,    0.80, "(B)", "  WT>WTmax");

  cFigCosThe->SaveAs("cFigCosThe.pdf");

}//FigCosThe


///////////////////////////////////////////////////////////////////////////////////
void FigNPhot()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNPhot =========================== "<<endl;
 //

  TH1D *hst_nPhAll     = (TH1D*)DiskFileA->Get("hst_nPhAll");
  TH1D *hst_nPhVis     = (TH1D*)DiskFileA->Get("hst_nPhVis");

  TH1D *hst_LnThPhAll  = (TH1D*)DiskFileA->Get("hst_LnThPhAll");
  TH1D *hst_LnThPhVis  = (TH1D*)DiskFileA->Get("hst_LnThPhVis");

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
  PlotSame2(hst_nPhAll, 0.20, ycapt, kBlue,  +2.0, "(a)", "All photons");
  PlotSame2(hst_nPhVis, 0.20, ycapt, kRed,   +1.0, "(b)", "Tagged photons");

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
  PlotSame2(hst_LnThPhAll,  0.20, ycapt, kBlue,    -3.0, "(a)", "All #gamma's");
  PlotSame2(hst_LnThPhVis,  0.20, ycapt, kRed,     -0.5, "(b)", "Tagged #gamma's");
  //

  //================================================
  if( g161GeVyes) cFigNPhot->SaveAs("cFigNPhot_161GeV.pdf");
  if( g125GeVyes) cFigNPhot->SaveAs("cFigNPhot_125GeV.pdf");
  if( g105GeVyes) cFigNPhot->SaveAs("cFigNPhot_105GeV.pdf");

}//FigNPhot


///////////////////////////////////////////////////////////////////////////////////
void FigVPhot2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVPhot2 =========================== "<<endl;
 //
  TH1D *hst_vtNuCeex2  = (TH1D*)DiskFileA->Get("hst_vtNuCeex2");  // Phot. untaged
  TH1D *hst_vaNuCeex2  = (TH1D*)DiskFileA->Get("hst_vaNuCeex2");  // Phot. tagged
 ////////////////////////////////////////
   TLatex *CaptT = new TLatex();
   CaptT->SetNDC(); // !!!
   CaptT->SetTextSize(0.04);
///////////////////////////////////////////////////////////////////////////////
    TCanvas *cFigVPhot2 = new TCanvas("cFigVPhot2","cFigVPhot2", gXcanv, gYcanv,    600, 600);
 //                                      Name    Title        xoff,yoff, WidPix,HeiPix
    gXcanv += 25, gYcanv += 25;
    cFigVPhot2->SetFillColor(10);
 ////////////////////////////////////////////////////////////////////////////////
    //====================plot1========================
    cFigVPhot2->cd();
    gPad->SetLogy(); // !!!!!!
    TH1D *Hst=hst_vtNuCeex2;
    Hst->SetStats(0);
    Hst->SetTitle(0);
    Hst->GetXaxis()->CenterTitle();
    Hst->GetXaxis()->SetTitleSize(0.04);
    Hst->GetXaxis()->SetTitle("v");

    Hst->SetMinimum(1e-5*Hst->GetMaximum());

    Hst->SetLineColor(kBlue);
    Hst->DrawCopy("h");

    CaptT->DrawLatex(0.05,0.94,"d#sigma/dv;     KKMC  e^{+}e^{-} -> #nu#bar{#nu}+n#gamma, #nu=#nu_{e}+#nu_{#mu}+#nu_{#tau}");
    double ycapt = 0.86;
    CaptT->DrawLatex(0.40, ycapt,gTextEne);
 //
    PlotSame2(hst_vtNuCeex2,  0.20, ycapt, kBlue, 0.4, "(a)", "#gamma's untagged, v=1-M^{2}_{#nu#bar{#nu}}/s");
    PlotSame2(hst_vaNuCeex2,  0.20, ycapt, kRed,  0.5, "(b)", "#gamma's tagged, v=E_{#gamma}/E_{beam}");
    //================================================
    if( g161GeVyes) cFigVPhot2->SaveAs("cFigVPhot2_161GeV.pdf");
    if( g125GeVyes) cFigVPhot2->SaveAs("cFigVPhot2_125GeV.pdf");
    if( g105GeVyes) cFigVPhot2->SaveAs("cFigVPhot2_105GeV.pdf");
}//FigVPhot2


///////////////////////////////////////////////////////////////////////////////////
void FigVPhot3()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVPhot3 =========================== "<<endl;
 //
  //TH1D *hst_vtNuCeex2  = (TH1D*)DiskFileA->Get("hst_vtNuCeex2");  // Phot. untaged
  //TH1D *hst_vaNuCeex2  = (TH1D*)DiskFileA->Get("hst_vaNuCeex2");  // Phot. tagged
  //
  TH1D *hst_vaNuMuCeex2  = (TH1D*)DiskFileA->Get("hst_vaNuMuCeex2");  // Phot. tagged
  TH1D *hst_vaNuElCeex2  = (TH1D*)DiskFileA->Get("hst_vaNuElCeex2");  // Phot. tagged


////////////////////////////////////////
   TLatex *CaptT = new TLatex();
   CaptT->SetNDC(); // !!!
   CaptT->SetTextSize(0.04);
///////////////////////////////////////////////////////////////////////////////
    TCanvas *cFigVPhot3 = new TCanvas("cFigVPhot3","cFigVPhot3", gXcanv, gYcanv,    600, 600);
 //                                      Name    Title        xoff,yoff, WidPix,HeiPix
    gXcanv += 25, gYcanv += 25;
    cFigVPhot3->SetFillColor(10);
 ////////////////////////////////////////////////////////////////////////////////
    //====================plot1========================
    cFigVPhot3->cd();
    gPad->SetLogy(); // !!!!!!
    TH1D *Hst=hst_vaNuElCeex2;
    Hst->SetStats(0);
    Hst->SetTitle(0);
    Hst->GetXaxis()->CenterTitle();
    Hst->GetXaxis()->SetTitleSize(0.04);
    Hst->GetXaxis()->SetTitle("v=E_{#gamma}/E_{beam}");

    Hst->SetMinimum(1e-4*Hst->GetMaximum());

    Hst->SetLineColor(kBlue);
    Hst->DrawCopy("h");

    CaptT->DrawLatex(0.05,0.94,"d#sigma/dv;     KKMC  e^{+}e^{-} -> #nu#bar{#nu}+n#gamma, #nu=#nu_{e},#nu_{#mu}");
    double ycapt = 0.85;
    CaptT->DrawLatex(0.40, ycapt,gTextEne);
 //
    PlotSame2(hst_vaNuElCeex2,  0.20, ycapt, kBlue, 0.4, "(a)", "#nu_{e}");
    PlotSame2(hst_vaNuMuCeex2,  0.20, ycapt, kRed,  0.5, "(b)", "#nu_{#mu}");
    //================================================
    if( g161GeVyes) cFigVPhot3->SaveAs("cFigVPhot3_161GeV.pdf");
    if( g125GeVyes) cFigVPhot3->SaveAs("cFigVPhot3_125GeV.pdf");
    if( g105GeVyes) cFigVPhot3->SaveAs("cFigVPhot3_105GeV.pdf");
}//FigVPhot3



///////////////////////////////////////////////////////////////////////////////////
void FigNuDif2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigNuDif2 =========================== "<<endl;
  ///
  TH1D *hst_vPhotNuel     = (TH1D*)DiskFileA->Get("hst_vPhotNuel");
  TH1D *hst_vPhotNumu     = (TH1D*)DiskFileA->Get("hst_vPhotNumu");
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
  RAT_NuelNumu->Divide(Hmu);
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
  Hst->SetMinimum(0);
  Hst->DrawCopy("h");
  ///
  double ycapt = 0.40;
  double vcapt = 1.0 - sqr(91.2/gCMSene);
  CaptT->DrawLatex(0.35, ycapt,gTextEne);
  PlotSame2(Hel,  0.20, ycapt, kBlue, vcapt-0.010, "(a)"," #nu = #nu_{el}");
  PlotSame2(Hmu,  0.20, ycapt, kRed,  vcapt+0.010, "(b)"," #nu = #nu_{#mu}");
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
  if( g125GeVyes) cNuDif2->SaveAs("cNuDif2_125GeV.pdf");
  if( g105GeVyes) cNuDif2->SaveAs("cNuDif2_105GeV.pdf");
}//FigNuDif2

///////////////////////////////////////////////////////////////////////////////////
void Table2(){
//------------------------------------------------------------------------
  cout<<" ========================= Table2 start=========================== "<<endl;
  TH1D *hst_vxNuElCeex2   = (TH1D*)DiskFileA->Get("hst_vxNuElCeex2");   //
  TH1D *hst_vxNuMuCeex2   = (TH1D*)DiskFileA->Get("hst_vxNuMuCeex2");   //
  TH1D *hst_vxNuTaCeex2   = (TH1D*)DiskFileA->Get("hst_vxNuTaCeex2");   //
  //
  TH1D  *hst_vxNuElcum   =HstCumul("hst_vxNuElcum", hst_vxNuElCeex2);
  TH1D  *hst_vxNuMucum   =HstCumul("hst_vxNuMucum", hst_vxNuMuCeex2);
  TH1D  *hst_vxNuTacum   =HstCumul("hst_vxNuTacum", hst_vxNuTaCeex2);

// Column captions
  int nPlt=3;   //
  Char_t *Capt[nPlt+1];
  for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
  strcpy(Capt[0],"{\\color{blue}$v_{\\max}$}");
  strcpy(Capt[1],"{\\color{blue} $\\nu_{e}$ }");
  strcpy(Capt[2],"{\\color{blue} $\\nu_{\\mu}$ }");
  strcpy(Capt[3],"{\\color{blue} $\\nu_{\\tau}$ }");

// multicolumn caption
  Char_t Mcapt[132];
  strcpy(Mcapt,"{\\color{red}$\\sigma(v_{\\max})$ [pb]}");

///************************************
  DiskFileTeX = fopen("Tab2.txp","w");
//************************************
// Initialization of the latex source file
  PlInitialize(DiskFileTeX, 2);

// pointers to histograms sigma(vmax)
  TH1D *iHst[nPlt+1];
  iHst[1]= hst_vxNuElcum;  //
  iHst[2]= hst_vxNuMucum;  //
  iHst[3]= hst_vxNuTacum;  //
  iHst[1]->Scale(1e3);    // nano- to pico-barns
  iHst[2]->Scale(1e3);    // nano- to pico-barns
  iHst[3]->Scale(1e3);    // nano- to pico-barns
  PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "B", 2, 2, 2); // for 100 bins
  PlTable2(-nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",10,90,20); // for 100 bins
  PlTable2(-nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "E",99,99, 2); // for 100 bins

// finalizing latex source file
  PlEnd(DiskFileTeX);
//************************************
  fclose(DiskFileTeX);
//************************************
  cout<<" ========================= Table2 end =========================== "<<endl;
}//Table2

///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

  DiskFileA = new TFile(FileA);

  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++

  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  int Nodes= HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
  if( fabs(gCMSene-161)<0.01) g161GeVyes = 1;
  if( fabs(gCMSene-125)<0.01) g125GeVyes = 1;
  if( fabs(gCMSene-105)<0.01) g105GeVyes = 1;
  //

  HistNormalize();     // Renormalization of MC histograms
  //========== PLOTTING ==========
  FigWtMain();
  FigCosThe();
  FigNPhot();
//
  FigVPhot2();
  FigVPhot3();

//  FigNuDif1();
  FigNuDif2();

  Table2();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA->ls();
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< "GeV"<<endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot. no. of events = "<<gNevTot<< endl;
//
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


